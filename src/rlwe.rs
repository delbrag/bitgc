use crate::{Fq, Poly};
use ark_ff::AdditiveGroup;
use ark_ff::PrimeField;
use rand::Rng;
use rand_chacha::ChaCha20Rng;
use rand_core::SeedableRng;
use rand_distr::{Distribution, Normal};
use std::ops::{Add, Mul, Sub};

/// RLWE encryption scheme over ring R_q = Z_q[x]/(x^n + 1)
/// Implements key generation, encryption, decryption, and bit encoding.
#[derive(Clone, Debug)]
pub struct RLWE {
    pub n: usize,         // Ring dimension (degree of modulus polynomial x^n + 1)
    pub q: u64,           // Prime modulus for coefficients
    pub std_dev: f64,     // Standard deviation of error distribution
    pub sk: Poly,         // Secret key polynomial s(x)
    pub pk: (Poly, Poly), // Public key (a(x), b(x) = a(x)*s(x) + e(x))
}

impl RLWE {
    /// Sample a small random polynomial with Gaussian-distributed coefficients
    /// Typically used for secret keys and error terms
    fn sample_small_poly(n: usize, std_dev: f64) -> Poly {
        let mut rng = ChaCha20Rng::seed_from_u64(42); // Use a fixed seed for reproducibility
        let normal = Normal::new(0.0, std_dev).unwrap();
        let coeffs = (0..n)
            .map(|_| {
                let x: f64 = normal.sample(&mut rng); // Sample from Gaussian
                let rounded = x.round() as i64; // Round to nearest integer
                Fq::from(rounded) // Map to finite field element
            })
            .collect();
        Poly { coeffs, n }
    }

    /// Sample a uniformly random polynomial in R_q
    /// Coefficients are sampled uniformly from [0, q)
    fn sample_uniform_poly(n: usize, q: u64) -> Poly {
        let mut rng = ChaCha20Rng::seed_from_u64(42); // Use a fixed seed for reproducibility
        let coeffs = (0..n).map(|_| Fq::from(rng.random_range(0..q))).collect();
        Poly { coeffs, n }
    }

    /// Generate a public/secret key pair for RLWE encryption
    /// Secret key: s ~ small noise
    /// Public key: a ← uniform, b = a*s + e
    pub fn keygen(n: usize, std_dev: f64) -> Self {
        let q = Fq::MODULUS.0[0] as u64; // get modulus from field
        let sk = Self::sample_small_poly(n, std_dev); // secret key s
        let a = Self::sample_uniform_poly(n, q); // random a(x)
        let e = Self::sample_small_poly(n, std_dev); // small error e
        let as_mul = a.mul(&sk); // compute a*s
        let b = as_mul.add(&e); // compute b = a*s + e
        RLWE {
            n,
            q,
            std_dev,
            sk: sk.clone(),
            pk: (a, b),
        }
    }

    /// Encode a bit vector {0,1}^n into a polynomial over R_q
    /// 0 → 0, 1 → q/2 (centered encoding)
    pub fn encode_bits(bits: &[u8], q: u64) -> Poly {
        let coeffs = bits
            .iter()
            .map(|b| if *b == 0 { Fq::ZERO } else { Fq::from(q / 2) })
            .collect();
        Poly {
            coeffs,
            n: bits.len(),
        }
    }

    /// Decode a polynomial back into a bit vector
    /// Each coefficient is interpreted as 0 or 1 depending on thresholding at q/4 and 3q/4
    pub fn decode_bits(poly: &Poly, q: u64) -> Vec<u8> {
        poly.coeffs
            .iter()
            .map(|c| {
                let v: u64 = (*c).into_bigint().0[0];
                if v > q / 4 && v < 3 * q / 4 { 1 } else { 0 }
            })
            .collect()
    }

    /// Encrypt a plaintext polynomial `m(x)` using the public key (a, b)
    /// Generates ephemeral secret r and error polynomials e1, e2
    /// Ciphertext: (c0 = a*r + e1, c1 = b*r + e2 + m)
    pub fn encrypt(&self, m: &Poly) -> (Poly, Poly) {
        let r = Self::sample_small_poly(self.n, self.std_dev); // ephemeral secret
        let e1 = Self::sample_small_poly(self.n, self.std_dev); // noise 1
        let e2 = Self::sample_small_poly(self.n, self.std_dev); // noise 2
        let a = &self.pk.0;
        let b = &self.pk.1;
        let c0 = a.mul(&r).add(&e1); // c0 = a*r + e1
        let mut c1 = b.mul(&r).add(&e2); // c1 = b*r + e2
        c1 = c1.add(m); // add message
        (c0, c1)
    }

    /// Decrypt a ciphertext (c0, c1) using the secret key s
    /// Computes c1 - s*c0 = m + noise
    pub fn decrypt(&self, c: &(Poly, Poly)) -> Poly {
        let (c0, c1) = c;
        let s_mul = c0.mul(&self.sk); // s * c0
        c1.sub(&s_mul) // c1 - s*c0 ≈ m
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;

    #[test]
    fn test_keygen_validity() {
        let rlwe = RLWE::keygen(8, 3.2);
        assert_eq!(rlwe.sk.coeffs.len(), 8);
        assert_eq!(rlwe.pk.0.coeffs.len(), 8);
        assert_eq!(rlwe.pk.1.coeffs.len(), 8);
    }

    #[test]
    fn test_encrypt_decrypt_random_bits() {
        let rlwe = RLWE::keygen(8, 3.2);
        let mut rng = ChaCha20Rng::seed_from_u64(42); // Use a fixed seed for reproducibility

        for i in 0..100 {
            // generate a 8-bit random 0-1 message
            let msg_bits: Vec<u8> = (0..8).map(|_| rng.random_range(0..2)).collect();

            let msg_poly = RLWE::encode_bits(&msg_bits, rlwe.q);

            let ct = rlwe.encrypt(&msg_poly);
            let decrypted_poly = rlwe.decrypt(&ct);
            let recovered_bits = RLWE::decode_bits(&decrypted_poly, rlwe.q);

            assert_eq!(msg_bits, recovered_bits, "Mismatch at iteration {}", i);
        }
    }

    #[test]
    fn test_encrypt_decrypt_all_zeros() {
        let rlwe = RLWE::keygen(8, 3.2);

        let msg_bits = vec![0; 8];
        let msg_poly = RLWE::encode_bits(&msg_bits, rlwe.q);

        let ct = rlwe.encrypt(&msg_poly);
        let decrypted_poly = rlwe.decrypt(&ct);
        let recovered_bits = RLWE::decode_bits(&decrypted_poly, rlwe.q);

        assert_eq!(msg_bits, recovered_bits);
    }

    #[test]
    fn test_encrypt_decrypt_all_ones() {
        let rlwe = RLWE::keygen(8, 3.2);

        let msg_bits = vec![1; 8];
        let msg_poly = RLWE::encode_bits(&msg_bits, rlwe.q);

        let ct = rlwe.encrypt(&msg_poly);
        let decrypted_poly = rlwe.decrypt(&ct);
        let recovered_bits = RLWE::decode_bits(&decrypted_poly, rlwe.q);

        assert_eq!(msg_bits, recovered_bits);
    }

    #[test]
    fn test_multiple_messages() {
        let rlwe = RLWE::keygen(8, 3.2);

        for i in 0..10 {
            let msg_bits = (0..8).map(|j| ((i + j) % 2) as u8).collect::<Vec<_>>();
            let msg_poly = RLWE::encode_bits(&msg_bits, rlwe.q);
            let ct = rlwe.encrypt(&msg_poly);
            let decrypted_poly = rlwe.decrypt(&ct);
            let recovered_bits = RLWE::decode_bits(&decrypted_poly, rlwe.q);
            assert_eq!(msg_bits, recovered_bits, "Mismatch at iteration {}", i);
        }
    }
}
