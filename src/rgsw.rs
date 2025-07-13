//! GSW13 protocol implementation as a Rust struct/class with unit tests.
//! This version uses i64 integers for clarity and testing. Adapt types as needed for your Fq/Poly definitions.
use crate::utils::{bit_decomp, bit_decomp_inv, flatten, powers_of_2};
use crate::{Fq, Poly, RLWE, gen_random_vector};
use ark_ff::AdditiveGroup;
use ark_ff::Zero;
use ark_ff::{Field, PrimeField};
use num_traits::pow;
use rand::Rng;
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use rand_distr::{Distribution, Normal, Uniform};
use std::ops::Neg;

/// Parameters for the GSW13 encryption scheme.
#[derive(Clone)]
pub struct GSW13Params {
    pub n: usize,     // Lattice dimension
    pub q: u64,       // Modulus
    pub l: usize,     // Bit-length for decomposition (log2(q))
    pub N: usize,     // (n+1)*l
    pub m: usize,     // Number of samples (m = O(n log q)), can be tuned
    pub std_dev: f64, // Noise standard deviation
}

/// RSGW struct: holds parameters, secret, and public key.
pub struct RGSW {
    pub params: GSW13Params,
    pub sk: Vec<Fq>,      // Secret key vector (1, -t1, ..., -tn)
    pub v: Vec<Fq>,       // PowersOf2(sk)
    pub pk: Vec<Vec<Fq>>, // Public key matrix A (m x (n+1))
}

#[derive(Clone)]
pub struct GSW13Ciphertext {
    pub C: Vec<Vec<Fq>>, // N x N matrix over Z_q
}

impl RGSW {
    /// Create parameters and generate secret/public keys
    pub fn new(n: usize, std_dev: f64) -> Self {
        let q = Fq::MODULUS.0[0] as u64; // Use Fq::MODULUS for the field modulus
        let l = (q as f64).log2().ceil() as usize;
        let N = (n + 1) * l;
        let m = n * l; // can be tuned; O(n log q)
        let params = GSW13Params {
            n,
            q,
            l,
            N,
            m,
            std_dev,
        };

        // Secret key generation: t ∈ Z_q^n
        let mut rng = ChaCha20Rng::seed_from_u64(42);
        let t: Vec<Fq> = (0..n).map(|_| Fq::from(rng.random_range(0..q))).collect();

        let mut sk = vec![Fq::ONE]; // s_0 = 1
        sk.extend(t.iter().map(|&ti| -ti)); // s_i = -t_i for i = 1..n

        let v = powers_of_2(&sk, l);

        // Public key generation
        let mut rng = ChaCha20Rng::seed_from_u64(43);
        // B ∈ Z_q^{m x n}
        let mut B = vec![vec![Fq::ZERO; n]; m];
        for i in 0..m {
            B[i] = gen_random_vector(n);
        }

        // e ∈ Z_q^m, sampled from discrete Gaussian
        let normal = Normal::new(0.0, std_dev).unwrap();
        let e: Vec<Fq> = (0..m)
            .map(|_| Fq::from(normal.sample(&mut rng).round() as i64))
            .collect();

        // t for multiplication (negated, as s = [1, -t])
        let t_plain: Vec<Fq> = sk[1..].iter().map(|x| -(*x)).collect();
        let b: Vec<Fq> = (0..m)
            .map(|i| {
                let mut acc = Fq::ZERO;
                for j in 0..n {
                    acc += B[i][j] * t_plain[j];
                }
                acc + e[i]
            })
            .collect();

        // A = [b | B] ∈ Z_q^{m x (n+1)}
        let mut A = vec![vec![Fq::ZERO; n + 1]; m];
        for i in 0..m {
            A[i][0] = b[i];
            for j in 0..n {
                A[i][j + 1] = B[i][j];
            }
        }

        Self {
            params,
            sk,
            v,
            pk: A,
        }
    }

    /// Encrypt a message μ ∈ Z_q as a RGSW ciphertext
    pub fn encrypt(&self, mu: Fq) -> GSW13Ciphertext {
        let p = &self.params;
        let mut rng = ChaCha20Rng::seed_from_u64(44);
        // R ∈ {0,1}^{N x m}
        let R: Vec<Vec<Fq>> = (0..p.N)
            .map(|_| {
                (0..p.m)
                    .map(|_| Fq::from(rng.random_range(0..=1)))
                    .collect()
            })
            .collect();

        // Matrix multiplication: R x A ∈ Z_q^{N x (n+1)}
        let mut RA = vec![vec![Fq::ZERO; self.pk[0].len()]; p.N];
        for i in 0..p.N {
            for j in 0..self.pk[0].len() {
                for k in 0..p.m {
                    RA[i][j] += R[i][k] * self.pk[k][j];
                }
            }
        }
        // Bit-decompose RA: N x (n+1) → N x ((n+1)*l)
        let bit_RA: Vec<Vec<Fq>> = RA.iter().map(|row| bit_decomp(row, p.l)).collect();

        // μ * I_N
        let mut MU_I_N = vec![vec![Fq::ZERO; p.N]; p.N];
        for i in 0..p.N {
            MU_I_N[i][i] = mu;
        }

        // C = Flatten(μ * I_N + BitDecomp(RA)) ∈ Z_q^{N x N}
        let mut C = vec![vec![Fq::ZERO; p.N]; p.N];
        for i in 0..p.N {
            for j in 0..p.N {
                C[i][j] = MU_I_N[i][j] + bit_RA[i][j];
            }
        }
        // Strong flatten to ensure boundedness
        for i in 0..p.N {
            C[i] = flatten(&C[i], p.l);
        }
        GSW13Ciphertext { C }
    }

    /// Decrypt a RGSW ciphertext using the secret key
    pub fn decrypt(&self, ct: &GSW13Ciphertext) -> Fq {
        let p = &self.params;
        // Compute C · v ∈ Z_q^N
        let mut x = vec![Fq::ZERO; p.N];
        for i in 0..p.N {
            for j in 0..p.N {
                x[i] += ct.C[i][j] * self.v[j];
            }
        }
        //let mut mu_est = x[0].into_bigint().0[0] % p.q;
        //if (mu_est as i64 - (p.q / 2) as i64).abs() < ((p.q / 4) as i64) {
        //    Fq::from(p.q / 2)
        //} else {
        //    Fq::ZERO
        //}
        // Recover μ by dividing first coefficient by v_0 = 1
        x[0] // In practice, may need rounding/modulo for noise
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_ff::UniformRand;

    #[test]
    fn test_gsw13_encrypt_decrypt() {
        // Use small n/q for test
        let n = 2;
        let std_dev = 1.2;
        let gsw = RGSW::new(n, std_dev);

        let mut rng = ark_std::rand::thread_rng();
        let mu = Fq::rand(&mut rng);
        let ct = gsw.encrypt(mu);
        let mu_rec = gsw.decrypt(&ct);

        // For small noise, should recover exactly
        assert_eq!(mu, mu_rec, "RGSW encryption/decryption failed for Fq");
    }

    #[test]
    fn test_bit_decomp_and_inv() {
        let q = 16;
        let l = 4;
        let v = vec![Fq::from(7u64), Fq::from(13u64)];
        let bits = bit_decomp(&v, l);
        let v2 = bit_decomp_inv(&bits, l);
        assert_eq!(v, v2, "BitDecomp/Inv failed");
    }

    #[test]
    fn test_flatten() {
        let q = 8;
        let l = 3;
        let v = vec![Fq::from(5u64), Fq::from(6u64)];
        let bits = bit_decomp(&v, l);
        let flat = flatten(&bits, l);
        assert_eq!(bits, flat, "Flatten failed");
    }

    #[test]
    fn test_powers_of_2() {
        let l = 3;
        let v = vec![Fq::from(2u64), Fq::from(3u64)];
        let pow2 = powers_of_2(&v, l);
        assert_eq!(
            pow2,
            vec![
                Fq::from(2),
                Fq::from(4),
                Fq::from(8),
                Fq::from(3),
                Fq::from(6),
                Fq::from(12)
            ]
        );
    }

    #[test]
    fn test_gsw13_encrypt_decrypt_binary() {
        let n = 2;
        let std_dev = 1.2;
        let gsw = RGSW::new(n, std_dev);

        for bit in vec![Fq::ZERO, Fq::ONE] {
            let ct = gsw.encrypt(bit);
            let bit_rec = gsw.decrypt(&ct);
            assert_eq!(bit, bit_rec, "GSW13 encrypt/decrypt failed (bit={})", bit);
        }
    }
}
