use crate::{Fq, Poly, RLWE};
use ark_ff::Zero;
use rand::Rng;
use ark_ff::AdditiveGroup;

/// GSW ciphertext is a matrix of RLWE encryptions
#[derive(Debug, Clone)]
pub struct GSW {
    pub n: usize,
    pub l: usize, // bit decomposition length
    pub q: u64,
    pub std_dev: f64,
    pub sk_rlwe: RLWE,
    pub ek: Vec<(Poly, Poly)>, // Encryption of powers of 2: Enc(s * 2^i)
}

impl GSW {
    /// Generate a GSW instance and public key encryptions for key switching
    pub fn keygen(n: usize, q: u64, std_dev: f64, l: usize, rng: &mut impl Rng) -> Self {
        let sk_rlwe = RLWE::keygen(n, q, std_dev, rng);

        let mut ek = vec![];
        for i in 0..l {
            let scalar = Fq::from(1u64 << i); // 2^i
            let scaled = Poly::new(sk_rlwe.sk.coeffs.iter().map(|c| *c * scalar).collect(), n);
            let enc = sk_rlwe.encrypt(&scaled, rng);
            ek.push(enc);
        }

        GSW {
            n,
            l,
            q,
            std_dev,
            sk_rlwe,
            ek,
        }
    }

    /// Encrypt a message m ∈ {0,1} under GSW as a binary-decomposed matrix
    pub fn encrypt(&self, m: u8, rng: &mut impl Rng) -> Vec<(Poly, Poly)> {
        let mut ct = vec![];
        for i in 0..self.l {
            let bit = if ((m >> i) & 1) == 1 { Fq::from(self.q / 2) } else { Fq::ZERO };
            let msg = Poly::new(vec![bit; self.n], self.n); // constant poly
            let c = self.sk_rlwe.encrypt(&msg, rng);
            ct.push(c);
        }
        ct
    }

    /// External product: computes GSW ciphertext ⊗ RLWE ciphertext
    pub fn external_product(&self, gsw_ct: &[ (Poly, Poly) ], rlwe_ct: &(Poly, Poly)) -> (Poly, Poly) {
        let mut acc0 = Poly::zero(self.n);
        let mut acc1 = Poly::zero(self.n);

        for i in 0..self.l {
            let (a_i, b_i) = &gsw_ct[i];
            let (c0, c1) = rlwe_ct;
            let scaled_c0 = Poly::new(c0.coeffs.iter().map(|x| *x * Fq::from(1 << i)).collect(), self.n);
            let scaled_c1 = Poly::new(c1.coeffs.iter().map(|x| *x * Fq::from(1 << i)).collect(), self.n);
            acc0 = &acc0 + &(&scaled_c0 * a_i);
            acc1 = &acc1 + &(&scaled_c1 * b_i);
        }

        (acc0, acc1)
    }

    /// Key switching using encrypted secret key bits
    pub fn key_switch(&self, ct: &(Poly, Poly)) -> (Poly, Poly) {
        let (c0, c1) = ct;
        let mut new_c0 = Poly::zero(self.n);
        let mut new_c1 = Poly::zero(self.n);

        for i in 0..self.l {
            let scalar = Fq::from(1u64 << i);
            let mut s_coeffs = c0.coeffs.iter().map(|c| *c * scalar).collect::<Vec<_>>();
            let scaled = Poly::new(s_coeffs, self.n);
            let (ek0, ek1) = &self.ek[i];
            new_c0 = &new_c0 + &(&scaled * ek0);
            new_c1 = &new_c1 + &(&scaled * ek1);
        }

        (new_c0, new_c1)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{RLWE, Poly};
    use rand_chacha::ChaCha20Rng;
    use rand::SeedableRng;

    #[test]
    fn test_gsw_keygen() {
        let mut rng = ChaCha20Rng::seed_from_u64(42);
        let gsw = GSW::keygen(8, 40961, 3.2, 8, &mut rng);
        assert_eq!(gsw.ek.len(), gsw.l);
        assert_eq!(gsw.sk_rlwe.n, 8);
    }

    #[test]
    fn test_gsw_encrypt_bit() {
        let mut rng = ChaCha20Rng::seed_from_u64(123);
        let gsw = GSW::keygen(8, 40961, 3.2, 8, &mut rng);
        let gsw_ct = gsw.encrypt(1, &mut rng);
        assert_eq!(gsw_ct.len(), gsw.l);
    }

    #[test]
    fn test_external_product_preserves_bit() {
        let mut rng = ChaCha20Rng::seed_from_u64(2024);
        let gsw = GSW::keygen(8, 40961, 3.2, 8, &mut rng);

        // Encrypt bit 1 with GSW
        let gsw_ct = gsw.encrypt(1, &mut rng);

        // Encrypt message under RLWE
        let msg = Poly::encode_msg(&[1, 0, 0, 0, 0, 0, 0, 0], 8, 40961);
        let rlwe_ct = gsw.sk_rlwe.encrypt(&msg, &mut rng);

        // Perform external product
        let result = gsw.external_product(&gsw_ct, &rlwe_ct);
        let decrypted = gsw.sk_rlwe.decrypt(&result);
        let bits = Poly::decode_msg(&decrypted, 40961);

        // Only the first bit should be 1
        assert_eq!(bits[0], 1);
        assert!(bits[1..].iter().all(|&b| b == 0));
    }

    #[test]
    fn test_key_switch_preserves_decryption() {
        let mut rng = ChaCha20Rng::seed_from_u64(8888);
        let gsw = GSW::keygen(8, 40961, 3.2, 8, &mut rng);

        let msg = Poly::encode_msg(&[1, 0, 0, 0, 0, 0, 0, 0], 8, 40961);
        let ct = gsw.sk_rlwe.encrypt(&msg, &mut rng);

        // Apply key switching to re-encrypt the ciphertext
        let switched = gsw.key_switch(&ct);

        // Decrypt using RLWE
        let decrypted = gsw.sk_rlwe.decrypt(&switched);
        let bits = Poly::decode_msg(&decrypted, 40961);
        assert_eq!(bits[0], 1);
        assert!(bits[1..].iter().all(|&b| b == 0));
    }
}
