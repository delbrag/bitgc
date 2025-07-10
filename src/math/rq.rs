use std::fmt;

use ark_ff::AdditiveGroup;
use ark_ff::{Fp, MontBackend, MontConfig};
use rand_chacha::ChaCha20Rng;
use rand_distr::{Distribution, Normal};
use std::ops::{Add, Mul, Sub};

#[derive(MontConfig)]
#[modulus = "40961"]
#[generator = "3"]
pub struct FqConfig;
pub type Fq = Fp<MontBackend<FqConfig, 1>, 1>;

#[derive(Clone, PartialEq, Eq)]
pub struct Poly {
    pub coeffs: Vec<Fq>,
    pub n: usize, // Modulus degree (for x^n + 1)
}

impl Poly {
    pub fn new(coeffs: Vec<Fq>, n: usize) -> Self {
        let mut c = coeffs;
        c.resize(n, Fq::ZERO);
        Poly { coeffs: c, n }
    }

    pub fn zero(n: usize) -> Self {
        Poly::new(vec![Fq::ZERO; n], n)
    }

    pub fn from_u64(vec: &[u64], n: usize) -> Self {
        let coeffs = vec.iter().map(|&x| Fq::from(x)).collect();
        Poly::new(coeffs, n)
    }

    pub fn degree(&self) -> usize {
        self.n
    }

    pub fn sample_noise(n: usize, std_dev: f64, rng: &mut ChaCha20Rng) -> Self {
        let normal = Normal::new(0.0, std_dev).unwrap();
        let coeffs = (0..n)
            .map(|_| {
                let x: f64 = normal.sample(rng);
                let rounded = x.round() as i64;
                if rounded >= 0 {
                    Fq::from(rounded as u64)
                } else {
                    -Fq::from((-rounded) as u64)
                }
            })
            .collect();
        Poly::new(coeffs, n)
    }

    pub fn encode_msg(msg: &[u8], n: usize, q: u64) -> Self {
        let delta = q / 2;
        let coeffs = msg
            .iter()
            .map(|&b| if b == 1 { Fq::from(delta) } else { Fq::ZERO })
            .collect();
        Poly::new(coeffs, n)
    }

    pub fn decode_msg(&self, q: u64) -> Vec<u8> {
        let half_q = Fq::from(q / 2);
        self.coeffs
            .iter()
            .map(|&c| if c > half_q { 1 } else { 0 })
            .collect()
    }
}

impl fmt::Debug for Poly {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let terms: Vec<String> = self
            .coeffs
            .iter()
            .enumerate()
            .filter(|&(_, &c)| c != Fq::ZERO)
            .map(|(i, c)| format!("{}x^{}", c, i))
            .collect();
        write!(f, "{}", terms.join(" + "))
    }
}

impl Add for &Poly {
    type Output = Poly;

    fn add(self, rhs: Self) -> Self::Output {
        let coeffs = self
            .coeffs
            .iter()
            .zip(rhs.coeffs.iter())
            .map(|(a, b)| *a + *b)
            .collect();
        Poly::new(coeffs, self.n)
    }
}

impl std::ops::Add<&Poly> for Poly {
    type Output = Poly;
    fn add(self, rhs: &Poly) -> Self::Output {
        let coeffs = self
            .coeffs
            .iter()
            .zip(rhs.coeffs.iter())
            .map(|(a, b)| *a + *b)
            .collect();
        Poly { coeffs, n: self.n }
    }
}

impl Sub for &Poly {
    type Output = Poly;

    fn sub(self, rhs: Self) -> Self::Output {
        let coeffs = self
            .coeffs
            .iter()
            .zip(rhs.coeffs.iter())
            .map(|(a, b)| *a - *b)
            .collect();
        Poly::new(coeffs, self.n)
    }
}

impl Mul for &Poly {
    type Output = Poly;

    fn mul(self, rhs: Self) -> Self::Output {
        let mut res = vec![Fq::ZERO; 2 * self.n - 1];
        for i in 0..self.n {
            for j in 0..self.n {
                res[i + j] += self.coeffs[i] * rhs.coeffs[j];
            }
        }
        // Reduce mod x^n + 1
        for i in self.n..2 * self.n - 1 {
            let tmp = res[i];
            res[i - self.n] -= tmp;
        }
        res.truncate(self.n);
        Poly::new(res, self.n)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_poly_ops() {
        let n = 4;
        let a = Poly::from_u64(&[1, 2, 3, 4], n);
        let b = Poly::from_u64(&[4, 3, 2, 1], n);

        let sum = &a + &b;
        let expected = Poly::from_u64(&[5, 5, 5, 5], n);
        assert_eq!(sum, expected);

        let diff = &a - &b;
        let expected = Poly::from_u64(&[40958, 40960, 1, 3], n); // mod 40961
        assert_eq!(diff, expected);

        let prod = &a * &b;
        println!("a * b = {:?}", prod); // just show the poly for now
    }
}
