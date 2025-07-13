use crate::math::Fq;
use ark_ff::AdditiveGroup;
use ark_ff::{BigInteger, PrimeField};
use ark_ff::{Field, Zero};

/// Bit-decompose a vector `v` ∈ Z_q^n into binary representation (least to most significant).
/// Each coefficient is decomposed into `l` bits (low to high).
pub fn bit_decomp(v: &[Fq], l: usize) -> Vec<Fq> {
    assert!(
        v.iter().all(|x| x.into_bigint().0[0] < (1u64 << l)),
        "element exceeds allowed bit width"
    );
    let mut result = Vec::with_capacity(v.len() * l);
    for &coeff in v {
        let bigint = coeff.into_bigint();
        let mut acc = bigint.to_bits_le(); // little-endian bit vector
        acc.resize(l, false); // ensure l bits
        for i in 0..l {
            result.push(if acc[i] { Fq::ONE } else { Fq::ZERO });
        }
    }
    result
}

/// Inverse of bit-decomposition: reconstruct vector of integers from bit-decomposed form.
pub fn bit_decomp_inv(bits: &[Fq], l: usize) -> Vec<Fq> {
    let n = bits.len() / l;
    let mut result = Vec::with_capacity(n);

    for i in 0..n {
        let mut acc = Fq::ZERO;
        let mut base = Fq::ONE;
        for j in 0..l {
            let bit = bits[i * l + j];
            acc += base * bit;
            base = base.double(); // base *= 2
        }
        result.push(acc);
    }

    result
}

/// Flatten: BitDecomp(BitDecompInv(v)) – ensures strong boundedness.
pub fn flatten(v: &[Fq], l: usize) -> Vec<Fq> {
    assert!(v.len() % l == 0, "Length must be a multiple of l");
    bit_decomp(&bit_decomp_inv(v, l), l)
}

/// Powersof2: expand [v_0, ..., v_{n-1}] into [v_0, 2v_0, ..., 2^{l-1}v_0, ..., v_{n-1}, ..., 2^{l-1}v_{n-1}] by scaling.
pub fn powers_of_2(v: &[Fq], l: usize) -> Vec<Fq> {
    let mut result = Vec::with_capacity(v.len() * l);
    for &coeff in v {
        let mut base = Fq::ONE;
        for i in 0..l {
            result.push(coeff * base);
            base = base.double(); // multiply by 2
        }
    }
    result
}

#[cfg(test)]
mod tests {
    use num_traits::pow;
    use rand_core::le;
    use sha2::digest::typenum::bit;

    use super::*;
    use crate::math::{Fq, gen_random_vector};

    #[test]
    fn test_bit_decomp_and_inverse() {
        let l = 32;
        let n = 4;
        let v: Vec<Fq> = gen_random_vector(n);

        let bits = bit_decomp(&v, l);
        assert_eq!(bits.len(), l * n, "Bit decomposition length incorrect");

        let recovered = bit_decomp_inv(&bits, l);
        assert_eq!(v.len(), recovered.len(), "Length mismatch");

        for (i, (orig, recon)) in v.iter().zip(recovered.iter()).enumerate() {
            assert_eq!(
                orig, recon,
                "Mismatch at index {}: {:?} != {:?}",
                i, orig, recon
            );
        }
    }

    #[test]
    fn test_bit_decomp_zero_vector() {
        // special case
        let l = 2;
        let n = 4;
        let v = vec![Fq::ZERO; n];

        let bits = bit_decomp(&v, l);
        assert!(
            bits.iter().all(|b| b.is_zero()),
            "All bits of zero should be zero"
        );

        let recovered = bit_decomp_inv(&bits, l);
        assert_eq!(v, recovered, "Inverse of zero vector failed");
    }

    #[test]
    fn test_bit_decomp_all_ones() {
        // special case
        let l = 3;
        let n = 2;
        let one_val = Fq::from((1u64 << l) - 1);
        let v = vec![one_val; n];

        let bits = bit_decomp(&v, l);
        assert!(
            bits.iter().all(|b| *b == Fq::ONE),
            "Bits of all-ones should be 1"
        );

        let recovered = bit_decomp_inv(&bits, l);
        for x in &recovered {
            assert_eq!(*x, one_val);
        }
    }

    // Verify equation <BitDecomp(a), Powersof2(b)> = <a, b>.
    #[test]
    fn test_bit_decomp_powers_of_2_inner_product() {
        let l = 32;
        let n = 2;
        let v = gen_random_vector(n);
        let u = gen_random_vector(n);
        let decomp_v = bit_decomp(&v, l);
        let pow_u = powers_of_2(&u, l);
        let lhs: Fq = decomp_v
            .iter()
            .zip(pow_u.iter())
            .map(|(a, b)| *a * *b)
            .sum();
        let rhs: Fq = v.iter().zip(u.iter()).map(|(a, b)| *a * *b).sum();
        assert_eq!(lhs, rhs, "Bit decomposition inner product mismatch");
    }

    // For N-dim a, <a, powers_of_2(b)> = <bit_decomp_inv(a), b> = <flatten(bit_decomp(a)), powers_of_2(b)>
    #[test]
    fn test_flatten_preserves_inner_product() {
        let l = 32;
        let n = 2;

        let v = gen_random_vector(n);
        let a = bit_decomp(&v, l);

        let b = gen_random_vector(n);

        let bit_decomp_inv = bit_decomp_inv(&a, l);
        let pow_u = powers_of_2(&b, l);

        let lhs: Fq = a.iter().zip(pow_u.iter()).map(|(a, b)| *a * *b).sum();
        let rhs: Fq = bit_decomp_inv
            .iter()
            .zip(b.iter())
            .map(|(a, b)| *a * *b)
            .sum();
        assert_eq!(lhs, rhs, "Flatten inner product mismatch");

        let flattened = flatten(&a, l);
        let lhs_flat: Fq = flattened
            .iter()
            .zip(pow_u.iter())
            .map(|(a, b)| *a * *b)
            .sum();
        assert_eq!(lhs_flat, rhs, "Flattened inner product mismatch");
    }

    #[test]
    fn test_powers_of_2_correctness() {
        let l = 4;
        let v = vec![Fq::from(3u64)];
        let expected = vec![Fq::from(3), Fq::from(6), Fq::from(12), Fq::from(24)];
        let result = powers_of_2(&v, l);
        assert_eq!(result, expected);
    }
}
