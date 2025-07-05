// Rust version of Ring-GSW13 (parameterized with CLI configuration via clap)
use rand::Rng;
use rand_distr::{Distribution, Normal};
use clap::Parser;

#[derive(Parser, Debug)]
#[command(author, version, about = "Ring-GSW13 Demo with CLI-Configurable Security Level", long_about = None)]
struct Args {
    #[arg(short, long, default_value_t = 128)]
    n: usize,

    #[arg(short, long, default_value_t = 1 << 14)]
    q: u64,

    #[arg(short, long, default_value_t = 3.2)]
    std_dev: f64,
}

struct Params {
    n: usize, // maximal polynomial degree, x^n+1 
    q: u64, // Z_q 
    std_dev: f64, // noise standard derives  
    delta: u64, // extension factor 
    l: usize, // l = log2(q)
}

fn compute_params(args: Args) -> Params {
    let delta = args.q / 4;
    let l = (64 - args.q.leading_zeros()) as usize;
    Params { n: args.n, q: args.q, std_dev: args.std_dev, delta, l }
}

fn poly_add(a: &[u64], b: &[u64], p: &Params) -> Vec<u64> {
    a.iter().zip(b).map(|(&x, &y)| (x + y) % p.q).collect()
}

fn poly_sub(a: &[u64], b: &[u64], p: &Params) -> Vec<u64> {
    a.iter().zip(b).map(|(&x, &y)| (x + p.q - y) % p.q).collect()
}

fn poly_mul(a: &[u64], b: &[u64], p: &Params) -> Vec<u64> {
    let mut res = vec![0u64; 2 * p.n - 1];
    for i in 0..p.n {
        for j in 0..p.n {
            res[i + j] = (res[i + j] + a[i] * b[j]) % p.q;
        }
    }
    for i in p.n..2 * p.n - 1 {
        res[i - p.n] = (res[i - p.n] + p.q - res[i]) % p.q;
    }
    res.truncate(p.n);
    res
}

// s
fn rlwe_keygen(p: &Params) -> Vec<u64> {
    (0..p.n).map(|_| rand::thread_rng().gen_range(0..2)).collect()
}

// b = ⟨a, s⟩ + mΔ + e (mod q) => (a, b) 
fn rlwe_encrypt(m: u64, s: &[u64], p: &Params) -> (Vec<u64>, Vec<u64>) {
    let mut rng = rand::thread_rng();
    let a: Vec<u64> = (0..p.n).map(|_| rng.gen_range(0..p.q)).collect();
    let normal = Normal::new(0.0, p.std_dev).unwrap();
    let e: Vec<i64> = (0..p.n).map(|_| normal.sample(&mut rng) as i64).collect();

    let mut m_poly = vec![0u64; p.n];
    m_poly[0] = (m * p.delta) % p.q;

    let as_mul = poly_mul(&a, s, p);
    let mut b = vec![0u64; p.n];
    for i in 0..p.n {
        let sum = as_mul[i] as i64 + m_poly[i] as i64 + e[i];
        b[i] = ((sum % p.q as i64 + p.q as i64) % p.q as i64) as u64;
    }

    (a, b)
}

fn rlwe_decrypt(ct: &(Vec<u64>, Vec<u64>), s: &[u64], p: &Params) -> u64 {
    let (a, b) = ct;
    let as_mul = poly_mul(&a, s, p);
    let mut phase = vec![0u64; p.n];
    for i in 0..p.n {
        phase[i] = (b[i] + p.q - as_mul[i]) % p.q;
    }
    ((phase[0] + p.delta / 2) / p.delta) % 2
}

fn rgsw_encrypt(m: u64, s: &[u64], p: &Params) -> Vec<(Vec<u64>, Vec<u64>)> {
    let mut rng = rand::thread_rng();
    let normal = Normal::new(0.0, p.std_dev).unwrap();
    let mut result = vec![];

    for i in 0..p.l {
        let mut m_poly = vec![0u64; p.n];
        m_poly[0] = (m * (1 << i)) % p.q;
        let a: Vec<u64> = (0..p.n).map(|_| rng.gen_range(0..p.q)).collect();
        let e: Vec<i64> = (0..p.n).map(|_| normal.sample(&mut rng) as i64).collect();
        let as_mul = poly_mul(&a, s, p);
        let mut b = vec![0u64; p.n];
        for j in 0..p.n {
            let sum = as_mul[j] as i64 + m_poly[j] as i64 + e[j];
            b[j] = ((sum % p.q as i64 + p.q as i64) % p.q as i64) as u64;
        }
        result.push((a, b));
    }
    result
}

fn ks_keygen(s_in: &[u64], s_out: &[u64], p: &Params) -> Vec<Vec<(Vec<u64>, Vec<u64>)>> {
    let mut result = vec![vec![]; p.n];
    let normal = Normal::new(0.0, p.std_dev).unwrap();
    let mut rng = rand::thread_rng();
    for i in 0..p.n {
        for j in 0..p.l {
            let val = (s_in[i] * (1 << j)) % p.q;
            let mut m_poly = vec![0u64; p.n];
            m_poly[0] = val;
            let a: Vec<u64> = (0..p.n).map(|_| rng.gen_range(0..p.q)).collect();
            let e: Vec<i64> = (0..p.n).map(|_| normal.sample(&mut rng) as i64).collect();
            let as_mul = poly_mul(&a, s_out, p);
            let mut b = vec![0u64; p.n];
            for k in 0..p.n {
                let sum = as_mul[k] as i64 + m_poly[k] as i64 + e[k];
                b[k] = ((sum % p.q as i64 + p.q as i64) % p.q as i64) as u64;
            }
            result[i].push((a, b));
        }
    }
    result
}

fn key_switch(ct: &(Vec<u64>, Vec<u64>), ksk: &[Vec<(Vec<u64>, Vec<u64>)>], p: &Params) -> (Vec<u64>, Vec<u64>) {
    let (a_in, b_in) = ct;
    let mut a_out = vec![0u64; p.n];
    let mut b_out = b_in.clone();
    for i in 0..p.n {
        for j in 0..p.l {
            let bit = (a_in[i] >> j) & 1;
            if bit == 1 {
                let (ref a_ksk, ref b_ksk) = ksk[i][j];
                a_out = poly_sub(&a_out, a_ksk, p);
                for k in 0..p.n {
                    b_out[k] = (b_out[k] + p.q - b_ksk[k]) % p.q;
                }
            }
        }
    }
    (a_out, b_out)
}

fn external_product(rgsw: &[(Vec<u64>, Vec<u64>)], lwe_ct: &(Vec<u64>, Vec<u64>), p: &Params) -> (Vec<u64>, Vec<u64>) {
    let (a_lwe, _b_lwe) = lwe_ct;
    let mut acc_a = vec![0u64; p.n];
    let mut acc_b = vec![0u64; p.n];
    for i in 0..p.l {
        let bit = (a_lwe[0] >> i) & 1;
        if bit == 1 {
            let (ref a_gsw, ref b_gsw) = rgsw[i];
            acc_a = poly_add(&acc_a, a_gsw, p);
            acc_b = poly_add(&acc_b, b_gsw, p);
        }
    }
    (acc_a, acc_b)
}

fn main() {
    let args = Args::parse();
    let params = compute_params(args);

    println!("--- Ring-GSW13 Demo with CLI-Configurable Security Level ---");
    println!("Using parameters: n = {}, q = {}, std_dev = {}", params.n, params.q, params.std_dev);

    let s_in = rlwe_keygen(&params);
    let s_out = rlwe_keygen(&params);
    let m = 1;

    let ct_in = rlwe_encrypt(m, &s_in, &params);
    println!("Decrypted under s_in: {}", rlwe_decrypt(&ct_in, &s_in, &params));

    let ksk = ks_keygen(&s_in, &s_out, &params);
    let ct_out = key_switch(&ct_in, &ksk, &params);
    println!("Decrypted after key-switch under s_out: {}", rlwe_decrypt(&ct_out, &s_out, &params));

    let gsw_ct = rgsw_encrypt(m, &s_in, &params);
    let lwe_ct = rlwe_encrypt(1, &s_in, &params);
    let ext = external_product(&gsw_ct, &lwe_ct, &params);
    println!("External product decryption: {}", rlwe_decrypt(&ext, &s_in, &params));
}