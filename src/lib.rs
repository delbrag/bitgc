// Rust version of Ring-GSW13 with full features: Cut-and-Choose, Key Switching, External Product, AND Gate Simulation, and Bootstrapping
use clap::Parser;

mod rlwe;
pub use rlwe::*;

mod rgsw;
use rgsw::*;

mod math;
pub use math::*;

#[derive(Parser, Debug)]
#[command(author, version, about = "Ring-GSW13 Full Demo", long_about = None)]
struct Args {
    #[arg(short, long, default_value_t = 128)]
    n: usize,
    #[arg(short, long, default_value_t = 1 << 14)]
    q: u64,
    #[arg(short, long, default_value_t = 3.2)]
    std_dev: f64,
    #[arg(long, default_value_t = 8)]
    cut_and_choose_rounds: usize,
}

//fn main() {
//    let args = Args::parse();
//    let p = compute_params(args);
//    println!("\n--- Ring-GSW13 Full Demo ---\nparams: n={}, q={}, std_dev={}, rounds={}\n", p.n, p.q, p.std_dev, p.cut_and_choose_rounds);
//
//    let s_in = rlwe_keygen(&p);
//    let s_out = rlwe_keygen(&p);
//
//    let ct = rlwe_encrypt(1, &s_in, &p);
//    println!("Decrypt before switch: {}", rlwe_decrypt(&ct, &s_in, &p));
//
//    let ksk = ks_keygen(&s_in, &s_out, &p);
//    let ct2 = key_switch(&ct, &ksk, &p);
//    println!("Decrypt after switch: {}", rlwe_decrypt(&ct2, &s_out, &p));
//
//    let rgsw = rgsw_encrypt(1, &s_in, &p);
//    let lwe = rlwe_encrypt(1, &s_in, &p);
//    let ext = external_product(&rgsw, &lwe, &p);
//    println!("External product: {}", rlwe_decrypt(&ext, &s_in, &p));
//
//    let and_ct = simulate_and_gate(&rgsw, &rgsw, &s_in, &p);
//    println!("Simulated AND gate: {}", rlwe_decrypt(&and_ct, &s_in, &p));
//
//    let ct3 = bootstrap(&ct, &ksk, &s_out, &p);
//    println!("Post bootstrap: {}", rlwe_decrypt(&ct3, &s_out, &p));
//
//    let cc_cts = cut_and_choose_encrypt(1, &s_in, &p);
//    if let Some(valid) = cut_and_choose_check(&cc_cts, &s_in, &p) {
//        println!("[Cut-and-Choose] Valid ciphertext: {}", rlwe_decrypt(&valid, &s_in, &p));
//    } else {
//        println!("[Cut-and-Choose] Verification failed.");
//    }
//}
