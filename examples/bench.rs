//! Rough timings for every operation (no external benchmark harness).
//!
//! Run with: `cargo run --release --example bench`

use fips_pqc::ml_dsa::{ml_dsa_44, ml_dsa_65, ml_dsa_87};
use fips_pqc::ml_kem::{ml_kem_1024, ml_kem_512, ml_kem_768};
use std::hint::black_box;
use std::time::Instant;

fn time<F: FnMut()>(label: &str, iters: u32, mut f: F) {
    let start = Instant::now();
    for _ in 0..iters {
        f();
    }
    let per = start.elapsed() / iters;
    println!("{label:<28} {:>9.1} µs", per.as_secs_f64() * 1e6);
}

macro_rules! bench_kem {
    ($m:ident, $name:literal) => {{
        let (ek, dk) = $m::generate_key_pair_from_seed(&[1u8; 64]);
        let (ct, _) = ek.encapsulate_internal(&[2u8; 32]);
        time(concat!($name, " keygen"), 2000, || {
            black_box($m::generate_key_pair_from_seed(black_box(&[1u8; 64])));
        });
        time(concat!($name, " encaps"), 2000, || {
            black_box(ek.encapsulate_internal(black_box(&[2u8; 32])));
        });
        time(concat!($name, " decaps"), 2000, || {
            black_box(dk.decapsulate(black_box(&ct)));
        });
    }};
}

macro_rules! bench_dsa {
    ($m:ident, $name:literal) => {{
        let (vk, sk) = $m::generate_key_pair_from_seed(&[1u8; 32]);
        let sig = sk.sign_deterministic(b"benchmark message", b"").unwrap();
        time(concat!($name, " keygen"), 500, || {
            black_box($m::generate_key_pair_from_seed(black_box(&[1u8; 32])));
        });
        let mut i = 0u32;
        time(concat!($name, " sign"), 500, || {
            i += 1;
            black_box(sk.sign_with_rnd(b"benchmark message", b"", &[i as u8; 32]).unwrap());
        });
        time(concat!($name, " verify"), 500, || {
            black_box(vk.verify(b"benchmark message", b"", black_box(&sig)).unwrap());
        });
    }};
}

fn main() {
    bench_kem!(ml_kem_512, "ML-KEM-512");
    bench_kem!(ml_kem_768, "ML-KEM-768");
    bench_kem!(ml_kem_1024, "ML-KEM-1024");
    bench_dsa!(ml_dsa_44, "ML-DSA-44");
    bench_dsa!(ml_dsa_65, "ML-DSA-65");
    bench_dsa!(ml_dsa_87, "ML-DSA-87");
}
