//! ML-KEM against the NIST ACVP-Server FIPS 203 vectors
//! (keyGen, encapsulation, decapsulation, encapsulation/decapsulation key checks).

mod common;
use common::*;
use fips_pqc::ml_kem::{ml_kem_1024, ml_kem_512, ml_kem_768};

/// Dispatches a block of code over the three parameter sets.
macro_rules! with_param_set {
    ($ps:expr, |$m:ident| $body:block) => {
        match $ps {
            "ML-KEM-512" => {
                use ml_kem_512 as $m;
                $body
            }
            "ML-KEM-768" => {
                use ml_kem_768 as $m;
                $body
            }
            "ML-KEM-1024" => {
                use ml_kem_1024 as $m;
                $body
            }
            other => panic!("unknown parameter set {other}"),
        }
    };
}

#[test]
fn acvp_keygen() {
    let recs = load("mlkem_keygen.txt");
    assert_eq!(recs.len(), 75);
    for r in &recs {
        let mut seed = [0u8; 64];
        seed[..32].copy_from_slice(&r.hex("d"));
        seed[32..].copy_from_slice(&r.hex("z"));
        with_param_set!(r.get("parameterSet"), |m| {
            let (ek, dk) = m::generate_key_pair_from_seed(&seed);
            assert_eq!(ek.as_bytes().to_vec(), r.hex("ek"), "tcId {}", r.id());
            assert_eq!(dk.as_bytes().to_vec(), r.hex("dk"), "tcId {}", r.id());
        });
    }
}

#[test]
fn acvp_encap_decap_and_key_checks() {
    let recs = load("mlkem_encap_decap.txt");
    assert_eq!(recs.len(), 165);
    let mut counts = std::collections::HashMap::new();
    for r in &recs {
        let func = r.get("function");
        *counts.entry(func.to_string()).or_insert(0) += 1;
        with_param_set!(r.get("parameterSet"), |m| {
            match func {
                "encapsulation" => {
                    let ek = m::EncapsulationKey::from_bytes(&r.hex("ek")).expect("valid ek");
                    let (c, k) = ek.encapsulate_internal(&arr(&r.hex("m")));
                    assert_eq!(c.as_bytes().to_vec(), r.hex("c"), "tcId {}", r.id());
                    assert_eq!(k.as_bytes().to_vec(), r.hex("k"), "tcId {}", r.id());
                    // The matching decapsulation key must recover the same secret.
                    let dk = m::DecapsulationKey::from_bytes(&r.hex("dk")).expect("valid dk");
                    assert_eq!(dk.decapsulate(&c), k);
                }
                "decapsulation" => {
                    let dk = m::DecapsulationKey::from_bytes(&r.hex("dk")).expect("valid dk");
                    let c = m::Ciphertext::from_bytes(&r.hex("c")).expect("valid ct length");
                    let k = dk.decapsulate(&c);
                    assert_eq!(k.as_bytes().to_vec(), r.hex("k"), "tcId {} ({})", r.id(), r.get("reason"));
                }
                "encapsulationKeyCheck" => {
                    let ok = m::EncapsulationKey::from_bytes(&r.hex("ek")).is_ok();
                    assert_eq!(ok, r.bool("testPassed"), "tcId {} ({})", r.id(), r.get("reason"));
                }
                "decapsulationKeyCheck" => {
                    let ok = m::DecapsulationKey::from_bytes(&r.hex("dk")).is_ok();
                    assert_eq!(ok, r.bool("testPassed"), "tcId {} ({})", r.id(), r.get("reason"));
                }
                other => panic!("unknown function {other}"),
            }
        });
    }
    assert_eq!(counts["encapsulation"], 75);
    assert_eq!(counts["decapsulation"], 30);
    assert_eq!(counts["encapsulationKeyCheck"], 30);
    assert_eq!(counts["decapsulationKeyCheck"], 30);
}
