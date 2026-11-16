//! ML-DSA against the NIST ACVP-Server FIPS 204 vectors: keyGen, sigGen and
//! sigVer, covering pure / pre-hash / internal / external-μ interfaces and
//! deterministic / hedged signing.

mod common;
use common::*;
use fips_pqc::ml_dsa::{ml_dsa_44, ml_dsa_65, ml_dsa_87, PreHashAlgorithm};
use fips_pqc::Error;

macro_rules! with_param_set {
    ($ps:expr, |$m:ident| $body:block) => {
        match $ps {
            "ML-DSA-44" => {
                use ml_dsa_44 as $m;
                $body
            }
            "ML-DSA-65" => {
                use ml_dsa_65 as $m;
                $body
            }
            "ML-DSA-87" => {
                use ml_dsa_87 as $m;
                $body
            }
            other => panic!("unknown parameter set {other}"),
        }
    };
}

fn hash_alg(name: &str) -> PreHashAlgorithm {
    match name {
        "SHA2-224" => PreHashAlgorithm::Sha2_224,
        "SHA2-256" => PreHashAlgorithm::Sha2_256,
        "SHA2-384" => PreHashAlgorithm::Sha2_384,
        "SHA2-512" => PreHashAlgorithm::Sha2_512,
        "SHA2-512/224" => PreHashAlgorithm::Sha2_512_224,
        "SHA2-512/256" => PreHashAlgorithm::Sha2_512_256,
        "SHA3-224" => PreHashAlgorithm::Sha3_224,
        "SHA3-256" => PreHashAlgorithm::Sha3_256,
        "SHA3-384" => PreHashAlgorithm::Sha3_384,
        "SHA3-512" => PreHashAlgorithm::Sha3_512,
        "SHAKE-128" => PreHashAlgorithm::Shake128,
        "SHAKE-256" => PreHashAlgorithm::Shake256,
        other => panic!("unknown hash {other}"),
    }
}

#[test]
fn acvp_keygen() {
    let recs = load("mldsa_keygen.txt");
    assert_eq!(recs.len(), 75);
    for r in &recs {
        let seed = arr::<32>(&r.hex("seed"));
        with_param_set!(r.get("parameterSet"), |m| {
            let (pk, sk) = m::generate_key_pair_from_seed(&seed);
            assert_eq!(pk.as_bytes().to_vec(), r.hex("pk"), "tcId {}", r.id());
            assert_eq!(sk.as_bytes().to_vec(), r.hex("sk"), "tcId {}", r.id());
            // Round-trips: parsing, and recomputing pk from sk.
            let sk2 = m::SigningKey::from_bytes(&r.hex("sk")).unwrap();
            assert_eq!(sk2.verifying_key(), pk);
        });
    }
}

#[test]
fn acvp_siggen() {
    let recs = load("mldsa_siggen.txt");
    assert_eq!(recs.len(), 360);
    let mut prehash_algs = std::collections::BTreeSet::new();
    for r in &recs {
        let rnd: [u8; 32] = if r.bool("deterministic") { [0; 32] } else { arr(&r.hex("rnd")) };
        let expected = r.hex("signature");
        with_param_set!(r.get("parameterSet"), |m| {
            let sk = m::SigningKey::from_bytes(&r.hex("sk")).unwrap();
            let vk = m::VerifyingKey::from_bytes(&r.hex("pk")).unwrap();
            let sig = match (r.get("signatureInterface"), r.get("preHash"), r.bool("externalMu")) {
                ("internal", _, true) => sk.sign_mu(&arr(&r.hex("mu")), &rnd),
                ("internal", _, false) => sk.sign_internal(&r.hex("message"), &rnd),
                ("external", "pure", _) => {
                    sk.sign_with_rnd(&r.hex("message"), &r.hex("context"), &rnd).unwrap()
                }
                ("external", "preHash", _) => {
                    let alg = hash_alg(r.get("hashAlg"));
                    prehash_algs.insert(r.get("hashAlg").to_string());
                    sk.sign_prehash_with_rnd(&r.hex("message"), &r.hex("context"), alg, &rnd).unwrap()
                }
                other => panic!("unexpected test group {other:?}"),
            };
            assert_eq!(sig.as_bytes().to_vec(), expected, "tcId {}", r.id());

            // Every generated signature must also verify through the matching interface.
            let ok = match (r.get("signatureInterface"), r.get("preHash"), r.bool("externalMu")) {
                ("internal", _, true) => vk.verify_mu(&arr(&r.hex("mu")), &sig),
                ("internal", _, false) => vk.verify_internal(&r.hex("message"), &sig),
                ("external", "pure", _) => vk.verify(&r.hex("message"), &r.hex("context"), &sig).is_ok(),
                _ => vk
                    .verify_prehash(&r.hex("message"), &r.hex("context"), hash_alg(r.get("hashAlg")), &sig)
                    .is_ok(),
            };
            assert!(ok, "self-verification failed for tcId {}", r.id());
        });
    }
    assert!(prehash_algs.len() >= 8, "pre-hash coverage too small: {prehash_algs:?}");
}

#[test]
fn acvp_sigver() {
    let recs = load("mldsa_sigver.txt");
    assert_eq!(recs.len(), 180);
    let (mut pass, mut fail) = (0, 0);
    for r in &recs {
        let expected = r.bool("testPassed");
        with_param_set!(r.get("parameterSet"), |m| {
            let vk = m::VerifyingKey::from_bytes(&r.hex("pk")).unwrap();
            let sig = m::Signature::from_bytes(&r.hex("signature")).unwrap();
            let ok = match (r.get("signatureInterface"), r.get("preHash"), r.bool("externalMu")) {
                ("internal", _, true) => vk.verify_mu(&arr(&r.hex("mu")), &sig),
                ("internal", _, false) => vk.verify_internal(&r.hex("message"), &sig),
                ("external", "pure", _) => {
                    match vk.verify(&r.hex("message"), &r.hex("context"), &sig) {
                        Ok(()) => true,
                        Err(Error::InvalidSignature) => false,
                        Err(e) => panic!("unexpected error {e:?}"),
                    }
                }
                ("external", "preHash", _) => vk
                    .verify_prehash(&r.hex("message"), &r.hex("context"), hash_alg(r.get("hashAlg")), &sig)
                    .is_ok(),
                other => panic!("unexpected test group {other:?}"),
            };
            assert_eq!(ok, expected, "tcId {} ({})", r.id(), r.get("reason"));
        });
        if expected {
            pass += 1
        } else {
            fail += 1
        }
    }
    assert!(pass > 20 && fail > 100, "pass={pass} fail={fail}");
}
