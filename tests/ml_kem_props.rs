//! Behavioral, negative and robustness tests for ML-KEM (all parameter sets).

mod common;
use common::*;
use fips_pqc::ml_kem::{ml_kem_1024, ml_kem_512, ml_kem_768, SharedSecret};
use fips_pqc::Error;
use fips_pqc::hash::Shake256;

macro_rules! kem_suite {
    ($suite:ident, $m:ident, k = $k:expr) => {
        mod $suite {
            use super::*;
            use $m as m;

            const K: usize = $k;

            #[test]
            fn roundtrip_many_keys() {
                let mut rng = TestRng::new(1);
                for _ in 0..64 {
                    let (ek, dk) = m::generate_key_pair(&mut rng).unwrap();
                    let (ct, ss) = ek.encapsulate(&mut rng).unwrap();
                    assert_eq!(dk.decapsulate(&ct), ss);
                    assert_ne!(ss.as_bytes(), &[0u8; 32]);
                }
            }

            #[test]
            fn keygen_is_deterministic_in_seed() {
                let seed = [0x42u8; 64];
                let (ek1, dk1) = m::generate_key_pair_from_seed(&seed);
                let (ek2, dk2) = m::generate_key_pair_from_seed(&seed);
                assert_eq!(ek1, ek2);
                assert_eq!(dk1, dk2);
                let mut other = seed;
                other[0] ^= 1;
                assert_ne!(m::generate_key_pair_from_seed(&other).0, ek1);
                // z does not influence ek but does influence dk.
                let mut z_changed = seed;
                z_changed[63] ^= 1;
                let (ek3, dk3) = m::generate_key_pair_from_seed(&z_changed);
                assert_eq!(ek3, ek1);
                assert_ne!(dk3, dk1);
            }

            #[test]
            fn key_layout_and_serialization() {
                let (ek, dk) = m::generate_key_pair(&mut TestRng::new(2)).unwrap();
                assert_eq!(ek.as_bytes().len(), m::ENCAPSULATION_KEY_SIZE);
                assert_eq!(dk.as_bytes().len(), m::DECAPSULATION_KEY_SIZE);
                assert_eq!(dk.encapsulation_key(), ek);
                // dk = dk_pke ‖ ek ‖ H(ek) ‖ z
                assert_eq!(&dk.as_bytes()[384 * K..768 * K + 32], &ek.as_bytes()[..]);
                assert_eq!(
                    &dk.as_bytes()[768 * K + 32..768 * K + 64],
                    &fips_pqc::hash::Sha3_256::digest(ek.as_bytes())[..]
                );
                let ek2 = m::EncapsulationKey::from_bytes(ek.as_bytes()).unwrap();
                let dk2 = m::DecapsulationKey::try_from(&dk.as_bytes()[..]).unwrap();
                assert_eq!(ek2, ek);
                assert_eq!(dk2, dk);
                let (ct, ss) = ek2.encapsulate(&mut TestRng::new(3)).unwrap();
                let ct2 = m::Ciphertext::from_bytes(ct.as_bytes()).unwrap();
                assert_eq!(dk2.decapsulate(&ct2), ss);
            }

            #[test]
            fn encapsulation_internal_is_deterministic() {
                let (ek, dk) = m::generate_key_pair(&mut TestRng::new(4)).unwrap();
                let msg = [9u8; 32];
                let (c1, k1) = ek.encapsulate_internal(&msg);
                let (c2, k2) = ek.encapsulate_internal(&msg);
                assert_eq!(c1, c2);
                assert_eq!(k1, k2);
                let (c3, k3) = ek.encapsulate_internal(&[10u8; 32]);
                assert_ne!(c1, c3);
                assert_ne!(k1, k3);
                assert_eq!(dk.decapsulate(&c3), k3);
            }

            #[test]
            fn implicit_rejection_on_every_bit_region() {
                let (ek, dk) = m::generate_key_pair(&mut TestRng::new(5)).unwrap();
                let (ct, ss) = ek.encapsulate(&mut TestRng::new(6)).unwrap();
                let z = &dk.as_bytes()[m::DECAPSULATION_KEY_SIZE - 32..];
                let mut rng = TestRng::new(7);
                let positions = [0, 1, m::CIPHERTEXT_SIZE / 2, m::CIPHERTEXT_SIZE - 1]
                    .into_iter()
                    .chain((0..60).map(|_| rng.below(m::CIPHERTEXT_SIZE)));
                for pos in positions {
                    let mut bad = *ct.as_bytes();
                    bad[pos] ^= 1 << (pos % 8);
                    let bad = m::Ciphertext::from(bad);
                    let k_bad = dk.decapsulate(&bad);
                    assert_ne!(k_bad, ss, "pos {pos}");
                    // K̄ = SHAKE256(z ‖ c, 32)
                    let mut expect = [0u8; 32];
                    let mut h = Shake256::new();
                    h.update(z);
                    h.update(bad.as_bytes());
                    h.finalize_xof().read(&mut expect);
                    assert_eq!(k_bad.as_bytes(), &expect, "pos {pos}");
                    // Deterministic.
                    assert_eq!(dk.decapsulate(&bad), k_bad);
                }
            }

            #[test]
            fn wrong_key_does_not_decapsulate() {
                let mut rng = TestRng::new(8);
                let (ek, _) = m::generate_key_pair(&mut rng).unwrap();
                let (_, dk_other) = m::generate_key_pair(&mut rng).unwrap();
                let (ct, ss) = ek.encapsulate(&mut rng).unwrap();
                assert_ne!(dk_other.decapsulate(&ct), ss);
            }

            #[test]
            fn modulus_check_rejects_out_of_range_coefficients() {
                let (ek, _) = m::generate_key_pair(&mut TestRng::new(9)).unwrap();
                for poly in 0..K {
                    for coeff in [0usize, 1, 255] {
                        let mut bad = *ek.as_bytes();
                        // Set coefficient `coeff` of polynomial `poly` to 3329 (= q).
                        let base = 384 * poly + 3 * (coeff / 2);
                        if coeff % 2 == 0 {
                            bad[base] = 0x01;
                            bad[base + 1] = (bad[base + 1] & 0xF0) | 0x0D;
                        } else {
                            bad[base + 1] = (bad[base + 1] & 0x0F) | 0x10;
                            bad[base + 2] = 0xD0;
                        }
                        assert_eq!(
                            m::EncapsulationKey::from_bytes(&bad),
                            Err(Error::InvalidEncapsulationKey),
                            "poly {poly} coeff {coeff}"
                        );
                    }
                }
                // Coefficient 4095 (all ones) is also rejected.
                let mut bad = *ek.as_bytes();
                bad[0] = 0xFF;
                bad[1] |= 0x0F;
                assert!(m::EncapsulationKey::from_bytes(&bad).is_err());
                // The ρ seed is unconstrained.
                let mut ok = *ek.as_bytes();
                ok[384 * K] ^= 0xFF;
                assert!(m::EncapsulationKey::from_bytes(&ok).is_ok());
            }

            #[test]
            fn hash_check_rejects_tampered_dk() {
                let (_, dk) = m::generate_key_pair(&mut TestRng::new(10)).unwrap();
                // Tamper with the embedded ek, and separately with H(ek).
                for pos in [384 * K, 768 * K, 768 * K + 32, 768 * K + 63] {
                    let mut bad = *dk.as_bytes();
                    bad[pos] ^= 0x80;
                    assert_eq!(
                        m::DecapsulationKey::from_bytes(&bad),
                        Err(Error::InvalidDecapsulationKey),
                        "pos {pos}"
                    );
                }
                // dk_pke and z are not covered by the hash check.
                let mut ok = *dk.as_bytes();
                ok[0] ^= 1;
                ok[m::DECAPSULATION_KEY_SIZE - 1] ^= 1;
                assert!(m::DecapsulationKey::from_bytes(&ok).is_ok());
            }

            #[test]
            fn length_checks() {
                let (ek, dk) = m::generate_key_pair(&mut TestRng::new(11)).unwrap();
                let ek_b = ek.as_bytes();
                let dk_b = dk.as_bytes();
                assert_eq!(m::EncapsulationKey::from_bytes(&ek_b[1..]), Err(Error::InvalidLength));
                assert_eq!(m::EncapsulationKey::from_bytes(&[ek_b.as_slice(), &[0]].concat()), Err(Error::InvalidLength));
                assert_eq!(m::DecapsulationKey::from_bytes(&dk_b[1..]).err(), Some(Error::InvalidLength));
                assert_eq!(m::Ciphertext::from_bytes(&[0u8; 10]).err(), Some(Error::InvalidLength));
                assert_eq!(m::Ciphertext::from_bytes(&vec![0u8; m::CIPHERTEXT_SIZE + 1]).err(), Some(Error::InvalidLength));
                assert!(m::Ciphertext::from_bytes(&vec![0u8; m::CIPHERTEXT_SIZE]).is_ok());
            }

            #[test]
            fn rng_failure_is_propagated() {
                assert_eq!(m::generate_key_pair(&mut FailingRng).err(), Some(Error::RandomnessFailure));
                let (ek, _) = m::generate_key_pair(&mut TestRng::new(12)).unwrap();
                assert_eq!(ek.encapsulate(&mut FailingRng).err(), Some(Error::RandomnessFailure));
            }

            #[test]
            fn all_zero_and_all_ones_ciphertexts_are_handled() {
                let (_, dk) = m::generate_key_pair(&mut TestRng::new(13)).unwrap();
                let a = dk.decapsulate(&m::Ciphertext::from([0u8; m::CIPHERTEXT_SIZE]));
                let b = dk.decapsulate(&m::Ciphertext::from([0xFFu8; m::CIPHERTEXT_SIZE]));
                assert_ne!(a, b);
            }

            #[test]
            fn debug_output_does_not_leak_secrets() {
                let (_, dk) = m::generate_key_pair(&mut TestRng::new(14)).unwrap();
                let s = format!("{dk:?}");
                assert!(s.contains("redacted"));
                let (ek, _) = m::generate_key_pair(&mut TestRng::new(15)).unwrap();
                let (_, ss) = ek.encapsulate(&mut TestRng::new(16)).unwrap();
                assert_eq!(format!("{ss:?}"), "SharedSecret(<redacted>)");
            }
        }
    };
}

kem_suite!(kem512, ml_kem_512, k = 2);
kem_suite!(kem768, ml_kem_768, k = 3);
kem_suite!(kem1024, ml_kem_1024, k = 4);

#[test]
fn shared_secret_equality_is_by_value() {
    let (ek, dk) = ml_kem_768::generate_key_pair(&mut TestRng::new(100)).unwrap();
    let (ct, a) = ek.encapsulate(&mut TestRng::new(101)).unwrap();
    let b: SharedSecret = dk.decapsulate(&ct);
    assert_eq!(a.clone(), b);
}

#[test]
fn parameter_sets_are_domain_separated() {
    // The same seed yields unrelated keys across parameter sets (k is hashed into G).
    let seed = [7u8; 64];
    let (ek512, _) = ml_kem_512::generate_key_pair_from_seed(&seed);
    let (ek768, _) = ml_kem_768::generate_key_pair_from_seed(&seed);
    assert_ne!(&ek512.as_bytes()[..32], &ek768.as_bytes()[..32]);
}

#[cfg(unix)]
#[test]
fn os_rng_roundtrip() {
    use fips_pqc::OsRng;
    let (ek, dk) = ml_kem_1024::generate_key_pair(&mut OsRng).unwrap();
    let (ct, ss) = ek.encapsulate(&mut OsRng).unwrap();
    assert_eq!(dk.decapsulate(&ct), ss);
}
