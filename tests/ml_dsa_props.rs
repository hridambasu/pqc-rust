//! Behavioral, negative and robustness tests for ML-DSA (all parameter sets).

mod common;
use common::*;
use fips_pqc::ml_dsa::{ml_dsa_44, ml_dsa_65, ml_dsa_87, PreHashAlgorithm, MAX_CONTEXT_SIZE};
use fips_pqc::Error;

macro_rules! dsa_suite {
    ($suite:ident, $m:ident, k = $k:expr, l = $l:expr, omega = $omega:expr, clen = $clen:expr,
     zbits = $zbits:expr, eta_bits = $eta_bits:expr) => {
        mod $suite {
            use super::*;
            use $m as m;

            const K: usize = $k;
            const L: usize = $l;
            const OMEGA: usize = $omega;
            const CLEN: usize = $clen;
            const ZB: usize = 32 * $zbits;

            fn keys(seed: u64) -> (m::VerifyingKey, m::SigningKey) {
                m::generate_key_pair(&mut TestRng::new(seed)).unwrap()
            }

            #[test]
            fn sign_verify_roundtrip_varied_inputs() {
                let mut rng = TestRng::new(1);
                for i in 0..24 {
                    let (vk, sk) = m::generate_key_pair(&mut rng).unwrap();
                    let msg_len = rng.below(2000);
                    let msg = rng.vec(msg_len);
                    let ctx_len = if i % 3 == 0 { 0 } else { rng.below(256) };
                    let ctx = rng.vec(ctx_len);
                    let sig = sk.sign(&msg, &ctx, &mut rng).unwrap();
                    assert_eq!(sig.as_bytes().len(), m::SIGNATURE_SIZE);
                    vk.verify(&msg, &ctx, &sig).unwrap();
                }
            }

            #[test]
            fn deterministic_and_hedged_modes() {
                let (vk, sk) = keys(2);
                let a = sk.sign_deterministic(b"msg", b"").unwrap();
                let b = sk.sign_deterministic(b"msg", b"").unwrap();
                assert_eq!(a, b);
                assert_eq!(sk.sign_with_rnd(b"msg", b"", &[0; 32]).unwrap(), a);
                let mut rng = TestRng::new(3);
                let h1 = sk.sign(b"msg", b"", &mut rng).unwrap();
                let h2 = sk.sign(b"msg", b"", &mut rng).unwrap();
                assert_ne!(h1, h2);
                assert_ne!(h1, a);
                for s in [&a, &h1, &h2] {
                    vk.verify(b"msg", b"", s).unwrap();
                }
            }

            #[test]
            fn tampering_is_detected() {
                let (vk, sk) = keys(4);
                let msg = b"the quick brown fox";
                let ctx = b"ctx";
                let sig = sk.sign(msg, ctx, &mut TestRng::new(5)).unwrap();
                assert!(vk.verify(msg, ctx, &sig).is_ok());

                assert_eq!(vk.verify(b"the quick brown fix", ctx, &sig), Err(Error::InvalidSignature));
                assert_eq!(vk.verify(msg, b"ctX", &sig), Err(Error::InvalidSignature));
                assert_eq!(vk.verify(msg, b"", &sig), Err(Error::InvalidSignature));

                let mut rng = TestRng::new(6);
                let positions = [0, CLEN - 1, CLEN, CLEN + L * ZB - 1, m::SIGNATURE_SIZE - 1]
                    .into_iter()
                    .chain((0..80).map(|_| rng.below(m::SIGNATURE_SIZE)));
                for pos in positions {
                    let mut bad = *sig.as_bytes();
                    bad[pos] ^= 1 << (pos % 8);
                    let bad = m::Signature::from(bad);
                    assert!(vk.verify(msg, ctx, &bad).is_err(), "bit flip at {pos} accepted");
                }

                let mut bad_pk = *vk.as_bytes();
                bad_pk[40] ^= 1;
                let vk2 = m::VerifyingKey::from_bytes(&bad_pk).unwrap();
                assert!(vk2.verify(msg, ctx, &sig).is_err());

                let (other_vk, _) = keys(7);
                assert!(other_vk.verify(msg, ctx, &sig).is_err());
            }

            #[test]
            fn context_length_limits() {
                let (vk, sk) = keys(8);
                let max = vec![0xAB; MAX_CONTEXT_SIZE];
                let sig = sk.sign_deterministic(b"m", &max).unwrap();
                vk.verify(b"m", &max, &sig).unwrap();
                let too_long = vec![0u8; MAX_CONTEXT_SIZE + 1];
                assert_eq!(sk.sign_deterministic(b"m", &too_long), Err(Error::ContextTooLong));
                assert_eq!(vk.verify(b"m", &too_long, &sig), Err(Error::ContextTooLong));
                assert_eq!(
                    sk.sign_prehash_deterministic(b"m", &too_long, PreHashAlgorithm::Sha2_512),
                    Err(Error::ContextTooLong)
                );
                assert_eq!(vk.compute_mu(b"m", &too_long), Err(Error::ContextTooLong));
            }

            #[test]
            fn prehash_all_algorithms() {
                let (vk, sk) = keys(9);
                let msg = b"pre-hashed message";
                let mut rng = TestRng::new(10);
                for alg in PreHashAlgorithm::ALL {
                    let sig = sk.sign_prehash(msg, b"c", alg, &mut rng).unwrap();
                    vk.verify_prehash(msg, b"c", alg, &sig).unwrap();
                    // Domain separation: a HashML-DSA signature is not a pure signature,
                    // and is bound to the hash algorithm.
                    assert!(vk.verify(msg, b"c", &sig).is_err());
                    let other = if alg == PreHashAlgorithm::Sha2_256 {
                        PreHashAlgorithm::Sha3_256
                    } else {
                        PreHashAlgorithm::Sha2_256
                    };
                    assert!(vk.verify_prehash(msg, b"c", other, &sig).is_err());
                    assert!(vk.verify_prehash(b"different", b"c", alg, &sig).is_err());
                }
                let pure = sk.sign_deterministic(msg, b"c").unwrap();
                assert!(vk.verify_prehash(msg, b"c", PreHashAlgorithm::Sha2_512, &pure).is_err());
            }

            #[test]
            fn prehash_digest_interface() {
                let (vk, sk) = keys(11);
                let msg = b"digest me";
                let digest = fips_pqc::hash::Sha512::digest(msg);
                let alg = PreHashAlgorithm::Sha2_512;
                let s1 = sk.sign_digest_with_rnd(&digest, b"", alg, &[0; 32]).unwrap();
                let s2 = sk.sign_prehash_deterministic(msg, b"", alg).unwrap();
                assert_eq!(s1, s2);
                vk.verify_digest(&digest, b"", alg, &s1).unwrap();
                assert_eq!(
                    sk.sign_digest_with_rnd(&digest[..32], b"", alg, &[0; 32]),
                    Err(Error::InvalidDigestLength)
                );
                assert_eq!(vk.verify_digest(&digest[..63], b"", alg, &s1), Err(Error::InvalidDigestLength));
            }

            #[test]
            fn external_mu_and_internal_interfaces_agree() {
                let (vk, sk) = keys(12);
                let msg = b"external mu";
                let ctx = b"app";
                let mu = vk.compute_mu(msg, ctx).unwrap();
                let rnd = [5u8; 32];
                let via_mu = sk.sign_mu(&mu, &rnd);
                let via_pure = sk.sign_with_rnd(msg, ctx, &rnd).unwrap();
                let m_prime = [&[0u8, ctx.len() as u8][..], ctx, msg].concat();
                let via_internal = sk.sign_internal(&m_prime, &rnd);
                assert_eq!(via_mu, via_pure);
                assert_eq!(via_mu, via_internal);
                assert!(vk.verify_mu(&mu, &via_mu));
                assert!(vk.verify_internal(&m_prime, &via_mu));
                let mut mu2 = mu;
                mu2[0] ^= 1;
                assert!(!vk.verify_mu(&mu2, &via_mu));
            }

            #[test]
            fn keygen_from_seed_and_key_parsing() {
                let seed = [0x11u8; 32];
                let (vk1, sk1) = m::generate_key_pair_from_seed(&seed);
                let (vk2, sk2) = m::generate_key_pair_from_seed(&seed);
                assert_eq!(vk1, vk2);
                assert_eq!(sk1, sk2);
                assert_eq!(sk1.verifying_key(), vk1);
                // sk = ρ ‖ K ‖ tr ‖ ...; pk = ρ ‖ t1
                assert_eq!(&sk1.as_bytes()[..32], &vk1.as_bytes()[..32]);
                let mut tr = [0u8; 64];
                fips_pqc::hash::Shake256::digest(vk1.as_bytes(), &mut tr);
                assert_eq!(&sk1.as_bytes()[64..128], &tr[..]);

                let sk3 = m::SigningKey::from_bytes(sk1.as_bytes()).unwrap();
                assert_eq!(sk3, sk1);
                assert_eq!(m::SigningKey::from_bytes(&sk1.as_bytes()[1..]).err(), Some(Error::InvalidLength));
                assert_eq!(m::VerifyingKey::from_bytes(&vk1.as_bytes()[1..]).err(), Some(Error::InvalidLength));
                assert_eq!(m::Signature::from_bytes(&[0u8; 3]).err(), Some(Error::InvalidLength));

                // Out-of-range s1 coefficient (all bits set ⇒ value 2^bits − 1 > 2η).
                let mut bad = *sk1.as_bytes();
                for b in &mut bad[128..128 + $eta_bits] {
                    *b = 0xFF;
                }
                assert_eq!(m::SigningKey::from_bytes(&bad).err(), Some(Error::InvalidSigningKey));
            }

            #[test]
            fn malformed_hints_are_rejected() {
                let (vk, sk) = keys(13);
                let sig = sk.sign_deterministic(b"hint", b"").unwrap();
                vk.verify(b"hint", b"", &sig).unwrap();
                let hstart = CLEN + L * ZB;
                let base = *sig.as_bytes();

                // Cumulative count larger than ω.
                let mut bad = base;
                bad[hstart + OMEGA + K - 1] = (OMEGA + 1) as u8;
                assert!(vk.verify(b"hint", b"", &m::Signature::from(bad)).is_err());

                // Decreasing cumulative counts.
                let mut bad = base;
                let total = usize::from(base[hstart + OMEGA + K - 1]);
                if total > 0 {
                    bad[hstart + OMEGA] = total as u8;
                    bad[hstart + OMEGA + 1] = 0;
                    assert!(vk.verify(b"hint", b"", &m::Signature::from(bad)).is_err());
                }

                // Non-zero padding byte after the last used hint index.
                let mut bad = base;
                if total < OMEGA {
                    bad[hstart + OMEGA - 1] = 1;
                    assert!(vk.verify(b"hint", b"", &m::Signature::from(bad)).is_err());
                }

                // Duplicate (non-strictly-increasing) indices.
                let mut bad = base;
                let first = usize::from(base[hstart + OMEGA]);
                if first >= 2 {
                    bad[hstart + 1] = bad[hstart];
                    assert!(vk.verify(b"hint", b"", &m::Signature::from(bad)).is_err());
                }
            }

            #[test]
            fn rejects_z_with_large_norm() {
                let (vk, sk) = keys(14);
                let sig = sk.sign_deterministic(b"z", b"").unwrap();
                let mut bad = *sig.as_bytes();
                // Encoded value 0 decodes to z = γ1, which violates ‖z‖∞ < γ1 − β.
                for b in &mut bad[CLEN..CLEN + 3] {
                    *b = 0;
                }
                assert!(vk.verify(b"z", b"", &m::Signature::from(bad)).is_err());
            }

            #[test]
            fn rng_failure_is_propagated() {
                assert_eq!(m::generate_key_pair(&mut FailingRng).err(), Some(Error::RandomnessFailure));
                let (_, sk) = keys(15);
                assert_eq!(sk.sign(b"m", b"", &mut FailingRng).err(), Some(Error::RandomnessFailure));
                assert_eq!(
                    sk.sign_prehash(b"m", b"", PreHashAlgorithm::Sha2_256, &mut FailingRng).err(),
                    Some(Error::RandomnessFailure)
                );
            }

            #[test]
            fn large_message() {
                let (vk, sk) = keys(16);
                let msg = vec![0x5Au8; 1 << 20];
                let sig = sk.sign(&msg, b"", &mut TestRng::new(17)).unwrap();
                vk.verify(&msg, b"", &sig).unwrap();
                let mut msg2 = msg.clone();
                *msg2.last_mut().unwrap() ^= 1;
                assert!(vk.verify(&msg2, b"", &sig).is_err());
            }

            #[test]
            fn debug_output_does_not_leak_secrets() {
                let (_, sk) = keys(18);
                assert!(format!("{sk:?}").contains("redacted"));
            }
        }
    };
}

dsa_suite!(dsa44, ml_dsa_44, k = 4, l = 4, omega = 80, clen = 32, zbits = 18, eta_bits = 3);
dsa_suite!(dsa65, ml_dsa_65, k = 6, l = 5, omega = 55, clen = 48, zbits = 20, eta_bits = 4);
dsa_suite!(dsa87, ml_dsa_87, k = 8, l = 7, omega = 75, clen = 64, zbits = 20, eta_bits = 3);

#[test]
fn parameter_sets_are_domain_separated() {
    let seed = [3u8; 32];
    let (a, _) = ml_dsa_44::generate_key_pair_from_seed(&seed);
    let (b, _) = ml_dsa_65::generate_key_pair_from_seed(&seed);
    assert_ne!(&a.as_bytes()[..32], &b.as_bytes()[..32]);
}

#[cfg(unix)]
#[test]
fn os_rng_roundtrip() {
    use fips_pqc::OsRng;
    let (vk, sk) = ml_dsa_87::generate_key_pair(&mut OsRng).unwrap();
    let sig = sk.sign(b"hello", b"", &mut OsRng).unwrap();
    vk.verify(b"hello", b"", &sig).unwrap();
}
