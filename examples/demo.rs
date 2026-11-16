//! End-to-end usage of ML-KEM-768 and ML-DSA-65.
//!
//! Run with: `cargo run --example demo`

use fips_pqc::ml_dsa::{ml_dsa_65, PreHashAlgorithm};
use fips_pqc::ml_kem::ml_kem_768;
use fips_pqc::OsRng;

fn main() -> Result<(), fips_pqc::Error> {
    let mut rng = OsRng;

    // --- ML-KEM: establish a shared secret -------------------------------
    let (ek, dk) = ml_kem_768::generate_key_pair(&mut rng)?;
    // The receiver publishes `ek`; the sender parses it (FIPS 203 §7.2 checks).
    let ek = ml_kem_768::EncapsulationKey::from_bytes(ek.as_bytes())?;
    let (ct, sender_secret) = ek.encapsulate(&mut rng)?;
    let ct = ml_kem_768::Ciphertext::from_bytes(ct.as_bytes())?;
    let receiver_secret = dk.decapsulate(&ct);
    assert_eq!(sender_secret, receiver_secret);
    println!("ML-KEM-768: ek {} B, ct {} B, shared secret agreed", ek.as_bytes().len(), ct.as_bytes().len());

    // --- ML-DSA: sign and verify -----------------------------------------
    let (vk, sk) = ml_dsa_65::generate_key_pair(&mut rng)?;
    let msg = b"post-quantum hello";
    let sig = sk.sign(msg, b"demo-context", &mut rng)?;
    vk.verify(msg, b"demo-context", &sig)?;
    println!("ML-DSA-65: signature {} B verified", sig.as_bytes().len());

    // HashML-DSA with SHA-512 pre-hashing.
    let sig = sk.sign_prehash(msg, b"", PreHashAlgorithm::Sha2_512, &mut rng)?;
    vk.verify_prehash(msg, b"", PreHashAlgorithm::Sha2_512, &sig)?;
    println!("HashML-DSA-65 (SHA-512): verified");

    // Tampering is detected.
    assert!(vk.verify(b"post-quantum hellO", b"demo-context", &sig).is_err());
    Ok(())
}
