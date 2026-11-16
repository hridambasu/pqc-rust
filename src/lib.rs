//! # fips-pqc
//!
//! A pure-Rust, dependency-free, `no_std`-compatible implementation of the NIST
//! post-quantum standards:
//!
//! * **ML-KEM** — FIPS 203, *Module-Lattice-Based Key-Encapsulation Mechanism*
//!   ([`ml_kem::ml_kem_512`], [`ml_kem::ml_kem_768`], [`ml_kem::ml_kem_1024`]).
//! * **ML-DSA** — FIPS 204, *Module-Lattice-Based Digital Signature Algorithm*
//!   ([`ml_dsa::ml_dsa_44`], [`ml_dsa::ml_dsa_65`], [`ml_dsa::ml_dsa_87`]),
//!   including HashML-DSA and the external-μ interface.
//!
//! The SHA-3/SHAKE (FIPS 202) and SHA-2 (FIPS 180-4) functions the standards
//! depend on are implemented in [`hash`].
//!
//! The implementation is validated against the NIST ACVP test vectors for
//! ML-KEM key generation, encapsulation, decapsulation and key checks, and
//! ML-DSA key generation, signature generation and signature verification.
//!
//! ## Side channels
//!
//! Operations on secret data avoid secret-dependent branches, memory indexing
//! and variable-time division (no `%`/`/` on secrets; Barrett reduction is used
//! throughout). Rejection-sampling loops leak only what FIPS 204 permits (the
//! number of iterations). Secret intermediates and key types are zeroized on
//! drop. This crate has **not** been independently audited or CMVP validated.
//!
//! ## Randomness
//!
//! Randomized functions take any [`CryptoRng`]. With the `std` feature on Unix
//! systems, [`OsRng`] reads from the operating system CSPRNG.

#![cfg_attr(not(feature = "std"), no_std)]
#![forbid(unsafe_code)]
#![warn(missing_docs, rust_2018_idioms)]
#![allow(clippy::needless_range_loop)]

mod common;
mod error;
pub mod hash;
pub mod ml_dsa;
pub mod ml_kem;
mod rng;

pub use error::Error;
#[cfg(all(feature = "std", unix))]
pub use rng::OsRng;
pub use rng::CryptoRng;
