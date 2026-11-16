//! Hash functions used by ML-KEM and ML-DSA.
//!
//! * SHA-3 / SHAKE (FIPS 202) — the primitives inside both algorithms.
//! * SHA-2 (FIPS 180-4) — approved pre-hash functions for HashML-DSA.

mod sha2;
mod sha3;

pub use sha2::{Sha224, Sha256, Sha384, Sha512, Sha512_224, Sha512_256};
pub use sha3::{Sha3_224, Sha3_256, Sha3_384, Sha3_512, Shake128, Shake256, XofReader};

pub(crate) use sha3::{sha3_256, sha3_512, shake256, Sponge};
