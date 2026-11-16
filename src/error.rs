//! Error type shared by every algorithm in this crate.

use core::fmt;

/// Errors returned by the ML-KEM and ML-DSA APIs.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum Error {
    /// An input byte string has the wrong length for the parameter set.
    InvalidLength,
    /// The ML-KEM encapsulation key failed the FIPS 203 §7.2 modulus check.
    InvalidEncapsulationKey,
    /// The ML-KEM decapsulation key failed the FIPS 203 §7.3 hash check.
    InvalidDecapsulationKey,
    /// The ML-DSA signing key is malformed (coefficients of `s1`/`s2` out of range).
    InvalidSigningKey,
    /// The ML-DSA context string is longer than 255 bytes (FIPS 204 Alg. 2/3/4/5).
    ContextTooLong,
    /// A pre-hashed digest does not have the length of the selected hash function.
    InvalidDigestLength,
    /// An ML-DSA signature did not verify.
    InvalidSignature,
    /// The random bit generator failed to produce output.
    RandomnessFailure,
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let s = match self {
            Error::InvalidLength => "input has an invalid length",
            Error::InvalidEncapsulationKey => "encapsulation key failed the modulus check",
            Error::InvalidDecapsulationKey => "decapsulation key failed the hash check",
            Error::InvalidSigningKey => "signing key is malformed",
            Error::ContextTooLong => "context string exceeds 255 bytes",
            Error::InvalidDigestLength => "digest length does not match the hash function",
            Error::InvalidSignature => "signature verification failed",
            Error::RandomnessFailure => "random number generator failure",
        };
        f.write_str(s)
    }
}

#[cfg(feature = "std")]
impl std::error::Error for Error {}
