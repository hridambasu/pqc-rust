//! Randomness source abstraction.
//!
//! FIPS 203 and FIPS 204 require randomness from an approved random bit
//! generator (RBG). This crate does not ship a DRBG; callers plug in their own
//! RBG through the [`CryptoRng`] trait. With the `std` feature on Unix-like
//! systems, [`OsRng`] reads from the operating system's CSPRNG.

use crate::Error;

/// A cryptographically secure random bit generator.
///
/// Implementations must return an error rather than produce weak output.
pub trait CryptoRng {
    /// Fills `dest` entirely with random bytes.
    fn fill_bytes(&mut self, dest: &mut [u8]) -> Result<(), Error>;
}

impl<R: CryptoRng + ?Sized> CryptoRng for &mut R {
    fn fill_bytes(&mut self, dest: &mut [u8]) -> Result<(), Error> {
        (**self).fill_bytes(dest)
    }
}

/// The operating system's CSPRNG (`/dev/urandom`).
#[cfg(all(feature = "std", unix))]
#[derive(Debug, Clone, Copy, Default)]
pub struct OsRng;

#[cfg(all(feature = "std", unix))]
impl CryptoRng for OsRng {
    fn fill_bytes(&mut self, dest: &mut [u8]) -> Result<(), Error> {
        use std::io::Read;
        std::fs::File::open("/dev/urandom")
            .and_then(|mut f| f.read_exact(dest))
            .map_err(|_| Error::RandomnessFailure)
    }
}
