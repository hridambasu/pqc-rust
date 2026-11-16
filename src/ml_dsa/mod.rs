//! ML-DSA, the Module-Lattice-Based Digital Signature Algorithm (FIPS 204).
//!
//! Three parameter sets are provided, each in its own module:
//!
//! | Module            | Security category | public key | private key | signature |
//! |-------------------|-------------------|------------|-------------|-----------|
//! | [`ml_dsa_44`]     | 2                 | 1312       | 2560        | 2420      |
//! | [`ml_dsa_65`]     | 3                 | 1952       | 4032        | 3309      |
//! | [`ml_dsa_87`]     | 5                 | 2592       | 4896        | 4627      |
//!
//! Supported interfaces:
//!
//! * **Pure ML-DSA** (Algorithms 2/3): [`ml_dsa_65::SigningKey::sign`],
//!   [`ml_dsa_65::VerifyingKey::verify`]; hedged (default) or deterministic.
//! * **HashML-DSA** (Algorithms 4/5) with every pre-hash in [`PreHashAlgorithm`].
//! * **Internal / external-μ** functions (Algorithms 7/8) for protocols and testing.
//!
//! ```
//! # #[cfg(all(feature = "std", unix))] {
//! use fips_pqc::{ml_dsa::ml_dsa_65, OsRng};
//!
//! let (vk, sk) = ml_dsa_65::generate_key_pair(&mut OsRng).unwrap();
//! let sig = sk.sign(b"message", b"context", &mut OsRng).unwrap();
//! assert!(vk.verify(b"message", b"context", &sig).is_ok());
//! # }
//! ```

mod arith;
mod sample;
mod scheme;

use crate::common::{ct_eq, zeroize_bytes};
use crate::hash::*;
use crate::Error;
use core::fmt;

/// Size in bytes of the seed ξ from which a key pair is derived.
pub const SEED_SIZE: usize = 32;
/// Size in bytes of the message representative μ.
pub const MU_SIZE: usize = 64;
/// Size in bytes of the per-signature randomness `rnd`.
pub const RND_SIZE: usize = 32;
/// Maximum length in bytes of a context string.
pub const MAX_CONTEXT_SIZE: usize = 255;

/// Approved pre-hash functions for HashML-DSA (FIPS 204 §5.4).
///
/// FIPS 204 recommends a hash whose collision strength is at least that of the
/// ML-DSA parameter set (e.g. SHA-512 / SHAKE256 for ML-DSA-87), but all of the
/// listed functions are permitted.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum PreHashAlgorithm {
    /// SHA-224
    Sha2_224,
    /// SHA-256
    Sha2_256,
    /// SHA-384
    Sha2_384,
    /// SHA-512
    Sha2_512,
    /// SHA-512/224
    Sha2_512_224,
    /// SHA-512/256
    Sha2_512_256,
    /// SHA3-224
    Sha3_224,
    /// SHA3-256
    Sha3_256,
    /// SHA3-384
    Sha3_384,
    /// SHA3-512
    Sha3_512,
    /// SHAKE128 with 256-bit output
    Shake128,
    /// SHAKE256 with 512-bit output
    Shake256,
}

impl PreHashAlgorithm {
    /// All supported pre-hash algorithms.
    pub const ALL: [PreHashAlgorithm; 12] = [
        Self::Sha2_224,
        Self::Sha2_256,
        Self::Sha2_384,
        Self::Sha2_512,
        Self::Sha2_512_224,
        Self::Sha2_512_256,
        Self::Sha3_224,
        Self::Sha3_256,
        Self::Sha3_384,
        Self::Sha3_512,
        Self::Shake128,
        Self::Shake256,
    ];

    /// DER encoding of the algorithm's OID (`2.16.840.1.101.3.4.2.x`).
    pub const fn oid(self) -> [u8; 11] {
        let last = match self {
            Self::Sha2_256 => 0x01,
            Self::Sha2_384 => 0x02,
            Self::Sha2_512 => 0x03,
            Self::Sha2_224 => 0x04,
            Self::Sha2_512_224 => 0x05,
            Self::Sha2_512_256 => 0x06,
            Self::Sha3_224 => 0x07,
            Self::Sha3_256 => 0x08,
            Self::Sha3_384 => 0x09,
            Self::Sha3_512 => 0x0A,
            Self::Shake128 => 0x0B,
            Self::Shake256 => 0x0C,
        };
        [0x06, 0x09, 0x60, 0x86, 0x48, 0x01, 0x65, 0x03, 0x04, 0x02, last]
    }

    /// Length in bytes of the digest `PH(M)`.
    pub const fn digest_len(self) -> usize {
        match self {
            Self::Sha2_224 | Self::Sha2_512_224 | Self::Sha3_224 => 28,
            Self::Sha2_256 | Self::Sha2_512_256 | Self::Sha3_256 | Self::Shake128 => 32,
            Self::Sha2_384 | Self::Sha3_384 => 48,
            Self::Sha2_512 | Self::Sha3_512 | Self::Shake256 => 64,
        }
    }

    /// Computes `PH(M)` into the first [`digest_len`](Self::digest_len) bytes of `out`.
    fn digest(self, msg: &[u8], out: &mut [u8; 64]) -> usize {
        let n = self.digest_len();
        match self {
            Self::Sha2_224 => out[..n].copy_from_slice(&Sha224::digest(msg)),
            Self::Sha2_256 => out[..n].copy_from_slice(&Sha256::digest(msg)),
            Self::Sha2_384 => out[..n].copy_from_slice(&Sha384::digest(msg)),
            Self::Sha2_512 => out[..n].copy_from_slice(&Sha512::digest(msg)),
            Self::Sha2_512_224 => out[..n].copy_from_slice(&Sha512_224::digest(msg)),
            Self::Sha2_512_256 => out[..n].copy_from_slice(&Sha512_256::digest(msg)),
            Self::Sha3_224 => out[..n].copy_from_slice(&Sha3_224::digest(msg)),
            Self::Sha3_256 => out[..n].copy_from_slice(&Sha3_256::digest(msg)),
            Self::Sha3_384 => out[..n].copy_from_slice(&Sha3_384::digest(msg)),
            Self::Sha3_512 => out[..n].copy_from_slice(&Sha3_512::digest(msg)),
            Self::Shake128 => Shake128::digest(msg, &mut out[..n]),
            Self::Shake256 => Shake256::digest(msg, &mut out[..n]),
        }
        n
    }
}

fn check_ctx(ctx: &[u8]) -> Result<u8, Error> {
    u8::try_from(ctx.len()).map_err(|_| Error::ContextTooLong)
}

/// μ for pure ML-DSA: M′ = 0 ‖ |ctx| ‖ ctx ‖ M (Algorithm 2, line 10).
fn mu_pure(tr: &[u8], msg: &[u8], ctx: &[u8]) -> Result<[u8; 64], Error> {
    let n = check_ctx(ctx)?;
    Ok(scheme::compute_mu(tr, &[&[0, n], ctx, msg]))
}

/// μ for HashML-DSA: M′ = 1 ‖ |ctx| ‖ ctx ‖ OID ‖ PH(M) (Algorithm 4, line 18).
fn mu_prehashed(
    tr: &[u8],
    digest: &[u8],
    ctx: &[u8],
    alg: PreHashAlgorithm,
) -> Result<[u8; 64], Error> {
    let n = check_ctx(ctx)?;
    if digest.len() != alg.digest_len() {
        return Err(Error::InvalidDigestLength);
    }
    Ok(scheme::compute_mu(tr, &[&[1, n], ctx, &alg.oid(), digest]))
}

macro_rules! ml_dsa_parameter_set {
    (
        $(#[$meta:meta])*
        $module:ident, $name:literal, k = $k:literal, l = $l:literal, $params:expr
    ) => {
        $(#[$meta])*
        pub mod $module {
            use super::*;
            use crate::CryptoRng;

            const K: usize = $k;
            const L: usize = $l;
            const PARAMS: scheme::Params = $params;

            #[doc = concat!("Size in bytes of an ", $name, " public (verifying) key.")]
            pub const VERIFYING_KEY_SIZE: usize = PARAMS.pk_len(K);
            #[doc = concat!("Size in bytes of an ", $name, " private (signing) key.")]
            pub const SIGNING_KEY_SIZE: usize = PARAMS.sk_len(K, L);
            #[doc = concat!("Size in bytes of an ", $name, " signature.")]
            pub const SIGNATURE_SIZE: usize = PARAMS.sig_len(K, L);

            #[doc = concat!("An ", $name, " public (verifying) key.")]
            #[derive(Clone, PartialEq, Eq)]
            pub struct VerifyingKey([u8; VERIFYING_KEY_SIZE]);

            #[doc = concat!("An ", $name, " private (signing) key. Zeroized on drop.")]
            #[derive(Clone)]
            pub struct SigningKey([u8; SIGNING_KEY_SIZE]);

            #[doc = concat!("An ", $name, " signature.")]
            #[derive(Clone, PartialEq, Eq)]
            pub struct Signature([u8; SIGNATURE_SIZE]);

            /// Algorithm 1: ML-DSA.KeyGen, using a 32-byte seed from `rng`.
            pub fn generate_key_pair<R: CryptoRng + ?Sized>(
                rng: &mut R,
            ) -> Result<(VerifyingKey, SigningKey), Error> {
                let mut xi = [0u8; SEED_SIZE];
                rng.fill_bytes(&mut xi)?;
                let kp = generate_key_pair_from_seed(&xi);
                zeroize_bytes(&mut xi);
                Ok(kp)
            }

            /// Algorithm 6: ML-DSA.KeyGen_internal(ξ). The seed must come from an
            /// approved RBG and may be stored as a compact private key.
            pub fn generate_key_pair_from_seed(xi: &[u8; SEED_SIZE]) -> (VerifyingKey, SigningKey) {
                let mut pk = [0u8; VERIFYING_KEY_SIZE];
                let mut sk = [0u8; SIGNING_KEY_SIZE];
                scheme::keygen::<K, L>(&PARAMS, xi, &mut pk, &mut sk);
                let pair = (VerifyingKey(pk), SigningKey(sk));
                zeroize_bytes(&mut sk);
                pair
            }

            impl SigningKey {
                /// Parses a signing key, rejecting encodings whose `s1`/`s2`
                /// coefficients fall outside `[−η, η]`.
                pub fn from_bytes(bytes: &[u8]) -> Result<Self, Error> {
                    let arr: [u8; SIGNING_KEY_SIZE] =
                        bytes.try_into().map_err(|_| Error::InvalidLength)?;
                    if !scheme::sk_is_well_formed::<K, L>(&PARAMS, &arr) {
                        return Err(Error::InvalidSigningKey);
                    }
                    Ok(Self(arr))
                }

                /// Returns the encoded key.
                pub fn as_bytes(&self) -> &[u8; SIGNING_KEY_SIZE] {
                    &self.0
                }

                /// Recomputes the matching verifying key.
                pub fn verifying_key(&self) -> VerifyingKey {
                    let mut pk = [0u8; VERIFYING_KEY_SIZE];
                    scheme::pk_from_sk::<K, L>(&PARAMS, &self.0, &mut pk);
                    VerifyingKey(pk)
                }

                fn tr(&self) -> &[u8] {
                    &self.0[64..128]
                }

                /// Algorithm 2: ML-DSA.Sign (hedged), with randomness from `rng`.
                pub fn sign<R: CryptoRng + ?Sized>(
                    &self,
                    msg: &[u8],
                    ctx: &[u8],
                    rng: &mut R,
                ) -> Result<Signature, Error> {
                    let mut rnd = [0u8; RND_SIZE];
                    rng.fill_bytes(&mut rnd)?;
                    let sig = self.sign_with_rnd(msg, ctx, &rnd);
                    zeroize_bytes(&mut rnd);
                    sig
                }

                /// Algorithm 2: ML-DSA.Sign, deterministic variant (`rnd = {0}^32`).
                pub fn sign_deterministic(&self, msg: &[u8], ctx: &[u8]) -> Result<Signature, Error> {
                    self.sign_with_rnd(msg, ctx, &[0u8; RND_SIZE])
                }

                /// Algorithm 2: ML-DSA.Sign with caller-supplied `rnd`.
                ///
                /// `rnd` must be fresh RBG output (hedged) or all zeros (deterministic).
                pub fn sign_with_rnd(
                    &self,
                    msg: &[u8],
                    ctx: &[u8],
                    rnd: &[u8; RND_SIZE],
                ) -> Result<Signature, Error> {
                    let mu = mu_pure(self.tr(), msg, ctx)?;
                    Ok(self.sign_mu(&mu, rnd))
                }

                /// Algorithm 4: HashML-DSA.Sign (hedged). `msg` is hashed with `alg`.
                pub fn sign_prehash<R: CryptoRng + ?Sized>(
                    &self,
                    msg: &[u8],
                    ctx: &[u8],
                    alg: PreHashAlgorithm,
                    rng: &mut R,
                ) -> Result<Signature, Error> {
                    let mut rnd = [0u8; RND_SIZE];
                    rng.fill_bytes(&mut rnd)?;
                    let sig = self.sign_prehash_with_rnd(msg, ctx, alg, &rnd);
                    zeroize_bytes(&mut rnd);
                    sig
                }

                /// Algorithm 4: HashML-DSA.Sign, deterministic variant.
                pub fn sign_prehash_deterministic(
                    &self,
                    msg: &[u8],
                    ctx: &[u8],
                    alg: PreHashAlgorithm,
                ) -> Result<Signature, Error> {
                    self.sign_prehash_with_rnd(msg, ctx, alg, &[0u8; RND_SIZE])
                }

                /// Algorithm 4: HashML-DSA.Sign with caller-supplied `rnd`.
                pub fn sign_prehash_with_rnd(
                    &self,
                    msg: &[u8],
                    ctx: &[u8],
                    alg: PreHashAlgorithm,
                    rnd: &[u8; RND_SIZE],
                ) -> Result<Signature, Error> {
                    let mut ph = [0u8; 64];
                    let n = alg.digest(msg, &mut ph);
                    self.sign_digest_with_rnd(&ph[..n], ctx, alg, rnd)
                }

                /// Algorithm 4 for a digest `PH(M)` the caller already computed with `alg`.
                pub fn sign_digest_with_rnd(
                    &self,
                    digest: &[u8],
                    ctx: &[u8],
                    alg: PreHashAlgorithm,
                    rnd: &[u8; RND_SIZE],
                ) -> Result<Signature, Error> {
                    let mu = mu_prehashed(self.tr(), digest, ctx, alg)?;
                    Ok(self.sign_mu(&mu, rnd))
                }

                /// Algorithm 7: ML-DSA.Sign_internal on a formatted message `M′`.
                pub fn sign_internal(&self, m_prime: &[u8], rnd: &[u8; RND_SIZE]) -> Signature {
                    let mu = scheme::compute_mu(self.tr(), &[m_prime]);
                    self.sign_mu(&mu, rnd)
                }

                /// Algorithm 7 starting from an externally computed μ
                /// (see [`VerifyingKey::compute_mu`]).
                pub fn sign_mu(&self, mu: &[u8; MU_SIZE], rnd: &[u8; RND_SIZE]) -> Signature {
                    let mut sig = [0u8; SIGNATURE_SIZE];
                    scheme::sign_mu::<K, L>(&PARAMS, &self.0, mu, rnd, &mut sig);
                    Signature(sig)
                }
            }

            impl TryFrom<&[u8]> for SigningKey {
                type Error = Error;
                fn try_from(b: &[u8]) -> Result<Self, Error> {
                    Self::from_bytes(b)
                }
            }

            impl PartialEq for SigningKey {
                fn eq(&self, other: &Self) -> bool {
                    ct_eq(&self.0, &other.0) == 1
                }
            }
            impl Eq for SigningKey {}

            impl fmt::Debug for SigningKey {
                fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                    f.write_str(concat!($name, "::SigningKey(<redacted>)"))
                }
            }

            impl Drop for SigningKey {
                fn drop(&mut self) {
                    zeroize_bytes(&mut self.0);
                }
            }

            impl VerifyingKey {
                /// Parses a verifying key (every byte string of the right length is valid).
                pub fn from_bytes(bytes: &[u8]) -> Result<Self, Error> {
                    bytes.try_into().map(Self).map_err(|_| Error::InvalidLength)
                }

                /// Returns the encoded key.
                pub fn as_bytes(&self) -> &[u8; VERIFYING_KEY_SIZE] {
                    &self.0
                }

                fn tr(&self) -> [u8; 64] {
                    let mut tr = [0u8; 64];
                    crate::hash::shake256(&[&self.0], &mut tr);
                    tr
                }

                /// Computes μ = H(tr ‖ 0 ‖ |ctx| ‖ ctx ‖ M) for external-μ signing.
                pub fn compute_mu(&self, msg: &[u8], ctx: &[u8]) -> Result<[u8; MU_SIZE], Error> {
                    mu_pure(&self.tr(), msg, ctx)
                }

                /// Algorithm 3: ML-DSA.Verify.
                pub fn verify(&self, msg: &[u8], ctx: &[u8], sig: &Signature) -> Result<(), Error> {
                    let mu = mu_pure(&self.tr(), msg, ctx)?;
                    self.verify_mu_result(&mu, sig)
                }

                /// Algorithm 5: HashML-DSA.Verify. `msg` is hashed with `alg`.
                pub fn verify_prehash(
                    &self,
                    msg: &[u8],
                    ctx: &[u8],
                    alg: PreHashAlgorithm,
                    sig: &Signature,
                ) -> Result<(), Error> {
                    let mut ph = [0u8; 64];
                    let n = alg.digest(msg, &mut ph);
                    self.verify_digest(&ph[..n], ctx, alg, sig)
                }

                /// Algorithm 5 for a digest `PH(M)` the caller already computed with `alg`.
                pub fn verify_digest(
                    &self,
                    digest: &[u8],
                    ctx: &[u8],
                    alg: PreHashAlgorithm,
                    sig: &Signature,
                ) -> Result<(), Error> {
                    let mu = mu_prehashed(&self.tr(), digest, ctx, alg)?;
                    self.verify_mu_result(&mu, sig)
                }

                /// Algorithm 8: ML-DSA.Verify_internal on a formatted message `M′`.
                pub fn verify_internal(&self, m_prime: &[u8], sig: &Signature) -> bool {
                    let mu = scheme::compute_mu(&self.tr(), &[m_prime]);
                    self.verify_mu(&mu, sig)
                }

                /// Algorithm 8 starting from an externally computed μ.
                pub fn verify_mu(&self, mu: &[u8; MU_SIZE], sig: &Signature) -> bool {
                    scheme::verify_mu::<K, L>(&PARAMS, &self.0, mu, &sig.0)
                }

                fn verify_mu_result(&self, mu: &[u8; MU_SIZE], sig: &Signature) -> Result<(), Error> {
                    if self.verify_mu(mu, sig) {
                        Ok(())
                    } else {
                        Err(Error::InvalidSignature)
                    }
                }
            }

            impl TryFrom<&[u8]> for VerifyingKey {
                type Error = Error;
                fn try_from(b: &[u8]) -> Result<Self, Error> {
                    Self::from_bytes(b)
                }
            }

            impl fmt::Debug for VerifyingKey {
                fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                    write!(f, concat!($name, "::VerifyingKey({:02x?}..)"), &self.0[..8])
                }
            }

            impl Signature {
                /// Parses a signature (length check only; structure is checked on verify).
                pub fn from_bytes(bytes: &[u8]) -> Result<Self, Error> {
                    bytes.try_into().map(Self).map_err(|_| Error::InvalidLength)
                }

                /// Returns the encoded signature.
                pub fn as_bytes(&self) -> &[u8; SIGNATURE_SIZE] {
                    &self.0
                }
            }

            impl From<[u8; SIGNATURE_SIZE]> for Signature {
                fn from(b: [u8; SIGNATURE_SIZE]) -> Self {
                    Self(b)
                }
            }

            impl TryFrom<&[u8]> for Signature {
                type Error = Error;
                fn try_from(b: &[u8]) -> Result<Self, Error> {
                    Self::from_bytes(b)
                }
            }

            impl fmt::Debug for Signature {
                fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                    write!(f, concat!($name, "::Signature({:02x?}..)"), &self.0[..8])
                }
            }
        }
    };
}

ml_dsa_parameter_set!(
    /// ML-DSA-44 (NIST security category 2).
    ml_dsa_44, "ML-DSA-44", k = 4, l = 4,
    scheme::Params { eta: 2, tau: 39, beta: 78, gamma1: 1 << 17, gamma2: (arith::Q - 1) / 88, omega: 80, lambda: 128 }
);
ml_dsa_parameter_set!(
    /// ML-DSA-65 (NIST security category 3).
    ml_dsa_65, "ML-DSA-65", k = 6, l = 5,
    scheme::Params { eta: 4, tau: 49, beta: 196, gamma1: 1 << 19, gamma2: (arith::Q - 1) / 32, omega: 55, lambda: 192 }
);
ml_dsa_parameter_set!(
    /// ML-DSA-87 (NIST security category 5).
    ml_dsa_87, "ML-DSA-87", k = 8, l = 7,
    scheme::Params { eta: 2, tau: 60, beta: 120, gamma1: 1 << 19, gamma2: (arith::Q - 1) / 32, omega: 75, lambda: 256 }
);

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sizes_match_fips204_table2() {
        assert_eq!(ml_dsa_44::VERIFYING_KEY_SIZE, 1312);
        assert_eq!(ml_dsa_44::SIGNING_KEY_SIZE, 2560);
        assert_eq!(ml_dsa_44::SIGNATURE_SIZE, 2420);
        assert_eq!(ml_dsa_65::VERIFYING_KEY_SIZE, 1952);
        assert_eq!(ml_dsa_65::SIGNING_KEY_SIZE, 4032);
        assert_eq!(ml_dsa_65::SIGNATURE_SIZE, 3309);
        assert_eq!(ml_dsa_87::VERIFYING_KEY_SIZE, 2592);
        assert_eq!(ml_dsa_87::SIGNING_KEY_SIZE, 4896);
        assert_eq!(ml_dsa_87::SIGNATURE_SIZE, 4627);
    }

    #[test]
    fn prehash_oids_and_lengths() {
        for alg in PreHashAlgorithm::ALL {
            let oid = alg.oid();
            assert_eq!(&oid[..10], &[0x06, 0x09, 0x60, 0x86, 0x48, 0x01, 0x65, 0x03, 0x04, 0x02]);
            let mut out = [0u8; 64];
            assert_eq!(alg.digest(b"abc", &mut out), alg.digest_len());
        }
        assert_eq!(PreHashAlgorithm::Sha2_256.oid()[10], 1);
        assert_eq!(PreHashAlgorithm::Shake256.oid()[10], 12);
    }
}
