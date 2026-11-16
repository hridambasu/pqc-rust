//! ML-KEM, the Module-Lattice-Based Key-Encapsulation Mechanism (FIPS 203).
//!
//! Three parameter sets are provided, each in its own module:
//!
//! | Module                     | Security category | ek    | dk    | ciphertext |
//! |----------------------------|-------------------|-------|-------|------------|
//! | [`ml_kem_512`]             | 1                 | 800   | 1632  | 768        |
//! | [`ml_kem_768`]             | 3                 | 1184  | 2400  | 1088       |
//! | [`ml_kem_1024`]            | 5                 | 1568  | 3168  | 1568       |
//!
//! ```
//! # #[cfg(all(feature = "std", unix))] {
//! use fips_pqc::{ml_kem::ml_kem_768, OsRng};
//!
//! let (ek, dk) = ml_kem_768::generate_key_pair(&mut OsRng).unwrap();
//! let (ct, ss_sender) = ek.encapsulate(&mut OsRng).unwrap();
//! let ss_receiver = dk.decapsulate(&ct);
//! assert_eq!(ss_sender, ss_receiver);
//! # }
//! ```
//!
//! Input checking follows FIPS 203 §7: [`ml_kem_768::EncapsulationKey::from_bytes`]
//! performs the type and modulus checks, [`ml_kem_768::DecapsulationKey::from_bytes`]
//! performs the type and hash checks, and the ciphertext type check is enforced by
//! [`ml_kem_768::Ciphertext`]'s fixed length.

mod arith;
mod kpke;

use crate::common::{ct_eq, ct_select, zeroize_bytes};
use crate::hash::{sha3_256, sha3_512, shake256};
use core::fmt;


/// Size in bytes of the shared secret key `K` (all parameter sets).
pub const SHARED_SECRET_SIZE: usize = 32;
/// Size in bytes of the 64-byte seed `d ‖ z` from which a key pair is derived.
pub const SEED_SIZE: usize = 64;


/// A 32-byte shared secret. Compared in constant time and zeroized on drop.
#[derive(Clone)]
pub struct SharedSecret([u8; SHARED_SECRET_SIZE]);

impl SharedSecret {
    /// Returns the secret bytes.
    pub fn as_bytes(&self) -> &[u8; SHARED_SECRET_SIZE] {
        &self.0
    }
}

impl PartialEq for SharedSecret {
    fn eq(&self, other: &Self) -> bool {
        ct_eq(&self.0, &other.0) == 1
    }
}
impl Eq for SharedSecret {}

impl fmt::Debug for SharedSecret {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("SharedSecret(<redacted>)")
    }
}

impl Drop for SharedSecret {
    fn drop(&mut self) {
        zeroize_bytes(&mut self.0);
    }
}

/// Largest ciphertext of any parameter set (ML-KEM-1024).
const MAX_CT: usize = 1568;
const _: () = assert!(MAX_CT >= ml_kem_1024::CIPHERTEXT_SIZE);

/// Algorithm 16: ML-KEM.KeyGen_internal(d, z).
fn keygen_internal<const K: usize>(
    p: &kpke::Params,
    d: &[u8; 32],
    z: &[u8; 32],
    ek: &mut [u8],
    dk: &mut [u8],
) {
    kpke::keygen::<K>(p, d, ek, &mut dk[..384 * K]);
    dk[384 * K..768 * K + 32].copy_from_slice(ek);
    dk[768 * K + 32..768 * K + 64].copy_from_slice(&sha3_256(&[ek]));
    dk[768 * K + 64..].copy_from_slice(z);
}

/// Algorithm 17: ML-KEM.Encaps_internal(ek, m). Writes `c`, returns `K`.
fn encaps_internal<const K: usize>(
    p: &kpke::Params,
    ek: &[u8],
    m: &[u8; 32],
    c: &mut [u8],
) -> SharedSecret {
    let h = sha3_256(&[ek]);
    let mut g = sha3_512(&[m, &h]);
    let mut k = [0u8; 32];
    let mut r = [0u8; 32];
    k.copy_from_slice(&g[..32]);
    r.copy_from_slice(&g[32..]);
    kpke::encrypt::<K>(p, ek, m, &r, c);
    zeroize_bytes(&mut g);
    zeroize_bytes(&mut r);
    let ss = SharedSecret(k);
    zeroize_bytes(&mut k);
    ss
}

/// Algorithm 18: ML-KEM.Decaps_internal(dk, c), with implicit rejection.
fn decaps_internal<const K: usize>(p: &kpke::Params, dk: &[u8], c: &[u8]) -> SharedSecret {
    let dk_pke = &dk[..384 * K];
    let ek_pke = &dk[384 * K..768 * K + 32];
    let h = &dk[768 * K + 32..768 * K + 64];
    let z = &dk[768 * K + 64..768 * K + 96];

    let mut m_prime = kpke::decrypt::<K>(p, dk_pke, c);
    let mut g = sha3_512(&[&m_prime, h]);
    let mut r_prime = [0u8; 32];
    r_prime.copy_from_slice(&g[32..]);

    let mut k_out = [0u8; 32];
    shake256(&[z, c], &mut k_out); // K̄ = J(z ‖ c)

    let mut c_prime_buf = [0u8; MAX_CT];
    let c_prime = &mut c_prime_buf[..c.len()];
    kpke::encrypt::<K>(p, ek_pke, &m_prime, &r_prime, c_prime);

    // If c == c′ return K′, otherwise K̄ — selected in constant time.
    let equal = ct_eq(c, c_prime);
    ct_select(&mut k_out, &g[..32], equal);

    zeroize_bytes(&mut m_prime);
    zeroize_bytes(&mut g);
    zeroize_bytes(&mut r_prime);
    let ss = SharedSecret(k_out);
    zeroize_bytes(&mut k_out);
    ss
}

/// FIPS 203 §7.2 modulus check: ByteEncode12(ByteDecode12(ek)) == ek.
fn ek_modulus_check<const K: usize>(ek: &[u8]) -> bool {
    let mut buf = [0u8; 384];
    for chunk in ek[..384 * K].chunks_exact(384) {
        let f = arith::byte_decode(12, chunk);
        arith::byte_encode(12, &f, &mut buf);
        if buf[..] != chunk[..] {
            return false;
        }
    }
    true
}

/// FIPS 203 §7.3 hash check: H(dk[384k : 768k+32]) == dk[768k+32 : 768k+64].
fn dk_hash_check<const K: usize>(dk: &[u8]) -> bool {
    let h = sha3_256(&[&dk[384 * K..768 * K + 32]]);
    ct_eq(&h, &dk[768 * K + 32..768 * K + 64]) == 1
}

macro_rules! ml_kem_parameter_set {
    (
        $(#[$meta:meta])*
        $module:ident, $name:literal, k = $k:literal, eta1 = $eta1:literal, eta2 = $eta2:literal,
        du = $du:literal, dv = $dv:literal
    ) => {
        $(#[$meta])*
        pub mod $module {
            use super::*;
            use crate::{CryptoRng, Error};

            const K: usize = $k;
            const PARAMS: kpke::Params = kpke::Params { eta1: $eta1, eta2: $eta2, du: $du, dv: $dv };

            #[doc = concat!("Size in bytes of an ", $name, " encapsulation key.")]
            pub const ENCAPSULATION_KEY_SIZE: usize = 384 * K + 32;
            #[doc = concat!("Size in bytes of an ", $name, " decapsulation key.")]
            pub const DECAPSULATION_KEY_SIZE: usize = 768 * K + 96;
            #[doc = concat!("Size in bytes of an ", $name, " ciphertext.")]
            pub const CIPHERTEXT_SIZE: usize = 32 * ($du * K + $dv);

            #[doc = concat!("An ", $name, " encapsulation (public) key that has passed the FIPS 203 §7.2 checks.")]
            #[derive(Clone, PartialEq, Eq)]
            pub struct EncapsulationKey([u8; ENCAPSULATION_KEY_SIZE]);

            #[doc = concat!("An ", $name, " decapsulation (private) key. Zeroized on drop.")]
            #[derive(Clone)]
            pub struct DecapsulationKey([u8; DECAPSULATION_KEY_SIZE]);

            #[doc = concat!("An ", $name, " ciphertext.")]
            #[derive(Clone, PartialEq, Eq)]
            pub struct Ciphertext([u8; CIPHERTEXT_SIZE]);

            /// Algorithm 19: ML-KEM.KeyGen, using randomness from `rng`.
            pub fn generate_key_pair<R: CryptoRng + ?Sized>(
                rng: &mut R,
            ) -> Result<(EncapsulationKey, DecapsulationKey), Error> {
                let mut seed = [0u8; SEED_SIZE];
                rng.fill_bytes(&mut seed)?;
                let kp = generate_key_pair_from_seed(&seed);
                zeroize_bytes(&mut seed);
                Ok(kp)
            }

            /// Derives a key pair deterministically from the 64-byte seed `d ‖ z`
            /// (ML-KEM.KeyGen_internal). The seed must come from an approved RBG
            /// and is itself a valid private-key representation.
            pub fn generate_key_pair_from_seed(seed: &[u8; SEED_SIZE]) -> (EncapsulationKey, DecapsulationKey) {
                let mut d = [0u8; 32];
                let mut z = [0u8; 32];
                d.copy_from_slice(&seed[..32]);
                z.copy_from_slice(&seed[32..]);
                let mut ek = [0u8; ENCAPSULATION_KEY_SIZE];
                let mut dk = [0u8; DECAPSULATION_KEY_SIZE];
                keygen_internal::<K>(&PARAMS, &d, &z, &mut ek, &mut dk);
                zeroize_bytes(&mut d);
                zeroize_bytes(&mut z);
                let pair = (EncapsulationKey(ek), DecapsulationKey(dk));
                zeroize_bytes(&mut dk);
                pair
            }

            impl EncapsulationKey {
                /// Parses an encapsulation key, performing the FIPS 203 §7.2 type
                /// (length) and modulus checks.
                pub fn from_bytes(bytes: &[u8]) -> Result<Self, Error> {
                    let arr: [u8; ENCAPSULATION_KEY_SIZE] =
                        bytes.try_into().map_err(|_| Error::InvalidLength)?;
                    if !ek_modulus_check::<K>(&arr) {
                        return Err(Error::InvalidEncapsulationKey);
                    }
                    Ok(Self(arr))
                }

                /// Returns the encoded key.
                pub fn as_bytes(&self) -> &[u8; ENCAPSULATION_KEY_SIZE] {
                    &self.0
                }

                /// Algorithm 20: ML-KEM.Encaps. Returns the ciphertext and shared secret.
                pub fn encapsulate<R: CryptoRng + ?Sized>(
                    &self,
                    rng: &mut R,
                ) -> Result<(Ciphertext, SharedSecret), Error> {
                    let mut m = [0u8; 32];
                    rng.fill_bytes(&mut m)?;
                    let out = self.encapsulate_internal(&m);
                    zeroize_bytes(&mut m);
                    Ok(out)
                }

                /// Algorithm 17: ML-KEM.Encaps_internal with caller-supplied randomness `m`.
                ///
                /// FIPS 203 only permits this for testing; `m` must otherwise be
                /// 32 fresh bytes from an approved RBG.
                pub fn encapsulate_internal(&self, m: &[u8; 32]) -> (Ciphertext, SharedSecret) {
                    let mut c = [0u8; CIPHERTEXT_SIZE];
                    let k = encaps_internal::<K>(&PARAMS, &self.0, m, &mut c);
                    (Ciphertext(c), k)
                }
            }

            impl TryFrom<&[u8]> for EncapsulationKey {
                type Error = Error;
                fn try_from(b: &[u8]) -> Result<Self, Error> {
                    Self::from_bytes(b)
                }
            }

            impl fmt::Debug for EncapsulationKey {
                fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                    write!(f, concat!($name, "::EncapsulationKey({:02x?}..)"), &self.0[..8])
                }
            }

            impl DecapsulationKey {
                /// Parses a decapsulation key, performing the FIPS 203 §7.3 type
                /// (length) and hash checks.
                pub fn from_bytes(bytes: &[u8]) -> Result<Self, Error> {
                    let arr: [u8; DECAPSULATION_KEY_SIZE] =
                        bytes.try_into().map_err(|_| Error::InvalidLength)?;
                    if !dk_hash_check::<K>(&arr) {
                        return Err(Error::InvalidDecapsulationKey);
                    }
                    Ok(Self(arr))
                }

                /// Returns the encoded key.
                pub fn as_bytes(&self) -> &[u8; DECAPSULATION_KEY_SIZE] {
                    &self.0
                }

                /// Returns the encapsulation key embedded in this decapsulation key.
                pub fn encapsulation_key(&self) -> EncapsulationKey {
                    let mut ek = [0u8; ENCAPSULATION_KEY_SIZE];
                    ek.copy_from_slice(&self.0[384 * K..768 * K + 32]);
                    EncapsulationKey(ek)
                }

                /// Algorithm 21: ML-KEM.Decaps. Never fails: an invalid ciphertext
                /// yields a pseudorandom key (implicit rejection).
                pub fn decapsulate(&self, c: &Ciphertext) -> SharedSecret {
                    decaps_internal::<K>(&PARAMS, &self.0, &c.0)
                }
            }

            impl TryFrom<&[u8]> for DecapsulationKey {
                type Error = Error;
                fn try_from(b: &[u8]) -> Result<Self, Error> {
                    Self::from_bytes(b)
                }
            }

            impl PartialEq for DecapsulationKey {
                fn eq(&self, other: &Self) -> bool {
                    ct_eq(&self.0, &other.0) == 1
                }
            }
            impl Eq for DecapsulationKey {}

            impl fmt::Debug for DecapsulationKey {
                fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                    f.write_str(concat!($name, "::DecapsulationKey(<redacted>)"))
                }
            }

            impl Drop for DecapsulationKey {
                fn drop(&mut self) {
                    zeroize_bytes(&mut self.0);
                }
            }

            impl Ciphertext {
                /// Parses a ciphertext (FIPS 203 §7.3 ciphertext type check).
                pub fn from_bytes(bytes: &[u8]) -> Result<Self, Error> {
                    bytes.try_into().map(Self).map_err(|_| Error::InvalidLength)
                }

                /// Returns the encoded ciphertext.
                pub fn as_bytes(&self) -> &[u8; CIPHERTEXT_SIZE] {
                    &self.0
                }
            }

            impl From<[u8; CIPHERTEXT_SIZE]> for Ciphertext {
                fn from(b: [u8; CIPHERTEXT_SIZE]) -> Self {
                    Self(b)
                }
            }

            impl TryFrom<&[u8]> for Ciphertext {
                type Error = Error;
                fn try_from(b: &[u8]) -> Result<Self, Error> {
                    Self::from_bytes(b)
                }
            }

            impl fmt::Debug for Ciphertext {
                fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                    write!(f, concat!($name, "::Ciphertext({:02x?}..)"), &self.0[..8])
                }
            }
        }
    };
}

ml_kem_parameter_set!(
    /// ML-KEM-512 (NIST security category 1).
    ml_kem_512, "ML-KEM-512", k = 2, eta1 = 3, eta2 = 2, du = 10, dv = 4
);
ml_kem_parameter_set!(
    /// ML-KEM-768 (NIST security category 3).
    ml_kem_768, "ML-KEM-768", k = 3, eta1 = 2, eta2 = 2, du = 10, dv = 4
);
ml_kem_parameter_set!(
    /// ML-KEM-1024 (NIST security category 5).
    ml_kem_1024, "ML-KEM-1024", k = 4, eta1 = 2, eta2 = 2, du = 11, dv = 5
);

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sizes_match_fips203_table3() {
        assert_eq!(ml_kem_512::ENCAPSULATION_KEY_SIZE, 800);
        assert_eq!(ml_kem_512::DECAPSULATION_KEY_SIZE, 1632);
        assert_eq!(ml_kem_512::CIPHERTEXT_SIZE, 768);
        assert_eq!(ml_kem_768::ENCAPSULATION_KEY_SIZE, 1184);
        assert_eq!(ml_kem_768::DECAPSULATION_KEY_SIZE, 2400);
        assert_eq!(ml_kem_768::CIPHERTEXT_SIZE, 1088);
        assert_eq!(ml_kem_1024::ENCAPSULATION_KEY_SIZE, 1568);
        assert_eq!(ml_kem_1024::DECAPSULATION_KEY_SIZE, 3168);
        assert_eq!(ml_kem_1024::CIPHERTEXT_SIZE, 1568);
    }
}
