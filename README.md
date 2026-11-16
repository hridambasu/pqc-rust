# fips-pqc

A pure-Rust implementation of the NIST post-quantum standards, with no dependencies and `no_std` support:

- **ML-KEM** ([FIPS 203](https://csrc.nist.gov/pubs/fips/203/final)): ML-KEM-512, ML-KEM-768, ML-KEM-1024
- **ML-DSA** ([FIPS 204](https://csrc.nist.gov/pubs/fips/204/final)): ML-DSA-44, ML-DSA-65, ML-DSA-87, covering pure ML-DSA, HashML-DSA (all 12 approved pre-hashes) and the internal / external-μ interfaces
- The SHA-3/SHAKE (FIPS 202) and SHA-2 (FIPS 180-4) functions these standards need are implemented in `fips_pqc::hash`

## Usage

```rust
use fips_pqc::{ml_dsa::ml_dsa_65, ml_kem::ml_kem_768, OsRng};

// Key encapsulation
let (ek, dk) = ml_kem_768::generate_key_pair(&mut OsRng)?;
let (ct, shared_sender) = ek.encapsulate(&mut OsRng)?;
assert_eq!(dk.decapsulate(&ct), shared_sender);

// Signatures
let (vk, sk) = ml_dsa_65::generate_key_pair(&mut OsRng)?;
let sig = sk.sign(b"message", b"context", &mut OsRng)?;
vk.verify(b"message", b"context", &sig)?;
# Ok::<(), fips_pqc::Error>(())
```

See `examples/demo.rs`. Run `cargo run --release --example bench` for timings.

## How it maps to the standards

| FIPS function | API |
|---|---|
| ML-KEM.KeyGen / KeyGen_internal | `generate_key_pair(rng)` / `generate_key_pair_from_seed(&[u8; 64])` (seed = `d ‖ z`) |
| ML-KEM.Encaps / Encaps_internal | `EncapsulationKey::encapsulate(rng)` / `encapsulate_internal(&m)` |
| ML-KEM.Decaps | `DecapsulationKey::decapsulate(&ct)`: implicit rejection, constant-time select |
| §7.2 encapsulation key check (type and modulus) | `EncapsulationKey::from_bytes` |
| §7.3 decapsulation input check (type and hash) | `DecapsulationKey::from_bytes`, `Ciphertext::from_bytes` |
| ML-DSA.KeyGen / KeyGen_internal | `generate_key_pair(rng)` / `generate_key_pair_from_seed(&ξ)` |
| ML-DSA.Sign (hedged / deterministic) | `sign(msg, ctx, rng)` / `sign_deterministic(msg, ctx)` / `sign_with_rnd` |
| HashML-DSA.Sign / Verify | `sign_prehash*`, `sign_digest_with_rnd`, `verify_prehash`, `verify_digest` |
| ML-DSA.Sign_internal / Verify_internal | `sign_internal(M′, rnd)` / `verify_internal(M′, sig)` |
| External μ | `VerifyingKey::compute_mu`, `SigningKey::sign_mu`, `VerifyingKey::verify_mu` |

## Testing

`cargo test` runs:

- **NIST ACVP vectors** (`tests/vectors/`, converted from [usnistgov/ACVP-Server](https://github.com/usnistgov/ACVP-Server) `internalProjection.json`):
  - ML-KEM: 75 keyGen, 75 encapsulation, 30 decapsulation (valid and modified ciphertexts), 30 encapsulation-key checks and 30 decapsulation-key checks
  - ML-DSA: 75 keyGen, 360 sigGen (pure, pre-hash, internal and external-μ, each deterministic and hedged) and 180 sigVer (valid and tampered)
- **Hash vectors:** 4,104 SHA-2/SHA-3/SHAKE digests from CPython `hashlib`, covering every input length from 0 to 299 and every rate boundary. Each one is also checked with chunked input and output.
- **Unit tests:** exhaustive checks of Decompose against the spec, of Barrett and division-by-q, and of Compress/Decompress against the rational definitions. Also NTT against schoolbook multiplication, the MakeHint/UseHint identity, bit-packing round trips and hint-encoding rejection rules.
- **Behavioral tests for every parameter set:** round trips; implicit rejection (checks K̄ = J(z ‖ c)); modulus and hash check failures; length errors; RNG failure propagation; context-length limits; domain separation between pure and pre-hash signing and between hash algorithms; bit-flip tampering across the whole signature; malformed hints; out-of-range `z`; 1 MiB messages.

## Security notes

- The code has no `unsafe`. Secret-dependent operations use no branches, no secret-indexed table lookups and no hardware division; all reduction is Barrett-based (this avoids KyberSlash-style issues).
- Secret keys, shared secrets and secret intermediates are zeroized on drop or after use.
- Rejection-sampling loops reveal only the iteration count, which FIPS 204 permits.
- **This code has not been independently audited or CMVP-validated.** Get a review before using it in production.

## License

MIT OR Apache-2.0
