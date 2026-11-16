//! SHA-3 hash functions and SHAKE extendable-output functions (FIPS 202).

use core::hint::black_box;
use core::sync::atomic::{compiler_fence, Ordering};

const RC: [u64; 24] = [
    0x0000_0000_0000_0001,
    0x0000_0000_0000_8082,
    0x8000_0000_0000_808A,
    0x8000_0000_8000_8000,
    0x0000_0000_0000_808B,
    0x0000_0000_8000_0001,
    0x8000_0000_8000_8081,
    0x8000_0000_0000_8009,
    0x0000_0000_0000_008A,
    0x0000_0000_0000_0088,
    0x0000_0000_8000_8009,
    0x0000_0000_8000_000A,
    0x0000_0000_8000_808B,
    0x8000_0000_0000_008B,
    0x8000_0000_0000_8089,
    0x8000_0000_0000_8003,
    0x8000_0000_0000_8002,
    0x8000_0000_0000_0080,
    0x0000_0000_0000_800A,
    0x8000_0000_8000_000A,
    0x8000_0000_8000_8081,
    0x8000_0000_0000_8080,
    0x0000_0000_8000_0001,
    0x8000_0000_8000_8008,
];

/// Rotation offsets, indexed by `x + 5*y`.
const RHO: [u32; 25] = [
    0, 1, 62, 28, 27, 36, 44, 6, 55, 20, 3, 10, 43, 25, 39, 41, 45, 15, 21, 8, 18, 2, 61, 56, 14,
];

/// The Keccak-p[1600, 24] permutation (FIPS 202 §3.3).
pub(crate) fn keccak_f1600(a: &mut [u64; 25]) {
    for &rc in RC.iter() {
        // θ
        let mut c = [0u64; 5];
        for x in 0..5 {
            c[x] = a[x] ^ a[x + 5] ^ a[x + 10] ^ a[x + 15] ^ a[x + 20];
        }
        for x in 0..5 {
            let d = c[(x + 4) % 5] ^ c[(x + 1) % 5].rotate_left(1);
            for y in 0..5 {
                a[x + 5 * y] ^= d;
            }
        }
        // ρ and π
        let mut b = [0u64; 25];
        for x in 0..5 {
            for y in 0..5 {
                b[y + 5 * ((2 * x + 3 * y) % 5)] = a[x + 5 * y].rotate_left(RHO[x + 5 * y]);
            }
        }
        // χ
        for y in 0..5 {
            for x in 0..5 {
                a[x + 5 * y] = b[x + 5 * y] ^ (!b[(x + 1) % 5 + 5 * y] & b[(x + 2) % 5 + 5 * y]);
            }
        }
        // ι
        a[0] ^= rc;
    }
}

/// A Keccak sponge with a byte-granular absorb/squeeze interface.
///
/// Absorbing after squeezing has started is a logic error.
#[derive(Clone)]
pub(crate) struct Sponge {
    state: [u64; 25],
    rate: usize,
    pos: usize,
    ds: u8,
    squeezing: bool,
}

impl Sponge {
    const fn new(rate: usize, ds: u8) -> Self {
        Self { state: [0; 25], rate, pos: 0, ds, squeezing: false }
    }
    pub(crate) const fn shake128() -> Self {
        Self::new(168, 0x1F)
    }
    pub(crate) const fn shake256() -> Self {
        Self::new(136, 0x1F)
    }
    pub(crate) const fn sha3_224() -> Self {
        Self::new(144, 0x06)
    }
    pub(crate) const fn sha3_256() -> Self {
        Self::new(136, 0x06)
    }
    pub(crate) const fn sha3_384() -> Self {
        Self::new(104, 0x06)
    }
    pub(crate) const fn sha3_512() -> Self {
        Self::new(72, 0x06)
    }

    pub(crate) fn absorb(&mut self, mut data: &[u8]) {
        assert!(!self.squeezing, "absorb called after squeeze");
        while !data.is_empty() {
            if self.pos == 0 && data.len() >= self.rate {
                for (lane, chunk) in self.state.iter_mut().zip(data[..self.rate].chunks_exact(8)) {
                    let mut w = [0u8; 8];
                    w.copy_from_slice(chunk);
                    *lane ^= u64::from_le_bytes(w);
                }
                keccak_f1600(&mut self.state);
                data = &data[self.rate..];
                continue;
            }
            let take = (self.rate - self.pos).min(data.len());
            for (k, &b) in data[..take].iter().enumerate() {
                let p = self.pos + k;
                self.state[p / 8] ^= u64::from(b) << (8 * (p % 8));
            }
            self.pos += take;
            data = &data[take..];
            if self.pos == self.rate {
                keccak_f1600(&mut self.state);
                self.pos = 0;
            }
        }
    }

    fn pad(&mut self) {
        let p = self.pos;
        self.state[p / 8] ^= u64::from(self.ds) << (8 * (p % 8));
        let last = self.rate - 1;
        self.state[last / 8] ^= 0x80u64 << (8 * (last % 8));
        keccak_f1600(&mut self.state);
        self.pos = 0;
        self.squeezing = true;
    }

    /// Squeezes output; the first call finalizes the absorb phase.
    pub(crate) fn squeeze(&mut self, out: &mut [u8]) {
        if !self.squeezing {
            self.pad();
        }
        for b in out.iter_mut() {
            if self.pos == self.rate {
                keccak_f1600(&mut self.state);
                self.pos = 0;
            }
            *b = (self.state[self.pos / 8] >> (8 * (self.pos % 8))) as u8;
            self.pos += 1;
        }
    }
}

impl Drop for Sponge {
    fn drop(&mut self) {
        for lane in self.state.iter_mut() {
            *lane = 0;
        }
        black_box(&self.state);
        compiler_fence(Ordering::SeqCst);
    }
}

macro_rules! sha3_fixed {
    ($name:ident, $ctor:ident, $out:expr, $doc:expr) => {
        #[doc = $doc]
        #[derive(Clone)]
        pub struct $name(Sponge);

        impl $name {
            /// Digest length in bytes.
            pub const OUTPUT_SIZE: usize = $out;

            /// Creates a fresh hash state.
            pub const fn new() -> Self {
                Self(Sponge::$ctor())
            }

            /// Absorbs more input.
            pub fn update(&mut self, data: &[u8]) {
                self.0.absorb(data);
            }

            /// Finishes the computation and returns the digest.
            pub fn finalize(mut self) -> [u8; $out] {
                let mut out = [0u8; $out];
                self.0.squeeze(&mut out);
                out
            }

            /// One-shot digest of `data`.
            pub fn digest(data: &[u8]) -> [u8; $out] {
                let mut h = Self::new();
                h.update(data);
                h.finalize()
            }
        }

        impl Default for $name {
            fn default() -> Self {
                Self::new()
            }
        }
    };
}

sha3_fixed!(Sha3_224, sha3_224, 28, "SHA3-224 (FIPS 202).");
sha3_fixed!(Sha3_256, sha3_256, 32, "SHA3-256 (FIPS 202).");
sha3_fixed!(Sha3_384, sha3_384, 48, "SHA3-384 (FIPS 202).");
sha3_fixed!(Sha3_512, sha3_512, 64, "SHA3-512 (FIPS 202).");

/// Output stream of a finalized SHAKE instance.
#[derive(Clone)]
pub struct XofReader(Sponge);

impl XofReader {
    /// Fills `out` with the next bytes of the output stream.
    pub fn read(&mut self, out: &mut [u8]) {
        self.0.squeeze(out);
    }
}

macro_rules! shake {
    ($name:ident, $ctor:ident, $doc:expr) => {
        #[doc = $doc]
        #[derive(Clone)]
        pub struct $name(Sponge);

        impl $name {
            /// Creates a fresh XOF state.
            pub const fn new() -> Self {
                Self(Sponge::$ctor())
            }

            /// Absorbs more input.
            pub fn update(&mut self, data: &[u8]) {
                self.0.absorb(data);
            }

            /// Finishes absorbing and returns a reader for the output stream.
            pub fn finalize_xof(self) -> XofReader {
                XofReader(self.0)
            }

            /// One-shot: absorbs `data` and fills `out` with output.
            pub fn digest(data: &[u8], out: &mut [u8]) {
                let mut h = Self::new();
                h.update(data);
                h.finalize_xof().read(out);
            }
        }

        impl Default for $name {
            fn default() -> Self {
                Self::new()
            }
        }
    };
}

shake!(Shake128, shake128, "SHAKE128 extendable-output function (FIPS 202).");
shake!(Shake256, shake256, "SHAKE256 extendable-output function (FIPS 202).");

/// SHA3-256 of the concatenation of `parts` (the function `H` of FIPS 203).
pub(crate) fn sha3_256(parts: &[&[u8]]) -> [u8; 32] {
    let mut s = Sponge::sha3_256();
    for p in parts {
        s.absorb(p);
    }
    let mut out = [0u8; 32];
    s.squeeze(&mut out);
    out
}

/// SHA3-512 of the concatenation of `parts` (the function `G` of FIPS 203).
pub(crate) fn sha3_512(parts: &[&[u8]]) -> [u8; 64] {
    let mut s = Sponge::sha3_512();
    for p in parts {
        s.absorb(p);
    }
    let mut out = [0u8; 64];
    s.squeeze(&mut out);
    out
}

/// SHAKE256 of the concatenation of `parts` into `out`.
pub(crate) fn shake256(parts: &[&[u8]], out: &mut [u8]) {
    let mut s = Sponge::shake256();
    for p in parts {
        s.absorb(p);
    }
    s.squeeze(out);
}

#[cfg(test)]
mod tests {
    use super::*;

    fn hex(b: &[u8]) -> String {
        use core::fmt::Write;
        let mut s = String::new();
        for x in b {
            write!(s, "{x:02x}").unwrap();
        }
        s
    }

    #[test]
    fn empty_string_digests() {
        assert_eq!(
            hex(&Sha3_256::digest(b"")),
            "a7ffc6f8bf1ed76651c14756a061d662f580ff4de43b49fa82d80a4b80f8434a"
        );
        assert_eq!(
            hex(&Sha3_512::digest(b"")),
            "a69f73cca23a9ac5c8b567dc185a756e97c982164fe25859e0d1dcc1475c80a6\
             15b2123af1f5f94c11e3e9402c3ac558f500199d95b6d3e301758586281dcd26"
        );
        let mut o = [0u8; 32];
        Shake128::digest(b"", &mut o);
        assert_eq!(hex(&o), "7f9c2ba4e88f827d616045507605853ed73b8093f6efbc88eb1a6eacfa66ef26");
        Shake256::digest(b"", &mut o);
        assert_eq!(hex(&o), "46b9dd2b0ba88d13233b3feb743eeb243fcd52ea62b81b82b50c27646ed5762f");
    }

    #[test]
    fn abc() {
        assert_eq!(
            hex(&Sha3_256::digest(b"abc")),
            "3a985da74fe225b2045c172d6bd390bd855f086e3e9d525b46bfe24511431532"
        );
    }

    #[test]
    fn incremental_squeeze_matches_one_shot() {
        let mut one = [0u8; 1000];
        Shake128::digest(b"split me", &mut one);
        let mut r = {
            let mut h = Shake128::new();
            h.update(b"split ");
            h.update(b"me");
            h.finalize_xof()
        };
        let mut many = [0u8; 1000];
        for chunk in many.chunks_mut(7) {
            r.read(chunk);
        }
        assert_eq!(one, many);
    }
}
