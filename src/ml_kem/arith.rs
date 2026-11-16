//! Arithmetic over R_q = Z_q[X]/(X^256 + 1) with q = 3329 (FIPS 203 §4.3),
//! plus the sampling, compression and encoding routines of §4.2.
//!
//! Coefficients are kept in canonical form `[0, q)` as `u32`. All operations
//! on secret data are branch-free and avoid variable-time division.

use crate::common::{bit_rev, pack, pow_mod, unpack};
use crate::hash::Sponge;

pub(crate) const N: usize = 256;
pub(crate) const Q: u32 = 3329;
pub(crate) type Poly = [u32; N];

/// `floor(2^32 / q)`, the Barrett constant.
const BARRETT_M: u64 = (1u64 << 32) / Q as u64;
/// `128^{-1} mod q`, the final scaling factor of NTT^{-1}.
const INV_128: u32 = 3303;

/// `ZETAS[i] = 17^{BitRev7(i)} mod q`.
pub(crate) const ZETAS: [u32; 128] = {
    let mut z = [0u32; 128];
    let mut i = 0;
    while i < 128 {
        z[i] = pow_mod(17, bit_rev(i, 7) as u64, Q as u64) as u32;
        i += 1;
    }
    z
};

/// `GAMMAS[i] = 17^{2·BitRev7(i)+1} mod q`.
pub(crate) const GAMMAS: [u32; 128] = {
    let mut g = [0u32; 128];
    let mut i = 0;
    while i < 128 {
        g[i] = pow_mod(17, 2 * bit_rev(i, 7) as u64 + 1, Q as u64) as u32;
        i += 1;
    }
    g
};

/// Conditionally subtracts q: maps `[0, 2q)` to `[0, q)` without branching.
#[inline(always)]
pub(crate) fn csub(a: u32) -> u32 {
    let t = a.wrapping_sub(Q);
    t.wrapping_add(Q & 0u32.wrapping_sub(t >> 31))
}

/// Barrett reduction of any `u32` to `[0, q)`.
#[inline(always)]
pub(crate) fn reduce(x: u32) -> u32 {
    let qe = ((u64::from(x) * BARRETT_M) >> 32) as u32;
    csub(x - qe * Q)
}

/// Exact `floor(n / q)` for `n < 2^31`, in constant time.
#[inline(always)]
fn div_q(n: u32) -> u32 {
    let qe = ((u64::from(n) * BARRETT_M) >> 32) as u32;
    let r = n - qe * Q; // r ∈ [0, 2q)
    qe + ((Q - 1).wrapping_sub(r) >> 31)
}

#[inline(always)]
pub(crate) fn add(a: u32, b: u32) -> u32 {
    csub(a + b)
}

#[inline(always)]
pub(crate) fn sub(a: u32, b: u32) -> u32 {
    csub(a + Q - b)
}

#[inline(always)]
pub(crate) fn mul(a: u32, b: u32) -> u32 {
    reduce(a * b)
}

/// Algorithm 9: NTT.
pub(crate) fn ntt(f: &mut Poly) {
    let mut i = 1;
    let mut len = 128;
    while len >= 2 {
        let mut start = 0;
        while start < N {
            let zeta = ZETAS[i];
            i += 1;
            for j in start..start + len {
                let t = mul(zeta, f[j + len]);
                f[j + len] = sub(f[j], t);
                f[j] = add(f[j], t);
            }
            start += 2 * len;
        }
        len /= 2;
    }
}

/// Algorithm 10: NTT^{-1}.
pub(crate) fn inv_ntt(f: &mut Poly) {
    let mut i = 127;
    let mut len = 2;
    while len <= 128 {
        let mut start = 0;
        while start < N {
            let zeta = ZETAS[i];
            i -= 1;
            for j in start..start + len {
                let t = f[j];
                f[j] = add(t, f[j + len]);
                f[j + len] = mul(zeta, sub(f[j + len], t));
            }
            start += 2 * len;
        }
        len *= 2;
    }
    for c in f.iter_mut() {
        *c = mul(*c, INV_128);
    }
}

/// Algorithms 11 and 12: `acc += f ×_T g` in the NTT domain.
pub(crate) fn mul_acc(acc: &mut Poly, f: &Poly, g: &Poly) {
    for i in 0..128 {
        let (a0, a1, b0, b1) = (f[2 * i], f[2 * i + 1], g[2 * i], g[2 * i + 1]);
        let c0 = add(mul(a0, b0), mul(mul(a1, b1), GAMMAS[i]));
        let c1 = add(mul(a0, b1), mul(a1, b0));
        acc[2 * i] = add(acc[2 * i], c0);
        acc[2 * i + 1] = add(acc[2 * i + 1], c1);
    }
}

pub(crate) fn poly_add(a: &mut Poly, b: &Poly) {
    for (x, y) in a.iter_mut().zip(b) {
        *x = add(*x, *y);
    }
}

/// Algorithm 7: SampleNTT. `seed` is ρ ‖ j ‖ i.
pub(crate) fn sample_ntt(seed: &[u8; 34]) -> Poly {
    let mut xof = Sponge::shake128();
    xof.absorb(seed);
    let mut a = [0u32; N];
    let mut j = 0;
    let mut buf = [0u8; 168];
    while j < N {
        xof.squeeze(&mut buf);
        for c in buf.chunks_exact(3) {
            let d1 = u32::from(c[0]) | ((u32::from(c[1]) & 0x0F) << 8);
            let d2 = (u32::from(c[1]) >> 4) | (u32::from(c[2]) << 4);
            if d1 < Q {
                a[j] = d1;
                j += 1;
            }
            if d2 < Q && j < N {
                a[j] = d2;
                j += 1;
            }
            if j == N {
                break;
            }
        }
    }
    a
}

/// Algorithm 8: SamplePolyCBD_η. `bytes` has length 64·η.
pub(crate) fn sample_poly_cbd(eta: usize, bytes: &[u8]) -> Poly {
    debug_assert_eq!(bytes.len(), 64 * eta);
    let bit = |k: usize| u32::from((bytes[k / 8] >> (k % 8)) & 1);
    let mut f = [0u32; N];
    for (i, c) in f.iter_mut().enumerate() {
        let mut x = 0;
        let mut y = 0;
        for j in 0..eta {
            x += bit(2 * i * eta + j);
            y += bit(2 * i * eta + eta + j);
        }
        *c = csub(x + Q - y);
    }
    f
}

/// Compress_d (§4.2.1): `round(2^d / q · x) mod 2^d`, constant time.
#[inline(always)]
pub(crate) fn compress(d: u32, x: u32) -> u32 {
    div_q((x << d) + Q / 2) & ((1 << d) - 1)
}

/// Decompress_d (§4.2.1): `round(q / 2^d · y)`.
#[inline(always)]
pub(crate) fn decompress(d: u32, y: u32) -> u32 {
    (y * Q + (1 << (d - 1))) >> d
}

/// ByteEncode_d (Algorithm 5) after optional compression.
pub(crate) fn byte_encode(d: u32, f: &Poly, out: &mut [u8]) {
    pack(f, d, out);
}

/// ByteDecode_d (Algorithm 6); for d = 12 the result is reduced mod q.
pub(crate) fn byte_decode(d: u32, input: &[u8]) -> Poly {
    let mut f = [0u32; N];
    unpack(input, d, &mut f);
    if d == 12 {
        for c in f.iter_mut() {
            *c = csub(*c);
        }
    }
    f
}

pub(crate) fn compress_poly(d: u32, f: &mut Poly) {
    for c in f.iter_mut() {
        *c = compress(d, *c);
    }
}

pub(crate) fn decompress_poly(d: u32, f: &mut Poly) {
    for c in f.iter_mut() {
        *c = decompress(d, *c);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn zeta_tables_match_fips203_appendix_a() {
        assert_eq!(&ZETAS[..8], &[1, 1729, 2580, 3289, 2642, 630, 1897, 848]);
        // 17^127 = -17^{-1} and 17^255 = 17^{-1} = 1175 (mod q).
        assert_eq!(ZETAS[127], Q - 1175);
        assert_eq!(&GAMMAS[..4], &[17, 3312, 2761, 568]);
        assert_eq!(GAMMAS[127], 1175);
    }

    #[test]
    fn reduce_is_exact_for_all_products() {
        for a in (0..Q).step_by(7) {
            for b in [0, 1, 2, 1664, 1665, 3327, 3328, a] {
                assert_eq!(mul(a, b), a * b % Q);
            }
        }
        for x in [0u32, 1, Q - 1, Q, 2 * Q, (Q - 1) * (Q - 1), u32::MAX, u32::MAX - 1] {
            assert_eq!(reduce(x), x % Q, "x = {x}");
        }
    }

    #[test]
    fn div_q_is_exact() {
        for n in 0..(3328u32 << 11) + Q {
            assert_eq!(div_q(n), n / Q);
        }
    }

    #[test]
    fn compress_matches_rational_definition() {
        for d in [1u32, 4, 5, 10, 11] {
            for x in 0..Q {
                // round(2^d * x / q) using exact rational arithmetic, round half up.
                let num = (u64::from(x) << d) * 2 + u64::from(Q);
                let expect = ((num / (2 * u64::from(Q))) % (1 << d)) as u32;
                assert_eq!(compress(d, x), expect, "d={d} x={x}");
            }
            for y in 0..(1u32 << d) {
                let num = u64::from(y) * u64::from(Q) * 2 + (1u64 << d);
                let expect = (num / (2u64 << d)) as u32;
                assert_eq!(decompress(d, y), expect);
            }
        }
    }

    #[test]
    fn compression_error_bound() {
        // FIPS 203 §4.2.1: |Decompress(Compress(x)) − x mod± q| ≤ round(q / 2^{d+1}).
        for d in [1u32, 4, 5, 10, 11] {
            let bound = (Q + (1 << d)) >> (d + 1);
            for x in 0..Q {
                let y = decompress(d, compress(d, x));
                let diff = (y + Q - x) % Q;
                let diff = diff.min(Q - diff);
                assert!(diff <= bound, "d={d} x={x} y={y}");
            }
        }
    }

    fn schoolbook(a: &Poly, b: &Poly) -> Poly {
        let mut c = [0u32; N];
        for i in 0..N {
            for j in 0..N {
                let p = a[i] * b[j] % Q;
                let k = i + j;
                if k < N {
                    c[k] = (c[k] + p) % Q;
                } else {
                    c[k - N] = (c[k - N] + Q - p) % Q;
                }
            }
        }
        c
    }

    fn pseudo_random_poly(seed: u8) -> Poly {
        let mut s = [seed; 34];
        s[33] = seed.wrapping_mul(7);
        sample_ntt(&s)
    }

    #[test]
    fn ntt_roundtrip_and_multiplication() {
        for seed in 0..8u8 {
            let a = pseudo_random_poly(seed);
            let b = pseudo_random_poly(seed.wrapping_add(100));
            let mut t = a;
            ntt(&mut t);
            let mut back = t;
            inv_ntt(&mut back);
            assert_eq!(back, a);

            let mut bh = b;
            ntt(&mut bh);
            let mut prod = [0u32; N];
            mul_acc(&mut prod, &t, &bh);
            inv_ntt(&mut prod);
            assert_eq!(prod, schoolbook(&a, &b));
        }
    }

    #[test]
    fn cbd_coefficients_are_bounded() {
        for eta in [2usize, 3] {
            let bytes: [u8; 192] = core::array::from_fn(|i| (i * 37 + 11) as u8);
            let f = sample_poly_cbd(eta, &bytes[..64 * eta]);
            for &c in &f {
                assert!(c <= eta as u32 || c >= Q - eta as u32);
            }
            let all_ones = sample_poly_cbd(eta, &[0xFF; 192][..64 * eta]);
            assert!(all_ones.iter().all(|&c| c == 0));
        }
    }

    #[test]
    fn byte_decode_12_reduces_mod_q() {
        let mut bytes = [0xFFu8; 384];
        bytes[0] = 0x01;
        bytes[1] = 0xF0;
        let f = byte_decode(12, &bytes);
        assert_eq!(f[0], 1);
        assert_eq!(f[1], 4095 - Q);
        assert!(f.iter().all(|&c| c < Q));
    }
}
