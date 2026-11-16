//! Arithmetic over R_q = Z_q[X]/(X^256 + 1) with q = 8380417 (FIPS 204 §7.5),
//! and the rounding functions of §7.4.
//!
//! Coefficients are stored in canonical form `[0, q)` as `u32`.

use crate::common::{bit_rev, pow_mod};

pub(crate) const N: usize = 256;
pub(crate) const Q: u32 = 8_380_417;
/// Number of dropped bits from t (FIPS 204 Table 1).
pub(crate) const D: u32 = 13;
pub(crate) type Poly = [u32; N];

/// `floor(2^64 / q)` for Barrett reduction of 64-bit products.
const BARRETT_M: u128 = (1u128 << 64) / Q as u128;
/// `256^{-1} mod q`.
const INV_256: u32 = 8_347_681;

/// `ZETAS[k] = 1753^{BitRev8(k)} mod q` (FIPS 204 Appendix B).
pub(crate) const ZETAS: [u32; 256] = {
    let mut z = [0u32; 256];
    let mut i = 0;
    while i < 256 {
        z[i] = pow_mod(1753, bit_rev(i, 8) as u64, Q as u64) as u32;
        i += 1;
    }
    z
};

#[inline(always)]
pub(crate) fn csub(a: u32) -> u32 {
    let t = a.wrapping_sub(Q);
    t.wrapping_add(Q & 0u32.wrapping_sub(t >> 31))
}

/// Barrett reduction of `x < 2^63` to `[0, q)`, constant time.
#[inline(always)]
pub(crate) fn reduce64(x: u64) -> u32 {
    let qe = ((u128::from(x) * BARRETT_M) >> 64) as u64;
    csub((x - qe * u64::from(Q)) as u32)
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
pub(crate) fn neg(a: u32) -> u32 {
    csub(Q - a)
}

#[inline(always)]
pub(crate) fn mul(a: u32, b: u32) -> u32 {
    reduce64(u64::from(a) * u64::from(b))
}

/// Algorithm 41: NTT.
pub(crate) fn ntt(w: &mut Poly) {
    let mut m = 0;
    let mut len = 128;
    while len > 0 {
        let mut start = 0;
        while start < N {
            m += 1;
            let z = ZETAS[m];
            for j in start..start + len {
                let t = mul(z, w[j + len]);
                w[j + len] = sub(w[j], t);
                w[j] = add(w[j], t);
            }
            start += 2 * len;
        }
        len /= 2;
    }
}

/// Algorithm 42: NTT^{-1}.
pub(crate) fn inv_ntt(w: &mut Poly) {
    let mut m = 256;
    let mut len = 1;
    while len < N {
        let mut start = 0;
        while start < N {
            m -= 1;
            let z = Q - ZETAS[m];
            for j in start..start + len {
                let t = w[j];
                w[j] = add(t, w[j + len]);
                w[j + len] = mul(z, sub(t, w[j + len]));
            }
            start += 2 * len;
        }
        len *= 2;
    }
    for c in w.iter_mut() {
        *c = mul(*c, INV_256);
    }
}

/// Algorithm 45: MultiplyNTT, accumulating: `acc += a ∘ b`.
pub(crate) fn mul_acc(acc: &mut Poly, a: &Poly, b: &Poly) {
    for i in 0..N {
        acc[i] = add(acc[i], mul(a[i], b[i]));
    }
}

/// Algorithm 45: MultiplyNTT.
pub(crate) fn pointwise(a: &Poly, b: &Poly) -> Poly {
    let mut c = [0u32; N];
    for i in 0..N {
        c[i] = mul(a[i], b[i]);
    }
    c
}

/// Returns the centered representative of `a` (in `[-(q-1)/2, (q-1)/2]`), branch-free.
#[inline(always)]
pub(crate) fn center(a: u32) -> i32 {
    let mask = (((Q - 1) / 2).wrapping_sub(a) as i32) >> 31; // −1 if a > (q−1)/2
    (a as i32) - (Q as i32 & mask)
}

/// `|a mod± q|`, branch-free.
#[inline(always)]
pub(crate) fn abs_centered(a: u32) -> u32 {
    let c = center(a);
    let m = c >> 31;
    ((c ^ m) - m) as u32
}

/// Infinity norm of a polynomial (FIPS 204 §2.3).
pub(crate) fn inf_norm(p: &Poly) -> u32 {
    p.iter().fold(0, |acc, &c| acc.max(abs_centered(c)))
}

/// Algorithm 35: Power2Round. Returns `(r1, r0 mod q)`.
#[inline(always)]
pub(crate) fn power2round(r: u32) -> (u32, u32) {
    let r1 = (r + (1 << (D - 1)) - 1) >> D;
    let r0 = r as i32 - (r1 << D) as i32; // r0 ∈ (−2^{d−1}, 2^{d−1}]
    (r1, csub((r0 + Q as i32) as u32))
}

/// Algorithm 36: Decompose. Returns `(r1, r0)` with `r0` centered. Constant time.
#[inline(always)]
pub(crate) fn decompose(gamma2: u32, r: u32) -> (u32, i32) {
    let a = r as i32;
    let mut a1 = (a + 127) >> 7;
    if gamma2 == (Q - 1) / 32 {
        a1 = (a1 * 1025 + (1 << 21)) >> 22;
        a1 &= 15;
    } else {
        debug_assert_eq!(gamma2, (Q - 1) / 88);
        a1 = (a1 * 11275 + (1 << 23)) >> 24;
        a1 ^= ((43 - a1) >> 31) & a1;
    }
    let mut a0 = a - a1 * 2 * gamma2 as i32;
    a0 -= (((Q as i32 - 1) / 2 - a0) >> 31) & Q as i32;
    (a1 as u32, a0)
}

/// Algorithm 37: HighBits.
#[inline(always)]
pub(crate) fn high_bits(gamma2: u32, r: u32) -> u32 {
    decompose(gamma2, r).0
}

/// Algorithm 38: LowBits (centered).
#[inline(always)]
pub(crate) fn low_bits(gamma2: u32, r: u32) -> i32 {
    decompose(gamma2, r).1
}

/// Algorithm 39: MakeHint. Returns 1 iff adding `z` to `r` changes the high bits.
#[inline(always)]
pub(crate) fn make_hint(gamma2: u32, z: u32, r: u32) -> u32 {
    let a = high_bits(gamma2, r);
    let b = high_bits(gamma2, add(r, z));
    u32::from(a != b)
}

/// Algorithm 40: UseHint.
#[inline(always)]
pub(crate) fn use_hint(gamma2: u32, h: u32, r: u32) -> u32 {
    let m = (Q - 1) / (2 * gamma2);
    let (r1, r0) = decompose(gamma2, r);
    if h == 1 {
        if r0 > 0 {
            (r1 + 1) % m
        } else {
            (r1 + m - 1) % m
        }
    } else {
        r1
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const GAMMAS2: [u32; 2] = [(Q - 1) / 88, (Q - 1) / 32];

    /// Reference Decompose, literally following FIPS 204 Algorithm 36.
    fn decompose_spec(gamma2: u32, r: u32) -> (u32, i32) {
        let alpha = 2 * gamma2 as i64;
        let rp = i64::from(r).rem_euclid(i64::from(Q));
        let mut r0 = rp.rem_euclid(alpha);
        if r0 > alpha / 2 {
            r0 -= alpha;
        }
        if rp - r0 == i64::from(Q) - 1 {
            (0, (r0 - 1) as i32)
        } else {
            (((rp - r0) / alpha) as u32, r0 as i32)
        }
    }

    #[test]
    fn zetas_match_fips204_appendix_b() {
        assert_eq!(&ZETAS[1..4], &[4_808_194, 3_765_607, 3_761_513]);
        assert_eq!(ZETAS[255], 7_648_983);
    }

    #[test]
    fn decompose_matches_spec_exhaustively() {
        for &g in &GAMMAS2 {
            for r in 0..Q {
                assert_eq!(decompose(g, r), decompose_spec(g, r), "gamma2={g} r={r}");
            }
        }
    }

    #[test]
    fn power2round_matches_spec() {
        for r in (0..Q).step_by(3).chain([0, 1, 4095, 4096, 4097, Q - 1]) {
            let (r1, r0) = power2round(r);
            let r0c = center(r0);
            assert!(r0c > -(1 << 12) && r0c <= 1 << 12, "r={r}");
            assert_eq!((r1 << D) as i64 + i64::from(r0c), i64::from(r));
        }
    }

    #[test]
    fn hint_roundtrip_property() {
        // FIPS 204 §7.4: UseHint(MakeHint(z, r), r) = HighBits(r + z) for ‖z‖∞ ≤ γ2.
        let mut s = 0x9E37_79B9_7F4A_7C15u64;
        let mut next = || {
            s ^= s << 13;
            s ^= s >> 7;
            s ^= s << 17;
            s
        };
        for &g in &GAMMAS2 {
            for _ in 0..200_000 {
                let r = (next() % u64::from(Q)) as u32;
                let zc = (next() % (2 * u64::from(g) + 1)) as i64 - i64::from(g);
                let z = (zc.rem_euclid(i64::from(Q))) as u32;
                let h = make_hint(g, z, r);
                assert_eq!(use_hint(g, h, r), high_bits(g, add(r, z)));
            }
        }
    }

    #[test]
    fn reduce_and_center() {
        for x in [0u64, 1, u64::from(Q), u64::from(Q) * u64::from(Q), (1 << 46) - 1, (1 << 62) + 12345] {
            assert_eq!(u64::from(reduce64(x)), x % u64::from(Q));
        }
        assert_eq!(center(0), 0);
        assert_eq!(center((Q - 1) / 2), ((Q - 1) / 2) as i32);
        assert_eq!(center((Q + 1) / 2), -(((Q - 1) / 2) as i32));
        assert_eq!(center(Q - 1), -1);
        assert_eq!(abs_centered(Q - 5), 5);
    }

    #[test]
    fn ntt_roundtrip_and_multiplication() {
        let mut a = [0u32; N];
        let mut b = [0u32; N];
        for i in 0..N {
            a[i] = (i as u32 * 7_654_321 + 13) % Q;
            b[i] = (i as u32 * 1_234_567 + 99) % Q;
        }
        let mut ah = a;
        ntt(&mut ah);
        let mut back = ah;
        inv_ntt(&mut back);
        assert_eq!(back, a);

        let mut bh = b;
        ntt(&mut bh);
        let mut c = pointwise(&ah, &bh);
        inv_ntt(&mut c);
        let mut expect = [0u64; N];
        for i in 0..N {
            for j in 0..N {
                let p = u64::from(a[i]) * u64::from(b[j]) % u64::from(Q);
                let k = (i + j) % N;
                if i + j < N {
                    expect[k] = (expect[k] + p) % u64::from(Q);
                } else {
                    expect[k] = (expect[k] + u64::from(Q) - p) % u64::from(Q);
                }
            }
        }
        for i in 0..N {
            assert_eq!(u64::from(c[i]), expect[i]);
        }
    }
}
