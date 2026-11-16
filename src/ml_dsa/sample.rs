//! Pseudorandom sampling (FIPS 204 §7.3) and encoding helpers (§7.1–7.2).

use super::arith::{csub, Poly, N, Q};
use crate::common::{pack, unpack, zeroize_bytes};
use crate::hash::Sponge;

/// Algorithm 29: SampleInBall. Returns a polynomial with τ coefficients in {−1, 1}.
pub(crate) fn sample_in_ball(tau: usize, rho: &[u8]) -> Poly {
    let mut h = Sponge::shake256();
    h.absorb(rho);
    let mut s = [0u8; 8];
    h.squeeze(&mut s);
    let signs = u64::from_le_bytes(s);
    let mut c = [0u32; N];
    for i in (N - tau)..N {
        let mut j = [0u8; 1];
        loop {
            h.squeeze(&mut j);
            if usize::from(j[0]) <= i {
                break;
            }
        }
        let j = usize::from(j[0]);
        c[i] = c[j];
        c[j] = if (signs >> (i + tau - N)) & 1 == 1 { Q - 1 } else { 1 };
    }
    c
}

/// Algorithm 30: RejNTTPoly. `seed` is ρ ‖ s ‖ r (34 bytes).
pub(crate) fn rej_ntt_poly(seed: &[u8; 34]) -> Poly {
    let mut g = Sponge::shake128();
    g.absorb(seed);
    let mut a = [0u32; N];
    let mut j = 0;
    let mut buf = [0u8; 168];
    while j < N {
        g.squeeze(&mut buf);
        for b in buf.chunks_exact(3) {
            // Algorithm 14: CoeffFromThreeBytes
            let z = u32::from(b[0]) | (u32::from(b[1]) << 8) | ((u32::from(b[2]) & 0x7F) << 16);
            if z < Q {
                a[j] = z;
                j += 1;
                if j == N {
                    break;
                }
            }
        }
    }
    a
}

/// Algorithm 15: CoeffFromHalfByte. Returns η − (b mod 5) / η − b as an element of Z_q.
#[inline(always)]
fn coeff_from_half_byte(eta: u32, b: u32) -> Option<u32> {
    if eta == 2 && b < 15 {
        Some(csub(2 + Q - b % 5))
    } else if eta == 4 && b < 9 {
        Some(csub(4 + Q - b))
    } else {
        None
    }
}

/// Algorithm 31: RejBoundedPoly. `seed` is ρ′ ‖ IntegerToBytes(r, 2) (66 bytes).
pub(crate) fn rej_bounded_poly(eta: u32, seed: &[u8; 66]) -> Poly {
    let mut h = Sponge::shake256();
    h.absorb(seed);
    let mut a = [0u32; N];
    let mut j = 0;
    let mut buf = [0u8; 136];
    'outer: loop {
        h.squeeze(&mut buf);
        for &z in buf.iter() {
            if let Some(z0) = coeff_from_half_byte(eta, u32::from(z & 0x0F)) {
                a[j] = z0;
                j += 1;
            }
            if j < N {
                if let Some(z1) = coeff_from_half_byte(eta, u32::from(z >> 4)) {
                    a[j] = z1;
                    j += 1;
                }
            }
            if j == N {
                break 'outer;
            }
        }
    }
    zeroize_bytes(&mut buf);
    a
}

/// Algorithm 34 (inner part): one polynomial of ExpandMask.
/// `seed` is ρ″ (64 bytes); `nonce` is κ + r.
pub(crate) fn expand_mask_poly(gamma1: u32, z_bits: u32, rho2: &[u8; 64], nonce: u16) -> Poly {
    let mut h = Sponge::shake256();
    h.absorb(rho2);
    h.absorb(&nonce.to_le_bytes());
    let mut buf = [0u8; 32 * 20];
    let v = &mut buf[..32 * z_bits as usize];
    h.squeeze(v);
    let y = bit_unpack(v, z_bits, gamma1);
    zeroize_bytes(&mut buf);
    y
}

/// Algorithm 17: BitPack(w, a, b) — packs `b − w_i` using `bits = bitlen(a + b)` bits.
pub(crate) fn bit_pack(w: &Poly, bits: u32, b: u32, out: &mut [u8]) {
    let mut t = [0u32; N];
    for (ti, &wi) in t.iter_mut().zip(w.iter()) {
        *ti = csub(b + Q - wi);
    }
    pack(&t, bits, out);
}

/// Algorithm 19: BitUnpack(v, a, b) — returns `b − z_i mod q`.
pub(crate) fn bit_unpack(v: &[u8], bits: u32, b: u32) -> Poly {
    let mut w = [0u32; N];
    unpack(v, bits, &mut w);
    for c in w.iter_mut() {
        *c = csub(b + Q - *c);
    }
    w
}

/// Algorithm 20: HintBitPack. `out` has ω + k bytes.
pub(crate) fn hint_bit_pack<const K: usize>(omega: usize, h: &[[u32; N]; K], out: &mut [u8]) {
    debug_assert_eq!(out.len(), omega + K);
    out.fill(0);
    let mut index = 0;
    for (i, hi) in h.iter().enumerate() {
        for (j, &bit) in hi.iter().enumerate() {
            if bit != 0 {
                out[index] = j as u8;
                index += 1;
            }
        }
        out[omega + i] = index as u8;
    }
}

/// Algorithm 21: HintBitUnpack. Returns `None` (⊥) on malformed input.
pub(crate) fn hint_bit_unpack<const K: usize>(omega: usize, y: &[u8]) -> Option<[[u32; N]; K]> {
    debug_assert_eq!(y.len(), omega + K);
    let mut h = [[0u32; N]; K];
    let mut index = 0usize;
    for (i, hi) in h.iter_mut().enumerate() {
        let limit = usize::from(y[omega + i]);
        if limit < index || limit > omega {
            return None;
        }
        let first = index;
        while index < limit {
            if index > first && y[index - 1] >= y[index] {
                return None;
            }
            hi[usize::from(y[index])] = 1;
            index += 1;
        }
    }
    if y[index..omega].iter().any(|&b| b != 0) {
        return None;
    }
    Some(h)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ml_dsa::arith::abs_centered;

    #[test]
    fn sample_in_ball_has_weight_tau() {
        for tau in [39usize, 49, 60] {
            for seed in 0..50u8 {
                let c = sample_in_ball(tau, &[seed; 48]);
                assert_eq!(c.iter().filter(|&&x| x != 0).count(), tau);
                assert!(c.iter().all(|&x| x == 0 || x == 1 || x == Q - 1));
            }
        }
    }

    #[test]
    fn rej_bounded_range() {
        for eta in [2u32, 4] {
            let a = rej_bounded_poly(eta, &[7u8; 66]);
            assert!(a.iter().all(|&c| abs_centered(c) <= eta));
        }
    }

    #[test]
    fn expand_mask_range() {
        for (g1, bits) in [(1u32 << 17, 18u32), (1 << 19, 20)] {
            let y = expand_mask_poly(g1, bits, &[3u8; 64], 5);
            for &c in &y {
                let v = crate::ml_dsa::arith::center(c);
                assert!(v > -(g1 as i32) && v <= g1 as i32);
            }
        }
    }

    #[test]
    fn bit_pack_roundtrip() {
        let mut w = [0u32; N];
        for (i, c) in w.iter_mut().enumerate() {
            let v = (i as i32 % 9) - 4; // [−4, 4]
            *c = ((v + Q as i32) as u32) % Q;
        }
        let mut out = [0u8; 128];
        bit_pack(&w, 4, 4, &mut out);
        assert_eq!(bit_unpack(&out, 4, 4), w);
    }

    #[test]
    fn hint_pack_roundtrip_and_rejections() {
        const K: usize = 4;
        let omega = 80;
        let mut h = [[0u32; N]; K];
        h[0][3] = 1;
        h[0][200] = 1;
        h[2][0] = 1;
        h[3][255] = 1;
        let mut out = [0u8; 84];
        hint_bit_pack::<K>(omega, &h, &mut out);
        assert_eq!(hint_bit_unpack::<K>(omega, &out), Some(h));

        // Non-increasing indices within a polynomial.
        let mut bad = out;
        bad.swap(0, 1);
        assert_eq!(hint_bit_unpack::<K>(omega, &bad), None);
        // Count exceeding ω.
        let mut bad = out;
        bad[omega + 3] = (omega + 1) as u8;
        assert_eq!(hint_bit_unpack::<K>(omega, &bad), None);
        // Decreasing cumulative counts.
        let mut bad = out;
        bad[omega + 1] = 1;
        assert_eq!(hint_bit_unpack::<K>(omega, &bad), None);
        // Non-zero padding after the last index.
        let mut bad = out;
        bad[omega - 1] = 1;
        assert_eq!(hint_bit_unpack::<K>(omega, &bad), None);
    }
}
