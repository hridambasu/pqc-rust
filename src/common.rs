//! Helpers shared by ML-KEM and ML-DSA: bit packing, constant-time utilities,
//! zeroization, and compile-time table generation.

use core::hint::black_box;
use core::sync::atomic::{compiler_fence, Ordering};

/// `base^exp mod m`, evaluated at compile time for the NTT tables.
pub(crate) const fn pow_mod(base: u64, mut exp: u64, m: u64) -> u64 {
    let mut result = 1u64;
    let mut b = base % m;
    while exp > 0 {
        if exp & 1 == 1 {
            result = result * b % m;
        }
        b = b * b % m;
        exp >>= 1;
    }
    result
}

/// Reverses the low `bits` bits of `x` (BitRev7 / BitRev8 in the standards).
pub(crate) const fn bit_rev(x: usize, bits: u32) -> usize {
    let mut r = 0;
    let mut i = 0;
    while i < bits {
        r |= ((x >> i) & 1) << (bits - 1 - i);
        i += 1;
    }
    r
}

/// Packs 256 coefficients of `bits` bits each, little-endian bit order.
///
/// Implements FIPS 203 `ByteEncode_d` and FIPS 204 `SimpleBitPack`.
/// Only the low `bits` bits of each coefficient are used.
pub(crate) fn pack(coeffs: &[u32; 256], bits: u32, out: &mut [u8]) {
    debug_assert!((1..=24).contains(&bits));
    debug_assert_eq!(out.len(), 32 * bits as usize);
    let mask = (1u64 << bits) - 1;
    let mut acc = 0u64;
    let mut n = 0u32;
    let mut idx = 0;
    for &c in coeffs {
        acc |= (u64::from(c) & mask) << n;
        n += bits;
        while n >= 8 {
            out[idx] = acc as u8;
            idx += 1;
            acc >>= 8;
            n -= 8;
        }
    }
}

/// Inverse of [`pack`] (FIPS 203 `ByteDecode_d` without the mod-q reduction,
/// FIPS 204 `SimpleBitUnpack`).
pub(crate) fn unpack(input: &[u8], bits: u32, out: &mut [u32; 256]) {
    debug_assert!((1..=24).contains(&bits));
    debug_assert_eq!(input.len(), 32 * bits as usize);
    let mask = (1u64 << bits) - 1;
    let mut acc = 0u64;
    let mut n = 0u32;
    let mut idx = 0;
    for c in out.iter_mut() {
        while n < bits {
            acc |= u64::from(input[idx]) << n;
            idx += 1;
            n += 8;
        }
        *c = (acc & mask) as u32;
        acc >>= bits;
        n -= bits;
    }
}

/// Constant-time equality: returns 1 if `a == b`, else 0.
///
/// Lengths are treated as public.
pub(crate) fn ct_eq(a: &[u8], b: &[u8]) -> u8 {
    if a.len() != b.len() {
        return 0;
    }
    let mut acc = 0u8;
    for (x, y) in a.iter().zip(b) {
        acc |= x ^ y;
    }
    let acc = u32::from(black_box(acc));
    (acc.wrapping_sub(1) >> 31) as u8
}

/// Constant-time conditional copy: `dst = src` if `choice == 1`, unchanged if 0.
pub(crate) fn ct_select(dst: &mut [u8], src: &[u8], choice: u8) {
    debug_assert_eq!(dst.len(), src.len());
    debug_assert!(choice <= 1);
    let mask = black_box(0u8.wrapping_sub(choice));
    for (d, s) in dst.iter_mut().zip(src) {
        *d ^= mask & (*d ^ *s);
    }
}

/// Overwrites a byte buffer with zeros in a way the optimizer cannot elide.
pub(crate) fn zeroize_bytes(buf: &mut [u8]) {
    for b in buf.iter_mut() {
        *b = 0;
    }
    black_box(&*buf);
    compiler_fence(Ordering::SeqCst);
}

/// Overwrites a coefficient buffer with zeros in a way the optimizer cannot elide.
pub(crate) fn zeroize_u32(buf: &mut [u32]) {
    for b in buf.iter_mut() {
        *b = 0;
    }
    black_box(&*buf);
    compiler_fence(Ordering::SeqCst);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pack_unpack_roundtrip_all_widths() {
        let mut state = 0x1234_5678_9abc_def0u64;
        for bits in 1..=24u32 {
            let mut coeffs = [0u32; 256];
            for c in coeffs.iter_mut() {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                *c = (state as u32) & ((1u32 << bits) - 1);
            }
            let mut buf = [0u8; 32 * 24];
            pack(&coeffs, bits, &mut buf[..32 * bits as usize]);
            let mut back = [0u32; 256];
            unpack(&buf[..32 * bits as usize], bits, &mut back);
            assert_eq!(coeffs, back, "bits = {bits}");
        }
    }

    #[test]
    fn pack_is_little_endian_bit_order() {
        let mut coeffs = [0u32; 256];
        coeffs[0] = 1;
        coeffs[1] = 1;
        let mut out = [0u8; 32];
        pack(&coeffs, 1, &mut out);
        assert_eq!(out[0], 0b11);
        let mut out = [0u8; 32 * 12];
        coeffs[0] = 0xABC;
        coeffs[1] = 0x123;
        pack(&coeffs, 12, &mut out);
        assert_eq!(&out[..3], &[0xBC, 0x3A, 0x12]);
    }

    #[test]
    fn ct_helpers() {
        assert_eq!(ct_eq(b"abc", b"abc"), 1);
        assert_eq!(ct_eq(b"abc", b"abd"), 0);
        assert_eq!(ct_eq(b"abc", b"ab"), 0);
        assert_eq!(ct_eq(&[], &[]), 1);
        for i in 0..=255u8 {
            assert_eq!(ct_eq(&[i], &[0]), u8::from(i == 0));
        }
        let mut d = [1u8, 2, 3];
        ct_select(&mut d, &[4, 5, 6], 0);
        assert_eq!(d, [1, 2, 3]);
        ct_select(&mut d, &[4, 5, 6], 1);
        assert_eq!(d, [4, 5, 6]);
        zeroize_bytes(&mut d);
        assert_eq!(d, [0, 0, 0]);
    }

    #[test]
    fn bit_rev_and_pow() {
        assert_eq!(bit_rev(1, 7), 64);
        assert_eq!(bit_rev(0b1010000, 7), 0b0000101);
        assert_eq!(bit_rev(1, 8), 128);
        assert_eq!(pow_mod(17, 128, 3329), 3328); // 17 is a primitive 256-th root of unity
        assert_eq!(pow_mod(1753, 512, 8_380_417), 1);
        assert_eq!(pow_mod(1753, 256, 8_380_417), 8_380_416); // 1753 is a primitive 512-th root
    }
}
