//! K-PKE, the IND-CPA public-key encryption scheme underlying ML-KEM
//! (FIPS 203 §5, Algorithms 13–15).

use super::arith::*;
use crate::common::{zeroize_bytes, zeroize_u32};
use crate::hash::Sponge;

/// Parameters that vary across ML-KEM parameter sets (FIPS 203 Table 2).
/// `k` is a const generic on every function.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Params {
    pub(crate) eta1: usize,
    pub(crate) eta2: usize,
    pub(crate) du: u32,
    pub(crate) dv: u32,
}

/// `SamplePolyCBD_η(PRF_η(s, b))`.
fn prf_cbd(eta: usize, s: &[u8], b: u8) -> Poly {
    let mut buf = [0u8; 64 * 3];
    let out = &mut buf[..64 * eta];
    let mut h = Sponge::shake256();
    h.absorb(s);
    h.absorb(&[b]);
    h.squeeze(out);
    let f = sample_poly_cbd(eta, out);
    zeroize_bytes(&mut buf);
    f
}

/// Generates Â (or Âᵀ when `transpose`), where Â[i][j] = SampleNTT(ρ ‖ j ‖ i).
fn sample_matrix<const K: usize>(rho: &[u8], transpose: bool) -> [[Poly; K]; K] {
    let mut a = [[[0u32; N]; K]; K];
    let mut seed = [0u8; 34];
    seed[..32].copy_from_slice(rho);
    for i in 0..K {
        for j in 0..K {
            seed[32] = j as u8;
            seed[33] = i as u8;
            let p = sample_ntt(&seed);
            if transpose {
                a[j][i] = p;
            } else {
                a[i][j] = p;
            }
        }
    }
    a
}

/// Algorithm 13: K-PKE.KeyGen(d). Writes `ek` (384k+32 bytes) and `dk` (384k bytes).
pub(crate) fn keygen<const K: usize>(p: &Params, d: &[u8; 32], ek: &mut [u8], dk: &mut [u8]) {
    debug_assert_eq!(ek.len(), 384 * K + 32);
    debug_assert_eq!(dk.len(), 384 * K);

    let mut g = [0u8; 64];
    let mut h = Sponge::sha3_512();
    h.absorb(d);
    h.absorb(&[K as u8]);
    h.squeeze(&mut g);
    let (rho, sigma) = g.split_at(32);

    let a_hat = sample_matrix::<K>(rho, false);
    let mut n = 0u8;
    let mut s = [[0u32; N]; K];
    let mut e = [[0u32; N]; K];
    for si in s.iter_mut() {
        *si = prf_cbd(p.eta1, sigma, n);
        n += 1;
    }
    for ei in e.iter_mut() {
        *ei = prf_cbd(p.eta1, sigma, n);
        n += 1;
    }
    for (si, ei) in s.iter_mut().zip(e.iter_mut()) {
        ntt(si);
        ntt(ei);
    }

    for i in 0..K {
        let mut t = e[i];
        for j in 0..K {
            mul_acc(&mut t, &a_hat[i][j], &s[j]);
        }
        byte_encode(12, &t, &mut ek[384 * i..384 * (i + 1)]);
        byte_encode(12, &s[i], &mut dk[384 * i..384 * (i + 1)]);
    }
    ek[384 * K..].copy_from_slice(rho);

    zeroize_bytes(&mut g);
    zeroize_u32(s.as_flattened_mut());
    zeroize_u32(e.as_flattened_mut());
}

/// Algorithm 14: K-PKE.Encrypt(ek, m, r). Writes the ciphertext into `c`.
pub(crate) fn encrypt<const K: usize>(
    p: &Params,
    ek: &[u8],
    m: &[u8; 32],
    r: &[u8; 32],
    c: &mut [u8],
) {
    let du_len = 32 * p.du as usize;
    let dv_len = 32 * p.dv as usize;
    debug_assert_eq!(ek.len(), 384 * K + 32);
    debug_assert_eq!(c.len(), K * du_len + dv_len);

    let mut t_hat = [[0u32; N]; K];
    for (i, t) in t_hat.iter_mut().enumerate() {
        *t = byte_decode(12, &ek[384 * i..384 * (i + 1)]);
    }
    let a_t = sample_matrix::<K>(&ek[384 * K..], true);

    let mut n = 0u8;
    let mut y = [[0u32; N]; K];
    let mut e1 = [[0u32; N]; K];
    for yi in y.iter_mut() {
        *yi = prf_cbd(p.eta1, r, n);
        n += 1;
    }
    for e in e1.iter_mut() {
        *e = prf_cbd(p.eta2, r, n);
        n += 1;
    }
    let mut e2 = prf_cbd(p.eta2, r, n);
    for yi in y.iter_mut() {
        ntt(yi);
    }

    // u = NTT^{-1}(Âᵀ ∘ ŷ) + e1
    for i in 0..K {
        let mut u = [0u32; N];
        for j in 0..K {
            mul_acc(&mut u, &a_t[i][j], &y[j]);
        }
        inv_ntt(&mut u);
        poly_add(&mut u, &e1[i]);
        compress_poly(p.du, &mut u);
        byte_encode(p.du, &u, &mut c[i * du_len..(i + 1) * du_len]);
    }

    // v = NTT^{-1}(t̂ᵀ ∘ ŷ) + e2 + μ
    let mut v = [0u32; N];
    for j in 0..K {
        mul_acc(&mut v, &t_hat[j], &y[j]);
    }
    inv_ntt(&mut v);
    poly_add(&mut v, &e2);
    let mut mu = byte_decode(1, m);
    decompress_poly(1, &mut mu);
    poly_add(&mut v, &mu);
    compress_poly(p.dv, &mut v);
    byte_encode(p.dv, &v, &mut c[K * du_len..]);

    zeroize_u32(y.as_flattened_mut());
    zeroize_u32(e1.as_flattened_mut());
    zeroize_u32(&mut e2);
    zeroize_u32(&mut mu);
    zeroize_u32(&mut v);
}

/// Algorithm 15: K-PKE.Decrypt(dk, c).
pub(crate) fn decrypt<const K: usize>(p: &Params, dk_pke: &[u8], c: &[u8]) -> [u8; 32] {
    let du_len = 32 * p.du as usize;
    debug_assert_eq!(dk_pke.len(), 384 * K);
    debug_assert_eq!(c.len(), K * du_len + 32 * p.dv as usize);

    // w = v' − NTT^{-1}(ŝᵀ ∘ NTT(u'))
    let mut acc = [0u32; N];
    for i in 0..K {
        let mut u = byte_decode(p.du, &c[i * du_len..(i + 1) * du_len]);
        decompress_poly(p.du, &mut u);
        ntt(&mut u);
        let mut s = byte_decode(12, &dk_pke[384 * i..384 * (i + 1)]);
        mul_acc(&mut acc, &s, &u);
        zeroize_u32(&mut s);
    }
    inv_ntt(&mut acc);
    let mut w = byte_decode(p.dv, &c[K * du_len..]);
    decompress_poly(p.dv, &mut w);
    for (wi, ai) in w.iter_mut().zip(acc.iter()) {
        *wi = sub(*wi, *ai);
    }
    compress_poly(1, &mut w);
    let mut m = [0u8; 32];
    byte_encode(1, &w, &mut m);
    zeroize_u32(&mut acc);
    zeroize_u32(&mut w);
    m
}
