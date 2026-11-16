//! ML-DSA internal algorithms (FIPS 204 §6): KeyGen_internal, Sign_internal,
//! Verify_internal, and the key/signature encodings of §7.2.

use super::arith::*;
use super::sample::*;
use crate::common::{ct_eq, pack, unpack, zeroize_bytes, zeroize_u32};
use crate::hash::{shake256, Sponge};

/// Parameters of an ML-DSA parameter set (FIPS 204 Table 1); `k` and `ℓ` are
/// const generics on every function.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Params {
    pub(crate) eta: u32,
    pub(crate) tau: usize,
    pub(crate) beta: u32,
    pub(crate) gamma1: u32,
    pub(crate) gamma2: u32,
    pub(crate) omega: usize,
    /// Collision strength λ in bits.
    pub(crate) lambda: usize,
}

impl Params {
    pub(crate) const fn eta_bits(&self) -> u32 {
        if self.eta == 2 {
            3
        } else {
            4
        }
    }
    pub(crate) const fn z_bits(&self) -> u32 {
        if self.gamma1 == 1 << 17 {
            18
        } else {
            20
        }
    }
    pub(crate) const fn w1_bits(&self) -> u32 {
        if self.gamma2 == (Q - 1) / 88 {
            6
        } else {
            4
        }
    }
    pub(crate) const fn ctilde_len(&self) -> usize {
        self.lambda / 4
    }
    pub(crate) const fn pk_len(&self, k: usize) -> usize {
        32 + 320 * k
    }
    pub(crate) const fn sk_len(&self, k: usize, l: usize) -> usize {
        128 + 32 * self.eta_bits() as usize * (k + l) + 32 * D as usize * k
    }
    pub(crate) const fn sig_len(&self, k: usize, l: usize) -> usize {
        self.ctilde_len() + 32 * self.z_bits() as usize * l + self.omega + k
    }
}

const T1_BITS: u32 = 10; // bitlen(q − 1) − d

/// Algorithm 32: ExpandA. Â[r][s] = RejNTTPoly(ρ ‖ s ‖ r).
fn expand_a<const K: usize, const L: usize>(rho: &[u8]) -> [[Poly; L]; K] {
    let mut a = [[[0u32; N]; L]; K];
    let mut seed = [0u8; 34];
    seed[..32].copy_from_slice(rho);
    for (r, row) in a.iter_mut().enumerate() {
        for (s, entry) in row.iter_mut().enumerate() {
            seed[32] = s as u8;
            seed[33] = r as u8;
            *entry = rej_ntt_poly(&seed);
        }
    }
    a
}

/// Algorithm 33: ExpandS.
fn expand_s<const K: usize, const L: usize>(eta: u32, rho1: &[u8]) -> ([Poly; L], [Poly; K]) {
    let mut seed = [0u8; 66];
    seed[..64].copy_from_slice(rho1);
    let mut s1 = [[0u32; N]; L];
    let mut s2 = [[0u32; N]; K];
    for (r, p) in s1.iter_mut().enumerate() {
        seed[64..].copy_from_slice(&(r as u16).to_le_bytes());
        *p = rej_bounded_poly(eta, &seed);
    }
    for (r, p) in s2.iter_mut().enumerate() {
        seed[64..].copy_from_slice(&((r + L) as u16).to_le_bytes());
        *p = rej_bounded_poly(eta, &seed);
    }
    zeroize_bytes(&mut seed);
    (s1, s2)
}

/// Â ∘ v̂ in the NTT domain.
fn mat_vec<const K: usize, const L: usize>(a: &[[Poly; L]; K], v: &[Poly; L]) -> [Poly; K] {
    let mut out = [[0u32; N]; K];
    for (o, row) in out.iter_mut().zip(a.iter()) {
        for (aij, vj) in row.iter().zip(v.iter()) {
            mul_acc(o, aij, vj);
        }
    }
    out
}

/// Computes t = NTT^{-1}(Â ∘ NTT(s1)) + s2 and returns (t1, t0).
fn compute_t<const K: usize, const L: usize>(
    a: &[[Poly; L]; K],
    s1: &[Poly; L],
    s2: &[Poly; K],
) -> ([Poly; K], [Poly; K]) {
    let mut s1_hat = *s1;
    for p in s1_hat.iter_mut() {
        ntt(p);
    }
    let mut t = mat_vec(a, &s1_hat);
    let mut t1 = [[0u32; N]; K];
    let mut t0 = [[0u32; N]; K];
    for i in 0..K {
        inv_ntt(&mut t[i]);
        for j in 0..N {
            let (r1, r0) = power2round(add(t[i][j], s2[i][j]));
            t1[i][j] = r1;
            t0[i][j] = r0;
        }
    }
    zeroize_u32(s1_hat.as_flattened_mut());
    zeroize_u32(t.as_flattened_mut());
    (t1, t0)
}

/// Algorithm 22: pkEncode.
fn pk_encode<const K: usize>(rho: &[u8], t1: &[Poly; K], pk: &mut [u8]) {
    pk[..32].copy_from_slice(rho);
    for (i, t) in t1.iter().enumerate() {
        pack(t, T1_BITS, &mut pk[32 + 320 * i..32 + 320 * (i + 1)]);
    }
}

/// Algorithm 6: ML-DSA.KeyGen_internal(ξ).
pub(crate) fn keygen<const K: usize, const L: usize>(
    p: &Params,
    xi: &[u8; 32],
    pk: &mut [u8],
    sk: &mut [u8],
) {
    debug_assert_eq!(pk.len(), p.pk_len(K));
    debug_assert_eq!(sk.len(), p.sk_len(K, L));

    let mut seeds = [0u8; 128];
    shake256(&[xi, &[K as u8, L as u8]], &mut seeds);
    let (rho, rest) = seeds.split_at(32);
    let (rho1, key) = rest.split_at(64);

    let a = expand_a::<K, L>(rho);
    let (mut s1, mut s2) = expand_s::<K, L>(p.eta, rho1);
    let (t1, mut t0) = compute_t(&a, &s1, &s2);

    pk_encode(rho, &t1, pk);
    let mut tr = [0u8; 64];
    shake256(&[pk], &mut tr);

    // Algorithm 24: skEncode.
    let eb = 32 * p.eta_bits() as usize;
    sk[..32].copy_from_slice(rho);
    sk[32..64].copy_from_slice(key);
    sk[64..128].copy_from_slice(&tr);
    let mut off = 128;
    for s in s1.iter().chain(s2.iter()) {
        bit_pack(s, p.eta_bits(), p.eta, &mut sk[off..off + eb]);
        off += eb;
    }
    for t in t0.iter() {
        bit_pack(t, D, 1 << (D - 1), &mut sk[off..off + 32 * D as usize]);
        off += 32 * D as usize;
    }
    debug_assert_eq!(off, sk.len());

    zeroize_bytes(&mut seeds);
    zeroize_u32(s1.as_flattened_mut());
    zeroize_u32(s2.as_flattened_mut());
    zeroize_u32(t0.as_flattened_mut());
}

/// Decoded signing key (Algorithm 25: skDecode).
struct SecretKey<const K: usize, const L: usize> {
    rho: [u8; 32],
    key: [u8; 32],
    s1: [Poly; L],
    s2: [Poly; K],
    t0: [Poly; K],
}

impl<const K: usize, const L: usize> Drop for SecretKey<K, L> {
    fn drop(&mut self) {
        zeroize_bytes(&mut self.key);
        zeroize_u32(self.s1.as_flattened_mut());
        zeroize_u32(self.s2.as_flattened_mut());
        zeroize_u32(self.t0.as_flattened_mut());
    }
}

fn sk_decode<const K: usize, const L: usize>(p: &Params, sk: &[u8]) -> SecretKey<K, L> {
    let eb = 32 * p.eta_bits() as usize;
    let mut out = SecretKey {
        rho: [0; 32],
        key: [0; 32],
        s1: [[0; N]; L],
        s2: [[0; N]; K],
        t0: [[0; N]; K],
    };
    out.rho.copy_from_slice(&sk[..32]);
    out.key.copy_from_slice(&sk[32..64]);
    let mut off = 128;
    for s in out.s1.iter_mut().chain(out.s2.iter_mut()) {
        *s = bit_unpack(&sk[off..off + eb], p.eta_bits(), p.eta);
        off += eb;
    }
    for t in out.t0.iter_mut() {
        *t = bit_unpack(&sk[off..off + 32 * D as usize], D, 1 << (D - 1));
        off += 32 * D as usize;
    }
    out
}

/// Checks that every packed coefficient of s1 and s2 lies in [−η, η], i.e.
/// that the key is a possible output of skEncode for these fields.
pub(crate) fn sk_is_well_formed<const K: usize, const L: usize>(p: &Params, sk: &[u8]) -> bool {
    let eb = 32 * p.eta_bits() as usize;
    let mut raw = [0u32; N];
    let mut ok = true;
    for chunk in sk[128..128 + eb * (K + L)].chunks_exact(eb) {
        unpack(chunk, p.eta_bits(), &mut raw);
        ok &= raw.iter().all(|&c| c <= 2 * p.eta);
    }
    zeroize_u32(&mut raw);
    ok
}

/// Recomputes the public key from a signing key.
pub(crate) fn pk_from_sk<const K: usize, const L: usize>(p: &Params, sk: &[u8], pk: &mut [u8]) {
    let s = sk_decode::<K, L>(p, sk);
    let a = expand_a::<K, L>(&s.rho);
    let (t1, mut t0) = compute_t(&a, &s.s1, &s.s2);
    pk_encode(&s.rho, &t1, pk);
    zeroize_u32(t0.as_flattened_mut());
}

/// μ = H(tr ‖ M′, 64), with M′ given as a list of parts.
pub(crate) fn compute_mu(tr: &[u8], m_prime: &[&[u8]]) -> [u8; 64] {
    let mut h = Sponge::shake256();
    h.absorb(tr);
    for part in m_prime {
        h.absorb(part);
    }
    let mut mu = [0u8; 64];
    h.squeeze(&mut mu);
    mu
}

/// Absorbs w1Encode(w1) (Algorithm 28) into `h`.
fn absorb_w1(h: &mut Sponge, w1: &Poly, bits: u32) {
    let mut buf = [0u8; 32 * 6];
    let out = &mut buf[..32 * bits as usize];
    pack(w1, bits, out);
    h.absorb(out);
}

/// Algorithm 7: ML-DSA.Sign_internal, starting from μ (lines 6+). Writes σ.
pub(crate) fn sign_mu<const K: usize, const L: usize>(
    p: &Params,
    sk: &[u8],
    mu: &[u8; 64],
    rnd: &[u8; 32],
    sig: &mut [u8],
) {
    debug_assert_eq!(sk.len(), p.sk_len(K, L));
    debug_assert_eq!(sig.len(), p.sig_len(K, L));

    let sk = sk_decode::<K, L>(p, sk);
    let mut s1_hat = sk.s1;
    let mut s2_hat = sk.s2;
    let mut t0_hat = sk.t0;
    for v in s1_hat.iter_mut().chain(s2_hat.iter_mut()).chain(t0_hat.iter_mut()) {
        ntt(v);
    }
    let a = expand_a::<K, L>(&sk.rho);

    let mut rho2 = [0u8; 64];
    shake256(&[&sk.key, rnd, mu], &mut rho2);

    let clen = p.ctilde_len();
    let gamma1 = p.gamma1;
    let gamma2 = p.gamma2;
    let beta = p.beta;
    let mut kappa: u16 = 0;

    let mut y = [[0u32; N]; L];
    let mut z = [[0u32; N]; L];
    let mut w: [Poly; K];
    let mut h = [[0u32; N]; K];
    let mut ctilde = [0u8; 64];

    loop {
        // y ← ExpandMask(ρ″, κ)
        for (r, yr) in y.iter_mut().enumerate() {
            *yr = expand_mask_poly(gamma1, p.z_bits(), &rho2, kappa.wrapping_add(r as u16));
        }
        kappa = kappa.wrapping_add(L as u16);

        // w ← NTT^{-1}(Â ∘ NTT(y))
        let mut y_hat = y;
        for v in y_hat.iter_mut() {
            ntt(v);
        }
        w = mat_vec(&a, &y_hat);
        zeroize_u32(y_hat.as_flattened_mut());

        // c̃ ← H(μ ‖ w1Encode(w1), λ/4)
        let mut hasher = Sponge::shake256();
        hasher.absorb(mu);
        for wi in w.iter_mut() {
            inv_ntt(wi);
            let mut w1 = [0u32; N];
            for (o, &c) in w1.iter_mut().zip(wi.iter()) {
                *o = high_bits(gamma2, c);
            }
            absorb_w1(&mut hasher, &w1, p.w1_bits());
        }
        hasher.squeeze(&mut ctilde[..clen]);

        let mut c_hat = sample_in_ball(p.tau, &ctilde[..clen]);
        ntt(&mut c_hat);

        // z ← y + ⟨⟨c·s1⟩⟩ ; reject if ‖z‖∞ ≥ γ1 − β
        let mut reject = false;
        for r in 0..L {
            let mut cs1 = pointwise(&c_hat, &s1_hat[r]);
            inv_ntt(&mut cs1);
            for j in 0..N {
                z[r][j] = add(y[r][j], cs1[j]);
            }
            reject |= inf_norm(&z[r]) >= gamma1 - beta;
        }
        if reject {
            continue;
        }

        // r0 ← LowBits(w − ⟨⟨c·s2⟩⟩) ; reject if ‖r0‖∞ ≥ γ2 − β
        for i in 0..K {
            let mut cs2 = pointwise(&c_hat, &s2_hat[i]);
            inv_ntt(&mut cs2);
            for j in 0..N {
                w[i][j] = sub(w[i][j], cs2[j]);
                let r0 = low_bits(gamma2, w[i][j]);
                reject |= r0.unsigned_abs() >= gamma2 - beta;
            }
        }
        if reject {
            continue;
        }

        // h ← MakeHint(−⟨⟨c·t0⟩⟩, w − ⟨⟨c·s2⟩⟩ + ⟨⟨c·t0⟩⟩)
        let mut ones = 0usize;
        for i in 0..K {
            let mut ct0 = pointwise(&c_hat, &t0_hat[i]);
            inv_ntt(&mut ct0);
            reject |= inf_norm(&ct0) >= gamma2;
            for j in 0..N {
                let hint = make_hint(gamma2, neg(ct0[j]), add(w[i][j], ct0[j]));
                h[i][j] = hint;
                ones += hint as usize;
            }
        }
        if reject || ones > p.omega {
            continue;
        }

        // Algorithm 26: sigEncode(c̃, z mod± q, h)
        let zb = 32 * p.z_bits() as usize;
        sig[..clen].copy_from_slice(&ctilde[..clen]);
        for (r, zr) in z.iter().enumerate() {
            bit_pack(zr, p.z_bits(), gamma1, &mut sig[clen + r * zb..clen + (r + 1) * zb]);
        }
        hint_bit_pack::<K>(p.omega, &h, &mut sig[clen + L * zb..]);
        break;
    }

    zeroize_bytes(&mut rho2);
    zeroize_u32(s1_hat.as_flattened_mut());
    zeroize_u32(s2_hat.as_flattened_mut());
    zeroize_u32(t0_hat.as_flattened_mut());
    zeroize_u32(y.as_flattened_mut());
    zeroize_u32(z.as_flattened_mut());
    zeroize_u32(w.as_flattened_mut());
}

/// Algorithm 8: ML-DSA.Verify_internal, starting from μ.
pub(crate) fn verify_mu<const K: usize, const L: usize>(
    p: &Params,
    pk: &[u8],
    mu: &[u8; 64],
    sig: &[u8],
) -> bool {
    if pk.len() != p.pk_len(K) || sig.len() != p.sig_len(K, L) {
        return false;
    }
    let clen = p.ctilde_len();
    let zb = 32 * p.z_bits() as usize;

    // Algorithm 27: sigDecode
    let ctilde = &sig[..clen];
    let mut z_hat = [[0u32; N]; L];
    for (r, zr) in z_hat.iter_mut().enumerate() {
        *zr = bit_unpack(&sig[clen + r * zb..clen + (r + 1) * zb], p.z_bits(), p.gamma1);
    }
    let h = match hint_bit_unpack::<K>(p.omega, &sig[clen + L * zb..]) {
        Some(h) => h,
        None => return false,
    };
    if z_hat.iter().any(|zr| inf_norm(zr) >= p.gamma1 - p.beta) {
        return false;
    }

    let a = expand_a::<K, L>(&pk[..32]);
    let mut c_hat = sample_in_ball(p.tau, ctilde);
    ntt(&mut c_hat);
    for zr in z_hat.iter_mut() {
        ntt(zr);
    }
    let az = mat_vec(&a, &z_hat);

    let mut hasher = Sponge::shake256();
    hasher.absorb(mu);
    for i in 0..K {
        // Algorithm 23: pkDecode, then t1·2^d (always < q since t1 < 2^10).
        let mut t1 = [0u32; N];
        unpack(&pk[32 + 320 * i..32 + 320 * (i + 1)], T1_BITS, &mut t1);
        for c in t1.iter_mut() {
            *c <<= D;
        }
        ntt(&mut t1);
        let mut w = [0u32; N];
        for j in 0..N {
            w[j] = sub(az[i][j], mul(c_hat[j], t1[j]));
        }
        inv_ntt(&mut w);
        for j in 0..N {
            w[j] = use_hint(p.gamma2, h[i][j], w[j]);
        }
        absorb_w1(&mut hasher, &w, p.w1_bits());
    }
    let mut ctilde2 = [0u8; 64];
    hasher.squeeze(&mut ctilde2[..clen]);
    ct_eq(ctilde, &ctilde2[..clen]) == 1
}
