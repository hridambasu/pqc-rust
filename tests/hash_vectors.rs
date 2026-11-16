//! SHA-2 / SHA-3 / SHAKE against reference digests produced by CPython's
//! hashlib (OpenSSL), covering every length 0..300 and every rate boundary.

mod common;
use common::*;
use fips_pqc::hash::*;

fn compute(alg: &str, msg: &[u8], out_len: usize) -> Vec<u8> {
    match alg {
        "sha224" => Sha224::digest(msg).to_vec(),
        "sha256" => Sha256::digest(msg).to_vec(),
        "sha384" => Sha384::digest(msg).to_vec(),
        "sha512" => Sha512::digest(msg).to_vec(),
        "sha512_224" => Sha512_224::digest(msg).to_vec(),
        "sha512_256" => Sha512_256::digest(msg).to_vec(),
        "sha3_224" => Sha3_224::digest(msg).to_vec(),
        "sha3_256" => Sha3_256::digest(msg).to_vec(),
        "sha3_384" => Sha3_384::digest(msg).to_vec(),
        "sha3_512" => Sha3_512::digest(msg).to_vec(),
        "shake_128" => {
            let mut o = vec![0; out_len];
            Shake128::digest(msg, &mut o);
            o
        }
        "shake_256" => {
            let mut o = vec![0; out_len];
            Shake256::digest(msg, &mut o);
            o
        }
        _ => panic!("unknown alg {alg}"),
    }
}

/// Same computation, but feeding input in irregular chunks and squeezing in pieces.
fn compute_chunked(alg: &str, msg: &[u8], out_len: usize, step: usize) -> Vec<u8> {
    macro_rules! fixed {
        ($t:ty) => {{
            let mut h = <$t>::new();
            for c in msg.chunks(step) {
                h.update(c);
            }
            h.finalize().to_vec()
        }};
    }
    macro_rules! xof {
        ($t:ty) => {{
            let mut h = <$t>::new();
            for c in msg.chunks(step) {
                h.update(c);
            }
            let mut r = h.finalize_xof();
            let mut o = vec![0; out_len];
            for c in o.chunks_mut(step) {
                r.read(c);
            }
            o
        }};
    }
    match alg {
        "sha224" => fixed!(Sha224),
        "sha256" => fixed!(Sha256),
        "sha384" => fixed!(Sha384),
        "sha512" => fixed!(Sha512),
        "sha512_224" => fixed!(Sha512_224),
        "sha512_256" => fixed!(Sha512_256),
        "sha3_224" => fixed!(Sha3_224),
        "sha3_256" => fixed!(Sha3_256),
        "sha3_384" => fixed!(Sha3_384),
        "sha3_512" => fixed!(Sha3_512),
        "shake_128" => xof!(Shake128),
        "shake_256" => xof!(Shake256),
        _ => panic!("unknown alg {alg}"),
    }
}

#[test]
fn reference_digests() {
    let recs = load("hashes.txt");
    let mut per_alg = std::collections::BTreeMap::new();
    for (i, r) in recs.iter().enumerate() {
        let alg = r.get("alg");
        let msg = r.hex("msg");
        let expect = r.hex("out");
        assert_eq!(compute(alg, &msg, expect.len()), expect, "{alg} len={}", msg.len());
        let step = 1 + (i % 97);
        assert_eq!(compute_chunked(alg, &msg, expect.len(), step), expect, "{alg} chunked step={step}");
        *per_alg.entry(alg.to_string()).or_insert(0) += 1;
    }
    assert_eq!(per_alg.len(), 12);
    assert!(per_alg.values().all(|&n| n > 300));
}

#[test]
fn nist_fips180_examples() {
    let msg = b"abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq";
    assert_eq!(
        Sha256::digest(msg).to_vec(),
        hex("248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1")
    );
    assert_eq!(
        Sha224::digest(b"abc").to_vec(),
        hex("23097d223405d8228642a477bda255b32aadbce4bda0b3f7e36c9da7")
    );
    assert_eq!(
        Sha384::digest(b"abc").to_vec(),
        hex("cb00753f45a35e8bb5a03d699ac65007272c32ab0eded1631a8b605a43ff5bed8086072ba1e7cc2358baeca134c825a7")
    );
}

#[test]
fn million_a() {
    let msg = vec![b'a'; 1_000_000];
    assert_eq!(
        Sha256::digest(&msg).to_vec(),
        hex("cdc76e5c9914fb9281a1c7e284d73e67f1809a48a497200e046d39ccc7112cd0")
    );
    assert_eq!(
        Sha3_256::digest(&msg).to_vec(),
        hex("5c8875ae474a3634ba4fd55ec85bffd661f32aca75c6d699d0cdcb6c115891c1")
    );
}
