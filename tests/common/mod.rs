//! Shared helpers for integration tests: a parser for the vector files in
//! `tests/vectors/`, hex decoding, and a deterministic test RNG.

#![allow(dead_code)]

use fips_pqc::hash::{Shake256, XofReader};
use fips_pqc::{CryptoRng, Error};
use std::collections::HashMap;
use std::path::PathBuf;

/// One test case: a map of `key = value` lines.
#[derive(Debug, Clone)]
pub struct Record(pub HashMap<String, String>);

impl Record {
    pub fn get(&self, k: &str) -> &str {
        self.0.get(k).unwrap_or_else(|| panic!("missing field {k:?} in {:?}", self.0.get("tcId")))
    }
    pub fn opt(&self, k: &str) -> Option<&str> {
        self.0.get(k).map(String::as_str)
    }
    pub fn hex(&self, k: &str) -> Vec<u8> {
        hex(self.get(k))
    }
    pub fn bool(&self, k: &str) -> bool {
        match self.get(k) {
            "true" => true,
            "false" => false,
            v => panic!("bad bool {v}"),
        }
    }
    pub fn id(&self) -> &str {
        self.opt("tcId").unwrap_or("?")
    }
}

/// Loads `tests/vectors/<name>`; records are separated by blank lines.
pub fn load(name: &str) -> Vec<Record> {
    let path: PathBuf = [env!("CARGO_MANIFEST_DIR"), "tests", "vectors", name].iter().collect();
    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{path:?}: {e}"));
    let mut out = Vec::new();
    let mut cur = HashMap::new();
    for line in text.lines() {
        let line = line.trim_end();
        if line.starts_with('#') {
            continue;
        }
        if line.is_empty() {
            if !cur.is_empty() {
                out.push(Record(std::mem::take(&mut cur)));
            }
            continue;
        }
        let (k, v) = line.split_once(" = ").or_else(|| line.strip_suffix(" =").map(|k| (k, ""))).expect(line);
        cur.insert(k.to_string(), v.to_string());
    }
    if !cur.is_empty() {
        out.push(Record(cur));
    }
    assert!(!out.is_empty(), "no records in {name}");
    out
}

pub fn hex(s: &str) -> Vec<u8> {
    assert!(s.len() % 2 == 0, "odd hex length");
    (0..s.len()).step_by(2).map(|i| u8::from_str_radix(&s[i..i + 2], 16).unwrap()).collect()
}

pub fn arr<const N: usize>(v: &[u8]) -> [u8; N] {
    v.try_into().unwrap_or_else(|_| panic!("expected {N} bytes, got {}", v.len()))
}

/// Deterministic SHAKE256-based RNG for reproducible tests.
pub struct TestRng(XofReader);

impl TestRng {
    pub fn new(seed: u64) -> Self {
        let mut h = Shake256::new();
        h.update(b"fips-pqc test rng");
        h.update(&seed.to_le_bytes());
        Self(h.finalize_xof())
    }
    pub fn bytes<const N: usize>(&mut self) -> [u8; N] {
        let mut b = [0u8; N];
        self.0.read(&mut b);
        b
    }
    pub fn vec(&mut self, n: usize) -> Vec<u8> {
        let mut b = vec![0u8; n];
        self.0.read(&mut b);
        b
    }
    pub fn below(&mut self, n: usize) -> usize {
        (u64::from_le_bytes(self.bytes()) % n as u64) as usize
    }
}

impl CryptoRng for TestRng {
    fn fill_bytes(&mut self, dest: &mut [u8]) -> Result<(), Error> {
        self.0.read(dest);
        Ok(())
    }
}

/// An RNG that always fails.
pub struct FailingRng;

impl CryptoRng for FailingRng {
    fn fill_bytes(&mut self, _: &mut [u8]) -> Result<(), Error> {
        Err(Error::RandomnessFailure)
    }
}
