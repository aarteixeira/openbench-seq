//! The similarity search's prefilter: how many of a query's distinct k-mers a candidate holds (on either strand for
//! DNA), how many chance alone would give, and the candidate window holding the most of them, so only that window is
//! aligned. No index is kept: one pass over the candidate.

use std::collections::{HashMap, HashSet};

/// k-mer length: 16 bases for DNA (A, C, G and T only), 5 residues for protein (the 20 standard amino acids).
pub const DNA_K: usize = 16;
pub const PROTEIN_K: usize = 5;
const AMINO: &[u8; 20] = b"ACDEFGHIKLMNPQRSTVWY";

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Screen {
    /// The query's distinct k-mers.
    pub query_kmers: usize,
    /// Distinct query k-mers found in the candidate (either strand for DNA).
    pub shared: usize,
    /// Distinct query k-mers chance alone would put in a candidate this long: query_kmers × positions × strands /
    /// alphabet^k.
    pub expected: f64,
    /// The window to align: start and length on the candidate (on a circle it may run past the end and wrap).
    pub window: (usize, usize),
    /// The strand most hits lie on (1 or −1; 1 when there are none).
    pub strand: i8,
}

fn coder(protein: bool) -> (usize, u32, fn(u8) -> Option<u64>) {
    if protein {
        (PROTEIN_K, 5, |x| AMINO.iter().position(|&c| c == x.to_ascii_uppercase()).map(|p| p as u64))
    } else {
        (DNA_K, 2, |x| b"ACGT".iter().position(|&c| c == x.to_ascii_uppercase()).map(|p| p as u64))
    }
}

/// Every k-mer of `s` with its start (windows holding a letter outside the alphabet are skipped).
fn kmers(s: &[u8], protein: bool, mut f: impl FnMut(u64, usize)) {
    let (k, bits, code) = coder(protein);
    let mask = (1u64 << (bits * k as u32)) - 1;
    let (mut v, mut run) = (0u64, 0usize);
    for (p, &x) in s.iter().enumerate() {
        match code(x) {
            Some(c) => {
                v = ((v << bits) | c) & mask;
                run += 1;
                if run >= k {
                    f(v, p + 1 - k);
                }
            }
            None => run = 0,
        }
    }
}

/// The reverse complement of a 2-bit DNA k-mer (A=0, C=1, G=2, T=3).
fn revcomp_kmer(mut v: u64, k: usize) -> u64 {
    let mut r = 0u64;
    for _ in 0..k {
        r = (r << 2) | (3 - (v & 3));
        v >>= 2;
    }
    r
}

/// The distinct k-mers of `query`.
pub fn query_kmers(query: &[u8], protein: bool) -> usize {
    let mut set = HashSet::new();
    kmers(query, protein, |v, _| {
        set.insert(v);
    });
    set.len()
}

/// Screen `subject` for `query`'s k-mers. Hit positions beyond the first million are counted but not placed.
pub fn screen(query: &[u8], subject: &[u8], protein: bool, circular: bool) -> Screen {
    let (k, _, _) = coder(protein);
    // Each query k-mer, and for DNA its reverse complement, maps to the query k-mer it stands for and its strand.
    let mut wanted: HashMap<u64, (u64, i8)> = HashMap::new();
    kmers(query, protein, |v, _| {
        wanted.insert(v, (v, 1));
    });
    let query_kmers = wanted.len();
    if !protein {
        let fwd: Vec<u64> = wanted.keys().copied().collect();
        for v in fwd {
            wanted.entry(revcomp_kmer(v, k)).or_insert((v, -1));
        }
    }
    let n = subject.len();
    let text: Vec<u8> = if circular && n >= k { [subject, &subject[..k - 1]].concat() } else { subject.to_vec() };
    let (mut found, mut at) = (HashSet::new(), [Vec::new(), Vec::new()]);
    kmers(&text, protein, |v, p| {
        if p >= n {
            return;
        }
        if let Some(&(q, strand)) = wanted.get(&v) {
            found.insert(q);
            let side = &mut at[(strand < 0) as usize];
            if side.len() < 1_000_000 {
                side.push(p);
            }
        }
    });
    let strands = if protein { 1.0 } else { 2.0 };
    let alphabet: f64 = if protein { 20.0 } else { 4.0 };
    let positions = n.saturating_sub(k - 1) as f64;
    let expected = query_kmers as f64 * positions * strands / alphabet.powi(k as i32);
    // The strand with more hits, then the span of 2 × query length on it holding the most, widened by the query length
    // on each side.
    let side = if at[1].len() > at[0].len() { 1 } else { 0 };
    let hits = &at[side];
    let q = query.len();
    let window = if hits.is_empty() {
        (0, 0)
    } else {
        let (mut best, mut lo, mut range) = (0usize, 0usize, (hits[0], hits[0]));
        for hi in 0..hits.len() {
            while hits[hi] - hits[lo] > 2 * q {
                lo += 1;
            }
            if hi - lo + 1 > best {
                best = hi - lo + 1;
                range = (hits[lo], hits[hi]);
            }
        }
        let (start, end) = (range.0 as isize - q as isize, (range.1 + k + q) as isize);
        if circular {
            let len = ((end - start) as usize).min(n);
            (start.rem_euclid(n as isize) as usize, len)
        } else {
            let (s, e) = (start.max(0) as usize, (end as usize).min(n));
            (s, e - s)
        }
    };
    Screen { query_kmers, shared: found.len(), expected, window, strand: if side == 1 { -1 } else { 1 } }
}

#[cfg(test)]
mod tests {
    use super::*;

    // 40 bases with no repeated 16-mer, built by hand from a fixed rule.
    const Q: &[u8] = b"ACGTTGCAAGCTTCGAGGATCCTTAAGCGCATATGCCGAT";

    fn rc(s: &[u8]) -> Vec<u8> {
        crate::alphabet::revcomp(s)
    }

    #[test]
    fn counts_both_strands_and_places_the_window() {
        assert_eq!(query_kmers(Q, false), Q.len() - DNA_K + 1);
        let flank = [b'T'; 200];
        let subject = [&flank[..], Q, &flank[..]].concat();
        let s = screen(Q, &subject, false, false);
        assert_eq!((s.shared, s.strand), (25, 1));
        assert_eq!(s.window, (200 - 40, 40 + 25 + 15 + 40));
        // The reverse strand: every query k-mer found, on strand −1.
        let r = [&flank[..], &rc(Q)[..], &flank[..]].concat();
        let s = screen(Q, &r, false, false);
        assert_eq!((s.shared, s.strand), (25, -1));
    }

    #[test]
    fn a_circle_wraps_its_window_and_finds_k_mers_across_the_origin() {
        // The query split across the origin of a 100-base circle.
        let mut subject = vec![b'T'; 100];
        subject[..20].copy_from_slice(&Q[20..]);
        subject[80..].copy_from_slice(&Q[..20]);
        let linear = screen(Q, &subject, false, false);
        let circle = screen(Q, &subject, false, true);
        assert!(linear.shared < 25);
        assert_eq!(circle.shared, 25);
        assert_eq!(circle.window.1, 100); // the window wraps and covers the circle
    }

    #[test]
    fn chance_scales_with_length_and_alphabet() {
        let s = screen(Q, &[b'A'; 1_000_015], false, false);
        // 25 query k-mers × 1,000,000 positions × 2 strands / 4^16.
        assert!((s.expected - 25.0 * 1_000_000.0 * 2.0 / 4f64.powi(16)).abs() < 1e-9);
        assert_eq!(s.shared, 0);
        assert_eq!(s.window, (0, 0));
    }
}
