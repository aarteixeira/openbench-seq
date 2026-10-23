//! Pairwise alignment of `a` (the reference) and `b` (the query) with Gotoh's affine gaps: global, local, or glocal
//! (all of `b`, free ends in `a`), in a band of diagonals for long inputs, and Mott trimming of base qualities.
//!
//! Scores: DNA match 2, mismatch −3 (a position matches when the bases one may be are a subset of the other's code, the
//! IUPAC rule of `alphabet`); protein BLOSUM62 (NCBI's table, `blosum62.txt`). A gap of length L costs open + extend × L
//! (BLAST's convention): DNA 5 and 2, protein 11 and 1.

use crate::alphabet::{mask, satisfies};
use std::sync::OnceLock;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Mode {
    Global,
    Local,
    Glocal,
}

#[derive(Clone, Copy, Debug)]
pub struct Scoring {
    pub protein: bool,
    pub open: i32,
    pub extend: i32,
}

impl Scoring {
    pub const DNA: Self = Self { protein: false, open: 5, extend: 2 };
    pub const PROTEIN: Self = Self { protein: true, open: 11, extend: 1 };

    pub fn pair(&self, x: u8, y: u8) -> i32 {
        if self.protein { blosum62(x, y) } else if same_base(x, y) { 2 } else { -3 }
    }

    /// Whether two letters are the same for identity: bases by the IUPAC rule, residues by letter.
    pub fn same(&self, x: u8, y: u8) -> bool {
        if self.protein { x.eq_ignore_ascii_case(&y) } else { same_base(x, y) }
    }

    /// The score of `b` aligned to itself in full: the most any alignment of `b` can reach.
    pub fn best(&self, b: &[u8]) -> i32 {
        b.iter().map(|&y| self.pair(y, y)).sum()
    }
}

fn same_base(x: u8, y: u8) -> bool {
    let (p, q) = (mask(x), mask(y));
    satisfies(p, q) || satisfies(q, p)
}

const LETTERS: &[u8; 24] = b"ARNDCQEGHILKMFPSTWYVBZX*";

fn matrix() -> &'static [[i8; 24]; 24] {
    static M: OnceLock<[[i8; 24]; 24]> = OnceLock::new();
    M.get_or_init(|| {
        let mut m = [[0i8; 24]; 24];
        let rows = include_str!("blosum62.txt").lines().filter(|l| !l.starts_with('#') && !l.starts_with(' '));
        for (r, line) in rows.enumerate() {
            for (c, v) in line.split_whitespace().skip(1).enumerate() {
                m[r][c] = v.parse().expect("blosum62.txt is a 24 by 24 table");
            }
        }
        m
    })
}

fn index(x: u8) -> usize {
    LETTERS.iter().position(|&l| l == x.to_ascii_uppercase()).unwrap_or(22) // anything else reads as X
}

pub fn blosum62(x: u8, y: u8) -> i32 {
    matrix()[index(x)][index(y)] as i32
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Op {
    /// a and b letters aligned, the same.
    Match,
    /// a and b letters aligned, different.
    Mismatch,
    /// A letter of b with no letter of a (inserted in b).
    Ins,
    /// A letter of a with no letter of b (deleted from b).
    Del,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Alignment {
    pub score: i32,
    /// Aligned ranges [start, end) of a and of b.
    pub a: (usize, usize),
    pub b: (usize, usize),
    /// Runs of operations, in order.
    pub ops: Vec<(Op, u32)>,
}

impl Alignment {
    /// Identical letters over the alignment's columns.
    pub fn identity(&self) -> f64 {
        let cols: u32 = self.ops.iter().map(|&(_, n)| n).sum();
        let same: u32 = self.ops.iter().filter(|(o, _)| *o == Op::Match).map(|&(_, n)| n).sum();
        if cols == 0 { 0.0 } else { same as f64 / cols as f64 }
    }
}

/// Inputs of at most this many cells are aligned in full; larger ones in a band found from shared k-mers.
pub const FULL_CELLS: usize = 16_000_000;
/// The traceback of a banded alignment is at most this many cells (one byte each).
pub const MAX_CELLS: usize = 400_000_000;

const NEG: i32 = i32::MIN / 4;

/// Align a and b. `band`, when given, is the range of diagonals d = j − i (j in b, i in a) the alignment may use; when
/// absent, inputs over FULL_CELLS get one from `diagonal` (refused when a and b share no k-mer).
pub fn align(a: &[u8], b: &[u8], mode: Mode, sc: &Scoring, band: Option<(isize, isize)>) -> Result<Alignment, String> {
    let (n, m) = (a.len() as isize, b.len() as isize);
    let (lo, hi) = match band {
        Some(b) => b,
        None if (n as usize).saturating_mul(m as usize + 1) <= FULL_CELLS => (-n, m),
        None => {
            let d = diagonal(a, b, sc.protein).ok_or("the sequences share no stretch to place the alignment on, and are too long to align in full")?;
            let pad = 100 + match mode { Mode::Global => (m - n).abs(), _ => m / 20 };
            match mode {
                // Global runs from (0, 0) to (n, m): the band must hold both corners.
                Mode::Global => ((d - pad).min(0).min(m - n), (d + pad).max(0).max(m - n)),
                _ => (d - pad, d + pad),
            }
        }
    };
    let (lo, hi) = (lo.max(-n), hi.min(m));
    if lo > hi {
        return Err("the band holds no cell".into());
    }
    let w = (hi - lo + 1) as usize;
    if (n as usize + 1).saturating_mul(w) > MAX_CELLS {
        return Err("too long to align: the band would be over 400 million cells".into());
    }
    dp(a, b, mode, sc, lo, w)
}

// Traceback byte: bits 0-1 the state before an H cell (0 H, 1 E, 2 F, 3 the start); bit 2 E opened (else extended),
// bit 3 E opened from F (else from H); bit 4 F opened, bit 5 F opened from E. E consumes b (Ins), F consumes a (Del).
const START: u8 = 3;

fn dp(a: &[u8], b: &[u8], mode: Mode, sc: &Scoring, lo: isize, w: usize) -> Result<Alignment, String> {
    let (n, m) = (a.len(), b.len());
    let (open, ext) = (sc.open + sc.extend, sc.extend);
    let mut trace = vec![0u8; (n + 1) * w];
    let (mut hp, mut ep, mut fp) = (vec![NEG; w], vec![NEG; w], vec![NEG; w]);
    let (mut hc, mut ec, mut fc) = (vec![NEG; w], vec![NEG; w], vec![NEG; w]);
    // The best end: (score, i, j, state).
    let mut end = (NEG, 0usize, 0usize, 0u8);
    for i in 0..=n {
        hc.fill(NEG);
        ec.fill(NEG);
        fc.fill(NEG);
        let jlo = (i as isize + lo).max(0) as usize;
        let jhi = ((i as isize + lo + w as isize - 1).min(m as isize)).max(-1);
        if jhi < jlo as isize {
            std::mem::swap(&mut hp, &mut hc);
            std::mem::swap(&mut ep, &mut ec);
            std::mem::swap(&mut fp, &mut fc);
            continue;
        }
        for j in jlo..=jhi as usize {
            let k = (j as isize - i as isize - lo) as usize;
            let t = &mut trace[i * w + k];
            if i == 0 && j == 0 {
                hc[k] = 0;
                *t = START;
                continue;
            }
            // H: a[i-1] against b[j-1].
            if i > 0 && j > 0 {
                let (h, e, f) = (hp[k], ep[k], fp[k]);
                let (mut best, mut from) = (h, 0u8);
                if e > best { (best, from) = (e, 1) }
                if f > best { (best, from) = (f, 2) }
                if mode == Mode::Local && best < 0 { (best, from) = (0, START) }
                if best > NEG {
                    hc[k] = best + sc.pair(a[i - 1], b[j - 1]);
                    *t = from;
                }
            } else if j == 0 && mode != Mode::Global {
                // Glocal and local start anywhere along a.
                hc[k] = 0;
                *t = START;
            }
            // E: b[j-1] against nothing (from the left).
            if j > 0 && k > 0 {
                let (h, e, f) = (hc[k - 1], ec[k - 1], fc[k - 1]);
                let ex = if e > NEG { e - ext } else { NEG };
                let (op, from) = if f > h { (f, 1u8) } else { (h, 0u8) };
                let op = if op > NEG { op - open } else { NEG };
                if op >= ex && op > NEG {
                    ec[k] = op;
                    *t |= 4 | (from << 3);
                } else {
                    ec[k] = ex;
                }
            }
            // F: a[i-1] against nothing (from above). Before b starts (j = 0) only a global alignment pays for it.
            if i > 0 && k + 1 < w && !(j == 0 && mode != Mode::Global) {
                let (h, e, f) = (hp[k + 1], ep[k + 1], fp[k + 1]);
                let ex = if f > NEG { f - ext } else { NEG };
                let (op, from) = if e > h { (e, 1u8) } else { (h, 0u8) };
                let op = if op > NEG { op - open } else { NEG };
                if op >= ex && op > NEG {
                    fc[k] = op;
                    *t |= 16 | (from << 5);
                } else {
                    fc[k] = ex;
                }
            }
            match mode {
                Mode::Local if hc[k] > end.0 => end = (hc[k], i, j, 0),
                Mode::Glocal if j == m => {
                    for (s, v) in [(0u8, hc[k]), (1, ec[k]), (2, fc[k])] {
                        if v > end.0 {
                            end = (v, i, j, s);
                        }
                    }
                }
                Mode::Global if i == n && j == m => {
                    for (s, v) in [(0u8, hc[k]), (1, ec[k]), (2, fc[k])] {
                        if v > end.0 {
                            end = (v, i, j, s);
                        }
                    }
                }
                _ => {}
            }
        }
        std::mem::swap(&mut hp, &mut hc);
        std::mem::swap(&mut ep, &mut ec);
        std::mem::swap(&mut fp, &mut fc);
    }
    if end.0 <= NEG {
        return Err("no alignment fits the band".into());
    }
    // Trace back from the end.
    let (score, mut i, mut j, mut state) = end;
    let (ai, bj) = (i, j);
    let mut ops: Vec<(Op, u32)> = Vec::new();
    let mut push = |o: Op| match ops.last_mut() {
        Some((p, c)) if *p == o => *c += 1,
        _ => ops.push((o, 1)),
    };
    loop {
        let k = (j as isize - i as isize - lo) as usize;
        let t = trace[i * w + k];
        match state {
            0 => {
                // An H cell on the first row or column is a start, not a column of the alignment; a local alignment's
                // first column is a cell whose predecessor is the start.
                if i == 0 || j == 0 {
                    break;
                }
                let from = t & 3;
                push(if sc.same(a[i - 1], b[j - 1]) { Op::Match } else { Op::Mismatch });
                i -= 1;
                j -= 1;
                if from == START {
                    break;
                }
                state = from;
            }
            1 => {
                push(Op::Ins);
                j -= 1;
                if t & 4 != 0 {
                    state = if t & 8 != 0 { 2 } else { 0 };
                }
            }
            _ => {
                push(Op::Del);
                i -= 1;
                if t & 16 != 0 {
                    state = if t & 32 != 0 { 1 } else { 0 };
                }
            }
        }
    }
    ops.reverse();
    Ok(Alignment { score, a: (i, ai), b: (j, bj), ops })
}

/// The diagonal (j − i) most shared k-mers of a and b lie on: 12 bases for DNA (A, C, G and T only), 5 residues for
/// protein. None when they share none.
pub fn diagonal(a: &[u8], b: &[u8], protein: bool) -> Option<isize> {
    let (k, bits) = if protein { (5usize, 5u32) } else { (12, 2) };
    let code = |x: u8| -> Option<u64> {
        let x = x.to_ascii_uppercase();
        if protein { (x.is_ascii_uppercase()).then(|| (x - b'A') as u64) } else { b"ACGT".iter().position(|&c| c == x).map(|p| p as u64) }
    };
    let mask = (1u64 << (bits * k as u32)) - 1;
    let kmers = |s: &[u8], mut f: Box<dyn FnMut(u64, usize) + '_>| {
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
    };
    let mut at: std::collections::HashMap<u64, Vec<usize>> = std::collections::HashMap::new();
    kmers(b, Box::new(|v, p| at.entry(v).or_default().push(p)));
    // A k-mer repeated many times in b says little about where b lies: it is left out.
    at.retain(|_, v| v.len() <= 8);
    let mut votes: std::collections::HashMap<isize, u32> = std::collections::HashMap::new();
    kmers(a, Box::new(|v, i| {
        if let Some(js) = at.get(&v) {
            for &j in js {
                *votes.entry(j as isize - i as isize).or_default() += 1;
            }
        }
    }));
    votes.into_iter().max_by_key(|&(d, c)| (c, std::cmp::Reverse(d.abs()))).map(|(d, _)| d)
}

/// Mott's trimming (as Phred's): with p the error probability of each base's Phred value, the stretch where the sum
/// of (limit − p) is largest. [start, end), empty when no base is better than the limit.
pub fn trim(quality: &[u8], limit: f64) -> (usize, usize) {
    let (mut best, mut range) = (0.0f64, (0, 0));
    let (mut sum, mut start) = (0.0f64, 0usize);
    for (i, &q) in quality.iter().enumerate() {
        sum += limit - 10f64.powf(-(q as f64) / 10.0);
        if sum <= 0.0 {
            sum = 0.0;
            start = i + 1;
        } else if sum > best {
            best = sum;
            range = (start, i + 1);
        }
    }
    range
}

#[cfg(test)]
mod tests {
    use super::*;
    use Op::*;

    fn run(a: &str, b: &str, mode: Mode) -> Alignment {
        align(a.as_bytes(), b.as_bytes(), mode, &Scoring::DNA, None).unwrap()
    }

    #[test]
    fn blosum62_is_ncbi_s_table() {
        // Symmetric, and known entries (diagonal and pairs) of NCBI's file.
        for &x in LETTERS {
            for &y in LETTERS {
                assert_eq!(blosum62(x, y), blosum62(y, x), "{} {}", x as char, y as char);
            }
        }
        let diag: Vec<i32> = b"ARNDCQEGHILKMFPSTWYV".iter().map(|&x| blosum62(x, x)).collect();
        assert_eq!(diag, [4, 5, 6, 6, 9, 5, 5, 6, 8, 4, 4, 5, 5, 6, 7, 4, 5, 11, 7, 4]);
        assert_eq!([blosum62(b'A', b'R'), blosum62(b'I', b'V'), blosum62(b'W', b'F'), blosum62(b'E', b'Z'), blosum62(b'*', b'*'), blosum62(b'J', b'A')], [-1, 3, 1, 4, 1, 0]);
    }

    #[test]
    fn global_alignments_by_hand() {
        // Identical: 8 matches, 16.
        assert_eq!(run("ACGTACGT", "ACGTACGT", Mode::Global), Alignment { score: 16, a: (0, 8), b: (0, 8), ops: vec![(Match, 8)] });
        // One substitution: 7 × 2 − 3 = 11.
        assert_eq!(run("ACGTACGT", "ACGAACGT", Mode::Global).ops, vec![(Match, 3), (Mismatch, 1), (Match, 4)]);
        assert_eq!(run("ACGTACGT", "ACGAACGT", Mode::Global).score, 11);
        // One base deleted from b (the C at 4 of a; no other placement scores as well): 7 × 2 − (5 + 2) = 7.
        let d = run("ACGTCACG", "ACGTACG", Mode::Global);
        assert_eq!((d.score, d.ops.clone()), (7, vec![(Match, 4), (Del, 1), (Match, 3)]));
        // One base inserted in b, the mirror.
        assert_eq!(run("ACGTACG", "ACGTCACG", Mode::Global).ops, vec![(Match, 4), (Ins, 1), (Match, 3)]);
        // Affine: one gap of three (5 + 6 = 11) rather than three of one (21): 9 × 2 − 11 = 7.
        let g = run("ACGTAAACGTTT", "ACGTCGTTT", Mode::Global);
        assert_eq!((g.score, g.ops.iter().filter(|(o, _)| *o == Del).count()), (7, 1));
    }

    #[test]
    fn glocal_and_local_by_hand() {
        // All of b inside a, a's ends free: 7 matches, 14, at 4-11 of a.
        assert_eq!(run("TTTTACGTACGTTTTT", "ACGTACG", Mode::Glocal), Alignment { score: 14, a: (4, 11), b: (0, 7), ops: vec![(Match, 7)] });
        // Local: the best common stretch, ACGT (8), at 3-7 of a and 2-6 of b.
        assert_eq!(run("GGGACGTGGG", "TTACGTTT", Mode::Local), Alignment { score: 8, a: (3, 7), b: (2, 6), ops: vec![(Match, 4)] });
        // IUPAC: N in b matches any base of a.
        assert_eq!(run("ACGTACGT", "ACNTACGT", Mode::Global).ops, vec![(Match, 8)]);
    }

    #[test]
    fn banded_equals_full_and_long_inputs_are_banded() {
        let mut r = 7u32;
        let mut bases = |n: usize| -> Vec<u8> { (0..n).map(|_| { r ^= r << 13; r ^= r >> 17; r ^= r << 5; b"ACGT"[(r >> 30) as usize] }).collect() };
        let a = bases(3000);
        let mut b = a[1000..1800].to_vec();
        b[100] = if b[100] == b'A' { b'C' } else { b'A' };
        b.remove(400);
        let full = align(&a, &b, Mode::Glocal, &Scoring::DNA, None).unwrap();
        let banded = align(&a, &b, Mode::Glocal, &Scoring::DNA, Some((-1000 - 150, -1000 + 150))).unwrap();
        assert_eq!((full.score, full.a, full.b), (banded.score, banded.a, banded.b));
        assert_eq!((full.a, full.score), ((1000, 1800), 798 * 2 - 3 - 7));
        // Over FULL_CELLS: the band comes from the shared 12-mers, with the same answer.
        let long = bases(30_000);
        let mut read = long[20_000..21_000].to_vec();
        read[500] = if read[500] == b'G' { b'T' } else { b'G' };
        let x = align(&long, &read, Mode::Glocal, &Scoring::DNA, None).unwrap();
        assert_eq!((x.a, x.ops.clone()), ((20_000, 21_000), vec![(Match, 500), (Mismatch, 1), (Match, 499)]));
        assert_eq!(diagonal(&long, &read, false), Some(-20_000));
        // Nothing shared and too long: refused, not filled.
        assert!(align(&vec![b'A'; 5000], &vec![b'C'; 5000], Mode::Global, &Scoring::DNA, None).is_err());
    }

    #[test]
    fn protein_alignment_uses_blosum62() {
        let p = align(b"MKVLAAGG", b"MKILAAGG", Mode::Global, &Scoring::PROTEIN, None).unwrap();
        // V/I is 3: 5 + 5 + 3 + 4 + 4 + 4 + 6 + 6.
        assert_eq!((p.score, p.ops.clone()), (37, vec![(Match, 2), (Mismatch, 1), (Match, 5)]));
    }

    #[test]
    fn mott_trimming_by_hand() {
        // Q10 is p 0.1 (0.05 − 0.1 < 0), Q30 is p 0.001 (+0.049): the Q30 stretch is kept.
        assert_eq!(trim(&[10, 10, 30, 30, 30, 30, 10, 10], 0.05), (2, 6));
        // A single Q10 between good stretches is kept when the stretch after it outweighs it.
        assert_eq!(trim(&[30, 30, 30, 10, 30, 30, 30, 10, 10], 0.05), (0, 7));
        assert_eq!(trim(&[5, 5, 5], 0.05), (0, 0));
    }
}
