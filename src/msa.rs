//! Progressive multiple alignment of 2-200 sequences of one alphabet, as supplied (no reverse-complement step).
//!
//! Distance: d = 1 − shared / min(|Ka|, |Kb|) over each sequence's distinct k-mers (12 bases for DNA, A/C/G/T only; 5
//! residues for protein, the 20 standard amino acids); no such k-mer gives d = 1. Guide tree: UPGMA, ties broken by
//! the lower indexes. Profiles are joined bottom-up in the scorer's own symbols (DNA: the 15 IUPAC codes with the
//! subset rule 2/−3; protein: the 24 BLOSUM62 symbols): a column pair scores the frequency-weighted sum of the pair
//! scores over the letters present (each letter's frequency is its count over all the profile's rows, so a column's
//! gaps lower its weight); a gap run inserted into a profile costs the affine penalties once (DNA 5/2, protein 11/1,
//! open + extend × length), whatever the number of rows; end gaps are free. No refinement pass.
//!
//! Bounds: the static ones (count, lengths, the guide tree's estimated work) are checked before any work; the ones
//! that depend on results (each join's actual profile widths, the work so far, the aligned width) before each join's
//! matrix and before the output is built.

use crate::alphabet::mask;
use std::cell::Cell;
use std::collections::HashSet;

pub struct Limits {
    pub count: usize,
    pub length: usize,
    pub total: usize,
    /// Cells of one join's matrix (one traceback byte each).
    pub join_cells: usize,
    /// Cells × symbol terms, summed over the joins.
    pub work: f64,
}

pub const LIMITS: Limits = Limits { count: 200, length: 10_000, total: 1_000_000, join_cells: 100_000_000, work: 1.5e10 };

thread_local! {
    /// Join matrices and outputs allocated (tests read it to see a refusal came before the allocation it guards).
    pub static ALLOCATIONS: Cell<usize> = const { Cell::new(0) };
}

const DNA_SYMBOLS: usize = 15;
const PROTEIN_SYMBOLS: usize = 24;
const GAP: u8 = b'-';

fn symbol(protein: bool, b: u8) -> Option<usize> {
    if protein {
        b"ARNDCQEGHILKMFPSTWYVBZX*".iter().position(|&l| l == b.to_ascii_uppercase())
    } else {
        match mask(b) {
            0 => None,
            m => Some(m as usize - 1),
        }
    }
}

fn pair_table(protein: bool) -> Vec<Vec<f64>> {
    let s = if protein { PROTEIN_SYMBOLS } else { DNA_SYMBOLS };
    let letters: Vec<u8> = if protein { b"ARNDCQEGHILKMFPSTWYVBZX*".to_vec() } else { (1..=15u8).map(crate::alphabet::letter).collect() };
    let sc = if protein { crate::align::Scoring::PROTEIN } else { crate::align::Scoring::DNA };
    (0..s).map(|a| (0..s).map(|b| sc.pair(letters[a], letters[b]) as f64).collect()).collect()
}

fn kmer_set(s: &[u8], protein: bool) -> HashSet<u64> {
    let (k, bits): (usize, u32) = if protein { (5, 5) } else { (12, 2) };
    let code = |x: u8| -> Option<u64> {
        let x = x.to_ascii_uppercase();
        if protein { b"ACDEFGHIKLMNPQRSTVWY".iter().position(|&c| c == x).map(|p| p as u64) } else { b"ACGT".iter().position(|&c| c == x).map(|p| p as u64) }
    };
    let maskv = (1u64 << (bits * k as u32)) - 1;
    let (mut v, mut run, mut set) = (0u64, 0usize, HashSet::new());
    for &x in s {
        match code(x) {
            Some(c) => {
                v = ((v << bits) | c) & maskv;
                run += 1;
                if run >= k {
                    set.insert(v);
                }
            }
            None => run = 0,
        }
    }
    set
}

/// The guide tree's joins, in order: each join names two clusters (an input index below n, or n + the join that made
/// it).
pub fn guide_tree(seqs: &[Vec<u8>], protein: bool) -> Vec<(usize, usize)> {
    let n = seqs.len();
    let sets: Vec<HashSet<u64>> = seqs.iter().map(|s| kmer_set(s, protein)).collect();
    let mut d = vec![vec![0f64; n]; n];
    for i in 0..n {
        for j in i + 1..n {
            let (a, b) = if sets[i].len() <= sets[j].len() { (&sets[i], &sets[j]) } else { (&sets[j], &sets[i]) };
            let shared = a.iter().filter(|v| b.contains(v)).count();
            let v = if a.is_empty() { 1.0 } else { 1.0 - shared as f64 / a.len() as f64 };
            d[i][j] = v;
            d[j][i] = v;
        }
    }
    // Live clusters: (label, size); UPGMA averages distances by size.
    let mut live: Vec<(usize, usize)> = (0..n).map(|i| (i, 1)).collect();
    let mut joins = Vec::new();
    while live.len() > 1 {
        let (mut bi, mut bj, mut bd) = (0, 1, f64::INFINITY);
        for i in 0..live.len() {
            for j in i + 1..live.len() {
                if d[i][j] < bd {
                    (bi, bj, bd) = (i, j, d[i][j]);
                }
            }
        }
        let (li, si) = live[bi];
        let (lj, sj) = live[bj];
        joins.push((li, lj));
        let row: Vec<f64> = (0..live.len()).map(|k| (d[bi][k] * si as f64 + d[bj][k] * sj as f64) / (si + sj) as f64).collect();
        // The new cluster takes bi's place; bj is removed.
        live[bi] = (n + joins.len() - 1, si + sj);
        for k in 0..live.len() {
            d[bi][k] = row[k];
            d[k][bi] = row[k];
        }
        d[bi][bi] = 0.0;
        live.remove(bj);
        d.remove(bj);
        for r in d.iter_mut() {
            r.remove(bj);
        }
    }
    joins
}

/// Align `seqs`: the rows in input order, gaps as `-`.
pub fn msa(seqs: &[Vec<u8>], protein: bool) -> Result<Vec<Vec<u8>>, String> {
    msa_with(seqs, protein, &LIMITS)
}

pub fn msa_with(seqs: &[Vec<u8>], protein: bool, lim: &Limits) -> Result<Vec<Vec<u8>>, String> {
    let n = seqs.len();
    if n < 2 || n > lim.count {
        return Err(format!("{n} sequences: 2 to {} are aligned at once", lim.count));
    }
    for (i, s) in seqs.iter().enumerate() {
        if s.is_empty() || s.len() > lim.length {
            return Err(format!("sequence {} has {} letters: 1 to {} are aligned", i + 1, s.len(), lim.length));
        }
        if let Some(p) = s.iter().position(|&b| symbol(protein, b).is_none()) {
            return Err(format!("sequence {} has a letter that is not {} at position {}", i + 1, if protein { "an amino acid" } else { "a DNA code" }, p + 1));
        }
    }
    let total: usize = seqs.iter().map(|s| s.len()).sum();
    if total > lim.total {
        return Err(format!("the sequences hold {total} letters: at most {} are aligned at once", lim.total));
    }
    let terms = if protein { PROTEIN_SYMBOLS } else { DNA_SYMBOLS } as f64;
    let joins = guide_tree(seqs, protein);
    // The estimate: each cluster at least as wide as its longest member.
    let mut widest: Vec<usize> = seqs.iter().map(|s| s.len()).collect();
    let mut estimate = 0f64;
    for &(a, b) in &joins {
        let cells = (widest[a] + 1) as f64 * (widest[b] + 1) as f64;
        if cells > lim.join_cells as f64 {
            return Err(format!("a join of the guide tree needs {cells} cells: at most {} a join", lim.join_cells));
        }
        estimate += cells * terms;
        widest.push(widest[a].max(widest[b]));
    }
    if estimate > lim.work {
        return Err(format!("the alignment needs about {estimate:.3e} steps: at most {:.3e}", lim.work));
    }
    let longest = seqs.iter().map(|s| s.len()).max().unwrap_or(0);
    let table = pair_table(protein);
    // Each cluster: its member indexes and their gapped rows.
    let mut clusters: Vec<(Vec<usize>, Vec<Vec<u8>>)> = seqs.iter().enumerate().map(|(i, s)| (vec![i], vec![s.to_vec()])).collect();
    let (mut work, mut done) = (0f64, 0usize);
    let total_cells = estimate / terms;
    for &(a, b) in &joins {
        let (wa, wb) = (clusters[a].1[0].len(), clusters[b].1[0].len());
        let cells = (wa + 1) * (wb + 1);
        if cells > lim.join_cells {
            return Err(format!("a join of the guide tree needs {cells} cells: at most {} a join", lim.join_cells));
        }
        work += cells as f64 * terms;
        if work > lim.work {
            return Err(format!("the alignment needs over {:.3e} steps: at most {:.3e}", work, lim.work));
        }
        let ops = join(&clusters[a].1, &clusters[b].1, protein, &table, &mut done, total_cells as usize);
        let width = ops.len();
        if width > 2 * longest {
            return Err(format!("the aligned width would be {width}: at most twice the longest sequence ({})", 2 * longest));
        }
        let mut rows = Vec::new();
        for (side, src) in [(0u8, &clusters[a].1), (1u8, &clusters[b].1)] {
            for r in src.iter() {
                let mut out = Vec::with_capacity(width);
                let mut p = 0;
                for &o in &ops {
                    // o: 0 both columns, 1 only b's column, 2 only a's column.
                    let takes = match (side, o) { (_, 0) => true, (0, 2) => true, (1, 1) => true, _ => false };
                    if takes {
                        out.push(r[p]);
                        p += 1;
                    } else {
                        out.push(GAP);
                    }
                }
                rows.push(out);
            }
        }
        let mut members = clusters[a].0.clone();
        members.extend(&clusters[b].0);
        clusters.push((members, rows));
    }
    ALLOCATIONS.with(|c| c.set(c.get() + 1));
    let (members, rows) = clusters.pop().expect("the tree has a root");
    let mut out = vec![Vec::new(); n];
    for (m, r) in members.into_iter().zip(rows) {
        out[m] = r;
    }
    Ok(out)
}

/// One profile's column frequencies (count over all rows) per symbol, kept sparse: (symbol, frequency) per column.
fn profile(rows: &[Vec<u8>], protein: bool) -> Vec<Vec<(usize, f64)>> {
    let width = rows[0].len();
    let s = if protein { PROTEIN_SYMBOLS } else { DNA_SYMBOLS };
    let r = rows.len() as f64;
    (0..width)
        .map(|c| {
            let mut counts = vec![0u32; s];
            for row in rows {
                if let Some(x) = symbol(protein, row[c]) {
                    counts[x] += 1;
                }
            }
            counts.iter().enumerate().filter(|&(_, &k)| k > 0).map(|(x, &k)| (x, k as f64 / r)).collect()
        })
        .collect()
}

/// Align two profiles (global, end gaps free): one op per output column, 0 both, 1 only b's column (a gap in a), 2
/// only a's column (a gap in b).
fn join(a: &[Vec<u8>], b: &[Vec<u8>], protein: bool, table: &[Vec<f64>], done: &mut usize, total: usize) -> Vec<u8> {
    let (pa, pb) = (profile(a, protein), profile(b, protein));
    let (n, m) = (pa.len(), pb.len());
    let s = table.len();
    // a's columns against every symbol, so a cell is a dot product over b's letters present.
    let va: Vec<Vec<f64>> = pa.iter().map(|col| (0..s).map(|y| col.iter().map(|&(x, f)| f * table[x][y]).sum()).collect()).collect();
    let sc = if protein { crate::align::Scoring::PROTEIN } else { crate::align::Scoring::DNA };
    let (open, ext) = ((sc.open + sc.extend) as f64, sc.extend as f64);
    let neg = f64::NEG_INFINITY;
    ALLOCATIONS.with(|c| c.set(c.get() + 1));
    let mut trace = vec![0u8; (n + 1) * (m + 1)];
    let (mut hp, mut ep, mut fp) = (vec![neg; m + 1], vec![0f64; m + 1], vec![neg; m + 1]);
    hp[0] = 0.0;
    ep[0] = neg;
    // Row 0: free leading gaps in a (E runs along it).
    for j in 1..=m {
        trace[j] = 1 | 4; // E, extended along the free row
    }
    let (mut best, mut end) = (neg, (0usize, 0usize, 0u8));
    let consider = |v: f64, st: u8, i: usize, j: usize, best: &mut f64, end: &mut (usize, usize, u8)| {
        if v > *best {
            *best = v;
            *end = (i, j, st);
        }
    };
    for j in 0..=m {
        let (v, st) = max3(hp[j], ep[j], fp[j]);
        if j == m {
            consider(v, st, 0, j, &mut best, &mut end);
        }
    }
    let mut told = *done;
    for i in 1..=n {
        let (mut hc, mut ec, mut fc) = (vec![neg; m + 1], vec![neg; m + 1], vec![neg; m + 1]);
        // Column 0: free leading gaps in b (F runs down it).
        fc[0] = 0.0;
        trace[i * (m + 1)] = 2 | 16;
        for j in 1..=m {
            let (pv, pst) = max3(hp[j - 1], ep[j - 1], fp[j - 1]);
            let cell: f64 = pb[j - 1].iter().map(|&(y, f)| va[i - 1][y] * f).sum();
            let h = pv + cell;
            let mut t = pst;
            // E: b's column j against a gap (from the left).
            let (ov, ofrom) = if fc[j - 1] > hc[j - 1] { (fc[j - 1], 1u8) } else { (hc[j - 1], 0u8) };
            let e = if ov - open >= ec[j - 1] - ext {
                t |= 4 | (ofrom << 3);
                ov - open
            } else {
                ec[j - 1] - ext
            };
            // F: a's column i against a gap (from above).
            let (ov, ofrom) = if ep[j] > hp[j] { (ep[j], 1u8) } else { (hp[j], 0u8) };
            let f = if ov - open >= fp[j] - ext {
                t |= 16 | (ofrom << 5);
                ov - open
            } else {
                fp[j] - ext
            };
            hc[j] = h;
            ec[j] = e;
            fc[j] = f;
            trace[i * (m + 1) + j] = t;
        }
        *done += m;
        if *done - told >= 1 << 20 {
            crate::progress(*done, total);
            told = *done;
        }
        // Ends: the last column (b used up) on every row, the last row (a used up) below.
        let (v, st) = max3(hc[m], ec[m], fc[m]);
        consider(v, st, i, m, &mut best, &mut end);
        if i == n {
            for j in 0..m {
                let (v, st) = max3(hc[j], ec[j], fc[j]);
                consider(v, st, n, j, &mut best, &mut end);
            }
        }
        (hp, ep, fp) = (hc, ec, fc);
    }
    // Trailing free gaps after the end cell, then the traceback to (0, 0).
    let (ie, je, mut state) = end;
    let mut ops: Vec<u8> = Vec::with_capacity(n + m);
    ops.extend(std::iter::repeat_n(1u8, m - je));
    ops.extend(std::iter::repeat_n(2u8, n - ie));
    let (mut i, mut j) = (ie, je);
    while i > 0 || j > 0 {
        if i == 0 {
            ops.push(1);
            j -= 1;
            continue;
        }
        if j == 0 {
            ops.push(2);
            i -= 1;
            continue;
        }
        let t = trace[i * (m + 1) + j];
        match state {
            0 => {
                ops.push(0);
                state = t & 3;
                i -= 1;
                j -= 1;
            }
            1 => {
                ops.push(1);
                if t & 4 != 0 {
                    state = if t & 8 != 0 { 2 } else { 0 };
                }
                j -= 1;
            }
            _ => {
                ops.push(2);
                if t & 16 != 0 {
                    state = if t & 32 != 0 { 1 } else { 0 };
                }
                i -= 1;
            }
        }
    }
    ops.reverse();
    ops
}

fn max3(h: f64, e: f64, f: f64) -> (f64, u8) {
    let (mut v, mut s) = (h, 0u8);
    if e > v {
        (v, s) = (e, 1);
    }
    if f > v {
        (v, s) = (f, 2);
    }
    (v, s)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rows(v: &[&str], protein: bool) -> Vec<String> {
        msa(&v.iter().map(|s| s.as_bytes().to_vec()).collect::<Vec<_>>(), protein).unwrap().into_iter().map(|r| String::from_utf8(r).unwrap()).collect()
    }

    #[test]
    fn identical_rows_stay_as_they_are() {
        assert_eq!(rows(&["ACGTACGT", "ACGTACGT", "ACGTACGT"], false), vec!["ACGTACGT"; 3]);
    }

    #[test]
    fn one_insertion_goes_where_it_scores_best() {
        // ACGT-TGCA against ACGTATGCA: 8 matches (16) minus one gap of 1 (5 + 2) = 9. A gap after the second T leaves T
        // against A (8 − 7 − 3 − 3 + 2…), lower; without a gap at most 4 matches line up. Worked by hand: unique.
        assert_eq!(rows(&["ACGTTGCA", "ACGTATGCA"], false), vec!["ACGT-TGCA", "ACGTATGCA"]);
    }

    #[test]
    fn the_closest_pair_joins_first_and_the_third_aligns_to_its_profile() {
        // Rows 0 and 1 are identical (distance 0): the guide tree joins them first, then row 2, which has T inserted
        // between A and C (neither neighbour a T, so the gap has one best place).
        let seqs = ["ACGTACGTACGTAC", "ACGTACGTACGTAC", "ACGTATCGTACGTAC"];
        let v: Vec<Vec<u8>> = seqs.iter().map(|s| s.as_bytes().to_vec()).collect();
        assert_eq!(guide_tree(&v, false)[0], (0, 1));
        assert_eq!(rows(&seqs, false), vec!["ACGTA-CGTACGTAC", "ACGTA-CGTACGTAC", "ACGTATCGTACGTAC"]);
    }

    #[test]
    fn protein_scoring_uses_blosum62() {
        // One insertion of W in a protein: MKTAYIAK against MKTAWYIAK.
        assert_eq!(rows(&["MKTAYIAK", "MKTAWYIAK"], true), vec!["MKTA-YIAK", "MKTAWYIAK"]);
    }

    #[test]
    fn end_gaps_are_free() {
        // A truncated copy lines up with its stretch instead of being spread with internal gaps.
        assert_eq!(rows(&["AAACCCGGGTTT", "CCCGGG"], false), vec!["AAACCCGGGTTT", "---CCCGGG---"]);
    }

    #[test]
    fn static_bounds_refuse_before_any_matrix() {
        ALLOCATIONS.with(|c| c.set(0));
        let one = vec![b"ACGT".to_vec()];
        assert!(msa(&one, false).unwrap_err().contains("2 to 200"));
        let long = vec![vec![b'A'; 10_001], b"ACGT".to_vec()];
        assert!(msa(&long, false).unwrap_err().contains("1 to 10000"));
        let bad = vec![b"ACGX".to_vec(), b"ACGT".to_vec()];
        assert!(msa(&bad, false).unwrap_err().contains("not a DNA code at position 4"));
        let tight = Limits { work: 100.0, ..LIMITS };
        assert!(msa_with(&[b"ACGTACGT".to_vec(), b"ACGTACGT".to_vec()], false, &tight).unwrap_err().contains("steps"));
        assert_eq!(ALLOCATIONS.with(|c| c.get()), 0);
    }

    #[test]
    fn a_result_dependent_bound_refuses_before_that_join_s_matrix() {
        // Three sequences of 9 letters: the first join (10 × 10 = 100 cells) widens its profile to 10 columns, so the
        // second join needs 11 × 10 = 110 cells where the estimate (each cluster as wide as its longest member) said
        // 10 × 10 = 100. With 105 cells a join, the estimate passes and the second join is refused after the first
        // matrix and before its own.
        ALLOCATIONS.with(|c| c.set(0));
        let lim = Limits { join_cells: 105, ..LIMITS };
        let seqs = vec![b"ACGTTGCAT".to_vec(), b"ACGTATGCA".to_vec(), b"ACGTAGCAT".to_vec()];
        let err = msa_with(&seqs, false, &lim).unwrap_err();
        assert!(err.contains("cells"), "{err}");
        assert_eq!(ALLOCATIONS.with(|c| c.get()), 1);
    }
}
