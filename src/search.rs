//! Motif search with up to k substitutions (no insertions or deletions), bit-parallel (shift-and with one state
//! word per allowed substitution count). Nucleotides follow the IUPAC subset policy in `alphabet`; proteins follow
//! their own ambiguity rules and are searched on one strand.

use crate::alphabet::{complement, mask, satisfies};

pub const MAX_PATTERN: usize = 64;
pub const MAX_SUBSTITUTIONS: usize = 3;

#[derive(Debug, PartialEq, Eq, Clone, Copy)]
pub struct Hit {
    /// 0-based start of the matched window on the top strand.
    pub pos: usize,
    pub strand: i8,
    pub substitutions: u8,
}

/// A protein residue `s` satisfies a query letter `q`: `X` in the query matches any residue; `B` matches D or N,
/// `Z` E or Q, `J` I or L; an ambiguous letter in the sequence matches only the same letter or a query `X`.
fn residue(s: u8, q: u8) -> bool {
    let (s, q) = (s.to_ascii_uppercase(), q.to_ascii_uppercase());
    s == q
        || q == b'X'
        || match q {
            b'B' => matches!(s, b'D' | b'N'),
            b'Z' => matches!(s, b'E' | b'Q'),
            b'J' => matches!(s, b'I' | b'L'),
            _ => false,
        }
}

// For every byte value, the query positions it satisfies (bit i = position i).
fn table(pattern: &[u8], ok: impl Fn(u8, u8) -> bool) -> Box<[u64; 256]> {
    let mut t = Box::new([0u64; 256]);
    for (b, slot) in t.iter_mut().enumerate() {
        for (i, &q) in pattern.iter().enumerate() {
            if ok(b as u8, q) {
                *slot |= 1 << i;
            }
        }
    }
    t
}

// Windows of `text` (starting before `limit`) within k substitutions of the pattern behind `t`, in position order:
// all are counted (each also in `bins`, equal parts of [0, limit)), the first `keep` are kept. Returns the count and
// the kept hits.
#[allow(clippy::too_many_arguments)] // one call per strand; a struct would only rename them
fn scan(text: &[u8], limit: usize, m: usize, k: usize, t: &[u64; 256], strand: i8, keep: usize, bins: &mut [u32]) -> (usize, Vec<Hit>) {
    let (mut count, mut out) = (0, Vec::new());
    let mut r = [0u64; MAX_SUBSTITUTIONS + 1];
    let top = 1u64 << (m - 1);
    for (j, &c) in text.iter().enumerate() {
        let b = t[c as usize];
        // Update from the highest count down so r[d - 1] is still the previous column's.
        for d in (0..=k).rev() {
            let shifted = (r[d] << 1) | 1;
            r[d] = (shifted & b) | if d > 0 { (r[d - 1] << 1) | 1 } else { 0 };
        }
        if j + 1 >= m {
            let start = j + 1 - m;
            if start >= limit {
                break;
            }
            if let Some(d) = (0..=k).find(|&d| r[d] & top != 0) {
                if count < keep {
                    out.push(Hit { pos: start, strand, substitutions: d as u8 });
                }
                if !bins.is_empty() {
                    let b = (start as u64 * bins.len() as u64 / limit as u64) as usize;
                    bins[b] = bins[b].saturating_add(1);
                }
                count += 1;
            }
        }
    }
    (count, out)
}

#[derive(Debug, PartialEq, Eq)]
pub enum SearchError {
    EmptyPattern,
    PatternTooLong,
    TooManySubstitutions,
    InvalidPattern,
}

/// Every window within `k` substitutions of `pattern`. Nucleotide search covers both strands (a palindromic
/// pattern is reported once, on the top strand); protein search covers the given strand only. On a circle, windows
/// may cross the origin. Hits are ordered by position, then strand.
pub fn search(seq: &[u8], pattern: &[u8], k: usize, protein: bool, circular: bool) -> Result<Vec<Hit>, SearchError> {
    search_first(seq, pattern, k, protein, circular, usize::MAX).map(|(_, h)| h)
}

/// As `search`, but only the first `keep` hits (in that order) are held in memory; returns the total count with them.
pub fn search_first(seq: &[u8], pattern: &[u8], k: usize, protein: bool, circular: bool, keep: usize) -> Result<(usize, Vec<Hit>), SearchError> {
    search_binned(seq, pattern, k, protein, circular, keep, &mut [])
}

/// As `search_first`, also counting every hit (kept or not) into `bins`: bin b counts hits starting in
/// [b * n / bins, (b + 1) * n / bins), so a map can show where all of them are.
pub fn search_binned(seq: &[u8], pattern: &[u8], k: usize, protein: bool, circular: bool, keep: usize, bins: &mut [u32]) -> Result<(usize, Vec<Hit>), SearchError> {
    let m = pattern.len();
    if m == 0 {
        return Err(SearchError::EmptyPattern);
    }
    if m > MAX_PATTERN {
        return Err(SearchError::PatternTooLong);
    }
    if k > MAX_SUBSTITUTIONS || k >= m {
        return Err(SearchError::TooManySubstitutions);
    }
    let n = seq.len();
    if n < m {
        return Ok((0, vec![]));
    }
    let text: Vec<u8> = if circular { [seq, &seq[..m - 1]].concat() } else { seq.to_vec() };
    if protein {
        return Ok(scan(&text, n, m, k, &table(pattern, residue), 1, keep, bins));
    }
    let fwd: Vec<u8> = pattern.iter().map(|&b| mask(b)).collect();
    if fwd.contains(&0) {
        return Err(SearchError::InvalidPattern);
    }
    let rev: Vec<u8> = fwd.iter().rev().map(|&q| complement(q)).collect();
    let by_mask = |pat: &[u8]| table(pat, |b, q| satisfies(mask(b), q));
    let (mut total, mut out) = scan(&text, n, m, k, &by_mask(&fwd), 1, keep, bins);
    if rev != fwd {
        // Each strand's first `keep` hits contain the first `keep` of both together.
        let (count, bottom) = scan(&text, n, m, k, &by_mask(&rev), -1, keep, bins);
        total += count;
        out.extend(bottom);
        out.sort_by_key(|h| (h.pos, -h.strand));
        out.truncate(keep);
    }
    Ok((total, out))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn hits(seq: &str, pat: &str, k: usize, protein: bool, circular: bool) -> Vec<(usize, i8, u8)> {
        search(seq.as_bytes(), pat.as_bytes(), k, protein, circular)
            .unwrap()
            .into_iter()
            .map(|h| (h.pos, h.strand, h.substitutions))
            .collect()
    }

    #[test]
    fn exact_palindrome_reported_once() {
        assert_eq!(hits("AAGGATCCAAGGATCC", "GGATCC", 0, false, false), vec![(2, 1, 0), (10, 1, 0)]);
    }

    #[test]
    fn reverse_strand_hit() {
        // GAATTC is palindromic; GGTACC... use a non-palindrome: ATGC on the top, GCAT read on the bottom.
        assert_eq!(hits("TTGCATTT", "ATGC", 0, false, false), vec![(2, -1, 0)]);
    }

    #[test]
    fn substitution_accepted_indel_refused() {
        // One substitution: GGATCC vs GGTTCC.
        assert_eq!(hits("AAGGTTCCAA", "GGATCC", 1, false, false), vec![(2, 1, 1)]);
        // GGATCC with one base deleted (GGTCC): one deletion, but at least two substitutions in every window.
        assert!(hits("AAAAGGTCCAAAAA", "GGATCC", 1, false, false).is_empty());
        assert!(hits("TTTTGGTCCTTTTT", "GGATCC", 1, false, false).is_empty());
    }

    #[test]
    fn iupac_subset_policy() {
        // Pattern R (A or G): A on the top strand; each C is a G on the bottom strand.
        assert_eq!(hits("CAC", "R", 0, false, false), vec![(0, -1, 0), (1, 1, 0), (2, -1, 0)]);
        // A sequence R or N is not known to be A.
        assert!(hits("CRC", "A", 0, false, false).is_empty());
        assert_eq!(hits("ANA", "A", 0, false, false), vec![(0, 1, 0), (2, 1, 0)]);
        assert_eq!(hits("CNC", "N", 0, false, false).len(), 3);
    }

    #[test]
    fn circular_window_crosses_origin() {
        assert_eq!(hits("TCCAAAAAGGA", "GGATCC", 0, false, true), vec![(8, 1, 0)]);
        assert!(hits("TCCAAAAAGGA", "GGATCC", 0, false, false).is_empty());
    }

    #[test]
    fn protein_rules_and_one_strand() {
        assert_eq!(hits("MKDNMR", "B", 0, true, false), vec![(2, 1, 0), (3, 1, 0)]);
        assert_eq!(hits("MKXR", "X", 0, true, false), vec![(0, 1, 0), (1, 1, 0), (2, 1, 0), (3, 1, 0)]);
        assert!(hits("MKDR", "Z", 0, true, false).is_empty());
        // Query X matches every residue, ambiguous ones included; a sequence B, Z or J matches only itself (or X).
        assert_eq!(hits("BZJ", "X", 0, true, false).len(), 3);
        assert!(hits("BZJ", "D", 0, true, false).is_empty() && hits("BZJ", "E", 0, true, false).is_empty());
        // M and R keep their residue meaning (no nucleotide reading of R as A/G).
        assert_eq!(hits("AGRM", "R", 0, true, false), vec![(2, 1, 0)]);
        assert_eq!(hits("MKAAR", "MKXXK", 1, true, false), vec![(0, 1, 1)]);
    }

    #[test]
    fn first_hits_and_the_total_without_holding_the_rest() {
        let all = search(b"ACACACACAC", b"AC", 0, false, false).unwrap(); // AC on top, GT read as AC on the bottom: none
        let (total, first) = search_first(b"ACACACACAC", b"AC", 0, false, false, 2).unwrap();
        assert_eq!((total, first.as_slice()), (all.len(), &all[..2]));
        // Both strands: the first three by position, then strand.
        let both = search(b"GTACGTAC", b"AC", 0, false, false).unwrap();
        let (t, f) = search_first(b"GTACGTAC", b"AC", 0, false, false, 3).unwrap();
        assert_eq!((t, f.as_slice()), (both.len(), &both[..3]));
    }

    #[test]
    fn bins_count_every_hit_including_those_not_kept() {
        // 40 bases, "AC" at 0, 2, ..., 38 (20 hits, the bottom strand's GT none); 4 bins of 10 bases; keep 3.
        let seq = "AC".repeat(20);
        let mut bins = [0u32; 4];
        let (total, kept) = search_binned(seq.as_bytes(), b"AC", 0, false, false, 3, &mut bins).unwrap();
        assert_eq!((total, kept.len()), (20, 3));
        assert_eq!(bins, [5, 5, 5, 5]);
        // Both strands are counted: "GATT" on the top strand at 0, its reverse complement "AATC" at 8 and 16.
        let mut bins = [0u32; 3];
        let (total, _) = search_binned(b"GATTCCCCAATCCCCCAATC", b"GATT", 0, false, false, 0, &mut bins).unwrap();
        assert_eq!(total, 3);
        assert_eq!(bins, [1, 1, 1]);
    }

    #[test]
    fn limits() {
        assert_eq!(search(b"ACGT", b"", 0, false, false), Err(SearchError::EmptyPattern));
        assert_eq!(search(b"ACGT", &[b'A'; 65], 0, false, false), Err(SearchError::PatternTooLong));
        assert_eq!(search(b"ACGT", b"ACGT", 4, false, false), Err(SearchError::TooManySubstitutions));
        assert_eq!(search(b"ACGT", b"AC!", 0, false, false), Err(SearchError::InvalidPattern));
    }

    #[test]
    fn matches_a_direct_count_on_a_long_random_sequence() {
        // Independent check: a plain window-by-window count with the same subset rule.
        let mut x = 12345u64;
        let seq: Vec<u8> = (0..5000)
            .map(|_| {
                x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                b"ACGT"[(x >> 33) as usize % 4]
            })
            .collect();
        let pat = b"GANTC";
        let direct = |p: &[u8]| {
            (0..=seq.len() - p.len())
                .filter(|&i| p.iter().enumerate().filter(|(j, q)| !satisfies(mask(seq[i + j]), mask(**q))).count() <= 1)
                .count()
        };
        let got = search(&seq, pat, 1, false, false).unwrap();
        let rc: Vec<u8> = crate::alphabet::revcomp(pat);
        assert_eq!(got.len(), direct(pat) + if rc == pat { 0 } else { direct(&rc) });
    }
}
