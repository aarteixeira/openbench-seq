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
/// `Z` E or Q, `J` I or L; an ambiguous letter in the sequence matches only the same letter.
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

// Windows of `text` (starting before `limit`) within k substitutions of the pattern behind `t`.
fn scan(text: &[u8], limit: usize, m: usize, k: usize, t: &[u64; 256], strand: i8, out: &mut Vec<Hit>) {
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
                out.push(Hit { pos: start, strand, substitutions: d as u8 });
            }
        }
    }
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
        return Ok(vec![]);
    }
    let text: Vec<u8> = if circular { [seq, &seq[..m - 1]].concat() } else { seq.to_vec() };
    let mut out = Vec::new();
    if protein {
        scan(&text, n, m, k, &table(pattern, residue), 1, &mut out);
        return Ok(out);
    }
    let fwd: Vec<u8> = pattern.iter().map(|&b| mask(b)).collect();
    if fwd.contains(&0) {
        return Err(SearchError::InvalidPattern);
    }
    let rev: Vec<u8> = fwd.iter().rev().map(|&q| complement(q)).collect();
    let by_mask = |pat: &[u8]| table(pat, |b, q| satisfies(mask(b), q));
    scan(&text, n, m, k, &by_mask(&fwd), 1, &mut out);
    if rev != fwd {
        scan(&text, n, m, k, &by_mask(&rev), -1, &mut out);
        out.sort_by_key(|h| (h.pos, -h.strand));
    }
    Ok(out)
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
        // M and R keep their residue meaning (no nucleotide reading of R as A/G).
        assert_eq!(hits("AGRM", "R", 0, true, false), vec![(2, 1, 0)]);
        assert_eq!(hits("MKAAR", "MKXXK", 1, true, false), vec![(0, 1, 1)]);
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
                .filter_map(|i| {
                    let d = p.iter().enumerate().filter(|(j, q)| !satisfies(mask(seq[i + j]), mask(**q))).count();
                    (d <= 1).then_some(i)
                })
                .count()
        };
        let got = search(&seq, pat, 1, false, false).unwrap();
        let rc: Vec<u8> = crate::alphabet::revcomp(pat);
        assert_eq!(got.len(), direct(pat) + if rc == pat { 0 } else { direct(&rc) });
    }
}
