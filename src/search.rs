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

// Each window within `k` substitutions, in the order `scan` finds them, passed to `f` with its start; `f` returns false
// to stop the scan.
fn scan_each(text: &[u8], limit: usize, m: usize, k: usize, t: &[u64; 256], f: &mut dyn FnMut(usize) -> bool) {
    let mut r = [0u64; MAX_SUBSTITUTIONS + 1];
    let top = 1u64 << (m - 1);
    for (j, &c) in text.iter().enumerate() {
        let b = t[c as usize];
        for d in (0..=k).rev() {
            let shifted = (r[d] << 1) | 1;
            r[d] = (shifted & b) | if d > 0 { (r[d - 1] << 1) | 1 } else { 0 };
        }
        if j + 1 >= m {
            let start = j + 1 - m;
            if start >= limit || (r[k] & top != 0 && !f(start)) {
                return;
            }
        }
    }
}

/// Substitutions a primer binding site may hold, and the 3′ bases that must pair exactly.
pub const BIND_SUBSTITUTIONS: usize = 2;
pub const BIND_CLAMP: usize = 5;

/// A site where a primer binds: its window's start on the top strand, the strand the oligo's sequence reads on (1:
/// the oligo equals the top strand there; −1: it equals the bottom strand, so it pairs with the top), and the oligo
/// positions (from its 5′ end, bit i) that do not pair.
#[derive(Debug, PartialEq, Eq, Clone, Copy)]
pub struct Bound {
    pub pos: usize,
    pub strand: i8,
    pub mismatches: u64,
}

/// Where `oligo` (5′→3′) binds `seq`: windows within 2 substitutions on either strand (across the origin of a circle)
/// whose oligo's last 5 bases pair exactly. Every raw hit is checked as it is found and only qualifying ones are
/// counted, so the count never comes from a truncated list. Counting stops at `stop_at` (2 decides "binds once");
/// at most `keep` sites are returned, in position order.
pub fn bind(seq: &[u8], oligo: &[u8], circular: bool, keep: usize, stop_at: usize) -> Result<(usize, Vec<Bound>), SearchError> {
    let m = oligo.len();
    if m == 0 {
        return Err(SearchError::EmptyPattern);
    }
    if m > MAX_PATTERN {
        return Err(SearchError::PatternTooLong);
    }
    if m <= BIND_CLAMP + BIND_SUBSTITUTIONS {
        return Err(SearchError::TooManySubstitutions);
    }
    let fwd: Vec<u8> = oligo.iter().map(|&b| mask(b)).collect();
    if fwd.contains(&0) {
        return Err(SearchError::InvalidPattern);
    }
    let n = seq.len();
    if n < m {
        return Ok((0, vec![]));
    }
    let text: Vec<u8> = if circular { [seq, &seq[..m - 1]].concat() } else { seq.to_vec() };
    let rev: Vec<u8> = fwd.iter().rev().map(|&q| complement(q)).collect();
    let by_mask = |pat: &[u8]| table(pat, |b, q| satisfies(mask(b), q));
    let (mut count, mut out) = (0usize, Vec::new());
    let strands: &[(i8, &[u8])] = if rev == fwd { &[(1, &fwd)] } else { &[(1, &fwd), (-1, &rev)] };
    for &(strand, pat) in strands {
        if count >= stop_at {
            break;
        }
        let t = by_mask(pat);
        scan_each(&text, n, m, BIND_SUBSTITUTIONS, &t, &mut |start| {
            // The oligo's base i against the template base it pairs with (literal: an N in the template does not pair).
            let mut mis = 0u64;
            for i in 0..m {
                let s = if strand == 1 { mask(text[start + i]) } else { complement(mask(text[start + m - 1 - i])) };
                if s != fwd[i] {
                    mis |= 1 << i;
                }
            }
            let clamp = ((1u64 << BIND_CLAMP) - 1) << (m - BIND_CLAMP);
            if mis & clamp == 0 {
                if out.len() < keep {
                    out.push(Bound { pos: start, strand, mismatches: mis });
                }
                count += 1;
            }
            count < stop_at
        });
    }
    out.sort_by_key(|b| (b.pos, -b.strand));
    out.truncate(keep);
    Ok((count, out))
}

/// The best short-query hit in one record: its window's start, strand, ungapped score and identical letters.
#[derive(Debug, PartialEq, Eq, Clone, Copy)]
pub struct Best {
    pub pos: usize,
    pub strand: i8,
    pub score: i32,
    pub identical: usize,
}

/// The best-scoring window of `seq` within `k` substitutions of `query` (both strands for DNA, across the origin of a
/// circle) with at least `min_identical` identical letters (the same letter in both: a subset match between different
/// letters does not count). Every raw hit is scored as it is found and only the best is kept (ties: more identical
/// letters, then the lower position, then the top strand; an ambiguity code scores as a match under the subset rule,
/// so a window of ambiguous letters can tie an exact copy), so a record full of near-matches costs time, not memory.
pub fn best_hit(seq: &[u8], query: &[u8], k: usize, protein: bool, circular: bool, min_identical: usize) -> Result<Option<Best>, SearchError> {
    let m = query.len();
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
        return Ok(None);
    }
    let text: Vec<u8> = if circular { [seq, &seq[..m - 1]].concat() } else { seq.to_vec() };
    let sc = if protein { crate::align::Scoring::PROTEIN } else { crate::align::Scoring::DNA };
    let mut best: Option<Best> = None;
    let mut consider = |start: usize, strand: i8| {
        let (mut score, mut identical) = (0i32, 0usize);
        for i in 0..m {
            let s = if strand == 1 { text[start + i] } else { crate::alphabet::letter(complement(mask(text[start + m - 1 - i]))) };
            score += sc.pair(query[i], s);
            if query[i].eq_ignore_ascii_case(&s) {
                identical += 1;
            }
        }
        let better = match best {
            None => true,
            Some(b) => (score, identical, std::cmp::Reverse((start, -strand))) > (b.score, b.identical, std::cmp::Reverse((b.pos, -b.strand))),
        };
        if identical >= min_identical && better {
            best = Some(Best { pos: start, strand, score, identical });
        }
    };
    if protein {
        scan_each(&text, n, m, k, &table(query, residue), &mut |start| {
            consider(start, 1);
            true
        });
        return Ok(best);
    }
    let fwd: Vec<u8> = query.iter().map(|&b| mask(b)).collect();
    if fwd.contains(&0) {
        return Err(SearchError::InvalidPattern);
    }
    let rev: Vec<u8> = fwd.iter().rev().map(|&q| complement(q)).collect();
    let by_mask = |pat: &[u8]| table(pat, |b, q| satisfies(mask(b), q));
    scan_each(&text, n, m, k, &by_mask(&fwd), &mut |start| {
        consider(start, 1);
        true
    });
    if rev != fwd {
        scan_each(&text, n, m, k, &by_mask(&rev), &mut |start| {
            consider(start, -1);
            true
        });
    }
    Ok(best)
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

    #[test]
    fn bind_counts_only_sites_whose_3_prime_end_pairs() {
        let oligo = b"ACGTTGCAAGCTTCGA"; // 16 bases, 3' end CTTCGA…
        // One exact site, and one with the oligo's last base mismatched (2 substitutions allowed overall).
        let mut seq = vec![b'T'; 30];
        seq.extend_from_slice(oligo);
        seq.extend(vec![b'T'; 30]);
        let mut decoy = oligo.to_vec();
        decoy[15] = b'C';
        seq.extend_from_slice(&decoy);
        seq.extend(vec![b'T'; 30]);
        let (count, sites) = bind(&seq, oligo, false, 10, usize::MAX).unwrap();
        assert_eq!((count, sites.len(), sites[0].pos, sites[0].strand, sites[0].mismatches), (1, 1, 30, 1, 0));
        // An internal mismatch is allowed and reported by position (from the oligo's 5' end).
        let mut inner = seq.clone();
        inner[30 + 3] = b'A';
        let (count, sites) = bind(&inner, oligo, false, 10, usize::MAX).unwrap();
        assert_eq!((count, sites[0].mismatches), (1, 1 << 3));
    }

    #[test]
    fn bind_counts_past_more_than_ten_decoys() {
        // Twelve decoys mismatched in the 3' clamp, then one qualifying site: a count read from the first ten raw hits
        // would miss it; with stop_at 2 the count is 1 (binds once).
        let oligo = b"ACGTTGCAAGCTTCGA";
        let mut decoy = oligo.to_vec();
        decoy[14] = b'T';
        let mut seq = Vec::new();
        for _ in 0..12 {
            seq.extend_from_slice(&decoy);
            seq.extend_from_slice(b"TTTTT");
        }
        seq.extend_from_slice(oligo);
        let (count, sites) = bind(&seq, oligo, false, 10, 2).unwrap();
        assert_eq!((count, sites.len(), sites[0].pos), (1, 1, 12 * 21));
    }

    #[test]
    fn bind_reads_a_reverse_primer_on_the_other_strand_and_across_the_origin() {
        let oligo = b"ACGTTGCAAGCTTCGA";
        let rc = crate::alphabet::revcomp(oligo);
        let mut seq = vec![b'T'; 20];
        seq.extend_from_slice(&rc);
        let (count, sites) = bind(&seq, oligo, false, 10, usize::MAX).unwrap();
        assert_eq!((count, sites[0].pos, sites[0].strand), (1, 20, -1));
        // Split across the origin of a circle.
        let mut circle = oligo[8..].to_vec();
        circle.extend(vec![b'T'; 20]);
        circle.extend_from_slice(&oligo[..8]);
        assert_eq!(bind(&circle, oligo, false, 10, usize::MAX).unwrap().0, 0);
        assert_eq!(bind(&circle, oligo, true, 10, usize::MAX).unwrap().0, 1);
        assert_eq!(bind(&circle, &[b'A'; 65], true, 10, 2), Err(SearchError::PatternTooLong));
    }

    #[test]
    fn best_hit_keeps_the_best_scoring_eligible_window() {
        // The antibody motif WVAGIWED: DVAGIWED (one substitution) scores 36 with BLOSUM62, WVSGVWED (two) scores 47; the
        // later, better-scoring one is kept (independent values from the planning review).
        let seq = b"MMMMDVAGIWEDMMMMMMWVSGVWEDMMMM";
        let b = best_hit(seq, b"WVAGIWED", 3, true, false, 3).unwrap().unwrap();
        assert_eq!((b.pos, b.score), (18, 47));
        let early = best_hit(b"MMMMDVAGIWEDMMMM", b"WVAGIWED", 3, true, false, 3).unwrap().unwrap();
        assert_eq!(early.score, 36);
    }

    #[test]
    fn best_hit_counts_literal_identity_before_choosing() {
        // Query N9A9: an earlier G9A9 satisfies N with no substitution and scores the same, but only 9 of 18 letters are
        // identical; with 11 required (60 %) the later exact N9A9 is the hit.
        let q = b"NNNNNNNNNAAAAAAAAA";
        let mut seq = b"GGGGGGGGGAAAAAAAAA".to_vec();
        seq.extend_from_slice(b"CCCCC");
        seq.extend_from_slice(q);
        let b = best_hit(&seq, q, 3, false, false, 11).unwrap().unwrap();
        assert_eq!((b.pos, b.identical), (23, 18));
    }

    #[test]
    fn best_hit_scans_a_record_of_near_matches_without_holding_them() {
        // 20,000 near-matches (one substitution each) before one exact copy: the exact copy wins, and nothing but the
        // best is held.
        let q = b"ACGTTGCAAGCTTCGAGG";
        let mut near = q.to_vec();
        near[5] = b'A';
        let mut seq = Vec::new();
        for _ in 0..20_000 {
            seq.extend_from_slice(&near);
        }
        seq.extend_from_slice(q);
        let b = best_hit(&seq, q, 3, false, false, 11).unwrap().unwrap();
        assert_eq!((b.pos, b.identical), (20_000 * 18, 18));
    }
}
