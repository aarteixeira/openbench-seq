//! Translation with NCBI genetic codes, and ATG-to-stop ORFs on both strands, across the origin of circles.

use crate::alphabet::{mask, revcomp};

/// NCBI translation tables (amino acids in TCAG codon order): 1 standard, 2 vertebrate mitochondrial,
/// 4 mold/protozoan mitochondrial and Mycoplasma, 11 bacterial (same amino acids as 1).
pub fn table(id: u32) -> Option<&'static [u8; 64]> {
    match id {
        1 | 11 => Some(b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
        2 => Some(b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG"),
        4 => Some(b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
        _ => None,
    }
}

// Index of a base in TCAG order, by mask bit (A=1, C=2, G=4, T=8).
const fn tcag(bit: u8) -> usize {
    match bit {
        8 => 0,
        2 => 1,
        1 => 2,
        _ => 3,
    }
}

/// One codon: the amino acid every reading of its ambiguity codes agrees on, else `X`.
pub fn codon(code: &[u8; 64], c: &[u8]) -> u8 {
    let [a, b, d] = [mask(c[0]), mask(c[1]), mask(c[2])];
    if a == 0 || b == 0 || d == 0 {
        return b'X';
    }
    let mut out = 0u8;
    for x in [1u8, 2, 4, 8].into_iter().filter(|x| a & x != 0) {
        for y in [1u8, 2, 4, 8].into_iter().filter(|y| b & y != 0) {
            for z in [1u8, 2, 4, 8].into_iter().filter(|z| d & z != 0) {
                let aa = code[16 * tcag(x) + 4 * tcag(y) + tcag(z)];
                if out != 0 && out != aa {
                    return b'X';
                }
                out = aa;
            }
        }
    }
    out
}

/// Translate whole codons from the start of `dna`; stops are `*`.
pub fn translate(code: &[u8; 64], dna: &[u8]) -> Vec<u8> {
    dna.chunks_exact(3).map(|c| codon(code, c)).collect()
}

#[derive(Debug, PartialEq, Eq)]
pub struct Orf {
    /// 0-based start on the top strand.
    pub start: usize,
    /// 0-based exclusive end on the top strand; `end < start` wraps the origin (circular only).
    pub end: usize,
    pub strand: i8,
    /// Amino acids, the stop excluded.
    pub aa: usize,
}

/// ATG..stop ORFs of at least `min_aa` amino acids (stop excluded) on both strands. On a circle an ORF may cross
/// the origin; for each stop only the longest ORF (its most upstream in-frame ATG) is reported.
pub fn orfs(code: &[u8; 64], seq: &[u8], circular: bool, min_aa: usize) -> Vec<Orf> {
    let n = seq.len();
    let mut out = Vec::new();
    for strand in [1i8, -1] {
        let s = if strand == 1 { seq.to_ascii_uppercase() } else { revcomp(seq) };
        for o in strand_orfs(code, &s, circular, min_aa) {
            // An ORF [a, e) on the reverse complement is [n - e, n - a) on the top strand (wrapping stays wrapping).
            let (start, end) = if strand == 1 { (o.0, o.1) } else { ((n - o.1) % n, n - o.0) };
            out.push(Orf { start, end, strand, aa: o.2 });
        }
    }
    out.sort_by_key(|o| (o.start, o.strand));
    out
}

// (start, end, amino acids) on one strand's text; end is in 1..=n, and end < start when the ORF wraps.
fn strand_orfs(code: &[u8; 64], s: &[u8], circular: bool, min_aa: usize) -> Vec<(usize, usize, usize)> {
    let n = s.len();
    if n < 6 {
        return vec![];
    }
    // A circle is scanned twice over, so every stop sees up to a whole turn of upstream context in its own frame.
    let text: Vec<u8> = if circular { [s, s].concat() } else { s.to_vec() };
    // By stop position mod n: (start, aa). Sparse, so a long sequence costs memory per ORF, not per base.
    let mut best: std::collections::HashMap<usize, (usize, usize)> = std::collections::HashMap::new();
    for frame in 0..3 {
        // Starts since the last stop in this frame, oldest first; each expires on its own once more than a turn away.
        let mut starts: std::collections::VecDeque<usize> = std::collections::VecDeque::new();
        let mut i = frame;
        while i + 3 <= text.len() {
            let c = &text[i..i + 3];
            while starts.front().is_some_and(|&a| i + 3 - a > n) {
                starts.pop_front();
            }
            if codon(code, c) == b'*' {
                // The oldest start within one turn gives the longest ORF ending here; on a circle the stop may be on the
                // second turn (its start then maps back to the first).
                if let Some(&a) = starts.front() {
                    let len = (i - a) / 3;
                    if len >= min_aa {
                        let slot = best.entry(i % n).or_insert((a % n, len));
                        if len > slot.1 {
                            *slot = (a % n, len);
                        }
                    }
                }
                starts.clear();
            } else if c == b"ATG" || c == b"AUG" {
                starts.push_back(i);
            }
            i += 3;
        }
    }
    let mut out: Vec<_> = best
        .into_iter()
        .map(|(stop, (a, len))| {
            let e = stop + 3;
            (a, if e > n { e - n } else { e }, len)
        })
        .collect();
    out.sort_unstable();
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn std() -> &'static [u8; 64] {
        table(1).unwrap()
    }

    #[test]
    fn standard_and_mitochondrial_codes() {
        assert_eq!(translate(std(), b"ATGTGAAGAATATAA"), b"M*RI*".to_vec());
        assert_eq!(translate(table(2).unwrap(), b"ATGTGAAGAATATAA"), b"MW*M*".to_vec());
        assert_eq!(translate(table(4).unwrap(), b"TGA"), b"W".to_vec());
        assert!(table(3).is_none());
    }

    #[test]
    fn ambiguity_translates_only_when_unanimous() {
        // GCN is always Ala; YTR is Leu (CTA CTG TTA TTG); RAY is Asn or Asp.
        assert_eq!(translate(std(), b"GCNYTRRAYNNN"), b"ALXX".to_vec());
        assert_eq!(translate(std(), b"ATGC"), b"M".to_vec()); // a trailing partial codon is ignored
    }

    #[test]
    fn linear_orfs_take_the_first_atg_per_stop() {
        // ATG AAA ATG CCC TAA: one ORF of 4 aa from the first ATG.
        let o = orfs(std(), b"ATGAAAATGCCCTAA", false, 1);
        assert_eq!(o, vec![Orf { start: 0, end: 15, strand: 1, aa: 4 }]);
        assert!(orfs(std(), b"ATGAAAATGCCCTAA", false, 5).is_empty());
    }

    #[test]
    fn reverse_strand_coordinates() {
        // Reverse complement of ATGAAATAG is CTATTTCAT; padded so the ORF sits at 3..12.
        let o = orfs(std(), b"GGGCTATTTCATGGG", false, 2);
        assert_eq!(o, vec![Orf { start: 3, end: 12, strand: -1, aa: 2 }]);
    }

    #[test]
    fn circular_orf_across_the_origin() {
        // Circle "AAATAGCC" + "CCATGAAA": ATG at 10 reads ATG AAA AAA TAG across the origin, stop ending at 6.
        let seq = b"AAATAGCCCCATGAAA";
        let o = orfs(std(), seq, true, 3);
        assert_eq!(o, vec![Orf { start: 10, end: 6, strand: 1, aa: 3 }]);
        assert!(orfs(std(), seq, false, 3).is_empty()); // linear: no stop after the ATG
    }

    #[test]
    fn circular_reverse_orf_across_the_origin() {
        // Reverse complement of the circle above: the same ORF on the minus strand, wrapping the origin.
        let seq = revcomp(b"AAATAGCCCCATGAAA"); // TTTCATGGGGCTATTT
        let o = orfs(std(), &seq, true, 3);
        assert_eq!(o, vec![Orf { start: 10, end: 6, strand: -1, aa: 3 }]);
    }

    #[test]
    fn rna_start_codon() {
        assert_eq!(orfs(std(), b"AUGAAAUAA", false, 1), vec![Orf { start: 0, end: 9, strand: 1, aa: 2 }]);
    }

    #[test]
    fn a_later_start_survives_an_earlier_one_that_expires() {
        // On this circle the ATG at 0 has no stop within a turn in its frame; the ATG at 9 reads ATG CAT GCC TAA
        // across the origin.
        assert_eq!(orfs(std(), b"ATGCCTAACATGC", true, 3), vec![Orf { start: 9, end: 8, strand: 1, aa: 3 }]);
        assert_eq!(orfs(std(), &revcomp(b"ATGCCTAACATGC"), true, 3), vec![Orf { start: 5, end: 4, strand: -1, aa: 3 }]);
        // From 6: ATG AAA TGA.
        assert_eq!(orfs(std(), b"ATGAAAATGAA", true, 2), vec![Orf { start: 6, end: 4, strand: 1, aa: 2 }]);
    }

    // Independent reference: from every ATG on the circle, read codons for at most one turn to the first stop; for each
    // stop keep the longest. Rotations of the circle must give the same ORFs, shifted.
    fn brute(seq: &[u8], min_aa: usize) -> Vec<(usize, usize, usize)> {
        let n = seq.len();
        let at = |p: usize| seq[p % n];
        let mut best: std::collections::BTreeMap<usize, (usize, usize)> = Default::default();
        for s in 0..n {
            if [at(s), at(s + 1), at(s + 2)] != *b"ATG" {
                continue;
            }
            let mut k = 0;
            while 3 * k + 3 <= n {
                let c = [at(s + 3 * k), at(s + 3 * k + 1), at(s + 3 * k + 2)];
                if codon(std(), &c) == b'*' {
                    if k >= min_aa {
                        let stop = (s + 3 * k) % n;
                        if best.get(&stop).is_none_or(|&(_, l)| k > l) {
                            best.insert(stop, (s, k));
                        }
                    }
                    break;
                }
                k += 1;
            }
        }
        let mut v: Vec<_> = best.into_iter().map(|(stop, (s, k))| (s, (stop + 3 - 1) % n + 1, k)).collect();
        v.sort_unstable();
        v
    }

    #[test]
    fn circular_orfs_match_an_independent_scan() {
        let mut x = 99u64;
        for _ in 0..3000 {
            x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            let n = 6 + (x >> 40) as usize % 40;
            // A, T and G rich so ATG and stops are frequent.
            let seq: Vec<u8> = (0..n).map(|i| b"ATGATGAC"[((x >> (i % 50)) as usize + i * 7) % 8]).collect();
            for min in [1, 2] {
                let got: Vec<_> = strand_orfs(std(), &seq, true, min);
                assert_eq!(got, brute(&seq, min), "{}", String::from_utf8_lossy(&seq));
            }
        }
    }

    #[test]
    fn circle_without_stop_has_no_orf() {
        assert!(orfs(std(), b"ATGAAAAAAAAA", true, 1).is_empty());
    }
}
