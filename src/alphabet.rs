//! Nucleotide codes as 4-bit masks (A=1, C=2, G=4, T/U=8) and their complements.

/// Mask of an IUPAC nucleotide code; 0 for anything else.
pub const fn mask(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
        b'A' => 1,
        b'C' => 2,
        b'G' => 4,
        b'T' | b'U' => 8,
        b'R' => 1 | 4,
        b'Y' => 2 | 8,
        b'S' => 2 | 4,
        b'W' => 1 | 8,
        b'K' => 4 | 8,
        b'M' => 1 | 2,
        b'B' => 2 | 4 | 8,
        b'D' => 1 | 4 | 8,
        b'H' => 1 | 2 | 8,
        b'V' => 1 | 2 | 4,
        b'N' => 15,
        _ => 0,
    }
}

/// Complement of a mask: A<->T, C<->G, bit by bit.
pub const fn complement(m: u8) -> u8 {
    (m & 1) << 3 | (m & 2) << 1 | (m & 4) >> 1 | (m & 8) >> 3
}

/// A sequence position `s` satisfies a pattern code `q` when every base `s` may be is one `q` allows.
/// `N` in the pattern matches anything; `N` in the sequence matches only `N`. Unknown characters never match.
pub const fn satisfies(s: u8, q: u8) -> bool {
    s != 0 && s & !q == 0
}

const CODES: &[u8; 16] = b"-ACMGRSVTWYHKDBN";

/// IUPAC letter of a mask.
pub fn letter(m: u8) -> u8 {
    CODES[(m & 15) as usize]
}

/// Reverse complement of IUPAC nucleotide text; unknown characters become `N`.
pub fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter().rev().map(|&b| match mask(b) { 0 => b'N', m => letter(complement(m)) }).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn complements_pair_bases_and_codes() {
        assert_eq!(revcomp(b"ACGTRYKMBVDHN"), b"NDHBVKMRYACGT".to_vec());
        assert_eq!(revcomp(b"acgu"), b"ACGT".to_vec());
    }

    #[test]
    fn subset_policy() {
        assert!(satisfies(mask(b'A'), mask(b'R')));
        assert!(!satisfies(mask(b'R'), mask(b'A'))); // the sequence may be G
        assert!(satisfies(mask(b'N'), mask(b'N')));
        assert!(!satisfies(mask(b'N'), mask(b'A')));
        assert!(!satisfies(mask(b'X'), mask(b'N')));
    }
}
