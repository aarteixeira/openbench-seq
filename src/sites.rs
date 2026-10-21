//! Restriction sites and cut positions, with a flag for sites that Dam, Dcm or CpG methylation may block.
//!
//! Cut offsets follow the table the page used before (fcut: top-strand cut from the site's start; rcut: the
//! bottom-strand cut, so a site read on the bottom strand at `i` cuts the top strand at `i + len - rcut`).
//! A site is flagged when the enzyme is sensitive to a methylation and a base that methylation modifies, on either
//! strand, lies inside the recognition sequence. That is the usual "may be blocked by overlapping methylation"
//! reading; whether a given overlap blocks or only slows an enzyme is not modelled.

use crate::alphabet::{complement, mask, satisfies};

pub const DAM: u8 = 1;
pub const DCM: u8 = 2;
pub const CPG: u8 = 4;

#[derive(Debug, Clone)]
pub struct Enzyme {
    pub site: Vec<u8>,
    pub fcut: i64,
    pub rcut: i64,
    /// DAM | DCM | CPG: the methylations the enzyme is sensitive to.
    pub sensitive: u8,
}

#[derive(Debug, PartialEq, Eq)]
pub struct Cut {
    pub enzyme: usize,
    /// Top-strand cut position, 0..n (circular) or 1..n-1 (linear).
    pub pos: usize,
    /// The methylations that may block the site behind this cut.
    pub blocked: u8,
}

// Methylated top-strand positions per kind, as prefix sums for window queries.
struct Methylation {
    prefix: [Vec<u32>; 3],
}

impl Methylation {
    fn new(text: &[u8], n: usize, circular: bool) -> Self {
        // (motif, offsets of methylated bases on the top strand, counting the base paired to a methylated
        // bottom-strand base as that position). Dam: GATC, A on both strands -> 1 and 2. Dcm: CCWGG, inner C on both
        // strands -> 1 and 3. CpG: CG, C on both strands -> 0 and 1.
        let kinds: [(&[u8], &[usize]); 3] = [(b"GATC", &[1, 2]), (b"CCWGG", &[1, 3]), (b"CG", &[0, 1])];
        let len = text.len();
        let prefix = kinds.map(|(motif, offs)| {
            let mut marked = vec![false; len];
            let last = if circular { n } else { len.saturating_sub(motif.len() - 1) };
            for i in 0..last {
                if i + motif.len() <= len && motif.iter().enumerate().all(|(j, &q)| satisfies(mask(text[i + j]), mask(q))) {
                    for &o in offs {
                        marked[(i + o) % if circular { n } else { len }] = true;
                    }
                }
            }
            if circular {
                // The text repeats the start of the circle after n; mirror the marks there.
                for i in n..len {
                    marked[i] = marked[i - n];
                }
            }
            let mut p = vec![0u32; len + 1];
            for i in 0..len {
                p[i + 1] = p[i] + marked[i] as u32;
            }
            p
        });
        Methylation { prefix }
    }

    fn blocked(&self, start: usize, len: usize, sensitive: u8) -> u8 {
        let mut b = 0;
        for (k, bit) in [DAM, DCM, CPG].into_iter().enumerate() {
            let p = &self.prefix[k];
            if sensitive & bit != 0 && p[start + len] > p[start] {
                b |= bit;
            }
        }
        b
    }
}

/// Cut positions of every enzyme, sorted by enzyme then position; each position once per enzyme.
pub fn cuts(seq: &[u8], circular: bool, enzymes: &[Enzyme]) -> Vec<Cut> {
    let n = seq.len();
    let longest = enzymes.iter().map(|e| e.site.len()).max().unwrap_or(1).max(5);
    // Circular text repeats enough of the start for sites and methylation motifs across the origin.
    let text: Vec<u8> = if circular { [seq, &seq[..(longest - 1).min(n)]].concat() } else { seq.to_vec() };
    let meth = Methylation::new(&text, n, circular);
    let mut out = Vec::new();
    for (ei, e) in enzymes.iter().enumerate() {
        let l = e.site.len();
        if l == 0 || l > n {
            continue;
        }
        let fwd: Vec<u8> = e.site.iter().map(|&b| mask(b)).collect();
        let rev: Vec<u8> = fwd.iter().rev().map(|&q| complement(q)).collect();
        let mut found: Vec<(usize, u8)> = Vec::new();
        for i in 0..n {
            if i + l > text.len() {
                break;
            }
            let w = &text[i..i + l];
            for (pat, cut) in [(&fwd, i as i64 + e.fcut), (&rev, (i + l) as i64 - e.rcut)] {
                if !w.iter().zip(pat.iter()).all(|(&b, &q)| satisfies(mask(b), q)) {
                    continue;
                }
                let pos = if circular {
                    cut.rem_euclid(n as i64) as usize
                } else if cut > 0 && cut < n as i64 {
                    cut as usize
                } else {
                    continue;
                };
                found.push((pos, meth.blocked(i, l, e.sensitive)));
            }
        }
        found.sort();
        // One cut per position: blocked only if every site cutting there may be blocked.
        let mut k = 0;
        while k < found.len() {
            let pos = found[k].0;
            let mut b = 0xff;
            while k < found.len() && found[k].0 == pos {
                b &= found[k].1;
                k += 1;
            }
            out.push(Cut { enzyme: ei, pos, blocked: b });
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn enz(site: &str, fcut: i64, rcut: i64, sensitive: u8) -> Enzyme {
        Enzyme { site: site.as_bytes().to_vec(), fcut, rcut, sensitive }
    }

    fn pos(c: &[Cut]) -> Vec<(usize, u8)> {
        c.iter().map(|c| (c.pos, c.blocked)).collect()
    }

    #[test]
    fn ecori_linear_and_circular() {
        // EcoRI G^AATTC: fcut 1, rcut 5.
        let e = [enz("GAATTC", 1, 5, 0)];
        assert_eq!(pos(&cuts(b"AAGAATTCAA", false, &e)), vec![(3, 0)]);
        // Across the origin of a circle: GAA at the end (8..11), TTC at the start; the cut follows the G at 8.
        assert_eq!(pos(&cuts(b"TTCAAAAAGAA", true, &e)), vec![(9, 0)]);
        assert!(cuts(b"TTCAAAAAGAA", false, &e).is_empty());
    }

    #[test]
    fn non_palindromic_site_on_the_bottom_strand() {
        // BsaI GGTCTC(1/5), fcut 7, rcut 11: the top site at 2 cuts the top strand at 2+7 = 9; GAGACC at 12 is the
        // site on the bottom strand, whose top-strand cut lies 5 nt to its left: 12+6-11 = 7.
        let e = [enz("GGTCTC", 7, 11, 0)];
        let c = cuts(b"AAGGTCTCAAAAGAGACCAAAA", false, &e);
        assert_eq!(pos(&c), vec![(7, 0), (9, 0)]);
    }

    #[test]
    fn dam_overlap_blocks_xbai() {
        // XbaI T^CTAGA is Dam-sensitive when the site is followed by TC (TCTAGATC contains GATC).
        let e = [enz("TCTAGA", 1, 5, DAM)];
        assert_eq!(pos(&cuts(b"AATCTAGATCAA", false, &e)), vec![(3, DAM)]);
        assert_eq!(pos(&cuts(b"AATCTAGAAAAA", false, &e)), vec![(3, 0)]);
    }

    #[test]
    fn dcm_and_cpg_flags() {
        // A site whose last bases start CCAGG: the inner C (offset 1 of CCAGG) falls inside the site.
        let e = [enz("AGGCCT", 3, 3, DCM)];
        assert_eq!(pos(&cuts(b"TTAGGCCTGGTT", false, &e)), vec![(5, DCM)]);
        // A CpG-sensitive site containing CG.
        let e = [enz("GCGC", 3, 1, CPG)];
        assert_eq!(pos(&cuts(b"AAGCGCAA", false, &e)), vec![(5, CPG)]);
        // Insensitive enzymes are never flagged.
        let e = [enz("GCGC", 3, 1, 0)];
        assert_eq!(pos(&cuts(b"AAGCGCAA", false, &e)), vec![(5, 0)]);
    }

    #[test]
    fn dam_motif_across_the_origin_of_a_circle() {
        // Circle ...TCTAGA|TC...: GATC spans the origin and overlaps the XbaI site.
        let e = [enz("TCTAGA", 1, 5, DAM)];
        assert_eq!(pos(&cuts(b"TCAAAATCTAGA", true, &e)), vec![(7, DAM)]);
    }

    #[test]
    fn iupac_sites() {
        // HincII GTY^RAC.
        let e = [enz("GTYRAC", 3, 3, 0)];
        assert_eq!(pos(&cuts(b"AGTCGACAAGTTAACA", false, &e)), vec![(4, 0), (12, 0)]);
    }
}
