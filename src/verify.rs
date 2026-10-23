//! Clone verification: reads (Sanger traces, a consensus) aligned to a reference, evidence per base of a target region,
//! and a verdict (OpenBench plan, Sprint 4, decisions 11 and 12).
//!
//! Each read is trimmed (Mott, `align::trim`) when it has qualities, then aligned glocally (all of the read, free ends
//! in the reference) on both strands; the better strand is kept. A circular reference is read as itself followed by
//! as many of its first bases as the read is long, positions taken modulo its length. A base, a deletion or an
//! insertion is *qualifying* when its Phred is at least `min_q` (a deletion and an insertion take the lower Phred of
//! the read bases around them). When no read has qualities, every observation stands in for a qualifying one and the
//! verdicts say "no quality data".
//!
//! Verdict: `fail` when no read aligns; `incomplete` when a region base has no qualifying observation or qualifying
//! observations disagree with each other (a base, a deletion or an insertion); otherwise `mutation` when they differ from
//! the reference somewhere, else `pass`. Without qualities, `pass` and `mutation` read `agrees (no quality data)` and
//! `mutation (no quality data)`.

use crate::align::{self, Mode, Op, Scoring};
use crate::alphabet::{complement, letter, mask, satisfies};
use crate::json;
use std::collections::BTreeMap;

/// A read aligns when its score is at least this share of the best its trimmed bases could reach.
pub const MIN_SHARE: f64 = 0.2;
/// Reads shorter than this after trimming are not aligned.
pub const MIN_BASES: usize = 20;

pub struct Job {
    reference: Vec<u8>,
    circular: bool,
    /// Region positions, in order (a region across the origin of a circle wraps).
    region: Vec<usize>,
    min_q: u8,
    limit: f64,
    any_quality: bool,
    reads: Vec<ReadResult>,
    /// Per reference position: (observed allele, how it counts), alleles a base letter or '-'.
    bases: BTreeMap<usize, Vec<(u8, Weight)>>,
    /// Per reference position: reads whose alignment spans the junction after it (both neighbours aligned), with what
    /// they insert there ("" for nothing) and how that observation counts.
    after: BTreeMap<usize, Vec<(String, Weight)>>,
}

/// How an observation counts: `strong` when its read has qualities and it is at least `min_q` (and not an N); `seen`
/// when it is not an N. With qualities anywhere only strong observations count; without, every seen one does.
#[derive(Clone, Copy)]
struct Weight {
    strong: bool,
    seen: bool,
}

pub struct ReadResult {
    pub name: String,
    pub aligned: bool,
    pub reason: String,
    /// 1 forward, −1 reverse complement.
    pub strand: i32,
    /// The trimmed range, in the read's own coordinates and orientation.
    pub used: (usize, usize),
    /// Reference start of the alignment (modulo the length on a circle) and its operations.
    pub ref_start: usize,
    pub ops: Vec<(Op, u32)>,
    pub score: i32,
    pub best: i32,
    pub identity: f64,
}

#[derive(Debug, PartialEq, Eq, Clone)]
pub enum Diff {
    /// A substituted base at a reference position.
    Sub { pos: usize, reference: u8, read: u8 },
    /// Reference positions deleted, in order (consecutive positions merged into one run).
    Del { positions: Vec<usize> },
    /// Bases inserted after a reference position.
    Ins { after: usize, bases: String },
}

pub struct Report {
    pub verdict: &'static str,
    pub quality: bool,
    pub differences: Vec<Diff>,
    /// Region positions with no qualifying observation, as runs of reference positions.
    pub uncovered: Vec<Vec<usize>>,
    /// Reference positions (and junctions after them) where qualifying observations disagree.
    pub conflicts: Vec<usize>,
}

impl Job {
    /// `region` is [start, end) on the reference, 0-based; on a circle end < start wraps, and start == end is the whole
    /// circle.
    pub fn new(reference: &[u8], circular: bool, region: (usize, usize), min_q: u8, limit: f64) -> Result<Job, String> {
        let n = reference.len();
        let (s, e) = region;
        if n == 0 || s > n || e > n {
            return Err("the region is outside the reference".into());
        }
        let region: Vec<usize> = if circular {
            let len = if e > s { e - s } else { n - s + e };
            (0..len).map(|t| (s + t) % n).collect()
        } else if e > s {
            (s..e).collect()
        } else {
            return Err("on a linear reference the region's end comes after its start".into());
        };
        Ok(Job { reference: reference.to_ascii_uppercase(), circular, region, min_q, limit, any_quality: false, reads: Vec::new(), bases: BTreeMap::new(), after: BTreeMap::new() })
    }

    /// Trim and align one read; its observations join the evidence.
    pub fn add(&mut self, name: &str, bases: &[u8], quality: Option<&[u8]>) {
        let n = self.reference.len();
        let mut r = ReadResult { name: name.to_string(), aligned: false, reason: String::new(), strand: 1, used: (0, bases.len()), ref_start: 0, ops: Vec::new(), score: 0, best: 0, identity: 0.0 };
        if let Some(q) = quality {
            self.any_quality = true;
            r.used = align::trim(q, self.limit);
        }
        let (s, e) = r.used;
        if e - s < MIN_BASES {
            r.reason = format!("fewer than {MIN_BASES} bases after trimming");
            self.reads.push(r);
            return;
        }
        let fwd: Vec<u8> = bases[s..e].to_ascii_uppercase();
        let rev: Vec<u8> = fwd.iter().rev().map(|&b| match mask(b) { 0 => b'N', m => letter(complement(m)) }).collect();
        let quals: Option<Vec<u8>> = quality.map(|q| q[s..e].to_vec());
        // The reference, and on a circle as many of its first bases again as the read is long.
        let ext: Vec<u8> = if self.circular { [&self.reference[..], &self.reference[..fwd.len().min(n)]].concat() } else { self.reference.clone() };
        let sc = Scoring::DNA;
        let mut best: Option<(i32, align::Alignment)> = None;
        for strand in [1, -1] {
            let read = if strand == 1 { &fwd } else { &rev };
            match align::align(&ext, read, Mode::Glocal, &sc, None) {
                Ok(a) if best.as_ref().is_none_or(|(_, b)| a.score > b.score) => best = Some((strand, a)),
                Ok(_) => {}
                Err(err) => r.reason = err,
            }
        }
        let Some((strand, a)) = best else {
            self.reads.push(r);
            return;
        };
        r.best = sc.best(&fwd);
        r.score = a.score;
        r.strand = strand;
        r.identity = a.identity();
        if (a.score as f64) < MIN_SHARE * r.best as f64 {
            r.reason = "it does not match the reference".into();
            self.reads.push(r);
            return;
        }
        r.aligned = true;
        r.reason.clear();
        r.ref_start = a.a.0 % n;
        // The read as aligned (reverse-complemented on the reverse strand), with its qualities in the same order.
        let read = if strand == 1 { &fwd } else { &rev };
        let q: Option<Vec<u8>> = quals.map(|q| if strand == 1 { q } else { q.into_iter().rev().collect() });
        let qual = |j: usize| -> u8 { q.as_ref().map_or(0, |q| q[j.min(q.len() - 1)]) };
        let weigh = |p: u8, base: u8| Weight { strong: q.is_some() && p >= self.min_q && mask(base) != 15, seen: mask(base) != 15 };
        let (mut i, mut j) = (a.a.0, a.b.0);
        // The bases inserted since the last aligned reference position, and their lowest Phred.
        let mut pending: Option<(String, u8)> = None;
        let mut last: Option<usize> = None;
        for &(op, len) in &a.ops {
            for _ in 0..len {
                match op {
                    Op::Match | Op::Mismatch | Op::Del => {
                        let pos = i % n;
                        // The junction after the previous reference position: spanned by this read, with what it inserts.
                        if let Some(prev) = last {
                            let (ins, iq) = pending.take().unwrap_or((String::new(), u8::MAX));
                            let around = qual(j.saturating_sub(1)).min(qual(j)).min(iq);
                            self.after.entry(prev).or_default().push((ins, weigh(around, b'A')));
                        }
                        if op == Op::Del {
                            let around = qual(j.saturating_sub(1)).min(qual(j));
                            self.bases.entry(pos).or_default().push((b'-', weigh(around, b'A')));
                        } else {
                            // An N in the read is never counted.
                            self.bases.entry(pos).or_default().push((read[j], weigh(qual(j), read[j])));
                            j += 1;
                        }
                        last = Some(pos);
                        i += 1;
                    }
                    Op::Ins => {
                        let (mut s, q) = pending.take().unwrap_or((String::new(), u8::MAX));
                        s.push(read[j] as char);
                        pending = Some((s, q.min(qual(j))));
                        j += 1;
                    }
                }
            }
        }
        r.ops = a.ops;
        self.reads.push(r);
    }

    pub fn reads(&self) -> &[ReadResult] {
        &self.reads
    }

    pub fn report(&self) -> Report {
        let quality = self.any_quality;
        let counts = |w: &Weight| if quality { w.strong } else { w.seen };
        let mut differences = Vec::new();
        let mut uncovered: Vec<Vec<usize>> = Vec::new();
        let mut conflicts = Vec::new();
        let mut deleted: Vec<usize> = Vec::new();
        let flush = |deleted: &mut Vec<usize>, differences: &mut Vec<Diff>| {
            if !deleted.is_empty() {
                differences.push(Diff::Del { positions: std::mem::take(deleted) });
            }
        };
        let none: Vec<(u8, Weight)> = Vec::new();
        for (t, &pos) in self.region.iter().enumerate() {
            let obs: Vec<u8> = self.bases.get(&pos).unwrap_or(&none).iter().filter(|(_, w)| counts(w)).map(|&(b, _)| b).collect();
            if obs.is_empty() {
                flush(&mut deleted, &mut differences);
                match uncovered.last_mut() {
                    Some(run) if run.last().is_some_and(|&p| self.region[t - 1] == p) => run.push(pos),
                    _ => uncovered.push(vec![pos]),
                }
                continue;
            }
            let reference = self.reference[pos];
            // Qualifying observations agree when they are the same letter (a deletion is '-').
            if obs.iter().any(|&b| b != obs[0]) {
                flush(&mut deleted, &mut differences);
                conflicts.push(pos);
                continue;
            }
            let b = obs[0];
            if b == b'-' {
                deleted.push(pos);
            } else {
                flush(&mut deleted, &mut differences);
                if !satisfies(mask(b), mask(reference)) {
                    differences.push(Diff::Sub { pos, reference, read: b });
                }
            }
            // The junction after this base, when the next region base follows it.
            if t + 1 < self.region.len() {
                let spans: Vec<&String> = self.after.get(&pos).map(|v| v.iter().filter(|(_, w)| counts(w)).map(|(s, _)| s).collect()).unwrap_or_default();
                if let Some(first) = spans.first() {
                    if spans.iter().any(|s| s != first) {
                        flush(&mut deleted, &mut differences);
                        conflicts.push(pos);
                    } else if !first.is_empty() {
                        flush(&mut deleted, &mut differences);
                        differences.push(Diff::Ins { after: pos, bases: first.to_string() });
                    }
                }
            }
        }
        flush(&mut deleted, &mut differences);
        let aligned = self.reads.iter().any(|r| r.aligned);
        let verdict = if !aligned {
            "fail"
        } else if !uncovered.is_empty() || !conflicts.is_empty() {
            "incomplete"
        } else if !differences.is_empty() {
            if quality { "mutation" } else { "mutation (no quality data)" }
        } else if quality {
            "pass"
        } else {
            "agrees (no quality data)"
        };
        Report { verdict, quality, differences, uncovered, conflicts }
    }
}

/// The report and the reads as JSON for the page.
pub fn to_json(job: &Job, rep: &Report) -> String {
    let mut o = String::from("{\"verdict\":");
    json::string(&mut o, rep.verdict);
    o.push_str(&format!(",\"quality\":{},\"differences\":[", rep.quality));
    for (k, d) in rep.differences.iter().enumerate() {
        if k > 0 {
            o.push(',');
        }
        match d {
            Diff::Sub { pos, reference, read } => o.push_str(&format!("{{\"kind\":\"sub\",\"pos\":{pos},\"ref\":\"{}\",\"alt\":\"{}\"}}", *reference as char, *read as char)),
            Diff::Del { positions } => {
                o.push_str("{\"kind\":\"del\",\"positions\":");
                json::numbers(&mut o, positions.iter());
                o.push('}');
            }
            Diff::Ins { after, bases } => {
                o.push_str(&format!("{{\"kind\":\"ins\",\"after\":{after},\"alt\":"));
                json::string(&mut o, bases);
                o.push('}');
            }
        }
    }
    o.push_str("],\"uncovered\":[");
    for (k, run) in rep.uncovered.iter().enumerate() {
        if k > 0 {
            o.push(',');
        }
        // A run as its first and last positions.
        json::numbers(&mut o, [run[0], run[run.len() - 1]]);
    }
    o.push_str("],\"conflicts\":");
    json::numbers(&mut o, rep.conflicts.iter().take(1000));
    o.push_str(",\"reads\":[");
    for (k, r) in job.reads.iter().enumerate() {
        if k > 0 {
            o.push(',');
        }
        o.push_str("{\"name\":");
        json::string(&mut o, &r.name);
        o.push_str(&format!(",\"aligned\":{},\"reason\":", r.aligned));
        json::string(&mut o, &r.reason);
        o.push_str(&format!(",\"strand\":{},\"used\":[{},{}],\"refStart\":{},\"score\":{},\"best\":{},\"identity\":{:.4},\"ops\":", r.strand, r.used.0, r.used.1, r.ref_start, r.score, r.best, r.identity));
        json::numbers(&mut o, r.ops.iter().flat_map(|&(op, n)| [op as u32, n]));
        o.push('}');
    }
    o.push_str("]}");
    o
}

#[cfg(test)]
mod tests {
    use super::*;

    fn bases(n: usize, seed: u32) -> Vec<u8> {
        let mut r = seed;
        (0..n).map(|_| { r ^= r << 13; r ^= r >> 17; r ^= r << 5; b"ACGT"[(r >> 30) as usize] }).collect()
    }
    fn q(n: usize, v: u8) -> Vec<u8> {
        vec![v; n]
    }
    fn rc(s: &[u8]) -> Vec<u8> {
        s.iter().rev().map(|&b| letter(complement(mask(b)))).collect()
    }

    #[test]
    fn a_read_that_matches_passes_and_one_with_a_substitution_is_a_mutation() {
        let r = bases(400, 11);
        let mut job = Job::new(&r, false, (50, 350), 20, 0.05).unwrap();
        job.add("f", &r[20..380], Some(&q(360, 40)));
        assert_eq!(job.report().verdict, "pass");
        let mut read = r[20..380].to_vec();
        read[180] = if read[180] == b'A' { b'C' } else { b'A' };
        let mut job = Job::new(&r, false, (50, 350), 20, 0.05).unwrap();
        job.add("f", &read, Some(&q(360, 40)));
        let rep = job.report();
        assert_eq!((rep.verdict, rep.differences.clone()), ("mutation", vec![Diff::Sub { pos: 200, reference: r[200], read: read[180] }]));
    }

    #[test]
    fn a_reverse_read_low_quality_ends_and_partial_coverage() {
        let r = bases(400, 12);
        // The reverse complement aligns on the reverse strand; Q5 ends are trimmed away.
        let mut read = rc(&r[40..360]);
        let mut qual = q(read.len(), 40);
        qual[..15].fill(5);
        read[3] = b'A';
        let mut job = Job::new(&r, false, (60, 340), 20, 0.05).unwrap();
        job.add("r", &read, Some(&qual));
        let rep = job.report();
        assert_eq!((rep.verdict, job.reads()[0].strand, job.reads()[0].used), ("pass", -1, (15, read.len())));
        // A region the read does not reach: incomplete, with the uncovered run.
        let mut job = Job::new(&r, false, (10, 340), 20, 0.05).unwrap();
        job.add("f", &r[40..380], Some(&q(340, 40)));
        let rep = job.report();
        assert_eq!((rep.verdict, rep.uncovered.len(), rep.uncovered[0][0], *rep.uncovered[0].last().unwrap()), ("incomplete", 1, 10, 39));
    }

    #[test]
    fn reads_that_disagree_are_incomplete_and_low_quality_does_not_count() {
        let r = bases(400, 13);
        let mut other = r[20..380].to_vec();
        other[100] = if other[100] == b'G' { b'T' } else { b'G' };
        let mut job = Job::new(&r, false, (50, 350), 20, 0.05).unwrap();
        job.add("a", &r[20..380], Some(&q(360, 40)));
        job.add("b", &other, Some(&q(360, 40)));
        let rep = job.report();
        assert_eq!((rep.verdict, rep.conflicts.clone()), ("incomplete", vec![120]));
        // The same disagreement at Phred 15 (below 20) does not count: the good read alone passes.
        let mut qual = q(360, 40);
        qual[100] = 15;
        let mut job = Job::new(&r, false, (50, 350), 20, 0.05).unwrap();
        job.add("a", &r[20..380], Some(&q(360, 40)));
        job.add("b", &other, Some(&qual));
        assert_eq!(job.report().verdict, "pass");
    }

    #[test]
    fn deletions_and_insertions_are_found_and_merged() {
        let r = bases(400, 14);
        let mut read = r[20..380].to_vec();
        read.drain(150..153); // reference 170-172 deleted
        read.splice(250..250, b"GGGTTT".iter().copied()); // after reference 272 (read 250 is reference 273 before the insert)
        let mut job = Job::new(&r, false, (50, 350), 20, 0.05).unwrap();
        let qual = q(read.len(), 40);
        job.add("f", &read, Some(&qual));
        let rep = job.report();
        assert_eq!(rep.verdict, "mutation");
        assert!(rep.differences.iter().any(|d| matches!(d, Diff::Del { positions } if positions.len() == 3)), "{:?}", rep.differences);
        assert!(rep.differences.iter().any(|d| matches!(d, Diff::Ins { bases, .. } if bases.len() == 6)), "{:?}", rep.differences);
    }

    #[test]
    fn a_circular_reference_and_a_region_across_its_origin() {
        let r = bases(300, 15);
        // A read over the origin: the last 80 bases, then the first 120.
        let read = [&r[220..], &r[..120]].concat();
        let mut job = Job::new(&r, true, (250, 100), 20, 0.05).unwrap();
        job.add("f", &read, Some(&q(read.len(), 40)));
        let rep = job.report();
        assert_eq!((rep.verdict, job.reads()[0].ref_start), ("pass", 220));
    }

    #[test]
    fn evidence_without_qualities_never_passes_and_a_foreign_read_fails() {
        let r = bases(400, 16);
        let mut job = Job::new(&r, false, (0, 400), 20, 0.05).unwrap();
        job.add("consensus", &r, None);
        assert_eq!(job.report().verdict, "agrees (no quality data)");
        let mut c = r.clone();
        c[10] = if c[10] == b'A' { b'T' } else { b'A' };
        let mut job = Job::new(&r, false, (0, 400), 20, 0.05).unwrap();
        job.add("consensus", &c, None);
        assert_eq!(job.report().verdict, "mutation (no quality data)");
        let mut job = Job::new(&r, false, (50, 350), 20, 0.05).unwrap();
        job.add("other", &bases(300, 99), Some(&q(300, 40)));
        assert_eq!((job.report().verdict, job.reads()[0].aligned), ("fail", false));
    }
}
