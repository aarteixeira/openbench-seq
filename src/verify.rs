//! Clone verification: reads (Sanger traces, a consensus) aligned to a reference, evidence per base of a target region,
//! and a verdict (OpenBench plan, Sprint 4, decisions 11 and 12; 4B-1 plan revision 2).
//!
//! Each read is trimmed (Mott, `align::trim`) when it has qualities, its bases and qualities reverse-complemented and
//! reversed together for the reverse strand, and aligned glocally (all of the read, free ends in the reference) on
//! both strands. A circular reference is read as itself followed by as many of its first bases as the read is long,
//! positions taken modulo its length; on a linear one, read bases past either end are left unaligned (overlap mode). A
//! read is *placed* when at least 50 of its bases align and its score reaches 0.4 per aligned read base (an ungapped
//! identity of 68%), and *ambiguous* (no evidence) when the other strand, or another placement on the same strand,
//! scores as high (the aligner's exact ties). A read running more than once around a circle is no evidence.
//!
//! Observations: per reference position a base or a deletion, per junction between two spanned positions what is
//! inserted there (possibly nothing). Quality: a base's Phred; a deletion's the lower of the two read bases around it;
//! an insertion's the lowest of its bases and the two around it. A read letter other than A, C, G, T is a no-call. An
//! observation *qualifies* when its read has qualities and its quality is at least `min_q` (when no read of the clone
//! has qualities, every call qualifies and the verdict says "no quality data").
//!
//! Per region position: uncovered without a qualifying call (or when the reference letter is not A, C, G or T); a
//! conflict when two calls of placed reads differ, whatever their quality. A junction inside the region with no
//! qualifying observation is unsupported. Verdict: `fail` when a read had at least 50 bases after trimming and no
//! read matches the reference (none placed or ambiguous); `incomplete` when anything is
//! uncovered, unsupported or in conflict; `mutation` when a call differs from the reference; else `pass`. Without
//! qualities `agrees (no quality data)` and `mutation (no quality data)` replace `pass` and `mutation`.

use crate::align::{self, Mode, Op, Scoring};
use crate::alphabet::{complement, letter, mask};
use crate::json;
use std::collections::BTreeMap;

/// A read is placed when its score is at least this much per aligned read base.
pub const MIN_PER_BASE: f64 = 0.4;
/// Reads with fewer bases than this after trimming, or aligned to the reference, are not placed.
pub const MIN_BASES: usize = 50;

pub struct Job {
    reference: Vec<u8>,
    circular: bool,
    /// Region positions, in order (a region across the origin of a circle wraps).
    region: Vec<usize>,
    min_q: u8,
    limit: f64,
    any_quality: bool,
    /// Reads with at least MIN_BASES after trimming.
    usable: usize,
    reads: Vec<ReadResult>,
    /// Per reference position: what each placed read shows there (a base letter or '-').
    bases: BTreeMap<usize, Vec<Obs<u8>>>,
    /// Per reference position: what each placed read spanning the junction after it inserts there ("" for nothing).
    after: BTreeMap<usize, Vec<Obs<String>>>,
}

/// One read's observation: which read, what it shows, its quality (None without qualities), whether that is a call
/// (A, C, G, T, '-'; no N or other code), and whether it qualifies (a call with a quality of at least `min_q`).
struct Obs<T> {
    read: usize,
    allele: T,
    q: Option<u8>,
    call: bool,
    strong: bool,
}

pub struct ReadResult {
    pub name: String,
    pub placed: bool,
    /// Matches the reference, but as well at a second place: no evidence.
    pub ambiguous: bool,
    pub reason: String,
    /// 1 forward, −1 reverse complement.
    pub strand: i32,
    /// The trimmed range, in the read's own coordinates and orientation.
    pub used: (usize, usize),
    /// Reference start of the alignment (modulo the length on a circle), where in the trimmed read (as aligned: reverse
    /// complemented on the reverse strand) its operations start, and the operations (unaligned ends left out).
    pub ref_start: usize,
    pub read_start: usize,
    pub ops: Vec<(Op, u32)>,
    pub score: i32,
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
    /// Region positions with no qualifying call, as runs of reference positions.
    pub uncovered: Vec<Vec<usize>>,
    /// Reference positions (and junctions after them) where calls disagree, with each read's call and quality.
    pub conflicts: Vec<Conflict>,
    /// Positions whose junction with the next region position has no qualifying observation.
    pub unsupported: Vec<usize>,
    /// Region positions with a call, and those whose call rests on exactly one qualifying observation.
    pub covered: usize,
    pub single: usize,
}

pub struct Conflict {
    pub pos: usize,
    /// The junction after `pos` rather than the base at it.
    pub junction: bool,
    /// (read index, its call: a base, '-' or the inserted bases, its quality).
    pub calls: Vec<(usize, String, Option<u8>)>,
}

fn is_call(b: u8) -> bool {
    matches!(b, b'A' | b'C' | b'G' | b'T' | b'-')
}

impl Job {
    /// `region` is [start, end) on the reference, 0-based; on a circle end < start wraps the origin (the whole circle is
    /// [0, n)). An empty region, or one outside the reference, is refused.
    pub fn new(reference: &[u8], circular: bool, region: (usize, usize), min_q: u8, limit: f64) -> Result<Job, String> {
        let n = reference.len();
        let (s, e) = region;
        if n == 0 || s >= n || e > n {
            return Err("the region is outside the reference".into());
        }
        if s == e {
            return Err("the region is empty".into());
        }
        let region: Vec<usize> = if circular {
            let len = if e > s { e - s } else { n - s + e };
            (0..len).map(|t| (s + t) % n).collect()
        } else if e > s {
            (s..e).collect()
        } else {
            return Err("on a linear reference the region's end comes after its start".into());
        };
        Ok(Job { reference: reference.to_ascii_uppercase(), circular, region, min_q, limit, any_quality: false, usable: 0, reads: Vec::new(), bases: BTreeMap::new(), after: BTreeMap::new() })
    }

    /// Trim, place and read one read; a placed read's observations join the evidence.
    pub fn add(&mut self, name: &str, bases: &[u8], quality: Option<&[u8]>) {
        let n = self.reference.len();
        let mut r = ReadResult { name: name.to_string(), placed: false, ambiguous: false, reason: String::new(), strand: 1, used: (0, bases.len()), ref_start: 0, read_start: 0, ops: Vec::new(), score: 0, identity: 0.0 };
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
        self.usable += 1;
        let fwd: Vec<u8> = bases[s..e].to_ascii_uppercase();
        let rev: Vec<u8> = fwd.iter().rev().map(|&b| match mask(b) { 0 => b'N', m => letter(complement(m)) }).collect();
        // The reference, and on a circle as many of its first bases again as the read is long.
        let ext: Vec<u8> = if self.circular { [&self.reference[..], &self.reference[..fwd.len().min(n)]].concat() } else { self.reference.clone() };
        let sc = Scoring::DNA;
        let mode = if self.circular { Mode::Glocal } else { Mode::Overlap };
        let mut tries: Vec<(i32, align::Alignment)> = Vec::new();
        for strand in [1, -1] {
            match align::align(&ext, if strand == 1 { &fwd } else { &rev }, mode, &sc, None) {
                Ok(a) => tries.push((strand, a)),
                Err(err) => r.reason = err,
            }
        }
        let Some(best) = tries.iter().max_by_key(|(_, a)| a.score).map(|(s, _)| *s) else {
            self.reads.push(r);
            return;
        };
        let (strand, a) = tries.iter().find(|(s, _)| *s == best).cloned().unwrap();
        let other = tries.iter().find(|(s, _)| *s != strand).map(|(_, a)| a.score);
        let read: &[u8] = if strand == 1 { &fwd } else { &rev };
        r.strand = strand;
        r.score = a.score;
        r.identity = a.identity();
        let aligned = a.b.1 - a.b.0;
        if aligned < MIN_BASES || (a.score as f64) < MIN_PER_BASE * aligned as f64 {
            r.reason = "it does not match the reference (below 68% identity over at least 50 bases)".into();
            self.reads.push(r);
            return;
        }
        // One traversal of a circle at most: a read running further (a consensus with an extra copy) is no evidence.
        if a.a.1 - a.a.0 > n {
            r.reason = "it runs more than once around the circular reference".into();
            self.reads.push(r);
            return;
        }
        // Another placement as good: the other strand, or another end on this strand (on a circle, an end and the same
        // end one turn later are one placement).
        let key = |(i, j): (usize, usize)| if self.circular { (i % n, j) } else { (i, j) };
        let here = key((a.a.1, a.b.1));
        if other.is_some_and(|o| o >= a.score) || a.ties.iter().any(|&t| key(t) != here) {
            r.ambiguous = true;
            r.reason = "it matches more than one place on the reference".into();
            self.reads.push(r);
            return;
        }
        r.placed = true;
        r.reason.clear();
        r.ref_start = a.a.0 % n;
        r.read_start = a.b.0;
        // Qualities in the read's aligned order (reversed with the bases on the reverse strand).
        let q: Option<Vec<u8>> = quality.map(|q| {
            let q = q[s..e].to_vec();
            if strand == 1 { q } else { q.into_iter().rev().collect() }
        });
        let qual = |j: usize| -> u8 { q.as_ref().map_or(0, |q| q[j.min(q.len() - 1)]) };
        let strong = |p: u8| q.is_some() && p >= self.min_q;
        let (mut i, mut j) = (a.a.0, a.b.0);
        let at = self.reads.len();
        // The bases inserted since the last aligned reference position, and their lowest Phred.
        let mut pending: Option<(String, u8)> = None;
        let (mut last, mut lead): (Option<usize>, Option<(String, u8)>) = (None, None);
        for &(op, len) in &a.ops {
            for _ in 0..len {
                match op {
                    Op::Match | Op::Mismatch | Op::Del => {
                        let pos = i % n;
                        // Bases before the first aligned position are no junction's (kept for a whole circle below).
                        if last.is_none() {
                            lead = pending.take();
                        }
                        // The junction after the previous reference position: spanned by this read, with what it inserts.
                        if let Some(prev) = last {
                            let (ins, iq) = pending.take().unwrap_or((String::new(), u8::MAX));
                            let around = qual(j.saturating_sub(1)).min(qual(j)).min(iq);
                            let call = ins.bytes().all(|b| b != b'-' && is_call(b));
                            self.after.entry(prev).or_default().push(Obs { read: at, allele: ins, q: q.as_ref().map(|_| around), call, strong: call && strong(around) });
                        }
                        if op == Op::Del {
                            let around = qual(j.saturating_sub(1)).min(qual(j));
                            self.bases.entry(pos).or_default().push(Obs { read: at, allele: b'-', q: q.as_ref().map(|_| around), call: true, strong: strong(around) });
                        } else {
                            let call = is_call(read[j]);
                            self.bases.entry(pos).or_default().push(Obs { read: at, allele: read[j], q: q.as_ref().map(|_| qual(j)), call, strong: call && strong(qual(j)) });
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
        // A read covering a circle exactly once is the circle (a whole-plasmid consensus): its end joins its start, so
        // it also spans the junction there, holding the bases it has after its last position and before its first.
        if let (true, Some(prev)) = (self.circular && a.a.1 - a.a.0 == n, last) {
            let (tail, tq) = pending.take().unwrap_or((String::new(), u8::MAX));
            let (head, hq) = lead.take().unwrap_or((String::new(), u8::MAX));
            let ins = tail + &head;
            let around = qual(read.len() - 1).min(qual(0)).min(tq).min(hq);
            let call = ins.bytes().all(|b| b != b'-' && is_call(b));
            self.after.entry(prev).or_default().push(Obs { read: at, allele: ins, q: q.as_ref().map(|_| around), call, strong: call && strong(around) });
        }
        r.ops = a.ops;
        self.reads.push(r);
    }

    pub fn reads(&self) -> &[ReadResult] {
        &self.reads
    }

    pub fn report(&self) -> Report {
        let quality = self.any_quality;
        // Without qualities anywhere, every call qualifies.
        let qualifies = |call: bool, strong: bool| if quality { strong } else { call };
        let mut differences = Vec::new();
        let (mut uncovered, mut conflicts, mut unsupported): (Vec<Vec<usize>>, Vec<Conflict>, Vec<usize>) = (Vec::new(), Vec::new(), Vec::new());
        let (mut covered, mut single) = (0usize, 0usize);
        let mut deleted: Vec<usize> = Vec::new();
        let flush = |deleted: &mut Vec<usize>, differences: &mut Vec<Diff>| {
            if !deleted.is_empty() {
                differences.push(Diff::Del { positions: std::mem::take(deleted) });
            }
        };
        for (t, &pos) in self.region.iter().enumerate() {
            let obs = self.bases.get(&pos).map(|v| &v[..]).unwrap_or(&[]);
            let reference = self.reference[pos];
            let calls: Vec<u8> = obs.iter().filter(|o| o.call).map(|o| o.allele).collect();
            let good = obs.iter().filter(|o| qualifies(o.call, o.strong)).count();
            if calls.iter().any(|&b| b != calls[0]) {
                flush(&mut deleted, &mut differences);
                let calls = obs.iter().filter(|o| o.call).map(|o| (o.read, (o.allele as char).to_string(), o.q)).collect();
                conflicts.push(Conflict { pos, junction: false, calls });
            } else if good == 0 || !matches!(reference, b'A' | b'C' | b'G' | b'T') {
                flush(&mut deleted, &mut differences);
                match uncovered.last_mut() {
                    Some(run) if t > 0 && run.last() == Some(&self.region[t - 1]) => run.push(pos),
                    _ => uncovered.push(vec![pos]),
                }
            } else {
                covered += 1;
                single += (good == 1) as usize;
                let b = calls[0];
                if b == b'-' {
                    deleted.push(pos);
                } else {
                    flush(&mut deleted, &mut differences);
                    if b != reference {
                        differences.push(Diff::Sub { pos, reference, read: b });
                    }
                }
            }
            // The junction after this base, when the next region base follows it.
            if t + 1 < self.region.len() {
                let obs = self.after.get(&pos).map(|v| &v[..]).unwrap_or(&[]);
                let calls: Vec<&String> = obs.iter().filter(|o| o.call).map(|o| &o.allele).collect();
                if calls.iter().any(|s| *s != calls[0]) {
                    flush(&mut deleted, &mut differences);
                    let calls = obs.iter().filter(|o| o.call).map(|o| (o.read, o.allele.clone(), o.q)).collect();
                    conflicts.push(Conflict { pos, junction: true, calls });
                } else if !obs.iter().any(|o| qualifies(o.call, o.strong)) {
                    unsupported.push(pos);
                } else if !calls[0].is_empty() {
                    flush(&mut deleted, &mut differences);
                    differences.push(Diff::Ins { after: pos, bases: calls[0].to_string() });
                }
            }
        }
        flush(&mut deleted, &mut differences);
        // Failed only on usable evidence: a clone whose reads all match more than one place, or have too few good
        // bases to judge, is incomplete.
        let matched = self.reads.iter().any(|r| r.placed || r.ambiguous);
        let verdict = if !matched && self.usable > 0 {
            "fail"
        } else if !uncovered.is_empty() || !conflicts.is_empty() || !unsupported.is_empty() {
            "incomplete"
        } else if !differences.is_empty() {
            if quality { "mutation" } else { "mutation (no quality data)" }
        } else if quality {
            "pass"
        } else {
            "agrees (no quality data)"
        };
        Report { verdict, quality, differences, uncovered, conflicts, unsupported, covered, single }
    }
}

/// The report and the reads as JSON for the page.
pub fn to_json(job: &Job, rep: &Report) -> String {
    let mut o = String::from("{\"verdict\":");
    json::string(&mut o, rep.verdict);
    o.push_str(&format!(",\"quality\":{},\"covered\":{},\"single\":{},\"differences\":[", rep.quality, rep.covered, rep.single));
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
    o.push_str("],\"conflicts\":[");
    for (k, c) in rep.conflicts.iter().take(1000).enumerate() {
        if k > 0 {
            o.push(',');
        }
        o.push_str(&format!("{{\"pos\":{},\"junction\":{},\"calls\":[", c.pos, c.junction));
        for (x, (read, allele, q)) in c.calls.iter().enumerate() {
            if x > 0 {
                o.push(',');
            }
            o.push_str(&format!("{{\"read\":{read},\"call\":"));
            json::string(&mut o, allele);
            o.push_str(&match q { Some(q) => format!(",\"q\":{q}}}"), None => "}".into() });
        }
        o.push_str("]}");
    }
    o.push_str("],\"unsupported\":");
    json::numbers(&mut o, rep.unsupported.iter().take(1000));
    o.push_str(",\"reads\":[");
    for (k, r) in job.reads.iter().enumerate() {
        if k > 0 {
            o.push(',');
        }
        o.push_str("{\"name\":");
        json::string(&mut o, &r.name);
        o.push_str(&format!(",\"placed\":{},\"ambiguous\":{},\"reason\":", r.placed, r.ambiguous));
        json::string(&mut o, &r.reason);
        o.push_str(&format!(",\"strand\":{},\"used\":[{},{}],\"refStart\":{},\"readStart\":{},\"score\":{},\"identity\":{:.4},\"ops\":", r.strand, r.used.0, r.used.1, r.ref_start, r.read_start, r.score, r.identity));
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
    fn other(b: u8) -> u8 {
        if b == b'A' { b'C' } else { b'A' }
    }
    fn job(r: &[u8], region: (usize, usize)) -> Job {
        Job::new(r, false, region, 20, 0.05).unwrap()
    }

    #[test]
    fn a_read_that_matches_passes_and_one_with_a_substitution_is_a_mutation() {
        let r = bases(400, 11);
        let mut j = job(&r, (50, 350));
        j.add("f", &r[20..380], Some(&q(360, 40)));
        let rep = j.report();
        assert_eq!((rep.verdict, rep.single), ("pass", 300));
        let mut read = r[20..380].to_vec();
        read[180] = other(read[180]);
        let mut j = job(&r, (50, 350));
        j.add("f", &read, Some(&q(360, 40)));
        let rep = j.report();
        assert_eq!((rep.verdict, rep.differences.clone()), ("mutation", vec![Diff::Sub { pos: 200, reference: r[200], read: read[180] }]));
    }

    #[test]
    fn a_reverse_read_keeps_each_quality_with_its_base() {
        let r = bases(400, 12);
        // The reverse complement aligns on the reverse strand; Q5 ends are trimmed away.
        let read = rc(&r[40..360]);
        let mut qual = q(read.len(), 40);
        qual[..15].fill(5);
        // Read base 100 (as sequenced) is reference 360 − 1 − 100 = 259: Q10 there leaves 259 uncovered.
        qual[100] = 10;
        let mut j = job(&r, (60, 340));
        j.add("r", &read, Some(&qual));
        let rep = j.report();
        assert_eq!((j.reads()[0].strand, j.reads()[0].used, rep.verdict, rep.uncovered.clone()), (-1, (15, read.len()), "incomplete", vec![vec![259]]));
        // A region the read does not reach: incomplete, with the uncovered run.
        let mut j = job(&r, (10, 340));
        j.add("f", &r[40..380], Some(&q(340, 40)));
        let rep = j.report();
        assert_eq!((rep.verdict, rep.uncovered.len(), rep.uncovered[0][0], *rep.uncovered[0].last().unwrap()), ("incomplete", 1, 10, 39));
    }

    #[test]
    fn reads_that_disagree_are_a_conflict_at_any_quality_and_an_n_is_no_call() {
        let r = bases(400, 13);
        let mut b = r[20..380].to_vec();
        b[100] = other(b[100]);
        // Q15 against Q30: still a conflict (decision 11).
        let mut qual = q(360, 40);
        qual[100] = 15;
        let mut j = job(&r, (50, 350));
        j.add("a", &r[20..380], Some(&q(360, 30)));
        j.add("b", &b, Some(&qual));
        let rep = j.report();
        assert_eq!(rep.verdict, "incomplete");
        // Each read's call with its Phred, so the conflict can be judged.
        assert_eq!((rep.conflicts.len(), rep.conflicts[0].pos, rep.conflicts[0].junction), (1, 120, false));
        assert_eq!(rep.conflicts[0].calls, vec![(0, (r[120] as char).to_string(), Some(30)), (1, (b[100] as char).to_string(), Some(15))]);
        // An N (and another IUPAC code) neither conflicts nor counts.
        let mut n = r[20..380].to_vec();
        (n[100], n[101]) = (b'N', b'R');
        let mut j = job(&r, (50, 350));
        j.add("a", &r[20..380], Some(&q(360, 40)));
        j.add("n", &n, Some(&q(360, 40)));
        let rep = j.report();
        assert_eq!((rep.verdict, rep.covered, rep.single), ("pass", 300, 2));
    }

    #[test]
    fn deletions_and_insertions_are_named_once_and_need_quality_around_them() {
        let r = bases(400, 14);
        let mut read = r[20..380].to_vec();
        read.drain(150..153); // reference 170-172 deleted
        read.splice(250..250, b"GGGTTT".iter().copied()); // after reference 272 (read 250 is reference 273 before the insert)
        let mut j = job(&r, (50, 350));
        j.add("f", &read, Some(&q(read.len(), 40)));
        let rep = j.report();
        assert_eq!(rep.verdict, "mutation");
        assert!(rep.differences.iter().any(|d| matches!(d, Diff::Del { positions } if positions.len() == 3)), "{:?}", rep.differences);
        assert!(rep.differences.iter().any(|d| matches!(d, Diff::Ins { bases, .. } if bases.len() == 6)), "{:?}", rep.differences);
        // A second read without the insertion: the junction is a conflict, whatever else agrees.
        let mut j = job(&r, (50, 350));
        j.add("f", &read, Some(&q(read.len(), 40)));
        j.add("g", &r[20..380], Some(&q(360, 40)));
        let rep = j.report();
        assert_eq!((rep.verdict, rep.conflicts.iter().filter(|c| c.junction).count()), ("incomplete", 1));
        // The insertion's bases at Q10: the junction has no qualifying evidence.
        let mut qual = q(read.len(), 40);
        qual[250..256].fill(10);
        let mut j = job(&r, (50, 350));
        j.add("f", &read, Some(&qual));
        let rep = j.report();
        assert_eq!((rep.verdict, rep.unsupported.len()), ("incomplete", 1));
    }

    #[test]
    fn a_deletion_in_a_run_is_the_same_deletion_in_every_read() {
        // CAAAAG in the reference, CAAAG in both reads: one deletion, placed on the last A, no conflict.
        let mut r = bases(400, 17);
        r[200..206].copy_from_slice(b"CAAAAG");
        let read: Vec<u8> = [&r[20..203], &r[204..380]].concat();
        let mut j = job(&r, (50, 350));
        j.add("f", &read, Some(&q(read.len(), 40)));
        j.add("r", &rc(&read), Some(&q(read.len(), 40)));
        let rep = j.report();
        assert_eq!((rep.verdict, rep.differences.clone()), ("mutation", vec![Diff::Del { positions: vec![204] }]));
    }

    #[test]
    fn a_deletion_is_placed_on_the_whole_reference_whatever_the_region() {
        // A run of four As across the region's end: the deletion is the run's last A (204), outside a region ending at
        // 203, so that region passes; the whole record names it.
        let mut r = bases(400, 26);
        r[200..206].copy_from_slice(b"CAAAAG");
        let read: Vec<u8> = [&r[20..202], &r[203..380]].concat();
        let mut j = job(&r, (50, 204));
        j.add("f", &read, Some(&q(read.len(), 40)));
        assert_eq!(j.report().verdict, "pass");
        let mut j = job(&r, (50, 350));
        j.add("f", &read, Some(&q(read.len(), 40)));
        assert_eq!(j.report().differences, vec![Diff::Del { positions: vec![204] }]);
    }

    #[test]
    fn a_circular_reference_and_a_region_across_its_origin() {
        let r = bases(300, 15);
        // A read over the origin: the last 80 bases, then the first 120.
        let read = [&r[220..], &r[..120]].concat();
        let mut j = Job::new(&r, true, (250, 100), 20, 0.05).unwrap();
        j.add("f", &read, Some(&q(read.len(), 40)));
        let rep = j.report();
        assert_eq!((rep.verdict, j.reads()[0].ref_start), ("pass", 220));
    }

    #[test]
    fn evidence_without_qualities_never_passes_and_a_foreign_read_fails() {
        let r = bases(400, 16);
        let mut j = job(&r, (0, 400));
        j.add("consensus", &r, None);
        assert_eq!(j.report().verdict, "agrees (no quality data)");
        let mut c = r.clone();
        c[10] = other(c[10]);
        let mut j = job(&r, (0, 400));
        j.add("consensus", &c, None);
        assert_eq!(j.report().verdict, "mutation (no quality data)");
        let mut j = job(&r, (50, 350));
        j.add("other", &bases(300, 99), Some(&q(300, 40)));
        assert_eq!((j.report().verdict, j.reads()[0].placed), ("fail", false));
    }

    #[test]
    fn placement_needs_68_percent_and_one_place() {
        let r = bases(1000, 18);
        // Every fourth base changed (75% identity): placed, with its differences.
        let mut read = r[100..500].to_vec();
        for k in (0..read.len()).step_by(4) {
            read[k] = other(read[k]);
        }
        let mut j = job(&r, (150, 450));
        j.add("f", &read, Some(&q(400, 40)));
        assert!(j.reads()[0].placed, "{}", j.reads()[0].reason);
        // A reference holding the same 300 bases twice: a read inside them matches both places, so it is no evidence,
        // and the clone is incomplete (not failed: the read does match the reference).
        let unit = bases(300, 19);
        let twice = [&bases(100, 20)[..], &unit, &bases(100, 21), &unit, &bases(100, 22)].concat();
        let mut j = job(&twice, (100, 400));
        j.add("f", &unit[20..280], Some(&q(260, 40)));
        assert_eq!((j.reads()[0].placed, j.reads()[0].ambiguous, j.reads()[0].reason.as_str(), j.report().verdict), (false, true, "it matches more than one place on the reference", "incomplete"));
        // One base different in the second copy: the read's own copy scores 5 higher, so it is placed there.
        let mut near = twice.clone();
        near[600] = other(near[600]);
        let mut j = job(&near, (150, 350));
        j.add("f", &unit[20..280], Some(&q(260, 40)));
        assert_eq!((j.reads()[0].placed, j.reads()[0].ref_start, j.report().verdict), (true, 120, "pass"));
    }

    #[test]
    fn tandem_repeats_two_turns_of_a_circle_and_unusable_reads_are_not_evidence() {
        // A read inside tandem copies matches at every copy.
        let unit = bases(60, 27);
        let tandem = [&bases(200, 28)[..], &unit, &unit, &unit, &unit, &unit, &bases(200, 29)].concat();
        let mut j = job(&tandem, (200, 500));
        j.add("f", &[&unit[..], &unit[..]].concat(), Some(&q(120, 40)));
        assert_eq!((j.reads()[0].ambiguous, j.report().verdict), (true, "incomplete"));
        // A consensus holding the plasmid twice runs two turns: no evidence, so no pass.
        let r = bases(300, 30);
        let mut j = Job::new(&r, true, (0, 300), 20, 0.05).unwrap();
        j.add("twice", &[&r[..], &r[..]].concat(), None);
        assert_eq!((j.reads()[0].placed, j.reads()[0].reason.as_str()), (false, "it runs more than once around the circular reference"));
        // The plasmid once, from another origin: one placement, not a tie with its copy one turn on.
        let mut j = Job::new(&r, true, (0, 300), 20, 0.05).unwrap();
        j.add("once", &[&r[100..], &r[..100]].concat(), None);
        assert_eq!((j.reads()[0].placed, j.report().verdict), (true, "agrees (no quality data)"));
        // An exact read at Q0 trims to nothing: incomplete (no usable evidence), not fail.
        let mut j = job(&r, (50, 250));
        j.add("f", &r[..300], Some(&q(300, 0)));
        assert_eq!(j.report().verdict, "incomplete");
    }

    #[test]
    fn a_read_past_a_linear_end_is_left_unaligned_and_the_reference_must_be_unambiguous() {
        let r = bases(400, 23);
        // 60 bases before the reference's first: not aligned and not evidence, the rest passes.
        let read = [&bases(60, 24)[..], &r[..300]].concat();
        let mut j = job(&r, (0, 250));
        j.add("f", &read, Some(&q(read.len(), 40)));
        assert_eq!((j.report().verdict, j.reads()[0].read_start, j.reads()[0].ref_start), ("pass", 60, 0));
        // An N in the reference region can't be verified.
        let mut nr = r.clone();
        nr[100] = b'N';
        let mut j = job(&nr, (50, 250));
        j.add("f", &r[..300], Some(&q(300, 40)));
        assert_eq!((j.report().verdict, j.report().uncovered), ("incomplete", vec![vec![100]]));
    }

    #[test]
    fn an_empty_region_or_one_outside_is_refused() {
        let r = bases(100, 25);
        assert!(Job::new(&r, false, (10, 10), 20, 0.05).is_err());
        assert!(Job::new(&r, false, (100, 100), 20, 0.05).is_err());
        assert!(Job::new(&r, true, (100, 5), 20, 0.05).is_err());
        assert!(Job::new(&r, true, (10, 10), 20, 0.05).is_err());
        // The whole circle is [0, n); a region across the origin wraps.
        assert_eq!((Job::new(&r, true, (0, 100), 20, 0.05).unwrap().region.len(), Job::new(&r, true, (90, 10), 20, 0.05).unwrap().region.len()), (100, 20));
    }
}
