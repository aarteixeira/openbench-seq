//! Clone verification: reads (Sanger traces, a consensus) aligned to a reference, evidence per base of a target region,
//! and a verdict (OpenBench plan, Sprint 4, decisions 11 and 12; 4B-1 plan revision 4).
//!
//! Each read is trimmed (Mott, `align::trim`) when it has qualities; a letter other than A, C, G, T becomes N (a no-call)
//! before its strands are made, and bases and qualities are reversed together for the reverse strand. Both strands
//! are aligned: glocal on a circle (the reference repeated past one turn and a half of the read, positions modulo its
//! length), overlap on a linear reference (bases past its ends are not aligned). A read is *placed* when at least 50 of
//! its bases align and its score reaches 0.4 per aligned base; *ambiguous* (no evidence) when another placement scores
//! as high: the other strand over another span, another end (the aligner's ties) or another start (the ties of the
//! reversed alignment); positions one turn apart are one placement. An alignment refused for its size refuses the
//! whole clone (`add` returns the reason): no verdict rests on a calculation that was not made.
//!
//! On a circle a read gives evidence for one turn at most: a read past a turn is no evidence unless the caller marks
//! it closed (a whole-plasmid consensus). A closed read is aligned again end to end against one turn of the reference
//! rotated to its placement, and its end joins its start (when trimming kept both), so it observes every junction.
//!
//! Observations: per reference position a base or a deletion; per junction between two consecutive positions a read
//! spans, what it inserts there (possibly nothing). Quality: a base's Phred; a deletion's, the lower Phred of the read
//! bases aligned on either side; a junction's, the lowest of its inserted bases and of those two flanking bases. A call
//! is A, C, G, T or a deletion (an insertion: only those letters, with called flanks). With qualities in any read of the
//! clone a call *qualifies* at Phred >= `min_q`; without, every call qualifies and the verdict says "no quality data".
//!
//! Per region position: a *conflict* when two calls of placed reads differ, at any quality; else *uncovered* without a
//! qualifying call or when the reference letter is not A, C, G or T; else its call. Per junction between consecutive
//! region positions (and the closing junction of a whole circle): a conflict, *unsupported* without a qualifying
//! observation, or its call. Verdict: `fail` when a read was usable and none is placed or ambiguous; `incomplete` when
//! anything is uncovered, unsupported or in conflict, or no read was usable; `mutation` when a call differs from the
//! reference; else `pass` (without qualities: `agrees (no quality data)`, `mutation (no quality data)`). A deletion
//! reaching the region's edge is reported whole, with its called positions outside the region.

use crate::align::{self, Mode, Op, Scoring};
use crate::alphabet::{complement, letter, mask};
use crate::json;
use std::collections::BTreeMap;

/// A read is placed when its score is at least this much per aligned read base.
pub const MIN_PER_BASE: f64 = 0.4;
/// Reads with fewer bases than this after trimming, or aligned to the reference, are not placed.
pub const MIN_BASES: usize = 50;
/// Reads one clone takes.
pub const MAX_READS: usize = 96;
/// Conflicts listed with each read's call, and differences listed, at most (the totals and runs are complete).
pub const MAX_LISTED: usize = 1000;
pub const MAX_DIFFERENCES: usize = 10_000;

pub struct Job {
    reference: Vec<u8>,
    circular: bool,
    /// Region positions, in order (a region across the origin of a circle wraps).
    region: Vec<usize>,
    min_q: u8,
    limit: f64,
    any_quality: bool,
    /// Reads with at least MIN_BASES after trimming, and those that match the reference but give no evidence (past one
    /// turn of a circle, not closed).
    usable: usize,
    matched_unused: usize,
    reads: Vec<ReadResult>,
    /// Per reference position: what each placed read shows there (a base letter or '-').
    bases: BTreeMap<usize, Vec<Obs<u8>>>,
    /// Per reference position: what each placed read spanning the junction after it inserts there ("" for nothing).
    after: BTreeMap<usize, Vec<Obs<String>>>,
}

/// One read's observation: which read, what it shows, its quality (None without qualities), whether that is a call,
/// and whether it qualifies (a call with a quality of at least `min_q`).
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
    /// Matches the reference, but as well at another place: no evidence.
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

fn is_base(b: u8) -> bool {
    matches!(b, b'A' | b'C' | b'G' | b'T')
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
        Ok(Job { reference: reference.to_ascii_uppercase(), circular, region, min_q, limit, any_quality: false, usable: 0, matched_unused: 0, reads: Vec::new(), bases: BTreeMap::new(), after: BTreeMap::new() })
    }

    /// Trim, place and read one read; a placed read's observations join the evidence. `closed`: the read is a
    /// whole-plasmid consensus, its end joined to its start. Err: an alignment was refused (too long), so the clone has
    /// no verdict.
    pub fn add(&mut self, name: &str, bases: &[u8], quality: Option<&[u8]>, closed: bool) -> Result<(), String> {
        if self.reads.len() >= MAX_READS {
            return Err(format!("a clone takes at most {MAX_READS} reads"));
        }
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
            return Ok(());
        }
        self.usable += 1;
        // Anything but A, C, G, T is N before the strands are made, so a letter is a no-call on both.
        let fwd: Vec<u8> = bases[s..e].iter().map(|&b| { let b = b.to_ascii_uppercase(); if is_base(b) { b } else { b'N' } }).collect();
        let rev: Vec<u8> = fwd.iter().rev().map(|&b| letter(complement(mask(b)))).collect();
        let len = fwd.len();
        // On a circle the reference repeated to one turn plus twice the read (a read's span is under twice its length at
        // the placement score, so every qualifying placement, a deletion across the origin included, fits).
        let ext: Vec<u8> = if self.circular {
            let want = n + 2 * len + 100;
            self.reference.iter().copied().cycle().take(want.max(n)).collect()
        } else {
            self.reference.clone()
        };
        let sc = Scoring::DNA;
        let mode = if self.circular { Mode::Glocal } else { Mode::Overlap };
        let mut tries: Vec<(i32, align::Alignment)> = Vec::new();
        for strand in [1, -1] {
            tries.push((strand, align::align(&ext, if strand == 1 { &fwd } else { &rev }, mode, &sc, None)?));
        }
        let best = if tries[1].1.score > tries[0].1.score { 1 } else { 0 };
        let (strand, mut a) = tries[best].clone();
        let other = &tries[1 - best].1;
        let read: &[u8] = if strand == 1 { &fwd } else { &rev };
        r.strand = strand;
        r.score = a.score;
        r.identity = a.identity();
        let aligned = a.b.1 - a.b.0;
        if aligned < MIN_BASES || (a.score as f64) < MIN_PER_BASE * aligned as f64 {
            r.reason = format!("it does not match the reference (fewer than {MIN_BASES} bases aligned, or a score below 0.4 per aligned base)");
            self.reads.push(r);
            return Ok(());
        }
        // Another placement as good: positions one turn apart are one placement. The other strand as good is another
        // reading of the evidence unless the read is its own reverse complement over the same span (the same bases).
        let key = |i: usize| if self.circular { i % n } else { i };
        let span = |x: &align::Alignment| (key(x.a.0), key(x.a.1));
        let palindrome = fwd == rev && span(other) == span(&a);
        let mut ambiguous = other.score >= a.score && !palindrome;
        ambiguous |= a.ties.iter().any(|&(i, j)| (key(i), j) != (key(a.a.1), a.b.1));
        if !ambiguous {
            // Other starts reaching the same end: the ends of the reversed alignment.
            let (ra, rb): (Vec<u8>, Vec<u8>) = (ext.iter().rev().copied().collect(), read.iter().rev().copied().collect());
            let (_, ends) = align::best_ends(&ra, &rb, mode, &sc)?;
            let starts: std::collections::BTreeSet<(usize, usize)> = ends.iter().map(|&(i, j)| (key(ext.len() - i), len - j)).collect();
            ambiguous = starts.len() > 1;
        }
        if ambiguous {
            r.ambiguous = true;
            r.reason = "it matches more than one place on the reference".into();
            self.reads.push(r);
            return Ok(());
        }
        // A read past one turn that is not closed (a concatemer, overlapping ends) is no evidence. A closed read (a
        // whole-plasmid consensus) is the circle: it is aligned again end to end against one turn of the reference,
        // rotated to where it was placed, so a difference where its ends meet (a deletion, an insertion) is seen too.
        if self.circular && a.a.1 - a.a.0 > n && !closed {
            self.matched_unused += 1;
            r.reason = "it covers the circle more than once".into();
            self.reads.push(r);
            return Ok(());
        }
        if self.circular && closed {
            let s = a.a.0 % n;
            let turn: Vec<u8> = [&self.reference[s..], &self.reference[..s]].concat();
            let mut g = align::align(&turn, read, Mode::Global, &sc, None)?;
            g.a = (g.a.0 + s, g.a.1 + s);
            a = g;
            r.identity = a.identity();
        }
        r.placed = true;
        r.ref_start = a.a.0 % n;
        r.read_start = a.b.0;
        // Qualities in the read's aligned order (reversed with the bases on the reverse strand).
        let q: Option<Vec<u8>> = quality.map(|q| {
            let q = q[s..e].to_vec();
            if strand == 1 { q } else { q.into_iter().rev().collect() }
        });
        let qual = |j: usize| -> u8 { q.as_ref().map_or(0, |q| q[j]) };
        let strong = |p: u8| q.is_some() && p >= self.min_q;
        let at = self.reads.len();
        // A junction's observation from its flanking read bases (`left` and `right`, when the read has them) and what
        // is inserted between them.
        let junction = |ins: String, iq: u8, left: Option<usize>, right: Option<usize>| -> Obs<String> {
            let flanks = [left, right].into_iter().flatten();
            let call = ins.bytes().all(is_base) && flanks.clone().all(|f| is_base(read[f]));
            let around = flanks.map(qual).fold(iq, u8::min);
            Obs { read: at, allele: ins, q: q.as_ref().map(|_| around), call, strong: call && strong(around) }
        };
        let (mut i, mut j) = (a.a.0, a.b.0);
        // A closed read whose ends trimming kept is a circle: past its last base comes its first, before its first its
        // last (the flanks of a deletion where its ends meet).
        let joined = closed && self.circular && r.used == (0, bases.len());
        // The bases inserted since the last aligned reference position, and their lowest Phred.
        let mut pending: Option<(String, u8)> = None;
        // The last reference position, the last read base aligned to one, the first such read base, and the bases
        // before the first position.
        let (mut last, mut left, mut first, mut lead): (Option<usize>, Option<usize>, Option<usize>, Option<(String, u8)>) = (None, None, None, None);
        for &(op, c) in &a.ops {
            for _ in 0..c {
                match op {
                    Op::Match | Op::Mismatch | Op::Del => {
                        let pos = i % n;
                        // The read base after the junction: this one, or for a deletion the next one the read has.
                        let right = if j < read.len() { Some(j) } else { joined.then_some(0) };
                        let before = if left.is_none() && joined { Some(read.len() - 1) } else { left };
                        match last {
                            None => lead = pending.take(),
                            Some(prev) => {
                                let (ins, iq) = pending.take().unwrap_or((String::new(), u8::MAX));
                                self.after.entry(prev).or_default().push(junction(ins, iq, before, right));
                            }
                        }
                        if op == Op::Del {
                            // A deletion is a call only between two called bases of the read.
                            let call = matches!((before, right), (Some(l), Some(r)) if is_base(read[l]) && is_base(read[r]));
                            let around = [before, right].into_iter().flatten().map(qual).fold(u8::MAX, u8::min);
                            let around = if q.is_some() { around } else { 0 };
                            self.bases.entry(pos).or_default().push(Obs { read: at, allele: b'-', q: q.as_ref().map(|_| around), call, strong: call && strong(around) });
                        } else {
                            let call = is_base(read[j]);
                            self.bases.entry(pos).or_default().push(Obs { read: at, allele: read[j], q: q.as_ref().map(|_| qual(j)), call, strong: call && strong(qual(j)) });
                            first.get_or_insert(j);
                            left = Some(j);
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
        // A closed read covering the circle once joins its end to its start: the junction between its last position
        // and its first holds what it has after the one and before the other.
        // Only when trimming kept both of the read's ends: otherwise its trimmed ends never met.
        if let (true, Some(prev)) = (joined, last) {
            let (tail, tq) = pending.take().unwrap_or((String::new(), u8::MAX));
            let (head, hq) = lead.take().unwrap_or((String::new(), u8::MAX));
            self.after.entry(prev).or_default().push(junction(tail + &head, tq.min(hq), left, first));
        }
        r.ops = a.ops;
        self.reads.push(r);
        Ok(())
    }

    pub fn reads(&self) -> &[ReadResult] {
        &self.reads
    }

    pub fn report(&self) -> Report {
        let n = self.reference.len();
        let quality = self.any_quality;
        // Without qualities anywhere, every call qualifies.
        let qualifies = |call: bool, strong: bool| if quality { strong } else { call };
        // A position's call, and a junction's (after a position): Some when the calls agree and one qualifies.
        let call_at = |pos: usize| -> Option<u8> {
            let obs = self.bases.get(&pos)?;
            let calls: Vec<u8> = obs.iter().filter(|o| o.call).map(|o| o.allele).collect();
            (calls.iter().all(|&b| b == calls[0]) && obs.iter().any(|o| qualifies(o.call, o.strong))).then(|| calls[0])
        };
        let junction_at = |pos: usize| -> Option<&String> {
            let obs = self.after.get(&pos)?;
            let calls: Vec<&String> = obs.iter().filter(|o| o.call).map(|o| &o.allele).collect();
            (calls.iter().all(|s| *s == calls[0]) && obs.iter().any(|o| qualifies(o.call, o.strong))).then(|| calls[0])
        };
        let (mut uncovered, mut conflicts, mut unsupported): (Vec<Vec<usize>>, Vec<Conflict>, Vec<usize>) = (Vec::new(), Vec::new(), Vec::new());
        let (mut covered, mut single) = (0usize, 0usize);
        let whole = self.circular && self.region.len() == n;
        // Region membership, and each position's place in the region (the order differences are listed in).
        let mut place = vec![usize::MAX; n];
        for (t, &pos) in self.region.iter().enumerate() {
            place[pos] = t;
        }
        for (t, &pos) in self.region.iter().enumerate() {
            let obs = self.bases.get(&pos).map(|v| &v[..]).unwrap_or(&[]);
            let calls: Vec<u8> = obs.iter().filter(|o| o.call).map(|o| o.allele).collect();
            let good = obs.iter().filter(|o| qualifies(o.call, o.strong)).count();
            if calls.iter().any(|&b| b != calls[0]) {
                let calls = obs.iter().filter(|o| o.call).map(|o| (o.read, (o.allele as char).to_string(), o.q)).collect();
                conflicts.push(Conflict { pos, junction: false, calls });
            } else if good == 0 || !is_base(self.reference[pos]) {
                match uncovered.last_mut() {
                    Some(run) if t > 0 && run.last() == Some(&self.region[t - 1]) => run.push(pos),
                    _ => uncovered.push(vec![pos]),
                }
            } else {
                covered += 1;
                single += (good == 1) as usize;
            }
            // The junction after this base: when the next region base follows it, and the closing one of a whole circle.
            if t + 1 < self.region.len() || whole {
                let obs = self.after.get(&pos).map(|v| &v[..]).unwrap_or(&[]);
                let calls: Vec<&String> = obs.iter().filter(|o| o.call).map(|o| &o.allele).collect();
                if calls.iter().any(|s| *s != calls[0]) {
                    let calls = obs.iter().filter(|o| o.call).map(|o| (o.read, o.allele.clone(), o.q)).collect();
                    conflicts.push(Conflict { pos, junction: true, calls });
                } else if !obs.iter().any(|o| qualifies(o.call, o.strong)) {
                    unsupported.push(pos);
                }
            }
        }
        // Differences are built once over the whole reference from the agreed calls (a deletion is one event however
        // the region cuts it, across the origin of a circle too), then kept when they touch the region.
        let mut differences: Vec<(usize, Diff)> = Vec::new();
        let positions: Vec<usize> = self.bases.keys().copied().collect();
        let mut runs: Vec<Vec<usize>> = Vec::new();
        for &pos in &positions {
            let reference = self.reference[pos];
            match call_at(pos) {
                Some(b'-') => match runs.last_mut() {
                    Some(run) if *run.last().unwrap() + 1 == pos => run.push(pos),
                    _ => runs.push(vec![pos]),
                },
                Some(b) if is_base(reference) && b != reference && place[pos] != usize::MAX => differences.push((place[pos], Diff::Sub { pos, reference, read: b })),
                _ => {}
            }
        }
        // On a circle a run ending at the last base and one starting at the first are one deletion.
        if self.circular && runs.len() > 1 && runs[0][0] == 0 && *runs.last().unwrap().last().unwrap() == n - 1 {
            let head = runs.remove(0);
            runs.last_mut().unwrap().extend(head);
        }
        for run in runs {
            if let Some(first) = run.iter().map(|&p| place[p]).filter(|&t| t != usize::MAX).min() {
                differences.push((first, Diff::Del { positions: run }));
            }
        }
        for (&pos, _) in &self.after {
            let next = (pos + 1) % n;
            let inside = place[pos] != usize::MAX && place[next] != usize::MAX && (place[next] == place[pos] + 1 || (whole && next == self.region[0]));
            if let (true, Some(bases)) = (inside, junction_at(pos)) {
                if !bases.is_empty() {
                    differences.push((place[pos], Diff::Ins { after: pos, bases: bases.clone() }));
                }
            }
        }
        // In region order; at one place a deletion or substitution before the insertion after it.
        differences.sort_by_key(|(t, d)| (*t, matches!(d, Diff::Ins { .. })));
        let differences: Vec<Diff> = differences.into_iter().map(|(_, d)| d).collect();
        // Failed only on usable evidence: a clone whose reads all match more than one place, or have too few good
        // bases to judge, is incomplete.
        let matched = self.matched_unused > 0 || self.reads.iter().any(|r| r.placed || r.ambiguous);
        let verdict = if !matched && self.usable > 0 {
            "fail"
        } else if !uncovered.is_empty() || !conflicts.is_empty() || !unsupported.is_empty() || self.usable == 0 {
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

/// Runs of consecutive positions (in the order given) as first and last.
fn runs(positions: impl Iterator<Item = usize>) -> Vec<(usize, usize)> {
    let mut out: Vec<(usize, usize)> = Vec::new();
    for p in positions {
        match out.last_mut() {
            Some((_, b)) if *b + 1 == p => *b = p,
            _ => out.push((p, p)),
        }
    }
    out
}

/// The report and the reads as JSON for the page: the differences (at most MAX_DIFFERENCES, with the total), uncovered
/// runs, unsupported junctions as runs, every conflict's position as runs with the first MAX_LISTED in detail.
pub fn to_json(job: &Job, rep: &Report) -> String {
    let mut o = String::from("{\"verdict\":");
    json::string(&mut o, rep.verdict);
    o.push_str(&format!(",\"quality\":{},\"covered\":{},\"single\":{},\"differenceTotal\":{},\"differences\":[", rep.quality, rep.covered, rep.single, rep.differences.len()));
    for (k, d) in rep.differences.iter().take(MAX_DIFFERENCES).enumerate() {
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
    let pairs = |o: &mut String, v: &[(usize, usize)]| {
        o.push('[');
        for (k, &(a, b)) in v.iter().enumerate() {
            if k > 0 {
                o.push(',');
            }
            json::numbers(o, [a, b]);
        }
        o.push(']');
    };
    o.push_str("],\"uncovered\":");
    pairs(&mut o, &rep.uncovered.iter().map(|r| (r[0], r[r.len() - 1])).collect::<Vec<_>>());
    o.push_str(",\"unsupported\":");
    pairs(&mut o, &runs(rep.unsupported.iter().copied()));
    o.push_str(&format!(",\"conflictTotal\":{},\"conflictRuns\":", rep.conflicts.len()));
    pairs(&mut o, &runs(rep.conflicts.iter().filter(|c| !c.junction).map(|c| c.pos)));
    o.push_str(",\"conflicts\":[");
    for (k, c) in rep.conflicts.iter().take(MAX_LISTED).enumerate() {
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
    o.push_str("],\"reads\":[");
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
    /// The reverse complement, keeping a U as U (a read as a sequencer could write it).
    fn rc_keep_u(s: &[u8]) -> Vec<u8> {
        s.iter().rev().map(|&b| if b == b'U' { b'U' } else { letter(complement(mask(b))) }).collect()
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
        j.add("f", &r[20..380], Some(&q(360, 40)), false).unwrap();
        let rep = j.report();
        assert_eq!((rep.verdict, rep.single), ("pass", 300));
        let mut read = r[20..380].to_vec();
        read[180] = other(read[180]);
        let mut j = job(&r, (50, 350));
        j.add("f", &read, Some(&q(360, 40)), false).unwrap();
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
        j.add("r", &read, Some(&qual), false).unwrap();
        let rep = j.report();
        assert_eq!((j.reads()[0].strand, j.reads()[0].used, rep.verdict, rep.uncovered.clone()), (-1, (15, read.len()), "incomplete", vec![vec![259]]));
        // A region the read does not reach: incomplete, with the uncovered run.
        let mut j = job(&r, (10, 340));
        j.add("f", &r[40..380], Some(&q(340, 40)), false).unwrap();
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
        j.add("a", &r[20..380], Some(&q(360, 30)), false).unwrap();
        j.add("b", &b, Some(&qual), false).unwrap();
        let rep = j.report();
        assert_eq!(rep.verdict, "incomplete");
        // Each read's call with its Phred, so the conflict can be judged.
        assert_eq!((rep.conflicts.len(), rep.conflicts[0].pos, rep.conflicts[0].junction), (1, 120, false));
        assert_eq!(rep.conflicts[0].calls, vec![(0, (r[120] as char).to_string(), Some(30)), (1, (b[100] as char).to_string(), Some(15))]);
        // An N (and another IUPAC code) neither conflicts nor counts.
        let mut n = r[20..380].to_vec();
        (n[100], n[101]) = (b'N', b'R');
        let mut j = job(&r, (50, 350));
        j.add("a", &r[20..380], Some(&q(360, 40)), false).unwrap();
        j.add("n", &n, Some(&q(360, 40)), false).unwrap();
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
        j.add("f", &read, Some(&q(read.len(), 40)), false).unwrap();
        let rep = j.report();
        assert_eq!(rep.verdict, "mutation");
        assert!(rep.differences.iter().any(|d| matches!(d, Diff::Del { positions } if positions.len() == 3)), "{:?}", rep.differences);
        assert!(rep.differences.iter().any(|d| matches!(d, Diff::Ins { bases, .. } if bases.len() == 6)), "{:?}", rep.differences);
        // A second read without the insertion: the junction is a conflict, whatever else agrees.
        let mut j = job(&r, (50, 350));
        j.add("f", &read, Some(&q(read.len(), 40)), false).unwrap();
        j.add("g", &r[20..380], Some(&q(360, 40)), false).unwrap();
        let rep = j.report();
        assert_eq!((rep.verdict, rep.conflicts.iter().filter(|c| c.junction).count()), ("incomplete", 1));
        // The insertion's bases at Q10: the junction has no qualifying evidence.
        let mut qual = q(read.len(), 40);
        qual[250..256].fill(10);
        let mut j = job(&r, (50, 350));
        j.add("f", &read, Some(&qual), false).unwrap();
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
        j.add("f", &read, Some(&q(read.len(), 40)), false).unwrap();
        j.add("r", &rc(&read), Some(&q(read.len(), 40)), false).unwrap();
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
        j.add("f", &read, Some(&q(read.len(), 40)), false).unwrap();
        assert_eq!(j.report().verdict, "pass");
        let mut j = job(&r, (50, 350));
        j.add("f", &read, Some(&q(read.len(), 40)), false).unwrap();
        assert_eq!(j.report().differences, vec![Diff::Del { positions: vec![204] }]);
    }

    #[test]
    fn a_circular_reference_and_a_region_across_its_origin() {
        let r = bases(300, 15);
        // A read over the origin: the last 80 bases, then the first 120.
        let read = [&r[220..], &r[..120]].concat();
        let mut j = Job::new(&r, true, (250, 100), 20, 0.05).unwrap();
        j.add("f", &read, Some(&q(read.len(), 40)), false).unwrap();
        let rep = j.report();
        assert_eq!((rep.verdict, j.reads()[0].ref_start), ("pass", 220));
    }

    #[test]
    fn evidence_without_qualities_never_passes_and_a_foreign_read_fails() {
        let r = bases(400, 16);
        let mut j = job(&r, (0, 400));
        j.add("consensus", &r, None, false).unwrap();
        assert_eq!(j.report().verdict, "agrees (no quality data)");
        let mut c = r.clone();
        c[10] = other(c[10]);
        let mut j = job(&r, (0, 400));
        j.add("consensus", &c, None, false).unwrap();
        assert_eq!(j.report().verdict, "mutation (no quality data)");
        let mut j = job(&r, (50, 350));
        j.add("other", &bases(300, 99), Some(&q(300, 40)), false).unwrap();
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
        j.add("f", &read, Some(&q(400, 40)), false).unwrap();
        assert!(j.reads()[0].placed, "{}", j.reads()[0].reason);
        // A reference holding the same 300 bases twice: a read inside them matches both places, so it is no evidence,
        // and the clone is incomplete (not failed: the read does match the reference).
        let unit = bases(300, 19);
        let twice = [&bases(100, 20)[..], &unit, &bases(100, 21), &unit, &bases(100, 22)].concat();
        let mut j = job(&twice, (100, 400));
        j.add("f", &unit[20..280], Some(&q(260, 40)), false).unwrap();
        assert_eq!((j.reads()[0].placed, j.reads()[0].ambiguous, j.reads()[0].reason.as_str(), j.report().verdict), (false, true, "it matches more than one place on the reference", "incomplete"));
        // One base different in the second copy: the read's own copy scores 5 higher, so it is placed there.
        let mut near = twice.clone();
        near[600] = other(near[600]);
        let mut j = job(&near, (150, 350));
        j.add("f", &unit[20..280], Some(&q(260, 40)), false).unwrap();
        assert_eq!((j.reads()[0].placed, j.reads()[0].ref_start, j.report().verdict), (true, 120, "pass"));
    }

    #[test]
    fn tandem_repeats_turns_of_a_circle_and_unusable_reads() {
        // A read inside tandem copies matches at every copy.
        let unit = bases(60, 27);
        let tandem = [&bases(200, 28)[..], &unit, &unit, &unit, &unit, &unit, &bases(200, 29)].concat();
        let mut j = job(&tandem, (200, 500));
        j.add("f", &[&unit[..], &unit[..]].concat(), Some(&q(120, 40)), false).unwrap();
        assert_eq!((j.reads()[0].ambiguous, j.report().verdict), (true, "incomplete"));
        // A closed consensus holding the plasmid twice: one turn, then a whole copy inserted where it closes.
        let r = bases(300, 30);
        let mut j = Job::new(&r, true, (0, 300), 20, 0.05).unwrap();
        j.add("twice", &[&r[..], &r[..]].concat(), None, true).unwrap();
        let rep = j.report();
        assert_eq!((j.reads()[0].placed, rep.verdict), (true, "mutation (no quality data)"));
        assert!(matches!(&rep.differences[..], [Diff::Ins { bases, .. }] if bases.len() == 300), "{:?}", rep.differences);
        // The plasmid once, closed, unrotated and rotated: one placement (not a tie with its copy one turn on), and
        // the junction where it closes is observed.
        for read in [r.clone(), [&r[100..], &r[..100]].concat()] {
            let mut j = Job::new(&r, true, (0, 300), 20, 0.05).unwrap();
            j.add("once", &read, None, true).unwrap();
            assert_eq!((j.reads()[0].placed, j.report().verdict), (true, "agrees (no quality data)"));
        }
        // Not marked closed, the same read does not observe that junction.
        let mut j = Job::new(&r, true, (0, 300), 20, 0.05).unwrap();
        j.add("once", &[&r[100..], &r[..100]].concat(), None, false).unwrap();
        assert_eq!((j.report().verdict, j.report().unsupported), ("incomplete", vec![99]));
        // An exact read at Q0 trims to nothing: incomplete (no usable evidence), not fail.
        let mut j = job(&r, (50, 250));
        j.add("f", &r[..300], Some(&q(300, 0)), false).unwrap();
        assert_eq!(j.report().verdict, "incomplete");
    }

    #[test]
    fn a_read_past_one_turn_that_is_not_closed_is_no_evidence() {
        // A linear read of a 300-base circle running 1.2 turns (a concatemer): no evidence, so incomplete (not failed);
        // the closed reading of the same bases is a 360-base molecule, a 60-base insertion.
        let r = bases(300, 40);
        let read = [&r[..], &r[..60]].concat();
        let mut j = Job::new(&r, true, (0, 300), 20, 0.05).unwrap();
        j.add("raw", &read, Some(&q(read.len(), 40)), false).unwrap();
        assert_eq!((j.report().verdict, j.reads()[0].reason.as_str()), ("incomplete", "it covers the circle more than once"));
        let mut j = Job::new(&r, true, (0, 300), 20, 0.05).unwrap();
        j.add("consensus", &read, Some(&q(read.len(), 40)), true).unwrap();
        assert!(matches!(&j.report().differences[..], [Diff::Ins { bases, .. }] if bases.len() == 60));
    }

    #[test]
    fn p4_cases_trimmed_closure_mutated_palindrome_deletion_across_origin_and_called_anchors() {
        // A closed consensus G + reference, the G at Q10: trimming removes it, so the read's ends no longer meet and the
        // closing junction is unsupported (not a pass).
        let r = bases(300, 41);
        let read = [&b"G"[..], &r[..]].concat();
        let mut qual = q(read.len(), 40);
        qual[0] = 10;
        let mut j = Job::new(&r, true, (0, 300), 20, 0.05).unwrap();
        j.add("c", &read, Some(&qual), true).unwrap();
        assert_eq!((j.report().verdict, j.report().unsupported), ("incomplete", vec![299]));
        // A reference palindrome and a read of it with one base changed: forward and reverse readings put the change
        // at mirrored positions, as good: ambiguous.
        let half = bases(60, 42);
        let pal = [&half[..], &rc(&half)].concat();
        let ref2 = [&bases(150, 43)[..], &pal, &bases(150, 44)].concat();
        let mut read = pal.clone();
        read[10] = other(read[10]);
        let mut j = job(&ref2, (150, 270));
        j.add("p", &read, Some(&q(120, 40)), false).unwrap();
        assert!(j.reads()[0].ambiguous);
        // One deletion of 298, 299, 0, 1 on a 300-base circle: one event whether the region is the whole circle or
        // stops short of it.
        let r = bases(300, 45);
        let read = [&r[2..298]].concat();
        for region in [(0, 300), (1, 299)] {
            let mut j = Job::new(&r, true, region, 20, 0.05).unwrap();
            j.add("c", &read, None, true).unwrap();
            assert_eq!(j.report().differences, vec![Diff::Del { positions: vec![298, 299, 0, 1] }], "{region:?}");
        }
        // A deletion beside an N in the read is no call.
        let r = bases(400, 46);
        let mut read: Vec<u8> = [&r[20..200], &r[201..380]].concat();
        read[179] = b'N'; // the base before the deletion of 200
        let mut j = job(&r, (200, 201));
        j.add("f", &read, Some(&q(read.len(), 40)), false).unwrap();
        assert_eq!((j.report().verdict, j.report().differences.len()), ("incomplete", 0));
    }

    #[test]
    fn the_closing_junction_of_a_whole_circle_is_checked() {
        // Reads covering every position and every inner junction, one with GG inserted after the last base.
        let r = bases(300, 31);
        let mut j = Job::new(&r, true, (0, 300), 20, 0.05).unwrap();
        let ins = [&r[200..], b"GG", &r[..150]].concat();
        j.add("a", &ins, Some(&q(ins.len(), 40)), false).unwrap();
        j.add("b", &r[100..250], Some(&q(150, 40)), false).unwrap();
        let rep = j.report();
        assert_eq!((rep.verdict, rep.differences.clone()), ("mutation", vec![Diff::Ins { after: 299, bases: "GG".into() }]));
    }

    #[test]
    fn different_starts_to_one_end_are_ambiguous_and_a_palindrome_is_one_placement() {
        // ACGTAAGCCA then 60 Ts, and ACGTA then 60 Ts: the read fits from 0 (a 5-base deletion) or from 5 equally.
        let mut r = b"ACGTAAGCCA".to_vec();
        r.extend([b'T'; 60]);
        r.extend(bases(200, 32));
        let mut read = b"ACGTA".to_vec();
        read.extend([b'T'; 60]);
        let mut j = job(&r, (0, 65));
        j.add("f", &read, None, false).unwrap();
        assert!(j.reads()[0].ambiguous, "{}", j.reads()[0].reason);
        // A read equal to its own reverse complement matches both strands over one span: placed.
        let half = bases(60, 33);
        let pal = [&half[..], &rc(&half)].concat();
        let ref2 = [&bases(150, 34)[..], &pal, &bases(150, 35)].concat();
        let mut j = job(&ref2, (150, 270));
        j.add("p", &pal, Some(&q(120, 40)), false).unwrap();
        assert_eq!((j.reads()[0].placed, j.report().verdict), (true, "pass"));
    }

    #[test]
    fn insertion_quality_includes_its_left_flank_and_u_is_no_call_on_both_strands() {
        let r = bases(400, 36);
        let mut read = r[20..380].to_vec();
        read.splice(180..180, b"AC".iter().copied()); // after reference 199
        let mut qual = q(read.len(), 40);
        qual[179] = 15; // the base before the insertion
        // A second read covers reference 199 at Q40 but ends there, so 199 is covered and only the junction after it
        // rests on the first read, whose left flank is Q15: unsupported, no insertion called.
        let mut j = job(&r, (50, 350));
        j.add("f", &read, Some(&qual), false).unwrap();
        j.add("g", &r[20..200], Some(&q(180, 40)), false).unwrap();
        let rep = j.report();
        assert_eq!((rep.verdict, rep.uncovered.len(), rep.unsupported.clone(), rep.differences.len()), ("incomplete", 0, vec![199], 0));
        // U in a read is a no-call whichever strand it lands on.
        let mut u = r[20..380].to_vec();
        u[100] = b'U';
        for read in [u.clone(), rc_keep_u(&u)] {
            let mut j = job(&r, (50, 350));
            j.add("u", &read, Some(&q(360, 40)), false).unwrap();
            assert_eq!(j.report().uncovered, vec![vec![120]]);
        }
    }

    #[test]
    fn a_refused_alignment_refuses_the_clone_and_a_deletion_at_the_region_edge_is_whole() {
        // 32,000 bases against 32,000: over a billion cells, refused, so no verdict.
        let big = bases(32_000, 37);
        let mut j = job(&big, (0, 100));
        assert_eq!(j.add("f", &big, None, false), Err("too long to align: over 1 billion cells".into()));
        // Five bases deleted (170-174), the region from 172: reported as the whole deletion.
        let r = bases(400, 38);
        let read: Vec<u8> = [&r[20..170], &r[175..380]].concat();
        let mut j = job(&r, (172, 350));
        j.add("f", &read, Some(&q(read.len(), 40)), false).unwrap();
        assert_eq!(j.report().differences, vec![Diff::Del { positions: vec![170, 171, 172, 173, 174] }]);
    }

    #[test]
    fn placement_boundary_is_a_score_of_0_4_per_aligned_base() {
        // 100 bases, every third changed from the start: 34 changes (score 66 × 2 − 34 × 3 = 30 < 40): not placed;
        // 32 changes (68 × 2 − 32 × 3 = 40): placed.
        let r = bases(600, 39);
        let changed = |k: usize| {
            let mut x = r[200..300].to_vec();
            for p in (0..100).step_by(3).take(k) {
                x[p] = other(x[p]);
            }
            x
        };
        for (k, placed) in [(34, false), (32, true)] {
            let mut j = job(&r, (150, 350));
            j.add("f", &changed(k), None, false).unwrap();
            assert_eq!(j.reads()[0].placed, placed, "{k} changes: score {}", j.reads()[0].score);
        }
    }

    #[test]
    fn a_read_past_a_linear_end_is_left_unaligned_and_the_reference_must_be_unambiguous() {
        let r = bases(400, 23);
        // 60 bases before the reference's first: not aligned and not evidence, the rest passes.
        let read = [&bases(60, 24)[..], &r[..300]].concat();
        let mut j = job(&r, (0, 250));
        j.add("f", &read, Some(&q(read.len(), 40)), false).unwrap();
        assert_eq!((j.report().verdict, j.reads()[0].read_start, j.reads()[0].ref_start), ("pass", 60, 0));
        // An N in the reference region can't be verified.
        let mut nr = r.clone();
        nr[100] = b'N';
        let mut j = job(&nr, (50, 250));
        j.add("f", &r[..300], Some(&q(300, 40)), false).unwrap();
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
