//! Pairwise alignment of `a` (the reference) and `b` (the query) with Gotoh's affine gaps: global, local, glocal
//! (all of `b`, free ends in `a`) or overlap (glocal, with `b`'s ends past `a`'s free too), in a band of diagonals for long inputs, and Mott trimming of base qualities.
//!
//! Scores: DNA match 2, mismatch −3 (a position matches when the bases one may be are a subset of the other's code, the
//! IUPAC rule of `alphabet`); protein BLOSUM62 (NCBI's table, `blosum62.txt`). A gap of length L costs open + extend × L
//! (BLAST's convention): DNA 5 and 2, protein 11 and 1.

use crate::alphabet::{mask, satisfies};
use std::sync::OnceLock;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Mode {
    Global,
    Local,
    Glocal,
    /// Glocal, and b may also run past either end of a for free (a read overhanging a linear reference).
    Overlap,
}

#[derive(Clone, Copy, Debug)]
pub struct Scoring {
    pub protein: bool,
    pub open: i32,
    pub extend: i32,
}

impl Scoring {
    pub const DNA: Self = Self { protein: false, open: 5, extend: 2 };
    pub const PROTEIN: Self = Self { protein: true, open: 11, extend: 1 };

    pub fn pair(&self, x: u8, y: u8) -> i32 {
        if self.protein { blosum62(x, y) } else if same_base(x, y) { 2 } else { -3 }
    }

    /// Whether two letters are identical (a match column; identity counts these): the same letter. Scoring is
    /// separate: A against R scores as a match but is a mismatch column.
    pub fn same(&self, x: u8, y: u8) -> bool {
        x.eq_ignore_ascii_case(&y)
    }

    /// The score of `b` aligned to itself in full: the most any alignment of `b` can reach.
    pub fn best(&self, b: &[u8]) -> i32 {
        b.iter().map(|&y| self.pair(y, y)).sum()
    }
}

fn same_base(x: u8, y: u8) -> bool {
    let (p, q) = (mask(x), mask(y));
    satisfies(p, q) || satisfies(q, p)
}

const LETTERS: &[u8; 24] = b"ARNDCQEGHILKMFPSTWYVBZX*";

fn matrix() -> &'static [[i8; 24]; 24] {
    static M: OnceLock<[[i8; 24]; 24]> = OnceLock::new();
    M.get_or_init(|| {
        let mut m = [[0i8; 24]; 24];
        let rows = include_str!("blosum62.txt").lines().filter(|l| !l.starts_with('#') && !l.starts_with(' '));
        for (r, line) in rows.enumerate() {
            for (c, v) in line.split_whitespace().skip(1).enumerate() {
                m[r][c] = v.parse().expect("blosum62.txt is a 24 by 24 table");
            }
        }
        m
    })
}

fn index(x: u8) -> usize {
    LETTERS.iter().position(|&l| l == x.to_ascii_uppercase()).unwrap_or(22) // anything else reads as X
}

pub fn blosum62(x: u8, y: u8) -> i32 {
    matrix()[index(x)][index(y)] as i32
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Op {
    /// a and b letters aligned, the same.
    Match,
    /// a and b letters aligned, different.
    Mismatch,
    /// A letter of b with no letter of a (inserted in b).
    Ins,
    /// A letter of a with no letter of b (deleted from b).
    Del,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Alignment {
    pub score: i32,
    /// Aligned ranges [start, end) of a and of b.
    pub a: (usize, usize),
    pub b: (usize, usize),
    /// Runs of operations, in order.
    pub ops: Vec<(Op, u32)>,
    /// Other end cells (i in a, j in b) where an alignment reaches the same score: other placements as good (not
    /// computed for local alignments).
    pub ties: Vec<(usize, usize)>,
}

impl Alignment {
    /// Identical letters over the alignment's columns.
    pub fn identity(&self) -> f64 {
        let cols: u32 = self.ops.iter().map(|&(_, n)| n).sum();
        let same: u32 = self.ops.iter().filter(|(o, _)| *o == Op::Match).map(|&(_, n)| n).sum();
        if cols == 0 { 0.0 } else { same as f64 / cols as f64 }
    }
}

/// Inputs of at most this many cells ((n + 1) × (m + 1)) are aligned with a full traceback: one byte a cell, so at most
/// 100 MB, plus six rows of scores.
pub const TRACE_CELLS: usize = 100_000_000;
/// Larger inputs, up to this many cells, first get the exact best score and its end from a pass that keeps scores only
/// (linear memory); past it they are refused.
pub const MAX_WORK: usize = 1_000_000_000;

const NEG: i32 = i32::MIN / 4;

/// Align a and b. `band`, when given (tests), is the range of diagonals d = j − i the alignment may use. Otherwise
/// inputs of at most TRACE_CELLS are aligned in full. Larger ones (up to MAX_WORK) get the exact best score and its end
/// from a score-only pass over every cell; the traceback then runs in a window of a ending there (for global, a band of
/// diagonals around the shared k-mers), doubled until it reaches that exact score, and the input is refused when the
/// window would pass TRACE_CELLS first. So an alignment below the best is never returned. Gaps are then placed rightmost
/// on a (`right_gaps`).
pub fn align(a: &[u8], b: &[u8], mode: Mode, sc: &Scoring, band: Option<(isize, isize)>) -> Result<Alignment, String> {
    let (n, m) = (a.len(), b.len());
    let full = (-(n as isize), m as isize);
    let traced = |a: &[u8], b: &[u8], mode: Mode, band: (isize, isize), report: bool| -> Result<Alignment, String> {
        let f = fill(a, b, mode, sc, band, true, report)?;
        Ok(traceback(a, b, sc, &f))
    };
    let cells = (n + 1).saturating_mul(m + 1);
    let found = if let Some(band) = band {
        traced(a, b, mode, band, true)?
    } else if cells <= TRACE_CELLS {
        traced(a, b, mode, full, true)?
    } else if cells > MAX_WORK {
        return Err("too long to align: over 1 billion cells".into());
    } else {
        let exact = fill(a, b, mode, sc, full, false, true)?;
        let (_, ie, je, _) = exact.end;
        // The traceback runs in a band of diagonals (pad each side) inside a window of a ending at the exact end, both
        // doubled until it reaches the exact score.
        let mut pad = 100 + (m / 20) as isize;
        loop {
            let (got, whole) = if mode == Mode::Global {
                let (dlo, dhi) = diagonals(a, b, sc.protein).unwrap_or((0, 0));
                let d = m as isize - n as isize;
                let (lo, hi) = ((dlo - pad).min(0).min(d), (dhi + pad).max(0).max(d));
                (traced(a, b, mode, (lo, hi), false), lo <= full.0 && hi >= full.1)
            } else {
                // Overlap's free start in b holds only at a's start; inside a it is glocal.
                let i0 = ie.saturating_sub(m + 2 * pad as usize);
                let wmode = if mode == Mode::Overlap && i0 > 0 { Mode::Glocal } else { mode };
                let (wa, wb) = (&a[i0..ie], &b[..je]);
                let end = wb.len() as isize - wa.len() as isize;
                let (lo, hi) = (end - pad, end + pad);
                let got = traced(wa, wb, wmode, (lo, hi), false).map(|mut x| {
                    x.a = (x.a.0 + i0, x.a.1 + i0);
                    x
                });
                (got, i0 == 0 && lo <= -(wa.len() as isize) && hi >= wb.len() as isize)
            };
            match got {
                Ok(x) if x.score == exact.score => break Alignment { ties: exact.ties.clone(), ..x },
                Ok(_) | Err(_) if !whole && got_fits(&got) => pad *= 2,
                _ => return Err("too long to trace the alignment back within 100 million cells".into()),
            }
        }
    };
    Ok(right_gaps(a, b, sc, found))
}

/// The best score of aligning a and b and every end cell reaching it (none listed for local), from a score-only pass.
pub fn best_ends(a: &[u8], b: &[u8], mode: Mode, sc: &Scoring) -> Result<(i32, Vec<(usize, usize)>), String> {
    if (a.len() + 1).saturating_mul(b.len() + 1) > MAX_WORK {
        return Err("too long to align: over 1 billion cells".into());
    }
    let f = fill(a, b, mode, sc, (-(a.len() as isize), b.len() as isize), false, true)?;
    let mut ends = f.ties;
    ends.push((f.end.1, f.end.2));
    Ok((f.score, ends))
}

/// Whether a traceback attempt may be widened: it ran, or found no path in its band (not one refused for size).
fn got_fits(got: &Result<Alignment, String>) -> bool {
    match got {
        Ok(_) => true,
        Err(e) => e == "no alignment fits the band",
    }
}

/// Every gap moved rightmost on a (3′ for DNA, C-terminal for protein): a gap shifts right past an aligned column while the letter it would give
/// up equals the one it takes, which leaves every column's letters, and so the score, unchanged. Gaps that meet merge;
/// end gaps stay at the ends.
pub fn right_gaps(a: &[u8], b: &[u8], sc: &Scoring, found: Alignment) -> Alignment {
    let mut cols: Vec<Op> = found.ops.iter().flat_map(|&(o, c)| std::iter::repeat_n(o, c as usize)).collect();
    let aligned = |o: Op| matches!(o, Op::Match | Op::Mismatch);
    loop {
        let mut moved = false;
        let (mut i, mut j, mut p) = (found.a.0, found.b.0, 0usize);
        while p < cols.len() {
            let o = cols[p];
            if aligned(o) {
                i += 1;
                j += 1;
                p += 1;
                continue;
            }
            let len = cols[p..].iter().take_while(|&&c| c == o).count();
            // Shift the run [p, p + len) right while the column after it is aligned and the letters allow; a gap at
            // the start of the alignment is an end gap and stays there.
            while p > 0 && p + len < cols.len() && aligned(cols[p + len]) {
                let same = if o == Op::Del { a[i].eq_ignore_ascii_case(&a[i + len]) } else { b[j].eq_ignore_ascii_case(&b[j + len]) };
                if !same {
                    break;
                }
                cols[p] = cols[p + len];
                cols[p + len] = o;
                i += 1;
                j += 1;
                p += 1;
                moved = true;
            }
            if o == Op::Del { i += len } else { j += len }
            p += len;
        }
        if !moved {
            break;
        }
    }
    // Columns again as runs, each aligned column's kind from its letters, and the score of those columns.
    let (mut i, mut j) = (found.a.0, found.b.0);
    let mut ops: Vec<(Op, u32)> = Vec::new();
    let mut score = 0i32;
    for &o in &cols {
        let o = match o {
            Op::Match | Op::Mismatch => {
                score += sc.pair(a[i], b[j]);
                let same = sc.same(a[i], b[j]);
                i += 1;
                j += 1;
                if same { Op::Match } else { Op::Mismatch }
            }
            Op::Del => {
                i += 1;
                o
            }
            Op::Ins => {
                j += 1;
                o
            }
        };
        match ops.last_mut() {
            Some((p, c)) if *p == o => *c += 1,
            _ => {
                if matches!(o, Op::Ins | Op::Del) {
                    score -= sc.open;
                }
                ops.push((o, 1));
            }
        }
        if matches!(o, Op::Ins | Op::Del) {
            score -= sc.extend;
        }
    }
    Alignment { score, a: found.a, b: found.b, ops, ties: found.ties }
}

// Traceback byte: bits 0-1 the state before an H cell (0 H, 1 E, 2 F, 3 the start); bit 2 E opened (else extended),
// bit 3 E opened from F (else from H); bit 4 F opened, bit 5 F opened from E. E consumes b (Ins), F consumes a (Del).
const START: u8 = 3;

/// A filled matrix in diagonals [lo, hi]: the best end (score, i, j, state), the other end cells reaching that score
/// (not for local), and, when kept, the traceback bytes row by row.
struct Filled {
    lo: isize,
    hi: isize,
    score: i32,
    end: (i32, usize, usize, u8),
    ties: Vec<(usize, usize)>,
    trace: Vec<u8>,
    off: Vec<usize>,
}

/// The columns of row i inside diagonals [lo, hi] of a matrix with m + 1 columns, or None.
fn row(i: usize, lo: isize, hi: isize, m: usize) -> Option<(usize, usize)> {
    let (l, h) = ((i as isize + lo).max(0), (i as isize + hi).min(m as isize));
    (h >= l).then_some((l as usize, h as usize))
}

/// `report`: tell the host how far the fill has got, every 2^20 cells (only the pass over the whole input reports).
fn fill(a: &[u8], b: &[u8], mode: Mode, sc: &Scoring, (lo, hi): (isize, isize), keep: bool, report: bool) -> Result<Filled, String> {
    let (n, m) = (a.len(), b.len());
    let (open, ext) = (sc.open + sc.extend, sc.extend);
    let mut off = Vec::new();
    let mut trace = Vec::new();
    if keep {
        let mut total = 0usize;
        for i in 0..=n {
            off.push(total);
            if let Some((l, h)) = row(i, lo, hi, m) {
                total += h - l + 1;
            }
        }
        if total > TRACE_CELLS {
            return Err("too long to align: the traceback would be over 100 million cells".into());
        }
        trace = vec![0u8; total];
    }
    // Letters as codes once (DNA masks, protein matrix indexes), so a cell's score is a lookup.
    let (ac, bc): (Vec<u8>, Vec<u8>) = if sc.protein {
        (a.iter().map(|&x| index(x) as u8).collect(), b.iter().map(|&x| index(x) as u8).collect())
    } else {
        (a.iter().map(|&x| mask(x)).collect(), b.iter().map(|&x| mask(x)).collect())
    };
    let mx = matrix();
    let pair = |x: u8, y: u8| -> i32 {
        if sc.protein {
            mx[x as usize][y as usize] as i32
        } else if x != 0 && y != 0 && (x & !y == 0 || y & !x == 0) {
            2
        } else {
            -3
        }
    };
    // Score rows indexed by column; cells outside a row's range hold NEG (reset before a buffer is reused).
    let (mut hp, mut ep, mut fp) = (vec![NEG; m + 1], vec![NEG; m + 1], vec![NEG; m + 1]);
    let (mut hc, mut ec, mut fc) = (vec![NEG; m + 1], vec![NEG; m + 1], vec![NEG; m + 1]);
    // The range each buffer last held, to reset it.
    let (mut held_p, mut held_c): (Option<(usize, usize)>, Option<(usize, usize)>) = (None, None);
    let mut ends: Vec<(i32, usize, usize, u8)> = Vec::new();
    // A local alignment's best starts as the empty one (score 0) at the first cell.
    let mut local = (0, 0usize, 0usize, 0u8);
    // Cells filled so far, of about `total`, and when progress was last told.
    let total = (n + 1).saturating_mul((hi - lo + 1).max(0) as usize).min((n + 1).saturating_mul(m + 1));
    let (mut cells, mut told) = (0usize, 0usize);
    for i in 0..=n {
        if report && cells - told >= 1 << 20 {
            crate::progress(cells, total);
            told = cells;
        }
        if let Some((l, h)) = held_c.take() {
            hc[l..=h].fill(NEG);
            ec[l..=h].fill(NEG);
            fc[l..=h].fill(NEG);
        }
        let Some((jl, jh)) = row(i, lo, hi, m) else {
            // An empty row: the next one reads nothing above it.
            if let Some((l, h)) = held_p.take() {
                hp[l..=h].fill(NEG);
                ep[l..=h].fill(NEG);
                fp[l..=h].fill(NEG);
            }
            continue;
        };
        cells += jh - jl + 1;
        let (global, overlap, is_local) = (mode == Mode::Global, mode == Mode::Overlap, mode == Mode::Local);
        let row_off = if keep { off[i] } else { 0 };
        // The cell left of the next one: (H, E, F); none at the band's left edge.
        let (mut lh, mut le, mut lf) = (NEG, NEG, NEG);
        let mut j = jl;
        if i == 0 || j == 0 {
            // The first row, or the first column of this row: starts and end gaps (see the mode rules).
            for jj in jl..=if i == 0 { jh } else { jl } {
                let mut t = 0u8;
                let h = if (i == 0 && jj == 0) || (jj == 0 && !global) || (i == 0 && overlap) {
                    t = START;
                    0
                } else {
                    NEG
                };
                // E: from the left (row 0 only: column 0 has no left).
                let e = if jj > jl {
                    let ex = le - ext;
                    let (op, from) = if lf > lh { (lf, 1u8) } else { (lh, 0u8) };
                    if op - open >= ex {
                        t |= 4 | (from << 3);
                        op - open
                    } else {
                        ex
                    }
                } else {
                    NEG
                };
                // F: from above, in column 0 of a global alignment only (before b starts the others pay nothing).
                let f = if i > 0 && global {
                    let (uh, ue, uf) = (hp[0], ep[0], fp[0]);
                    let ex = uf - ext;
                    let (op, from) = if ue > uh { (ue, 1u8) } else { (uh, 0u8) };
                    if op - open >= ex {
                        t |= 16 | (from << 5);
                        op - open
                    } else {
                        ex
                    }
                } else {
                    NEG
                };
                (hc[jj], ec[jj], fc[jj]) = (h, e, f);
                (lh, le, lf) = (h, e, f);
                if keep {
                    trace[row_off + jj - jl] = t;
                }
            }
            j = if i == 0 { jh + 1 } else { jl + 1 };
        }
        // Inner cells: scores below NEG / 2 are unreachable, so no cell checks for NEG before adding.
        let ai = if i > 0 { ac[i - 1] } else { 0 };
        while j <= jh {
            let (ph, pe, pf) = (hp[j - 1], ep[j - 1], fp[j - 1]);
            let (mut best, mut from) = (ph, 0u8);
            if pe > best {
                (best, from) = (pe, 1);
            }
            if pf > best {
                (best, from) = (pf, 2);
            }
            if is_local && best < 0 {
                (best, from) = (0, START);
            }
            let h = best + pair(ai, bc[j - 1]);
            let mut t = from;
            let ex = le - ext;
            let (op, ef) = if lf > lh { (lf, 1u8) } else { (lh, 0u8) };
            let e = if op - open >= ex {
                t |= 4 | (ef << 3);
                op - open
            } else {
                ex
            };
            let (uh, ue, uf) = (hp[j], ep[j], fp[j]);
            let fx = uf - ext;
            let (op, ff) = if ue > uh { (ue, 1u8) } else { (uh, 0u8) };
            let f = if op - open >= fx {
                t |= 16 | (ff << 5);
                op - open
            } else {
                fx
            };
            (hc[j], ec[j], fc[j]) = (h, e, f);
            (lh, le, lf) = (h, e, f);
            if keep {
                trace[row_off + j - jl] = t;
            }
            if is_local && h > local.0 {
                local = (h, i, j, 0);
            }
            j += 1;
        }
        // Ends: glocal at b's end; overlap also anywhere on a's last row; global at both ends.
        let mut end_at = |j: usize| {
            let (h, e, f) = (hc[j], ec[j], fc[j]);
            let (mut v, mut s) = (h, 0u8);
            if e > v {
                (v, s) = (e, 1);
            }
            if f > v {
                (v, s) = (f, 2);
            }
            if v > NEG / 2 {
                ends.push((v, i, j, s));
            }
        };
        match mode {
            Mode::Overlap if i == n => (jl..=jh).for_each(&mut end_at),
            Mode::Glocal | Mode::Overlap if jh == m => end_at(m),
            Mode::Global if i == n && jh == m => end_at(m),
            _ => {}
        }
        std::mem::swap(&mut hp, &mut hc);
        std::mem::swap(&mut ep, &mut ec);
        std::mem::swap(&mut fp, &mut fc);
        (held_c, held_p) = (held_p, Some((jl, jh)));
    }
    let end = if mode == Mode::Local {
        local
    } else {
        // The first end reaching the best score, in row order; the others reaching it are ties.
        let best = ends.iter().map(|e| e.0).max().unwrap_or(NEG);
        *ends.iter().find(|e| e.0 == best).unwrap_or(&(NEG, 0, 0, 0))
    };
    if end.0 <= NEG / 2 {
        return Err("no alignment fits the band".into());
    }
    let ties = ends.iter().filter(|e| e.0 == end.0 && (e.1, e.2) != (end.1, end.2)).map(|e| (e.1, e.2)).collect();
    // A reporting pass ends at 1, however few cells it had, so the host can tell one pass from the next.
    if report {
        crate::progress(total, total);
    }
    Ok(Filled { lo, hi, score: end.0, end, ties, trace, off })
}

fn traceback(a: &[u8], b: &[u8], sc: &Scoring, f: &Filled) -> Alignment {
    let m = b.len();
    let (score, mut i, mut j, mut state) = f.end;
    let (ai, bj) = (i, j);
    let mut ops: Vec<(Op, u32)> = Vec::new();
    let mut push = |o: Op| match ops.last_mut() {
        Some((p, c)) if *p == o => *c += 1,
        _ => ops.push((o, 1)),
    };
    loop {
        let (jl, _) = row(i, f.lo, f.hi, m).expect("the traceback stays in the band");
        let t = f.trace[f.off[i] + j - jl];
        match state {
            0 => {
                // An H cell on the first row or column is a start, not a column of the alignment; a local alignment's
                // first column is a cell whose predecessor is the start.
                if i == 0 || j == 0 {
                    break;
                }
                let from = t & 3;
                push(if sc.same(a[i - 1], b[j - 1]) { Op::Match } else { Op::Mismatch });
                i -= 1;
                j -= 1;
                if from == START {
                    break;
                }
                state = from;
            }
            1 => {
                push(Op::Ins);
                j -= 1;
                if t & 4 != 0 {
                    state = if t & 8 != 0 { 2 } else { 0 };
                }
            }
            _ => {
                push(Op::Del);
                i -= 1;
                if t & 16 != 0 {
                    state = if t & 32 != 0 { 1 } else { 0 };
                }
            }
        }
    }
    ops.reverse();
    Alignment { score, a: (i, ai), b: (j, bj), ops, ties: f.ties.clone() }
}

/// The range of diagonals (j − i) shared k-mers of a and b lie on, 12 bases for DNA (A, C, G and T only), 5 residues
/// for protein: every diagonal with at least a tenth of the best one's k-mers (and at least 2, when the best has more),
/// so an indel between two well-matched stretches lies inside. None when they share none.
pub fn diagonals(a: &[u8], b: &[u8], protein: bool) -> Option<(isize, isize)> {
    let (k, bits) = if protein { (5usize, 5u32) } else { (12, 2) };
    let code = |x: u8| -> Option<u64> {
        let x = x.to_ascii_uppercase();
        if protein { (x.is_ascii_uppercase()).then(|| (x - b'A') as u64) } else { b"ACGT".iter().position(|&c| c == x).map(|p| p as u64) }
    };
    let mask = (1u64 << (bits * k as u32)) - 1;
    let kmers = |s: &[u8], mut f: Box<dyn FnMut(u64, usize) + '_>| {
        let (mut v, mut run) = (0u64, 0usize);
        for (p, &x) in s.iter().enumerate() {
            match code(x) {
                Some(c) => {
                    v = ((v << bits) | c) & mask;
                    run += 1;
                    if run >= k {
                        f(v, p + 1 - k);
                    }
                }
                None => run = 0,
            }
        }
    };
    let mut at: std::collections::HashMap<u64, Vec<usize>> = std::collections::HashMap::new();
    kmers(b, Box::new(|v, p| at.entry(v).or_default().push(p)));
    // A k-mer repeated many times in b says little about where b lies: it is left out.
    at.retain(|_, v| v.len() <= 8);
    let mut votes: std::collections::HashMap<isize, u32> = std::collections::HashMap::new();
    kmers(a, Box::new(|v, i| {
        if let Some(js) = at.get(&v) {
            for &j in js {
                *votes.entry(j as isize - i as isize).or_default() += 1;
            }
        }
    }));
    let best = *votes.values().max()?;
    let floor = if best >= 2 { (best / 10).max(2) } else { 1 };
    let kept = votes.into_iter().filter(|&(_, c)| c >= floor).map(|(d, _)| d);
    kept.fold(None, |r: Option<(isize, isize)>, d| Some(r.map_or((d, d), |(lo, hi)| (lo.min(d), hi.max(d)))))
}

/// Mott's trimming (as Phred's): with p the error probability of each base's Phred value, the stretch where the sum
/// of (limit − p) is largest. [start, end), empty when no base is better than the limit.
pub fn trim(quality: &[u8], limit: f64) -> (usize, usize) {
    let (mut best, mut range) = (0.0f64, (0, 0));
    let (mut sum, mut start) = (0.0f64, 0usize);
    for (i, &q) in quality.iter().enumerate() {
        sum += limit - 10f64.powf(-(q as f64) / 10.0);
        if sum <= 0.0 {
            sum = 0.0;
            start = i + 1;
        } else if sum > best {
            best = sum;
            range = (start, i + 1);
        }
    }
    range
}

#[cfg(test)]
mod tests {
    use super::*;
    use Op::*;

    fn run(a: &str, b: &str, mode: Mode) -> Alignment {
        align(a.as_bytes(), b.as_bytes(), mode, &Scoring::DNA, None).unwrap()
    }

    #[test]
    fn blosum62_is_ncbi_s_table() {
        // Symmetric, and known entries (diagonal and pairs) of NCBI's file.
        for &x in LETTERS {
            for &y in LETTERS {
                assert_eq!(blosum62(x, y), blosum62(y, x), "{} {}", x as char, y as char);
            }
        }
        let diag: Vec<i32> = b"ARNDCQEGHILKMFPSTWYV".iter().map(|&x| blosum62(x, x)).collect();
        assert_eq!(diag, [4, 5, 6, 6, 9, 5, 5, 6, 8, 4, 4, 5, 5, 6, 7, 4, 5, 11, 7, 4]);
        assert_eq!([blosum62(b'A', b'R'), blosum62(b'I', b'V'), blosum62(b'W', b'F'), blosum62(b'E', b'Z'), blosum62(b'*', b'*'), blosum62(b'J', b'A')], [-1, 3, 1, 4, 1, 0]);
    }

    #[test]
    fn global_alignments_by_hand() {
        // Identical: 8 matches, 16.
        assert_eq!(run("ACGTACGT", "ACGTACGT", Mode::Global), Alignment { score: 16, a: (0, 8), b: (0, 8), ops: vec![(Match, 8)], ties: vec![] });
        // One substitution: 7 × 2 − 3 = 11.
        assert_eq!(run("ACGTACGT", "ACGAACGT", Mode::Global).ops, vec![(Match, 3), (Mismatch, 1), (Match, 4)]);
        assert_eq!(run("ACGTACGT", "ACGAACGT", Mode::Global).score, 11);
        // One base deleted from b (the C at 4 of a; no other placement scores as well): 7 × 2 − (5 + 2) = 7.
        let d = run("ACGTCACG", "ACGTACG", Mode::Global);
        assert_eq!((d.score, d.ops.clone()), (7, vec![(Match, 4), (Del, 1), (Match, 3)]));
        // One base inserted in b, the mirror.
        assert_eq!(run("ACGTACG", "ACGTCACG", Mode::Global).ops, vec![(Match, 4), (Ins, 1), (Match, 3)]);
        // Affine: one gap of three (5 + 6 = 11) rather than three of one (21): 9 × 2 − 11 = 7.
        let g = run("ACGTAAACGTTT", "ACGTCGTTT", Mode::Global);
        assert_eq!((g.score, g.ops.iter().filter(|(o, _)| *o == Del).count()), (7, 1));
    }

    #[test]
    fn glocal_and_local_by_hand() {
        // All of b inside a, a's ends free: 7 matches, 14, at 4-11 of a.
        assert_eq!(run("TTTTACGTACGTTTTT", "ACGTACG", Mode::Glocal), Alignment { score: 14, a: (4, 11), b: (0, 7), ops: vec![(Match, 7)], ties: vec![] });
        // Local: the best common stretch, ACGT (8), at 3-7 of a and 2-6 of b.
        assert_eq!(run("GGGACGTGGG", "TTACGTTT", Mode::Local), Alignment { score: 8, a: (3, 7), b: (2, 6), ops: vec![(Match, 4)], ties: vec![] });
        // IUPAC: N in b scores as a match against any base of a (16), as a mismatch column.
        let x = run("ACGTACGT", "ACNTACGT", Mode::Global);
        assert_eq!((x.score, x.ops), (16, vec![(Match, 2), (Mismatch, 1), (Match, 5)]));
    }

    #[test]
    fn a_band_given_equals_full_and_too_long_is_refused() {
        let mut r = 7u32;
        let mut bases = |n: usize| -> Vec<u8> { (0..n).map(|_| { r ^= r << 13; r ^= r >> 17; r ^= r << 5; b"ACGT"[(r >> 30) as usize] }).collect() };
        let a = bases(3000);
        let mut b = a[1000..1800].to_vec();
        b[100] = if b[100] == b'A' { b'C' } else { b'A' };
        b.remove(400);
        let full = align(&a, &b, Mode::Glocal, &Scoring::DNA, None).unwrap();
        let banded = align(&a, &b, Mode::Glocal, &Scoring::DNA, Some((-1000 - 150, -1000 + 150))).unwrap();
        assert_eq!((full.score, full.a, full.b), (banded.score, banded.a, banded.b));
        assert_eq!((full.a, full.score), ((1000, 1800), 798 * 2 - 3 - 7));
        // A long input; and the diagonals its shared 12-mers lie on.
        let long = bases(30_000);
        let mut read = long[20_000..21_000].to_vec();
        read[500] = if read[500] == b'G' { b'T' } else { b'G' };
        let x = align(&long, &read, Mode::Glocal, &Scoring::DNA, None).unwrap();
        assert_eq!((x.a, x.ops.clone()), ((20_000, 21_000), vec![(Match, 500), (Mismatch, 1), (Match, 499)]));
        assert_eq!(diagonals(&long, &read, false), Some((-20_000, -20_000)));
        // Past MAX_WORK: refused before any cell is filled.
        assert_eq!(align(&vec![b'A'; 2_000_000], &vec![b'C'; 600], Mode::Global, &Scoring::DNA, None), Err("too long to align: over 1 billion cells".into()));
    }

    #[test]
    fn iupac_scores_the_same_in_both_orders_and_gaps_go_3_prime() {
        let (x, y) = (run("ACGAACGT", "ACGRACGT", Mode::Global), run("ACGRACGT", "ACGAACGT", Mode::Global));
        // A against R scores as a match (2) but is not identical: a mismatch column, identity 7/8.
        assert_eq!((x.score, y.score, x.ops.clone(), x.identity()), (16, 16, vec![(Match, 3), (Mismatch, 1), (Match, 4)], 0.875));
        // AAAA against AAA: the deletion is the last A, whichever A the traceback took.
        let d = run("CAAAAG", "CAAAG", Mode::Global);
        assert_eq!((d.score, d.ops.clone()), (10 - 7, vec![(Match, 4), (Del, 1), (Match, 1)]));
        let i = run("CAAAG", "CAAAAG", Mode::Global);
        assert_eq!(i.ops, vec![(Match, 4), (Ins, 1), (Match, 1)]);
        // A repeated unit: CAG CAG against one CAG deleted, placed after the last whole repeat it can reach.
        let r = run("TTCAGCAGTT", "TTCAGTT", Mode::Global);
        assert_eq!(r.ops, vec![(Match, 5), (Del, 3), (Match, 2)]);
    }

    #[test]
    fn a_long_deletion_lies_inside_the_band() {
        let mut r = 9u32;
        let mut bases = |n: usize| -> Vec<u8> { (0..n).map(|_| { r ^= r << 13; r ^= r >> 17; r ^= r << 5; b"ACGT"[(r >> 30) as usize] }).collect() };
        // 20 kb reference; a 4 kb read with 600 bases deleted in its middle: both stretches' diagonals are in the band.
        let a = bases(20_000);
        let read = [&a[10_000..12_000], &a[12_600..14_600]].concat();
        let x = align(&a, &read, Mode::Glocal, &Scoring::DNA, None).unwrap();
        // The deletion is placed 3′-most: here one base later, as reference 12,000 and 12,600 are the same base.
        assert_eq!(a[12_000], a[12_600]);
        assert_eq!((x.a, x.score, x.ops.clone()), ((10_000, 14_600), 8000 - 5 - 1200, vec![(Match, 2001), (Del, 600), (Match, 1999)]));
    }

    #[test]
    fn past_the_traceback_limit_the_exact_best_is_found_and_ties_are_reported() {
        let mut r = 31u32;
        let mut bases = |n: usize| -> Vec<u8> { (0..n).map(|_| { r ^= r << 13; r ^= r >> 17; r ^= r << 5; b"ACGT"[(r >> 30) as usize] }).collect() };
        // 200 kb reference (over TRACE_CELLS with a 1 kb read). The read's true place, at 50,000, has every twelfth
        // base changed, so it shares no 12-mer with the read; a decoy at 150,000 holds the read's first 400 bases
        // exactly (many shared 12-mers) and nothing after.
        let mut a = bases(200_000);
        let read = bases(1000);
        let mut copy = read.clone();
        for k in (0..1000).step_by(12) {
            copy[k] = if copy[k] == b'A' { b'C' } else { b'A' };
        }
        a[50_000..51_000].copy_from_slice(&copy);
        a[150_000..150_400].copy_from_slice(&read[..400]);
        assert!((a.len() + 1) * (read.len() + 1) > TRACE_CELLS);
        let x = align(&a, &read, Mode::Glocal, &Scoring::DNA, None).unwrap();
        let changed = 1000usize.div_ceil(12) as i32;
        assert_eq!((x.a, x.score, x.ties.len()), ((50_000, 51_000), (1000 - changed) * 2 - 3 * changed, 0));
        // The same read placed twice exactly: the second end is a tie.
        let mut b2 = a.clone();
        b2[150_000..151_000].copy_from_slice(&copy);
        let y = align(&b2, &read, Mode::Glocal, &Scoring::DNA, None).unwrap();
        assert_eq!((y.a, y.ties), ((50_000, 51_000), vec![(151_000, 1000)]));
        // A 12 kb consensus against itself (144 million cells): traced in a band inside the window, exactly.
        let c = bases(12_000);
        let w = align(&c, &c, Mode::Glocal, &Scoring::DNA, None).unwrap();
        assert_eq!((w.a, w.score, w.ops), ((0, 12_000), 24_000, vec![(Match, 12_000)]));
        // Tandem copies in a short reference: every copy's end is a tie.
        let z = align(b"AACGAACGAACGAACGAACGAACG", b"AACGAACG", Mode::Glocal, &Scoring::DNA, None).unwrap();
        assert_eq!(z.ties.len(), 4);
    }

    #[test]
    fn overlap_leaves_b_past_a_s_ends_unaligned() {
        // b runs 3 bases before a's start and 2 past its end: those are not aligned and cost nothing.
        let x = run("ACGTACGTAC", "TTTACGTACGTACGG", Mode::Overlap);
        assert_eq!((x.score, x.a, x.b, x.ops.clone()), (20, (0, 10), (3, 13), vec![(Match, 10)]));
        // Inside a, overlap is glocal.
        assert_eq!(run("TTTTACGTACGTTTTT", "ACGTACG", Mode::Overlap), run("TTTTACGTACGTTTTT", "ACGTACG", Mode::Glocal));
    }

    #[test]
    fn local_with_nothing_in_common_is_empty() {
        assert_eq!(run("A", "C", Mode::Local), Alignment { score: 0, a: (0, 0), b: (0, 0), ops: vec![], ties: vec![] });
        assert_eq!(run("", "ACGT", Mode::Local).score, 0);
    }

    #[test]
    fn protein_alignment_uses_blosum62() {
        let p = align(b"MKVLAAGG", b"MKILAAGG", Mode::Global, &Scoring::PROTEIN, None).unwrap();
        // V/I is 3: 5 + 5 + 3 + 4 + 4 + 4 + 6 + 6.
        assert_eq!((p.score, p.ops.clone()), (37, vec![(Match, 2), (Mismatch, 1), (Match, 5)]));
    }

    #[test]
    fn mott_trimming_by_hand() {
        // Q10 is p 0.1 (0.05 − 0.1 < 0), Q30 is p 0.001 (+0.049): the Q30 stretch is kept.
        assert_eq!(trim(&[10, 10, 30, 30, 30, 30, 10, 10], 0.05), (2, 6));
        // A single Q10 between good stretches is kept when the stretch after it outweighs it.
        assert_eq!(trim(&[30, 30, 30, 10, 30, 30, 30, 10, 10], 0.05), (0, 7));
        assert_eq!(trim(&[5, 5, 5], 0.05), (0, 0));
    }
}
