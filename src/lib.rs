//! OpenBench's sequence engine, compiled to WebAssembly and used by the browser's workers and the server.
//!
//! Calling convention (no unsafe code on this side): the host asks `ob_input(slot, len)` for a buffer and writes its
//! bytes there, then calls an operation. An operation returns a count, or a negative error code with a message in
//! the byte output (`-1` invalid input; its length from `ob_out_bytes_len`). Results are read from `ob_out_i32_ptr` (i32 values) or
//! `ob_out_bytes_ptr` (UTF-8 bytes) until the next call. The host re-reads memory views after every call, since any
//! call may grow memory.
//!
//! No unsafe code: the only `unsafe` is the `#[unsafe(no_mangle)]` export attribute (Rust 2024), allowed on each
//! export below and nowhere else.
#![deny(unsafe_code)]

pub mod abif;
pub mod align;
pub mod alphabet;
mod json;
pub mod search;
pub mod sites;
pub mod snapgene;
pub mod translate;
pub mod verify;

use std::cell::RefCell;

/// Largest input accepted in one slot: 10 Mb sequences (decision 7) and trace files fit well within it.
pub const MAX_INPUT: usize = 64 << 20;
const SLOTS: usize = 3;

thread_local! {
    static INPUT: RefCell<[Vec<u8>; SLOTS]> = RefCell::new(Default::default());
    static OUT_I32: RefCell<Vec<i32>> = const { RefCell::new(Vec::new()) };
    static OUT_BYTES: RefCell<Vec<u8>> = const { RefCell::new(Vec::new()) };
}

fn input(slot: usize) -> Vec<u8> {
    INPUT.with(|i| std::mem::take(&mut i.borrow_mut()[slot]))
}

fn ints(v: Vec<i32>) -> i32 {
    let n = v.len() as i32;
    OUT_I32.with(|o| *o.borrow_mut() = v);
    n
}

fn bytes(v: Vec<u8>) -> i32 {
    let n = v.len() as i32;
    OUT_BYTES.with(|o| *o.borrow_mut() = v);
    n
}

fn fail(code: i32, msg: impl Into<String>) -> i32 {
    bytes(msg.into().into_bytes());
    code
}

/// A zeroed buffer of `len` bytes for input `slot`; its address, or 0 when the slot or length is refused.
#[allow(unsafe_code)]
#[unsafe(no_mangle)]
pub extern "C" fn ob_input(slot: u32, len: u32) -> u32 {
    let (slot, len) = (slot as usize, len as usize);
    if slot >= SLOTS || len > MAX_INPUT {
        return 0;
    }
    INPUT.with(|i| {
        let mut i = i.borrow_mut();
        i[slot] = vec![0; len];
        i[slot].as_mut_ptr() as usize as u32
    })
}

#[allow(unsafe_code)]
#[unsafe(no_mangle)]
pub extern "C" fn ob_out_i32_ptr() -> u32 {
    OUT_I32.with(|o| o.borrow().as_ptr() as usize as u32)
}

#[allow(unsafe_code)]
#[unsafe(no_mangle)]
pub extern "C" fn ob_out_bytes_ptr() -> u32 {
    OUT_BYTES.with(|o| o.borrow().as_ptr() as usize as u32)
}

/// Length of the byte output: the result of the last byte operation, or the message of the last error.
#[allow(unsafe_code)]
#[unsafe(no_mangle)]
pub extern "C" fn ob_out_bytes_len() -> u32 {
    OUT_BYTES.with(|o| o.borrow().len() as u32)
}

/// Translate slot 0 with NCBI code `table`: the protein, as bytes.
#[allow(unsafe_code)]
#[unsafe(no_mangle)]
pub extern "C" fn ob_translate(table: u32) -> i32 {
    let seq = input(0);
    match translate::table(table) {
        Some(code) => bytes(translate::translate(code, &seq)),
        None => fail(-1, format!("unknown genetic code {table}")),
    }
}

/// ORFs of slot 0: (start, end, strand, amino acids) per ORF.
#[allow(unsafe_code)]
#[unsafe(no_mangle)]
pub extern "C" fn ob_orfs(circular: u32, min_aa: u32, table: u32) -> i32 {
    let seq = input(0);
    let Some(code) = translate::table(table) else { return fail(-1, format!("unknown genetic code {table}")) };
    let orfs = translate::orfs(code, &seq, circular != 0, min_aa as usize);
    ints(orfs.iter().flat_map(|o| [o.start as i32, o.end as i32, o.strand as i32, o.aa as i32]).collect())
}

/// Search slot 0 for the pattern in slot 1: the total hit count, then `bins` counts of all hits over equal parts of
/// the sequence, then (position, strand, substitutions) for the first `max_hits` hits.
#[allow(unsafe_code)]
#[unsafe(no_mangle)]
pub extern "C" fn ob_search(k: u32, protein: u32, circular: u32, max_hits: u32, bins: u32) -> i32 {
    let (seq, pat) = (input(0), input(1));
    let mut counts = vec![0u32; (bins as usize).min(1 << 16)];
    match search::search_binned(&seq, &pat, k as usize, protein != 0, circular != 0, max_hits as usize, &mut counts) {
        Ok((total, h)) => {
            let mut out = vec![total.min(i32::MAX as usize) as i32];
            out.extend(counts.iter().map(|&c| c.min(i32::MAX as u32) as i32));
            out.extend(h.iter().flat_map(|h| [h.pos as i32, h.strand as i32, h.substitutions as i32]));
            ints(out)
        }
        Err(e) => fail(-1, format!("{e:?}")),
    }
}

/// Cut positions in slot 0 for the enzymes in slot 1 (one per line: site, fcut, rcut, sensitivity letters d/c/g,
/// tab-separated): (enzyme line index, position, blocked mask) per cut.
#[allow(unsafe_code)]
#[unsafe(no_mangle)]
pub extern "C" fn ob_cuts(circular: u32) -> i32 {
    let (seq, table) = (input(0), input(1));
    let mut enzymes = Vec::new();
    for (i, line) in String::from_utf8_lossy(&table).lines().enumerate() {
        let f: Vec<&str> = line.split('\t').collect();
        let parsed = match f.as_slice() {
            [site, fcut, rcut, meth] => fcut.parse().ok().zip(rcut.parse().ok()).map(|(fcut, rcut)| sites::Enzyme {
                site: site.as_bytes().to_vec(),
                fcut,
                rcut,
                sensitive: meth.bytes().fold(0, |m, c| m | match c { b'd' => sites::DAM, b'c' => sites::DCM, b'g' => sites::CPG, _ => 0 }),
            }),
            _ => None,
        };
        match parsed {
            Some(e) if !e.site.is_empty() && e.site.len() <= 64 && e.site.iter().all(|&b| alphabet::mask(b) != 0) => enzymes.push(e),
            _ => return fail(-1, format!("enzyme line {} is not site, fcut, rcut, sensitivity", i + 1)),
        }
    }
    let cuts = sites::cuts(&seq, circular != 0, &enzymes);
    ints(cuts.iter().flat_map(|c| [c.enzyme as i32, c.pos as i32, c.blocked as i32]).collect())
}

/// Parse a SnapGene `.dna` file in slot 0 (sequences up to `max_len`): JSON bytes.
#[allow(unsafe_code)]
#[unsafe(no_mangle)]
pub extern "C" fn ob_parse_snapgene(max_len: u32) -> i32 {
    match snapgene::parse(&input(0), max_len as usize) {
        Ok(d) => bytes(snapgene::to_json(&d).into_bytes()),
        Err(e) => fail(-1, e),
    }
}

/// Parse an ABIF `.ab1` trace in slot 0: JSON bytes.
#[allow(unsafe_code)]
#[unsafe(no_mangle)]
pub extern "C" fn ob_parse_abif() -> i32 {
    match abif::parse(&input(0)) {
        Ok(t) => bytes(abif::to_json(&t).into_bytes()),
        Err(e) => fail(-1, e),
    }
}

/// Align slot 1 (b) to slot 0 (a): mode 0 global, 1 local, 2 glocal; protein scoring when `protein`. Ints: score,
/// a start, a end, b start, b end, then (operation, length) per run (0 match, 1 mismatch, 2 insertion, 3 deletion).
#[allow(unsafe_code)]
#[unsafe(no_mangle)]
pub extern "C" fn ob_align(mode: u32, protein: u32) -> i32 {
    let (a, b) = (input(0), input(1));
    let mode = match mode { 0 => align::Mode::Global, 1 => align::Mode::Local, 2 => align::Mode::Glocal, _ => return fail(-1, "unknown alignment mode") };
    let sc = if protein != 0 { align::Scoring::PROTEIN } else { align::Scoring::DNA };
    match align::align(&a, &b, mode, &sc, None) {
        Ok(x) => {
            let mut out = vec![x.score, x.a.0 as i32, x.a.1 as i32, x.b.0 as i32, x.b.1 as i32];
            out.extend(x.ops.iter().flat_map(|&(o, n)| [o as i32, n as i32]));
            ints(out)
        }
        Err(e) => fail(-1, e),
    }
}

/// Mott trimming of the Phred values in slot 0 at `limit_milli` / 1000: ints start, end.
#[allow(unsafe_code)]
#[unsafe(no_mangle)]
pub extern "C" fn ob_trim(limit_milli: u32) -> i32 {
    let (s, e) = align::trim(&input(0), limit_milli as f64 / 1000.0);
    ints(vec![s as i32, e as i32])
}

thread_local! {
    static VERIFY: RefCell<Option<verify::Job>> = const { RefCell::new(None) };
}

/// Start a verification against the reference in slot 0: region [start, end), Phred threshold, trimming limit in
/// thousandths. Reads follow one by one (`ob_verify_add`), then `ob_verify_end`.
#[allow(unsafe_code)]
#[unsafe(no_mangle)]
pub extern "C" fn ob_verify_begin(circular: u32, start: u32, end: u32, min_q: u32, limit_milli: u32) -> i32 {
    match verify::Job::new(&input(0), circular != 0, (start as usize, end as usize), min_q.min(255) as u8, limit_milli as f64 / 1000.0) {
        Ok(j) => {
            VERIFY.with(|v| *v.borrow_mut() = Some(j));
            0
        }
        Err(e) => fail(-1, e),
    }
}

/// Add a read: bases in slot 0, Phred values in slot 1 (one per base) when `has_quality`, its name in slot 2.
#[allow(unsafe_code)]
#[unsafe(no_mangle)]
pub extern "C" fn ob_verify_add(has_quality: u32) -> i32 {
    let (bases, quality, name) = (input(0), input(1), input(2));
    if has_quality != 0 && quality.len() != bases.len() {
        return fail(-1, "a read needs one quality value per base");
    }
    VERIFY.with(|v| match v.borrow_mut().as_mut() {
        Some(j) => {
            j.add(&String::from_utf8_lossy(&name), &bases, (has_quality != 0).then_some(&quality[..]));
            0
        }
        None => fail(-1, "no verification was started"),
    })
}

/// The verification's report as JSON bytes; the job ends.
#[allow(unsafe_code)]
#[unsafe(no_mangle)]
pub extern "C" fn ob_verify_end() -> i32 {
    match VERIFY.with(|v| v.borrow_mut().take()) {
        Some(j) => bytes(verify::to_json(&j, &j.report()).into_bytes()),
        None => fail(-1, "no verification was started"),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn put(slot: u32, data: &[u8]) {
        assert_ne!(ob_input(slot, data.len() as u32), 0);
        INPUT.with(|i| i.borrow_mut()[slot as usize].copy_from_slice(data));
    }

    #[test]
    fn calls_through_the_abi() {
        put(0, b"ATGAAATAG");
        assert_eq!(ob_translate(1), 3);
        assert_eq!(OUT_BYTES.with(|o| o.borrow().clone()), b"MK*");
        put(0, b"AAGGATCCAA");
        put(1, b"GGATCC");
        // Total 1, two bins of five bases (the hit at 2 in the first), then the hit.
        assert_eq!(ob_search(0, 0, 0, 10, 2), 6);
        assert_eq!(OUT_I32.with(|o| o.borrow().clone()), vec![1, 1, 0, 2, 1, 0]);
        put(0, b"AAGAATTCAA");
        put(1, b"GAATTC\t1\t5\t\nGATC\t0\t4\td");
        assert_eq!(ob_cuts(0), 3);
        assert_eq!(OUT_I32.with(|o| o.borrow().clone()), vec![0, 3, 0]);
    }

    #[test]
    fn aligns_trims_and_verifies_through_the_abi() {
        put(0, b"ACGTCACG");
        put(1, b"ACGTACG");
        assert_eq!(ob_align(0, 0), 11);
        assert_eq!(OUT_I32.with(|o| o.borrow().clone()), vec![7, 0, 8, 0, 7, 0, 4, 3, 1, 0, 3]);
        put(0, &[10, 10, 30, 30, 30, 30, 10, 10]);
        assert_eq!(ob_trim(50), 2);
        assert_eq!(OUT_I32.with(|o| o.borrow().clone()), vec![2, 6]);
        let reference: Vec<u8> = (0..120).map(|i| b"ACGTTGCAAC"[(i * 7 + i / 3) % 10]).collect();
        put(0, &reference);
        assert_eq!(ob_verify_begin(0, 10, 110, 20, 50), 0);
        put(0, &reference[5..115]);
        put(1, &[40u8; 110]);
        put(2, b"clone1_F");
        assert_eq!(ob_verify_add(1), 0);
        assert!(ob_verify_end() > 0);
        let json = OUT_BYTES.with(|o| String::from_utf8(o.borrow().clone()).unwrap());
        assert!(json.starts_with(r#"{"verdict":"pass","quality":true,"covered":100,"single":100,"differences":[]"#), "{json}");
        assert_eq!(ob_verify_end(), -1);
    }

    #[test]
    fn errors_carry_a_message() {
        put(0, b"ACGT");
        assert_eq!(ob_translate(99), -1);
        assert_eq!(ob_out_bytes_len(), 23);
        assert_eq!(OUT_BYTES.with(|o| String::from_utf8(o.borrow().clone()).unwrap()), "unknown genetic code 99");
        put(0, b"ACGT");
        put(1, b"GAATTC\tx\t5\t");
        assert_eq!(ob_cuts(0), -1);
        assert_eq!(ob_input(0, (MAX_INPUT + 1) as u32), 0);
        assert_eq!(ob_input(SLOTS as u32, 1), 0);
    }
}
