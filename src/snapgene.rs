//! SnapGene `.dna` files: a sequence of packets (type byte, big-endian u32 length, data). Read here: the cookie
//! (0x09, "SnapGene", sequence type 1 = DNA), the sequence (0x00: a flags byte, bit 0 = circular, then ASCII bases)
//! and the features (0x0A: XML). Other packets are skipped by length.

use crate::json;

pub const MAX_FEATURES: usize = 10_000;

#[derive(Debug, PartialEq)]
pub struct Feature {
    pub name: String,
    pub kind: String,
    /// 0-based start; `end` exclusive; `end < start` wraps the origin.
    pub start: usize,
    pub end: usize,
    pub direction: i8,
    pub color: Option<String>,
}

#[derive(Debug, PartialEq)]
pub struct Dna {
    pub seq: String,
    pub circular: bool,
    pub features: Vec<Feature>,
    /// Features left out because none of their segments lies within the sequence; the caller says so.
    pub skipped: usize,
}

pub fn parse(bytes: &[u8], max_len: usize) -> Result<Dna, String> {
    let mut at = 0;
    let mut cookie = false;
    let mut dna: Option<(String, bool)> = None;
    let mut xml: Option<&[u8]> = None;
    while at < bytes.len() {
        if bytes.len() - at < 5 {
            return Err("the file ends inside a packet header".into());
        }
        let kind = bytes[at];
        let len = u32::from_be_bytes(bytes[at + 1..at + 5].try_into().unwrap()) as usize;
        let body = bytes.get(at + 5..at + 5 + len).ok_or("a packet is longer than the file")?;
        at += 5 + len;
        match kind {
            0x09 => {
                if body.len() < 10 || &body[..8] != b"SnapGene" {
                    return Err("not a SnapGene file".into());
                }
                if u16::from_be_bytes([body[8], body[9]]) != 1 {
                    return Err("not a DNA file (SnapGene protein and RNA files are not supported)".into());
                }
                cookie = true;
            }
            _ if !cookie => return Err("not a SnapGene file".into()),
            0x00 => {
                let (&flags, seq) = body.split_first().ok_or("the sequence packet is empty")?;
                if seq.len() > max_len {
                    return Err(format!("the sequence is longer than {max_len} bases"));
                }
                if !seq.iter().all(|b| b.is_ascii_alphabetic()) {
                    return Err("the sequence holds characters that are not bases".into());
                }
                dna = Some((String::from_utf8_lossy(seq).to_ascii_uppercase(), flags & 1 == 1));
            }
            0x0A => xml = Some(body),
            _ => {}
        }
    }
    let (seq, circular) = dna.ok_or("the file has no sequence")?;
    let (features, skipped) = match xml {
        Some(x) => features(&String::from_utf8_lossy(x), seq.len(), circular)?,
        None => (vec![], 0),
    };
    Ok(Dna { seq, circular, features, skipped })
}

// Attribute `name` of the tag text starting after `<Tag`.
fn attr(tag: &str, name: &str) -> Option<String> {
    let mut rest = tag;
    while let Some(i) = rest.find(name) {
        let before_ok = i == 0 || rest[..i].ends_with(|c: char| c.is_whitespace());
        let after = rest[i + name.len()..].trim_start();
        if before_ok && after.starts_with('=') {
            let v = after[1..].trim_start();
            let q = v.chars().next()?;
            if q == '"' || q == '\'' {
                let end = v[1..].find(q)?;
                return Some(unescape(&v[1..1 + end]));
            }
            return None;
        }
        rest = &rest[i + name.len()..];
    }
    None
}

fn unescape(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    let mut rest = s;
    while let Some(i) = rest.find('&') {
        out.push_str(&rest[..i]);
        rest = &rest[i..];
        // Only a short window is searched, so a run of '&' costs linear time.
        let Some(semi) = rest.bytes().take(12).position(|b| b == b';') else {
            out.push('&');
            rest = &rest[1..];
            continue;
        };
        let ent = &rest[1..semi];
        let ch = match ent {
            "amp" => Some('&'),
            "lt" => Some('<'),
            "gt" => Some('>'),
            "quot" => Some('"'),
            "apos" => Some('\''),
            _ if ent.starts_with("#x") => u32::from_str_radix(&ent[2..], 16).ok().and_then(char::from_u32),
            _ if ent.starts_with('#') => ent[1..].parse().ok().and_then(char::from_u32),
            _ => None,
        };
        match ch {
            Some(c) => {
                out.push(c);
                rest = &rest[semi + 1..];
            }
            None => {
                out.push('&');
                rest = &rest[1..];
            }
        }
    }
    out.push_str(rest);
    out
}

// 1-based inclusive "a-b" (or a single "a") to a 0-based [start, end), checked against the sequence length. On a circle
// a > b wraps the origin (end < start), and b = a - 1 is the whole circle, [0, n); on a linear sequence a > b is
// invalid.
fn range(r: &str, n: usize, circular: bool) -> Option<(usize, usize)> {
    let (a, b) = r.split_once('-').unwrap_or((r, r));
    let (a, b): (usize, usize) = (a.trim().parse().ok()?, b.trim().parse().ok()?);
    if a < 1 || b < 1 || a > n || b > n {
        return None;
    }
    match (a <= b, circular) {
        (true, _) => Some((a - 1, b)),
        (false, false) => None,
        (false, true) if b == a - 1 => Some((0, n)),
        (false, true) => Some((a - 1, b)),
    }
}

fn features(xml: &str, n: usize, circular: bool) -> Result<(Vec<Feature>, usize), String> {
    let mut out = Vec::new();
    let mut skipped = 0;
    let mut rest = xml;
    while let Some(i) = rest.find("<Feature ") {
        rest = &rest[i + 9..];
        let close = rest.find('>').ok_or("an unterminated feature tag")?;
        let tag = &rest[..close];
        let body_end = if tag.ends_with('/') { 0 } else { rest.find("</Feature>").ok_or("an unterminated feature")? };
        let body = &rest[close..body_end.max(close)];
        // Every segment must lie on the sequence, or the feature is skipped (and counted) rather than shortened.
        // Pieces are (start, end, bases covered): the count, not the endpoints, says when a join has gone round.
        let mut segments: Vec<(usize, usize, usize)> = Vec::new();
        let mut color = None;
        let mut valid = true;
        let mut b = body;
        while let Some(j) = b.find("<Segment ") {
            b = &b[j + 9..];
            let end = b.find('>').ok_or("an unterminated segment tag")?;
            let seg = &b[..end];
            match attr(seg, "range").and_then(|r| range(&r, n, circular)) {
                Some(r) => {
                    let len = if r.1 > r.0 { r.1 - r.0 } else { n - r.0 + r.1 };
                    // A segment that continues the previous one (also across the origin of a circle) extends it; a
                    // piece that has covered a whole turn is the whole circle and stays so.
                    match segments.last_mut() {
                        Some(last) if last.2 >= n => {}
                        Some(last) if last.1 == r.0 || (circular && last.1 == n && r.0 == 0) => {
                            *last = if last.2 + len >= n { (0, n, n) } else { (last.0, r.1, last.2 + len) };
                        }
                        _ => segments.push((r.0, r.1, len)),
                    }
                    if color.is_none() {
                        color = attr(seg, "color").filter(|c| c.len() <= 32);
                    }
                }
                None => valid = false,
            }
            b = &b[end..];
        }
        rest = &rest[body_end.max(close)..];
        if !valid || segments.is_empty() {
            skipped += 1;
            continue;
        }
        if out.len() + segments.len() > MAX_FEATURES {
            return Err(format!("more than {MAX_FEATURES} features"));
        }
        let direction = match attr(tag, "directionality").as_deref() {
            Some("1") => 1,
            Some("2") => -1,
            _ => 0,
        };
        let name: String = attr(tag, "name").unwrap_or_else(|| "feature".into()).chars().take(200).collect();
        let kind: String = attr(tag, "type").unwrap_or_else(|| "misc_feature".into()).chars().take(40).collect();
        // Separate segments stay separate (one annotation each, numbered), so the gaps between them stay uncovered.
        let parts = segments.len();
        for (i, (start, end, _)) in segments.into_iter().enumerate() {
            let name = if parts == 1 { name.clone() } else { format!("{name} (part {} of {parts})", i + 1) };
            out.push(Feature { name, kind: kind.clone(), start, end, direction, color: color.clone() });
        }
    }
    Ok((out, skipped))
}

pub fn to_json(d: &Dna) -> String {
    let mut s = String::with_capacity(d.seq.len() + 64 * d.features.len() + 64);
    s.push_str("{\"seq\":");
    json::string(&mut s, &d.seq);
    s.push_str(if d.circular { ",\"circular\":true" } else { ",\"circular\":false" });
    s.push_str(&format!(",\"skipped\":{}", d.skipped));
    s.push_str(",\"features\":[");
    for (i, f) in d.features.iter().enumerate() {
        if i > 0 {
            s.push(',');
        }
        s.push_str("{\"name\":");
        json::string(&mut s, &f.name);
        s.push_str(",\"type\":");
        json::string(&mut s, &f.kind);
        s.push_str(&format!(",\"start\":{},\"end\":{},\"direction\":{}", f.start, f.end, f.direction));
        if let Some(c) = &f.color {
            s.push_str(",\"color\":");
            json::string(&mut s, c);
        }
        s.push('}');
    }
    s.push_str("]}");
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    fn packet(kind: u8, body: &[u8]) -> Vec<u8> {
        let mut p = vec![kind];
        p.extend((body.len() as u32).to_be_bytes());
        p.extend(body);
        p
    }

    fn cookie(seq_type: u16) -> Vec<u8> {
        let mut b = b"SnapGene".to_vec();
        b.extend(seq_type.to_be_bytes());
        b.extend([0, 15, 0, 19]);
        packet(0x09, &b)
    }

    fn file(seq: &[u8], flags: u8, xml: &str) -> Vec<u8> {
        let mut f = cookie(1);
        f.extend(packet(0x00, &[&[flags][..], seq].concat()));
        f.extend(packet(0x08, b"<AdditionalSequenceProperties/>")); // skipped
        f.extend(packet(0x0A, xml.as_bytes()));
        f
    }

    const XML: &str = r##"<?xml version="1.0"?><Features nextValidID="3">
<Feature recentID="0" name="lac &amp; promoter" directionality="1" type="promoter" swappedSegmentNumbering="1">
  <Segment range="2-5" color="#ffffff" type="standard"/><Q name="note"><V text="x"/></Q></Feature>
<Feature name='rev' directionality="2" type="CDS"><Segment range="9-3" color="#993366"/></Feature>
<Feature name="split" type="misc_feature"><Segment range="1-2"/><Segment range="4-6"/></Feature>
<Feature name="outside" type="misc_feature"><Segment range="5-99"/></Feature>
</Features>"##;

    #[test]
    fn reads_sequence_topology_and_features() {
        let d = parse(&file(b"acgtacgtac", 0x03, XML), 1000).unwrap();
        assert_eq!(d.seq, "ACGTACGTAC");
        assert!(d.circular);
        let f: Vec<_> = d.features.iter().map(|f| (f.name.as_str(), f.kind.as_str(), f.start, f.end, f.direction)).collect();
        assert_eq!(f, vec![
            ("lac & promoter", "promoter", 1, 5, 1),
            ("rev", "CDS", 8, 3, -1), // wraps the origin
            ("split (part 1 of 2)", "misc_feature", 0, 2, 0),
            ("split (part 2 of 2)", "misc_feature", 3, 6, 0),
        ]);
        assert_eq!(d.features[0].color.as_deref(), Some("#ffffff"));
        assert_eq!(d.skipped, 1); // "outside" lies past the end
        assert_eq!(to_json(&d).matches("\"name\"").count(), 4);
    }

    #[test]
    fn coverage_is_kept_exactly() {
        let feats = |xml: &str, flags: u8| {
            let d = parse(&file(b"ACGTACGTAC", flags, xml), 100).unwrap();
            (d.features.iter().map(|f| (f.name.clone(), f.start, f.end)).collect::<Vec<_>>(), d.skipped)
        };
        let one = |seg: &str| format!(r#"<Features><Feature name="f" type="misc_feature">{seg}</Feature></Features>"#);
        // The whole circle written from base 9 covers all ten bases.
        assert_eq!(feats(&one(r#"<Segment range="9-8"/>"#), 1), (vec![("f".into(), 0, 10)], 0));
        // Contiguous segments join, also across the origin.
        assert_eq!(feats(&one(r#"<Segment range="1-3"/><Segment range="4-6"/>"#), 1), (vec![("f".into(), 0, 6)], 0));
        assert_eq!(feats(&one(r#"<Segment range="9-10"/><Segment range="1-2"/>"#), 1), (vec![("f".into(), 8, 2)], 0));
        // Joined segments that go all the way round are the whole circle, never an empty range.
        assert_eq!(feats(&one(r#"<Segment range="9-10"/><Segment range="1-2"/><Segment range="3-8"/>"#), 1), (vec![("f".into(), 0, 10)], 0));
        assert_eq!(feats(&one(r#"<Segment range="5-8"/><Segment range="9-4"/>"#), 1), (vec![("f".into(), 0, 10)], 0));
        // ... also when they go past their start.
        for segs in [r#"<Segment range="9-10"/><Segment range="1-9"/>"#, r#"<Segment range="5-8"/><Segment range="9-6"/>"#,
            r#"<Segment range="1-10"/><Segment range="1-2"/>"#, r#"<Segment range="9-2"/><Segment range="3-9"/>"#] {
            assert_eq!(feats(&one(segs), 1), (vec![("f".into(), 0, 10)], 0), "{segs}");
        }
        // One segment off the sequence skips the whole feature; a wrapping range on a linear sequence too.
        assert_eq!(feats(&one(r#"<Segment range="1-2"/><Segment range="9-12"/>"#), 1), (vec![], 1));
        assert_eq!(feats(&one(r#"<Segment range="9-2"/>"#), 0), (vec![], 1));
    }

    // The bases a feature's annotations cover equal the bases its segments cover, for any segments on a circle.
    #[test]
    fn annotations_cover_exactly_the_segments() {
        let covered = |n: usize, ranges: &[(usize, usize)]| {
            let mut on = vec![false; n];
            for &(s, e) in ranges {
                let len = if e > s { e - s } else { n - s + e };
                for k in 0..len {
                    on[(s + k) % n] = true;
                }
            }
            on
        };
        let mut x = 5u64;
        let mut next = |m: u64| { x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407); (x >> 33) % m };
        for _ in 0..5000 {
            let n = 3 + next(10) as usize;
            let segs: Vec<(usize, usize)> = (0..1 + next(4)).map(|_| (1 + next(n as u64) as usize, 1 + next(n as u64) as usize)).collect();
            let xml = format!(r#"<Features><Feature name="f" type="misc_feature">{}</Feature></Features>"#,
                segs.iter().map(|(a, b)| format!(r#"<Segment range="{a}-{b}"/>"#)).collect::<String>());
            let d = parse(&file(&b"ACGTACGTACGTA"[..n], 1, &xml), 100).unwrap();
            // Expected: the bases of the 1-based inclusive segments (a > b wraps; b = a - 1 is the whole circle).
            let want = covered(n, &segs.iter().map(|&(a, b)| if a <= b { (a - 1, b) } else if b == a - 1 { (0, n) } else { (a - 1, b) }).collect::<Vec<_>>());
            let got = covered(n, &d.features.iter().map(|f| (f.start, f.end)).collect::<Vec<_>>());
            assert_eq!((got, d.skipped), (want, 0), "n={n} {segs:?}");
        }
    }

    #[test]
    fn a_run_of_ampersands_is_read_in_linear_time() {
        let name = "&".repeat(1 << 20);
        let xml = format!(r#"<Features><Feature name="{name}" type="misc_feature"><Segment range="1-2"/></Feature></Features>"#);
        let t = std::time::Instant::now();
        let d = parse(&file(b"ACGTACGTAC", 0, &xml), 100).unwrap();
        assert_eq!(d.features[0].name.len(), 200);
        assert!(t.elapsed().as_secs() < 2, "{:?}", t.elapsed());
    }

    #[test]
    fn linear_without_features() {
        let mut f = cookie(1);
        f.extend(packet(0x00, b"\x00ACGT"));
        let d = parse(&f, 10).unwrap();
        assert!(!d.circular && d.features.is_empty());
        assert_eq!(to_json(&d), r#"{"seq":"ACGT","circular":false,"skipped":0,"features":[]}"#);
    }

    #[test]
    fn refuses_bad_files() {
        let good = file(b"ACGT", 0, XML);
        assert_eq!(parse(&good[..good.len() - 3], 10).unwrap_err(), "a packet is longer than the file");
        assert_eq!(parse(&good[..3], 10).unwrap_err(), "the file ends inside a packet header");
        assert_eq!(parse(b"\x00\x00\x00\x00\x02\x00A", 10).unwrap_err(), "not a SnapGene file");
        assert!(parse(&cookie(2), 10).unwrap_err().starts_with("not a DNA file"));
        assert_eq!(parse(&cookie(1), 10).unwrap_err(), "the file has no sequence");
        assert!(parse(&file(b"ACGTACGTACGT", 0, ""), 10).unwrap_err().starts_with("the sequence is longer"));
        // A length field past the end of the file.
        let mut lie = cookie(1);
        lie.extend([0x00, 0xff, 0xff, 0xff, 0xff, 0x00]);
        assert_eq!(parse(&lie, 10).unwrap_err(), "a packet is longer than the file");
        assert!(parse(&file(b"AC<T", 0, ""), 10).unwrap_err().contains("not bases"));
    }

    #[test]
    fn mutated_files_give_an_error_or_a_result_never_a_panic() {
        let f = file(b"ACGTACGTAC", 0x03, XML);
        let mut x = 11u64;
        for _ in 0..20_000 {
            x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            let mut g = f.clone();
            let at = (x >> 33) as usize % g.len();
            g[at] = (x >> 17) as u8;
            let cut = if x & 1 == 0 { g.len() } else { (x >> 40) as usize % (g.len() + 1) };
            let _ = parse(&g[..cut], 1000);
        }
    }
}
