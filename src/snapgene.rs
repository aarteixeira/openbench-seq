//! SnapGene `.dna` files: a sequence of packets (type byte, big-endian u32 length, data). Read here: the cookie
//! (0x09, "SnapGene", sequence type 1 = DNA), the sequence (0x00: a flags byte, bit 0 = circular, then ASCII bases),
//! the features (0x0A: XML) and the primers (0x05: XML, their binding sites as `primer_bind` features). Other packets
//! are skipped by length. Checked against Biopython's SnapGene test files and expectations (tests below).

use crate::json;

/// The most features (parts and primer sites included) a file brings; more is refused while reading, so a file of tiny
/// features cannot build millions of them. OpenBench's text readers use the same limit (src/lib/seq.ts MAX_FEATURES).
pub const MAX_FEATURES: usize = 100_000;

fn too_many() -> String {
    format!("the file has more than {MAX_FEATURES} features; at most {MAX_FEATURES} are imported")
}

#[derive(Debug, PartialEq)]
pub struct Feature {
    pub name: String,
    pub kind: String,
    /// 0-based start; `end` exclusive; `end < start` wraps the origin.
    pub start: usize,
    pub end: usize,
    pub direction: i8,
    pub color: Option<String>,
    /// The parts of one feature share a group (written as "g<n>"), so a reader can join them again: a CDS in parts is
    /// translated across its parts, not part by part.
    pub group: Option<usize>,
}

#[derive(Debug, PartialEq)]
pub struct Dna {
    pub seq: String,
    pub circular: bool,
    pub features: Vec<Feature>,
    /// Features left out because none of their segments lies within the sequence; the caller says so.
    pub skipped: usize,
    /// Names (200 characters) and types (40) cut to OpenBench's stored lengths; the caller says so.
    pub cut: usize,
}

pub fn parse(bytes: &[u8], max_len: usize) -> Result<Dna, String> {
    let mut at = 0;
    let mut cookie = false;
    let mut dna: Option<(String, bool)> = None;
    let mut xml: Option<&[u8]> = None;
    let mut primer_xml: Option<&[u8]> = None;
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
            0x05 => primer_xml = Some(body),
            _ => {}
        }
    }
    let (seq, circular) = dna.ok_or("the file has no sequence")?;
    let mut cut = 0;
    let (mut features, mut skipped) = match xml {
        Some(x) => features(std::str::from_utf8(x).map_err(|_| "the feature list is not UTF-8 text")?, seq.len(), circular, &mut cut)?,
        None => (vec![], 0),
    };
    if let Some(x) = primer_xml {
        let (found, left_out) = primers(std::str::from_utf8(x).map_err(|_| "the primer list is not UTF-8 text")?, seq.len(), circular, MAX_FEATURES - features.len(), &mut cut)?;
        features.extend(found);
        skipped += left_out;
    }
    Ok(Dna { seq, circular, features, skipped, cut })
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
// a > b wraps the origin (end < start), and b = a - 1 is the whole circle from a: (a - 1, a - 1), its start kept (a CDS
// reads from it); on a linear sequence a > b is invalid.
fn range(r: &str, n: usize, circular: bool) -> Option<(usize, usize)> {
    let (a, b) = r.split_once('-').unwrap_or((r, r));
    let (a, b): (usize, usize) = (a.trim().parse().ok()?, b.trim().parse().ok()?);
    if a < 1 || b < 1 || a > n || b > n {
        return None;
    }
    match (a <= b, circular) {
        (true, _) => Some((a - 1, b)),
        (false, false) => None,
        (false, true) if b == a - 1 => Some((a - 1, a - 1)),
        (false, true) => Some((a - 1, b)),
    }
}

/// `s` cut to `max` characters, counted in `cut` when it was longer.
fn bounded(s: String, max: usize, cut: &mut usize) -> String {
    match s.char_indices().nth(max) {
        Some((at, _)) => {
            *cut += 1;
            s[..at].to_string()
        }
        None => s,
    }
}

fn features(xml: &str, n: usize, circular: bool, cut: &mut usize) -> Result<(Vec<Feature>, usize), String> {
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
            // A gap segment marks the stretch between two parts: not covered by the feature.
            if attr(seg, "type").as_deref() == Some("gap") {
                b = &b[end..];
                continue;
            }
            match attr(seg, "range").and_then(|r| range(&r, n, circular)) {
                Some(r) => {
                    let len = if r.1 > r.0 { r.1 - r.0 } else { n - r.0 + r.1 };
                    // A segment that continues the previous one (also across the origin of a circle) extends it; a
                    // piece that has covered a whole turn is the whole circle and stays so.
                    match segments.last_mut() {
                        Some(last) if last.2 >= n => {}
                        Some(last) if last.1 == r.0 || (circular && last.1 == n && r.0 == 0) => {
                            *last = if last.2 + len >= n { (last.0, last.0, n) } else { (last.0, r.1, last.2 + len) };
                        }
                        _ => {
                            if out.len() + segments.len() >= MAX_FEATURES {
                                return Err(too_many());
                            }
                            segments.push((r.0, r.1, len))
                        }
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
        let direction = match attr(tag, "directionality").as_deref() {
            Some("1") => 1,
            Some("2") => -1,
            _ => 0,
        };
        let name = bounded(attr(tag, "name").unwrap_or_else(|| "feature".into()), 200, cut);
        let kind = bounded(attr(tag, "type").unwrap_or_else(|| "misc_feature".into()), 40, cut);
        // A piece round the whole circle is [0, n), except a CDS's that starts elsewhere: its reading starts there, so it
        // stays two parts, from its start to the end and from the origin back to its start (as the GenBank reader keeps it).
        let mut segments: Vec<(usize, usize, usize)> = segments
            .into_iter()
            .flat_map(|(start, end, len)| match (len >= n, start != 0 && kind == "CDS") {
                (true, true) => vec![(start, n, n - start), (0, start, start)],
                (true, false) => vec![(0, n, n)],
                _ => vec![(start, end, len)],
            })
            .collect();
        if out.len() + segments.len() > MAX_FEATURES {
            return Err(too_many());
        }
        // Separate segments stay separate (one annotation each, numbered), so the gaps between them stay uncovered. Parts
        // are listed and numbered in reading order, part 1 first: from the last segment on the reverse strand.
        let parts = segments.len();
        let group = (parts > 1).then_some(out.len());
        if direction == -1 {
            segments.reverse();
        }
        for (i, (start, end, _)) in segments.into_iter().enumerate() {
            let k = i + 1;
            let name = if parts == 1 { name.clone() } else { format!("{name} (part {k} of {parts})") };
            out.push(Feature { name, kind: kind.clone(), start, end, direction, color: color.clone(), group });
        }
    }
    Ok((out, skipped))
}

// Primer binding sites as `primer_bind` features, as SnapGene shows them: a site whose annealed bases or melting
// temperature fall short of the file's hybridization parameters is not shown, nor a "simplified" site repeating one
// already listed. Binding-site locations are 0-based and inclusive ("751-775" is [751, 776)), unlike segments; a
// site bound to the bottom strand (boundStrand 1) is a reverse primer. A site outside the sequence is left out.
fn primers(xml: &str, n: usize, circular: bool, room: usize, cut: &mut usize) -> Result<(Vec<Feature>, usize), String> {
    let num = |s: Option<String>| s.and_then(|v| v.trim().parse::<i64>().ok());
    let (mut min_len, mut min_tm) = (0, 0);
    if let Some(i) = xml.find("<HybridizationParams") {
        let tag = &xml[i + 20..xml[i..].find('>').map_or(xml.len(), |e| i + e)];
        min_len = num(attr(tag, "minContinuousMatchLen")).unwrap_or(0);
        min_tm = num(attr(tag, "minMeltingTemperature")).unwrap_or(0);
    }
    let mut out = Vec::new();
    let mut skipped = 0;
    let mut rest = xml;
    // An unfinished tag or primer is a damaged file, refused as the feature list's are: nothing is dropped silently.
    while let Some(i) = rest.find("<Primer ") {
        rest = &rest[i + 8..];
        let close = rest.find('>').ok_or("an unterminated primer tag")?;
        let tag = &rest[..close];
        let body_end = if tag.ends_with('/') { close } else { rest.find("</Primer>").ok_or("an unterminated primer")? };
        let body = &rest[close..body_end.max(close)];
        let name = bounded(attr(tag, "name").unwrap_or_else(|| "primer".into()), 200, cut);
        // Sites listed so far, looked up in constant time: a primer with many sites stays linear.
        let mut seen: std::collections::HashSet<(usize, usize, i8)> = std::collections::HashSet::new();
        let mut b = body;
        while let Some(j) = b.find("<BindingSite ") {
            b = &b[j + 13..];
            let end = b.find('>').ok_or("an unterminated binding site tag")?;
            let site = &b[..end];
            b = &b[end..];
            // A site whose location is not on the sequence is left out and counted, as a feature's is.
            let Some((s0, s1)) = attr(site, "location").and_then(|l| {
                let (a, z) = l.split_once('-')?;
                Some((a.trim().parse::<usize>().ok()?, z.trim().parse::<usize>().ok()?))
            }) else {
                skipped += 1;
                continue;
            };
            if s0 >= n || s1 >= n || (s0 > s1 && !circular) {
                skipped += 1;
                continue;
            }
            // Across the origin of a circle when s0 > s1; one base short of a full turn reads as the whole circle.
            let (start, end) = if s0 <= s1 || s1 + 1 != s0 { (s0, s1 + 1) } else { (0, n) };
            let direction = match attr(site, "boundStrand").as_deref() {
                Some("1") => -1,
                None | Some("0") => 1,
                Some(_) => return Err("a binding site's strand is neither 0 nor 1".into()),
            };
            if attr(site, "simplified").as_deref() == Some("1") && seen.contains(&(start, end, direction)) {
                continue;
            }
            if attr(site, "annealedBases").is_some_and(|a| (a.len() as i64) < min_len)
                || num(attr(site, "meltingTemperature")).is_some_and(|t| t < min_tm)
            {
                continue;
            }
            if out.len() == room {
                return Err(too_many());
            }
            seen.insert((start, end, direction));
            out.push(Feature { name: name.clone(), kind: "primer_bind".into(), start, end, direction, color: None, group: None });
        }
        rest = &rest[body_end.max(close)..];
    }
    Ok((out, skipped))
}

pub fn to_json(d: &Dna) -> String {
    let mut s = String::with_capacity(d.seq.len() + 64 * d.features.len() + 64);
    s.push_str("{\"seq\":");
    json::string(&mut s, &d.seq);
    s.push_str(if d.circular { ",\"circular\":true" } else { ",\"circular\":false" });
    s.push_str(&format!(",\"skipped\":{},\"cut\":{}", d.skipped, d.cut));
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
        if let Some(g) = f.group {
            s.push_str(&format!(",\"group\":\"g{g}\""));
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
        assert_eq!((d.features[0].name.len(), d.cut), (200, 1));
        assert!(t.elapsed().as_secs() < 2, "{:?}", t.elapsed());
    }

    #[test]
    fn linear_without_features() {
        let mut f = cookie(1);
        f.extend(packet(0x00, b"\x00ACGT"));
        let d = parse(&f, 10).unwrap();
        assert!(!d.circular && d.features.is_empty());
        assert_eq!(to_json(&d), r#"{"seq":"ACGT","circular":false,"skipped":0,"cut":0,"features":[]}"#);
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

    #[test]
    fn gap_segments_are_not_covered_and_reverse_parts_count_from_the_end() {
        let xml = r#"<Features><Feature name="t" directionality="2" type="terminator"><Segment range="3-4" type="standard"/><Segment range="5-6" type="gap"/><Segment range="7-8" type="standard"/><Segment range="9-10" type="standard"/></Feature></Features>"#;
        let d = parse(&file(b"ACGTACGTACGT", 0, xml), 100).unwrap();
        let got: Vec<_> = d.features.iter().map(|f| (f.name.as_str(), f.start, f.end, f.direction, f.group)).collect();
        assert_eq!(got, [("t (part 1 of 2)", 6, 10, -1, Some(0)), ("t (part 2 of 2)", 2, 4, -1, Some(0))]);
        assert!(to_json(&d).contains(r#""direction":-1,"group":"g0"}"#));
    }

    #[test]
    fn primer_sites_are_zero_based_filtered_and_stranded() {
        let primers = r#"<Primers><HybridizationParams minContinuousMatchLen="10" minMeltingTemperature="40"/>
<Primer name="P1"><BindingSite location="2-5" boundStrand="0" annealedBases="ACGTACGTACGT" meltingTemperature="52"/>
<BindingSite simplified="1" location="2-5" boundStrand="0" annealedBases="ACGTACGTACGT" meltingTemperature="52"/>
<BindingSite location="8-1" boundStrand="1" annealedBases="ACGTACGTACGT" meltingTemperature="60"/>
<BindingSite location="3-6" boundStrand="0" annealedBases="ACG" meltingTemperature="60"/>
<BindingSite location="3-7" boundStrand="0" annealedBases="ACGTACGTACGT" meltingTemperature="20"/>
<BindingSite location="3-70" boundStrand="0" annealedBases="ACGTACGTACGT" meltingTemperature="60"/></Primer></Primers>"#;
        let mut f = file(b"ACGTACGTACGT", 1, "<Features/>");
        f.extend(packet(0x05, primers.as_bytes()));
        let d = parse(&f, 100).unwrap();
        let got: Vec<_> = d.features.iter().map(|f| (f.kind.as_str(), f.start, f.end, f.direction)).collect();
        // [2, 6) forward once (the simplified repeat dropped); [8, 2) across the origin, reverse; the short, cold and
        // out-of-range sites left out (the out-of-range one counted as skipped; the short and cold ones are filtered, as
        // SnapGene shows them).
        assert_eq!(got, [("primer_bind", 2, 6, 1), ("primer_bind", 8, 2, -1)]);
        assert_eq!(d.skipped, 1);
    }

    #[test]
    fn damaged_primer_lists_are_refused() {
        let with = |primers: &str| {
            let mut f = file(b"ACGTACGTACGT", 1, "<Features/>");
            f.extend(packet(0x05, primers.as_bytes()));
            parse(&f, 100)
        };
        // Framed packets whose XML stops inside a tag or before </Primer>, and a strand that is neither 0 nor 1.
        assert_eq!(with(r#"<Primers><Primer name="p"><BindingSite location="1-2""#).unwrap_err(), "an unterminated primer");
        assert_eq!(with(r#"<Primers><Primer name="p"><BindingSite location="1-2"</Primer></Primers>"#).unwrap_err(), "an unterminated binding site tag");
        assert_eq!(with(r#"<Primers><Primer name="p"><BindingSite location="1-2"/>"#).unwrap_err(), "an unterminated primer");
        assert_eq!(with(r#"<Primers><Primer name="p""#).unwrap_err(), "an unterminated primer tag");
        assert_eq!(with(r#"<Primers><Primer name="p"><BindingSite location="1-2" boundStrand="banana"/></Primer></Primers>"#).unwrap_err(),
            "a binding site's strand is neither 0 nor 1");
    }

    #[test]
    fn a_cds_round_the_whole_circle_keeps_its_start() {
        // AATGAAACCCTA, circular, a CDS 2-1: it reads from base 2 (ATG AAA CCC TAA), as two parts in reading order; another
        // type round the circle is [0, n).
        let cds = |dir: &str, kind: &str| format!(r#"<Features><Feature name="c" directionality="{dir}" type="{kind}"><Segment range="2-1"/></Feature></Features>"#);
        let got = |xml: String| parse(&file(b"AATGAAACCCTA", 1, &xml), 100).unwrap().features.iter().map(|f| (f.name.clone(), f.start, f.end, f.group)).collect::<Vec<_>>();
        assert_eq!(got(cds("1", "CDS")), [("c (part 1 of 2)".into(), 1, 12, Some(0)), ("c (part 2 of 2)".into(), 0, 1, Some(0))]);
        // On the reverse strand it reads from base 1 down: the part ending at its start comes first.
        assert_eq!(got(cds("2", "CDS")), [("c (part 1 of 2)".into(), 0, 1, Some(0)), ("c (part 2 of 2)".into(), 1, 12, Some(0))]);
        assert_eq!(got(cds("1", "gene")), [("c".into(), 0, 12, None)]);
    }

    #[test]
    fn features_parts_and_primer_sites_stop_at_the_limit() {
        let seg = r#"<Segment range="1-2"/>"#;
        // MAX_FEATURES one-part features are read; one more feature, or one more part, is refused.
        let feats = |k: usize, parts: usize| format!("<Features>{}</Features>", format!(r#"<Feature name="f" type="misc_feature">{}</Feature>"#,
            [seg, r#"<Segment range="4-5"/>"#][..parts].concat()).repeat(k));
        assert_eq!(parse(&file(b"ACGTACGTAC", 0, &feats(MAX_FEATURES, 1)), 100).unwrap().features.len(), MAX_FEATURES);
        assert_eq!(parse(&file(b"ACGTACGTAC", 0, &feats(MAX_FEATURES + 1, 1)), 100).unwrap_err(), too_many());
        assert_eq!(parse(&file(b"ACGTACGTAC", 0, &feats(MAX_FEATURES / 2 + 1, 2)), 100).unwrap_err(), too_many());
        // Primer sites count with the features.
        let site = r#"<Primer name="p"><BindingSite location="1-3" boundStrand="0"/></Primer>"#;
        let with_primers = |k: usize| { let mut f = file(b"ACGTACGTAC", 0, &feats(MAX_FEATURES - 1, 1)); f.extend(packet(0x05, format!("<Primers>{}</Primers>", site.repeat(k)).as_bytes())); f };
        assert_eq!(parse(&with_primers(1), 100).unwrap().features.len(), MAX_FEATURES);
        assert_eq!(parse(&with_primers(2), 100).unwrap_err(), too_many());
    }

    #[test]
    fn xml_that_is_not_utf8_is_refused_and_many_sites_stay_linear() {
        let mut bad = file(b"ACGTACGTAC", 0, "<Features><Feature name=\"x\" type=\"misc_feature\"><Segment range=\"1-2\"/></Feature></Features>");
        let at = bad.windows(4).position(|w| w == b"\"x\" ").unwrap() + 1;
        bad[at] = 0xff;
        assert_eq!(parse(&bad, 100).unwrap_err(), "the feature list is not UTF-8 text");
        // 40,000 distinct sites of one primer, each also listed as its simplified repeat: one pass, not one scan per site.
        let sites: String = (0..40_000).map(|i| format!(r#"<BindingSite location="{i}-{}" boundStrand="0"/><BindingSite simplified="1" location="{i}-{}" boundStrand="0"/>"#, i + 1, i + 1)).collect();
        let mut f = file(&b"ACGT".repeat(10_001), 0, "<Features/>");
        f.extend(packet(0x05, format!(r#"<Primers><Primer name="p">{sites}</Primer></Primers>"#).as_bytes()));
        let t = std::time::Instant::now();
        assert_eq!(parse(&f, 100_000).unwrap().features.len(), 40_000);
        assert!(t.elapsed().as_secs() < 2, "{:?}", t.elapsed());
    }
}
