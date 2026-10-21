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
    let features = match xml {
        Some(x) => features(&String::from_utf8_lossy(x), seq.len())?,
        None => vec![],
    };
    Ok(Dna { seq, circular, features })
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
        let Some(semi) = rest.find(';').filter(|&j| j <= 10) else {
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

// 1-based inclusive "a-b" (or a single "a") to a 0-based [start, end), checked against the sequence length.
fn range(r: &str, n: usize) -> Option<(usize, usize)> {
    let (a, b) = r.split_once('-').unwrap_or((r, r));
    let (a, b): (usize, usize) = (a.trim().parse().ok()?, b.trim().parse().ok()?);
    (a >= 1 && b >= 1 && a <= n && b <= n).then_some((a - 1, b))
}

fn features(xml: &str, n: usize) -> Result<Vec<Feature>, String> {
    let mut out = Vec::new();
    let mut rest = xml;
    while let Some(i) = rest.find("<Feature ") {
        rest = &rest[i + 9..];
        let close = rest.find('>').ok_or("an unterminated feature tag")?;
        let tag = &rest[..close];
        let body_end = if tag.ends_with('/') { 0 } else { rest.find("</Feature>").ok_or("an unterminated feature")? };
        let body = &rest[close..body_end.max(close)];
        let mut segments = Vec::new();
        let mut color = None;
        let mut b = body;
        while let Some(j) = b.find("<Segment ") {
            b = &b[j + 9..];
            let end = b.find('>').ok_or("an unterminated segment tag")?;
            let seg = &b[..end];
            if let Some(r) = attr(seg, "range").and_then(|r| range(&r, n)) {
                segments.push(r);
                if color.is_none() {
                    color = attr(seg, "color").filter(|c| c.len() <= 32);
                }
            }
            b = &b[end..];
        }
        rest = &rest[body_end.max(close)..];
        let (Some(&(start, _)), Some(&(_, end))) = (segments.first(), segments.last()) else { continue };
        if out.len() == MAX_FEATURES {
            return Err(format!("more than {MAX_FEATURES} features"));
        }
        let direction = match attr(tag, "directionality").as_deref() {
            Some("1") => 1,
            Some("2") => -1,
            _ => 0,
        };
        let name: String = attr(tag, "name").unwrap_or_else(|| "feature".into()).chars().take(200).collect();
        let kind: String = attr(tag, "type").unwrap_or_else(|| "misc_feature".into()).chars().take(40).collect();
        // A whole-circle feature ends where it starts; keep it as start..n rather than an empty wrap.
        let end = if end == start { n } else { end };
        out.push(Feature { name, kind, start, end, direction, color });
    }
    Ok(out)
}

pub fn to_json(d: &Dna) -> String {
    let mut s = String::with_capacity(d.seq.len() + 64 * d.features.len() + 64);
    s.push_str("{\"seq\":");
    json::string(&mut s, &d.seq);
    s.push_str(if d.circular { ",\"circular\":true" } else { ",\"circular\":false" });
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
            ("split", "misc_feature", 0, 6, 0),
        ]);
        assert_eq!(d.features[0].color.as_deref(), Some("#ffffff"));
        assert_eq!(to_json(&d).matches("\"name\"").count(), 3);
    }

    #[test]
    fn linear_without_features() {
        let mut f = cookie(1);
        f.extend(packet(0x00, b"\x00ACGT"));
        let d = parse(&f, 10).unwrap();
        assert!(!d.circular && d.features.is_empty());
        assert_eq!(to_json(&d), r#"{"seq":"ACGT","circular":false,"features":[]}"#);
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
}
