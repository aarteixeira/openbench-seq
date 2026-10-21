//! ABIF (`.ab1`) Sanger trace files: "ABIF", a version, then a 28-byte root directory entry pointing at the
//! directory. Each entry: name (4), number (i32), element type (i16), element size (i16), element count (i32),
//! data size (i32), data offset (i32), handle (i32), big-endian; data of 4 bytes or fewer sits in the offset field.
//! Read: base calls (PBAS 2, else 1), qualities (PCON 2, else 1), peak positions (PLOC 2, else 1), the four
//! channels (DATA 9-12, in the base order of FWO_ 1) and the sample name (SMPL 1).

use crate::json;

pub const MAX_BASES: usize = 100_000;
pub const MAX_POINTS: usize = 1_000_000;

#[derive(Debug, PartialEq)]
pub struct Trace {
    pub name: String,
    pub bases: String,
    pub quality: Vec<u8>,
    pub peaks: Vec<i32>,
    /// Channels in A, C, G, T order.
    pub channels: [Vec<i16>; 4],
}

struct Entry<'a> {
    name: [u8; 4],
    number: i32,
    count: usize,
    elsize: usize,
    data: &'a [u8],
}

fn be32(b: &[u8], at: usize) -> Result<i32, String> {
    b.get(at..at + 4).map(|x| i32::from_be_bytes(x.try_into().unwrap())).ok_or_else(|| "the file is cut short".to_string())
}

fn entry<'a>(f: &'a [u8], at: usize) -> Result<Entry<'a>, String> {
    let e = f.get(at..at + 28).ok_or("the directory runs past the end of the file")?;
    let count = i32::from_be_bytes(e[12..16].try_into().unwrap());
    let elsize = i16::from_be_bytes([e[10], e[11]]);
    let size = i32::from_be_bytes(e[16..20].try_into().unwrap());
    if count < 0 || elsize < 0 || size < 0 {
        return Err("a directory entry has a negative size".into());
    }
    let size = size as usize;
    let data = if size <= 4 {
        &e[20..20 + size]
    } else {
        let off = i32::from_be_bytes(e[20..24].try_into().unwrap());
        let off = usize::try_from(off).map_err(|_| "a directory entry has a negative offset")?;
        f.get(off..off.checked_add(size).ok_or("a directory entry overflows")?).ok_or("a directory entry points past the end of the file")?
    };
    let (count, elsize) = (count as usize, elsize as usize);
    if count.checked_mul(elsize).is_none_or(|n| n > data.len()) {
        return Err("a directory entry is larger than its data".into());
    }
    Ok(Entry { name: e[..4].try_into().unwrap(), number: i32::from_be_bytes(e[4..8].try_into().unwrap()), count, elsize, data })
}

pub fn parse(f: &[u8]) -> Result<Trace, String> {
    if f.len() < 34 || &f[..4] != b"ABIF" {
        return Err("not an ABIF trace file".into());
    }
    let root = entry(f, 6)?;
    let dir = usize::try_from(be32(f, 6 + 20)?).map_err(|_| "the directory offset is negative")?;
    if root.count > 100_000 {
        return Err("the directory is too large".into());
    }
    let entries: Vec<Entry> = (0..root.count).map(|i| entry(f, dir + 28 * i)).collect::<Result<_, _>>()?;
    let find = |name: &[u8; 4], numbers: &[i32]| numbers.iter().find_map(|&n| entries.iter().find(|e| &e.name == name && e.number == n));

    let bases = find(b"PBAS", &[2, 1]).ok_or("the file has no base calls")?;
    if bases.count > MAX_BASES {
        return Err(format!("more than {MAX_BASES} base calls"));
    }
    let bases_text = &bases.data[..bases.count];
    if !bases_text.iter().all(|b| b.is_ascii_alphabetic() || *b == b'-') {
        return Err("the base calls hold characters that are not bases".into());
    }
    let quality = find(b"PCON", &[2, 1]).map(|e| e.data[..e.count].to_vec()).unwrap_or_default();
    let peaks = match find(b"PLOC", &[2, 1]) {
        Some(e) if e.elsize == 2 => e.data[..2 * e.count].chunks_exact(2).map(|c| i16::from_be_bytes([c[0], c[1]]) as i32).collect(),
        Some(_) => return Err("peak positions have an unexpected element size".into()),
        None => vec![],
    };
    let order = find(b"FWO_", &[1]).map(|e| e.data[..e.count.min(4)].to_vec()).unwrap_or_else(|| b"GATC".to_vec());
    let mut channels: [Vec<i16>; 4] = Default::default();
    for (i, base) in order.iter().enumerate() {
        let slot = match base.to_ascii_uppercase() {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => return Err("the channel order (FWO_) is not a permutation of ACGT".into()),
        };
        let Some(e) = find(b"DATA", &[9 + i as i32]) else { continue };
        if e.elsize != 2 || e.count > MAX_POINTS {
            return Err("a trace channel has an unexpected size".into());
        }
        if !channels[slot].is_empty() {
            return Err("the channel order (FWO_) is not a permutation of ACGT".into());
        }
        channels[slot] = e.data[..2 * e.count].chunks_exact(2).map(|c| i16::from_be_bytes([c[0], c[1]])).collect();
    }
    let name = find(b"SMPL", &[1])
        .and_then(|e| e.data.split_first().map(|(&l, rest)| String::from_utf8_lossy(&rest[..(l as usize).min(rest.len())]).into_owned()))
        .unwrap_or_default();
    Ok(Trace { name, bases: String::from_utf8_lossy(bases_text).to_ascii_uppercase(), quality, peaks, channels })
}

pub fn to_json(t: &Trace) -> String {
    let mut s = String::with_capacity(t.bases.len() * 12 + t.channels.iter().map(Vec::len).sum::<usize>() * 6 + 64);
    s.push_str("{\"name\":");
    json::string(&mut s, &t.name);
    s.push_str(",\"bases\":");
    json::string(&mut s, &t.bases);
    s.push_str(",\"quality\":");
    json::numbers(&mut s, t.quality.iter());
    s.push_str(",\"peaks\":");
    json::numbers(&mut s, t.peaks.iter());
    s.push_str(",\"channels\":{");
    for (i, (b, c)) in ["A", "C", "G", "T"].iter().zip(&t.channels).enumerate() {
        if i > 0 {
            s.push(',');
        }
        s.push_str(&format!("\"{b}\":"));
        json::numbers(&mut s, c.iter());
    }
    s.push_str("}}");
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    // A hand-built ABIF file: entries (name, number, element type, element size, count, data).
    fn abif(tags: &[(&[u8; 4], i32, i16, i16, &[u8])]) -> Vec<u8> {
        let mut data: Vec<u8> = Vec::new();
        let data_start = 34usize;
        let mut dir = Vec::new();
        for (name, number, kind, elsize, bytes) in tags {
            let count = bytes.len() / *elsize as usize;
            dir.extend(*name);
            dir.extend(number.to_be_bytes());
            dir.extend(kind.to_be_bytes());
            dir.extend(elsize.to_be_bytes());
            dir.extend((count as i32).to_be_bytes());
            dir.extend((bytes.len() as i32).to_be_bytes());
            if bytes.len() <= 4 {
                let mut inline = bytes.to_vec();
                inline.resize(4, 0);
                dir.extend(inline);
            } else {
                dir.extend(((data_start + data.len()) as i32).to_be_bytes());
                data.extend(*bytes);
            }
            dir.extend([0; 4]);
        }
        let dir_at = data_start + data.len();
        let mut f = b"ABIF".to_vec();
        f.extend(101i16.to_be_bytes());
        f.extend(b"tdir");
        f.extend(1i32.to_be_bytes());
        f.extend(1023i16.to_be_bytes());
        f.extend(28i16.to_be_bytes());
        f.extend((tags.len() as i32).to_be_bytes());
        f.extend(((28 * tags.len()) as i32).to_be_bytes());
        f.extend((dir_at as i32).to_be_bytes());
        f.extend([0; 4]);
        f.extend(data);
        f.extend(dir);
        f
    }

    fn shorts(xs: &[i16]) -> Vec<u8> {
        xs.iter().flat_map(|x| x.to_be_bytes()).collect()
    }

    fn sample() -> Vec<u8> {
        let (g, a, t, c) = (shorts(&[1, 2, 3]), shorts(&[4, 5, 6]), shorts(&[7, 8, 9]), shorts(&[10, 11, 12]));
        abif(&[
            (b"PBAS", 1, 2, 1, b"NNNNNN"), // superseded by PBAS 2
            (b"PBAS", 2, 2, 1, b"acgtn"),
            (b"PCON", 2, 2, 1, &[40, 38, 20, 5, 0]),
            (b"PLOC", 2, 4, 2, &shorts(&[10, 22, 35, 47, 60])),
            (b"FWO_", 1, 2, 1, b"GATC"),
            (b"DATA", 9, 4, 2, &g),
            (b"DATA", 10, 4, 2, &a),
            (b"DATA", 11, 4, 2, &t),
            (b"DATA", 12, 4, 2, &c),
            (b"SMPL", 1, 18, 1, b"\x07clone12"),
        ])
    }

    #[test]
    fn reads_calls_qualities_peaks_and_channels() {
        let t = parse(&sample()).unwrap();
        assert_eq!(t.name, "clone12");
        assert_eq!(t.bases, "ACGTN");
        assert_eq!(t.quality, vec![40, 38, 20, 5, 0]);
        assert_eq!(t.peaks, vec![10, 22, 35, 47, 60]);
        // FWO_ GATC: DATA9 is G, 10 A, 11 T, 12 C; channels come back in A, C, G, T order.
        assert_eq!(t.channels, [vec![4, 5, 6], vec![10, 11, 12], vec![1, 2, 3], vec![7, 8, 9]]);
        assert!(to_json(&t).starts_with(r#"{"name":"clone12","bases":"ACGTN","quality":[40,38,20,5,0]"#));
    }

    #[test]
    fn refuses_bad_files() {
        assert_eq!(parse(b"ABIF").unwrap_err(), "not an ABIF trace file");
        assert_eq!(parse(&[0u8; 64]).unwrap_err(), "not an ABIF trace file");
        let f = sample();
        assert!(parse(&f[..f.len() - 10]).is_err()); // directory cut short
        // A data offset past the end of the file.
        let mut lie = f.clone();
        let dir_at = i32::from_be_bytes(lie[26..30].try_into().unwrap()) as usize;
        lie[dir_at + 28 + 20..dir_at + 28 + 24].copy_from_slice(&i32::MAX.to_be_bytes());
        assert!(parse(&lie).unwrap_err().contains("past the end"));
        assert_eq!(parse(&abif(&[(b"PCON", 2, 2, 1, &[1])])).unwrap_err(), "the file has no base calls");
    }
}
