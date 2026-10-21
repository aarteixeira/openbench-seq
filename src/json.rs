//! Just enough JSON writing for parser results.

pub fn string(out: &mut String, s: &str) {
    out.push('"');
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            c if (c as u32) < 0x20 || c == '\u{2028}' || c == '\u{2029}' => out.push_str(&format!("\\u{:04x}", c as u32)),
            c => out.push(c),
        }
    }
    out.push('"');
}

pub fn numbers<T: std::fmt::Display>(out: &mut String, xs: impl IntoIterator<Item = T>) {
    out.push('[');
    for (i, x) in xs.into_iter().enumerate() {
        if i > 0 {
            out.push(',');
        }
        out.push_str(&x.to_string());
    }
    out.push(']');
}

#[cfg(test)]
mod tests {
    #[test]
    fn escapes() {
        let mut s = String::new();
        super::string(&mut s, "a\"b\\c\n\u{1}");
        assert_eq!(s, r#""a\"b\\c\u000a\u0001""#);
    }
}
