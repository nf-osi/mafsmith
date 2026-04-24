use std::{
    fs,
    io::{Read, Seek, SeekFrom},
    path::{Path, PathBuf},
};

/// Fetch a single base from an indexed FASTA file.
///
/// Reads the .fai sidecar to locate the byte offset, then seeks directly
/// into the FASTA. Returns None if the chromosome or position is not found.
pub fn fasta_fetch_base(fasta_path: &Path, chrom: &str, pos: u64) -> Option<char> {
    let fai_path = PathBuf::from(format!("{}.fai", fasta_path.display()));
    let fai = fs::read_to_string(&fai_path).ok()?;

    let lookup = |name: &str| -> Option<(u64, u64, u64)> {
        for line in fai.lines() {
            let cols: Vec<&str> = line.split('\t').collect();
            if cols.len() >= 5 && cols[0] == name {
                let offset: u64 = cols[2].parse().ok()?;
                let bases_per_line: u64 = cols[3].parse().ok()?;
                let bytes_per_line: u64 = cols[4].parse().ok()?;
                return Some((offset, bases_per_line, bytes_per_line));
            }
        }
        None
    };

    let (offset, bases_per_line, bytes_per_line) = lookup(chrom).or_else(|| {
        if let Some(stripped) = chrom.strip_prefix("chr") {
            lookup(stripped)
        } else {
            lookup(&format!("chr{}", chrom))
        }
    })?;

    if bases_per_line == 0 || pos == 0 {
        return None;
    }
    let line_num = (pos - 1) / bases_per_line;
    let col = (pos - 1) % bases_per_line;
    let byte_pos = offset + line_num * bytes_per_line + col;

    let mut f = fs::File::open(fasta_path).ok()?;
    f.seek(SeekFrom::Start(byte_pos)).ok()?;
    let mut buf = [0u8; 1];
    f.read_exact(&mut buf).ok()?;
    let c = buf[0] as char;
    if c.is_ascii_alphabetic() {
        Some(c.to_ascii_uppercase())
    } else {
        None
    }
}
