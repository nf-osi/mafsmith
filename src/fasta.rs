use anyhow::{bail, Context, Result};
use memmap2::Mmap;
use std::{
    collections::HashMap,
    fs,
    path::{Path, PathBuf},
};

/// One line of a `.fai` FASTA index.
struct FaiEntry {
    /// Number of bases in the contig.
    length: u64,
    /// Byte offset of the first base of the contig in the FASTA.
    offset: u64,
    /// Bases per line (excluding the line terminator).
    bases_per_line: u64,
    /// Bytes per line (including the line terminator).
    bytes_per_line: u64,
}

/// Random-access reader for an uncompressed, `.fai`-indexed FASTA.
///
/// Opened once and shared across threads: the FASTA is memory-mapped, so a
/// lookup is a slice of already-mapped bytes with no syscall and no shared file
/// cursor, making `&self` safe to use concurrently. This matters because
/// `flanking_bps` needs one lookup per variant — a read syscall per variant
/// measurably slows whole-genome conversion.
pub struct IndexedFasta {
    mmap: Mmap,
    index: HashMap<String, FaiEntry>,
}

impl IndexedFasta {
    /// Open an indexed FASTA. Requires an uncompressed FASTA plus its `.fai`
    /// sidecar, because lookups address raw byte offsets.
    pub fn open(path: &Path) -> Result<Self> {
        let file = fs::File::open(path)
            .with_context(|| format!("Cannot open reference FASTA: {}", path.display()))?;

        // SAFETY: a concurrent external truncation of the reference FASTA could fault
        // on access. Reference genomes are static files, and samtools/htslib rely on
        // the same assumption.
        let mmap = unsafe { Mmap::map(&file) }
            .with_context(|| format!("Cannot memory-map reference FASTA: {}", path.display()))?;

        // Byte offsets in a .fai refer to the *uncompressed* stream, so raw reads
        // against a gzipped FASTA would return garbage. Reject it explicitly.
        if mmap.len() >= 2 && mmap[..2] == [0x1f, 0x8b] {
            bail!(
                "reference FASTA {} is gzip-compressed; an uncompressed FASTA is required \
                 for random access",
                path.display()
            );
        }

        let fai_path = PathBuf::from(format!("{}.fai", path.display()));
        let fai = fs::read_to_string(&fai_path)
            .with_context(|| format!("Cannot read FASTA index: {}", fai_path.display()))?;

        let mut index = HashMap::new();
        for line in fai.lines() {
            let cols: Vec<&str> = line.split('\t').collect();
            if cols.len() < 5 {
                continue;
            }
            let (Ok(length), Ok(offset), Ok(bases_per_line), Ok(bytes_per_line)) = (
                cols[1].parse::<u64>(),
                cols[2].parse::<u64>(),
                cols[3].parse::<u64>(),
                cols[4].parse::<u64>(),
            ) else {
                continue;
            };
            if bases_per_line == 0 {
                continue;
            }
            index.insert(
                cols[0].to_owned(),
                FaiEntry {
                    length,
                    offset,
                    bases_per_line,
                    bytes_per_line,
                },
            );
        }

        if index.is_empty() {
            bail!(
                "FASTA index {} contained no usable entries",
                fai_path.display()
            );
        }

        Ok(Self { mmap, index })
    }

    /// Look up a contig, tolerating a `chr` prefix mismatch between VCF and FASTA.
    fn entry(&self, chrom: &str) -> Option<&FaiEntry> {
        if let Some(e) = self.index.get(chrom) {
            return Some(e);
        }
        match chrom.strip_prefix("chr") {
            Some(stripped) => self.index.get(stripped),
            None => self.index.get(&format!("chr{}", chrom)),
        }
    }

    /// Fetch the reference bases in `[start, end]` (1-based, inclusive), uppercased.
    ///
    /// Coordinates are clamped to the contig, matching `samtools faidx`, which
    /// silently returns a shorter sequence at contig edges.
    pub fn fetch_range(&self, chrom: &str, start: u64, end: u64) -> Option<String> {
        let e = self.entry(chrom)?;
        let start = start.max(1);
        let end = end.min(e.length);
        if start > end {
            return None;
        }

        // Byte offset of a base = contig offset + full lines before it + column within its line.
        let byte_of = |pos: u64| -> u64 {
            let zero_based = pos - 1;
            e.offset
                + (zero_based / e.bases_per_line) * e.bytes_per_line
                + (zero_based % e.bases_per_line)
        };
        let first = byte_of(start) as usize;
        let last = byte_of(end) as usize;
        let raw = self.mmap.get(first..=last)?;

        // The span may straddle line breaks — keep bases only.
        let mut seq = String::with_capacity((end - start + 1) as usize);
        for b in raw {
            if b.is_ascii_alphabetic() {
                seq.push(b.to_ascii_uppercase() as char);
            }
        }
        if seq.is_empty() {
            None
        } else {
            Some(seq)
        }
    }

    /// Fetch a single reference base (uppercased).
    pub fn fetch_base(&self, chrom: &str, pos: u64) -> Option<char> {
        if pos == 0 {
            return None;
        }
        self.fetch_range(chrom, pos, pos)?.chars().next()
    }
}

/// Fetch a single base from an indexed FASTA file, opening it for this call only.
///
/// For repeated lookups build an [`IndexedFasta`] once instead.
pub fn fasta_fetch_base(fasta_path: &Path, chrom: &str, pos: u64) -> Option<char> {
    IndexedFasta::open(fasta_path).ok()?.fetch_base(chrom, pos)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    /// Write a tiny 2-contig FASTA (10 bases per line) plus its .fai.
    fn write_test_fasta(dir: &Path) -> PathBuf {
        let fa = dir.join("t.fa");
        let mut f = fs::File::create(&fa).unwrap();
        // chr1: 25 bases across 3 lines; chr2: 5 bases on 1 line.
        write!(f, ">chr1\nACGTACGTAC\nGGGGCCCCTT\naaaaa\n>chr2\nTTTTT\n").unwrap();
        let mut fai = fs::File::create(dir.join("t.fa.fai")).unwrap();
        // name, length, offset of first base, bases/line, bytes/line
        writeln!(fai, "chr1\t25\t6\t10\t11").unwrap();
        writeln!(fai, "chr2\t5\t40\t10\t11").unwrap();
        fa
    }

    #[test]
    fn fetch_range_spans_lines() {
        let dir = tempfile::tempdir().unwrap();
        let fa = write_test_fasta(dir.path());
        let f = IndexedFasta::open(&fa).unwrap();
        assert_eq!(f.fetch_range("chr1", 1, 10).unwrap(), "ACGTACGTAC");
        // Crosses a newline
        assert_eq!(f.fetch_range("chr1", 9, 12).unwrap(), "ACGG");
        // Lowercase bases are uppercased, like samtools output is normalized by vcf2maf.pl
        assert_eq!(f.fetch_range("chr1", 21, 25).unwrap(), "AAAAA");
        assert_eq!(f.fetch_base("chr1", 11), Some('G'));
    }

    #[test]
    fn fetch_range_clamps_to_contig() {
        let dir = tempfile::tempdir().unwrap();
        let fa = write_test_fasta(dir.path());
        let f = IndexedFasta::open(&fa).unwrap();
        // start=0 clamps to 1 (variant at position 1 has no upstream flank)
        assert_eq!(f.fetch_range("chr2", 0, 2).unwrap(), "TT");
        // end past the contig clamps to its length
        assert_eq!(f.fetch_range("chr2", 4, 99).unwrap(), "TT");
        assert_eq!(f.fetch_range("chr2", 99, 100), None);
    }

    #[test]
    fn chr_prefix_mismatch_tolerated() {
        let dir = tempfile::tempdir().unwrap();
        let fa = write_test_fasta(dir.path());
        let f = IndexedFasta::open(&fa).unwrap();
        // VCF says "1", FASTA says "chr1"
        assert_eq!(f.fetch_base("1", 1), Some('A'));
        assert_eq!(f.fetch_base("nonexistent", 1), None);
    }

    #[test]
    fn gzipped_fasta_rejected() {
        let dir = tempfile::tempdir().unwrap();
        let fa = dir.path().join("g.fa.gz");
        fs::write(&fa, [0x1f, 0x8b, 0x08, 0x00]).unwrap();
        assert!(IndexedFasta::open(&fa).is_err());
    }
}
