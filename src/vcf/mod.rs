pub mod normalization;

use std::sync::Arc;

/// A minimally parsed VCF record — we parse the raw VCF lines ourselves
/// rather than relying on noodles here, so we can handle the heterogeneous
/// FORMAT fields from multiple callers without schema enforcement.
#[derive(Debug, Clone)]
pub struct VcfRecord {
    pub chrom: String,
    pub pos: u64, // 1-based
    pub id: String,
    pub ref_allele: String,
    pub alt_allele: String,    // first ALT only (used throughout pipeline)
    pub all_alts: Vec<String>, // all ALT alleles for multi-allelic GT handling
    pub qual: String,
    pub filter: String,
    pub info: String,
    /// Shared when FORMAT column is unchanged across records — clone is an Arc refcount bump.
    pub format_keys: Arc<Vec<String>>,
    /// Raw sample column strings ("GT:AD:DP" form), one per sample, split on demand.
    pub samples_raw: Vec<String>,
    /// Shared across all records from the same VCF — clone is an Arc refcount bump.
    pub sample_names: Arc<Vec<String>>,
}

impl VcfRecord {
    /// Return split fields for the named sample column. Splits on demand (no eager parse).
    pub fn sample_values<'a>(&'a self, name: &str) -> Option<Vec<&'a str>> {
        let idx = self.sample_names.iter().position(|s| s == name)?;
        Some(self.samples_raw.get(idx)?.split(':').collect())
    }

    /// Return split fields for a sample by pre-computed index (avoids name linear scan).
    pub fn sample_fields_by_idx(&self, idx: usize) -> Option<Vec<&str>> {
        Some(self.samples_raw.get(idx)?.split(':').collect())
    }

    /// Look up a single INFO field value (handles both FLAG and KEY=VALUE).
    pub fn info_field(&self, key: &str) -> Option<&str> {
        for part in self.info.split(';') {
            if part == key {
                return Some("");
            }
            if let Some(rest) = part.strip_prefix(key) {
                if let Some(val) = rest.strip_prefix('=') {
                    return Some(val);
                }
            }
        }
        None
    }

    /// Return the CSQ INFO value (may contain multiple comma-separated entries).
    pub fn csq_value(&self) -> Option<&str> {
        self.info_field("CSQ")
    }
}

/// Iterate over VCF records in a text stream, yielding header lines and records separately.
pub struct VcfReader<R: std::io::BufRead> {
    reader: R,
    pub header_lines: Vec<String>,
    pub sample_names: Vec<String>,
    sample_names_arc: Arc<Vec<String>>,
    initialized: bool,
    /// Cache last-seen FORMAT string and its parsed Arc so identical FORMAT rows
    /// share one Arc (refcount bump) instead of allocating a fresh Vec<String> per record.
    last_format_str: String,
    last_format_arc: Arc<Vec<String>>,
}

impl<R: std::io::BufRead> VcfReader<R> {
    pub fn new(reader: R) -> Self {
        Self {
            reader,
            header_lines: Vec::new(),
            sample_names: Vec::new(),
            sample_names_arc: Arc::new(Vec::new()),
            initialized: false,
            last_format_str: String::new(),
            last_format_arc: Arc::new(Vec::new()),
        }
    }

    fn init(&mut self) -> anyhow::Result<()> {
        let mut line = String::new();
        loop {
            line.clear();
            let n = self.reader.read_line(&mut line)?;
            if n == 0 {
                break;
            }
            let trimmed = line.trim_end().to_owned();
            if trimmed.starts_with("##") {
                self.header_lines.push(trimmed);
            } else if trimmed.starts_with('#') {
                // #CHROM line
                let cols: Vec<&str> = trimmed.split('\t').collect();
                if cols.len() > 9 {
                    self.sample_names = cols[9..].iter().map(|s| s.to_string()).collect();
                }
                self.sample_names_arc = Arc::new(self.sample_names.clone());
                self.header_lines.push(trimmed);
                break;
            }
        }
        self.initialized = true;
        Ok(())
    }

    /// Read the next VCF record. Returns None at EOF.
    pub fn next_record(&mut self) -> anyhow::Result<Option<VcfRecord>> {
        if !self.initialized {
            self.init()?;
        }
        let mut line = String::new();
        loop {
            line.clear();
            let n = self.reader.read_line(&mut line)?;
            if n == 0 {
                return Ok(None);
            }
            let trimmed = line.trim_end();
            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }
            let cols: Vec<&str> = trimmed.split('\t').collect();
            if cols.len() < 8 {
                continue;
            }
            let chrom = cols[0].to_owned();
            let pos: u64 = cols[1].parse().unwrap_or(0);
            let id = cols[2].to_owned();
            let ref_allele = cols[3].to_owned();
            let all_alts: Vec<String> = cols[4].split(',').map(|s| s.to_owned()).collect();
            let alt_allele = all_alts[0].clone();
            let qual = cols[5].to_owned();
            let filter = cols[6].to_owned();
            let info = cols[7].to_owned();

            let (format_keys, samples_raw) = if cols.len() > 8 {
                let fmt_str = cols[8];
                // Reuse Arc when FORMAT is identical to last record (common case: same caller/file).
                let keys = if fmt_str == self.last_format_str {
                    Arc::clone(&self.last_format_arc)
                } else {
                    let arc =
                        Arc::new(fmt_str.split(':').map(|s| s.to_owned()).collect::<Vec<_>>());
                    self.last_format_str.clear();
                    self.last_format_str.push_str(fmt_str);
                    self.last_format_arc = Arc::clone(&arc);
                    arc
                };
                // Store each sample column as a single raw string; split on demand in sample_fields_by_idx.
                let raw: Vec<String> = cols[9..].iter().map(|s| s.to_string()).collect();
                (keys, raw)
            } else {
                (Arc::new(vec![]), vec![])
            };

            return Ok(Some(VcfRecord {
                chrom,
                pos,
                id,
                ref_allele,
                alt_allele,
                all_alts,
                qual,
                filter,
                info,
                format_keys,
                samples_raw,
                // Arc clone is a single atomic refcount increment — no heap allocation.
                sample_names: Arc::clone(&self.sample_names_arc),
            }));
        }
    }
}
