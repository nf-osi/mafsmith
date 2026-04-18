pub mod normalization;

/// A minimally parsed VCF record — we parse the raw VCF lines ourselves
/// rather than relying on noodles here, so we can handle the heterogeneous
/// FORMAT fields from multiple callers without schema enforcement.
#[derive(Debug, Clone)]
pub struct VcfRecord {
    pub chrom: String,
    pub pos: u64,    // 1-based
    pub id: String,
    pub ref_allele: String,
    pub alt_allele: String, // first ALT only
    pub qual: String,
    pub filter: String,
    pub info: String,
    pub format_keys: Vec<String>,
    pub samples: Vec<Vec<String>>,
    pub sample_names: Vec<String>,
}

impl VcfRecord {
    pub fn sample_values(&self, name: &str) -> Option<Vec<&str>> {
        let idx = self.sample_names.iter().position(|s| s == name)?;
        Some(
            self.samples
                .get(idx)?
                .iter()
                .map(|s| s.as_str())
                .collect(),
        )
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
    initialized: bool,
}

impl<R: std::io::BufRead> VcfReader<R> {
    pub fn new(reader: R) -> Self {
        Self {
            reader,
            header_lines: Vec::new(),
            sample_names: Vec::new(),
            initialized: false,
        }
    }

    fn init(&mut self) -> anyhow::Result<()> {
        use std::io::BufRead;
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
                let cols: Vec<&str> = trimmed.splitn(10, '\t').collect();
                if cols.len() > 9 {
                    self.sample_names = cols[9..].iter().map(|s| s.to_string()).collect();
                }
                self.header_lines.push(trimmed);
                break;
            }
        }
        self.initialized = true;
        Ok(())
    }

    /// Read the next VCF record. Returns None at EOF.
    pub fn next_record(&mut self) -> anyhow::Result<Option<VcfRecord>> {
        use std::io::BufRead;
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
            let cols: Vec<&str> = trimmed.splitn(10, '\t').collect();
            if cols.len() < 8 {
                continue;
            }
            let chrom = cols[0].to_owned();
            let pos: u64 = cols[1].parse().unwrap_or(0);
            let id = cols[2].to_owned();
            let ref_allele = cols[3].to_owned();
            // Take only first ALT allele
            let alt_allele = cols[4].split(',').next().unwrap_or(".").to_owned();
            let qual = cols[5].to_owned();
            let filter = cols[6].to_owned();
            let info = cols[7].to_owned();

            let (format_keys, samples) = if cols.len() > 8 {
                let keys: Vec<String> = cols[8].split(':').map(|s| s.to_owned()).collect();
                let samp: Vec<Vec<String>> = cols[9..]
                    .iter()
                    .map(|s| s.split(':').map(|v| v.to_owned()).collect())
                    .collect();
                (keys, samp)
            } else {
                (vec![], vec![])
            };

            return Ok(Some(VcfRecord {
                chrom,
                pos,
                id,
                ref_allele,
                alt_allele,
                qual,
                filter,
                info,
                format_keys,
                samples,
                sample_names: self.sample_names.clone(),
            }));
        }
    }
}
