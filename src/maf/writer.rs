use super::record::{MafRecord, STANDARD_COLUMNS};
use anyhow::Result;
use std::io::Write;

pub struct MafWriter<W: Write> {
    inner: W,
    extra_columns: Vec<String>,
}

impl<W: Write> MafWriter<W> {
    pub fn new(mut inner: W, extra_columns: Vec<String>) -> Result<Self> {
        writeln!(inner, "#version 2.4")?;
        let mut headers: Vec<&str> = STANDARD_COLUMNS.to_vec();
        let extra_refs: Vec<&str> = extra_columns.iter().map(|s| s.as_str()).collect();
        headers.extend_from_slice(&extra_refs);
        writeln!(inner, "{}", headers.join("\t"))?;
        Ok(Self {
            inner,
            extra_columns,
        })
    }

    pub fn write_record(&mut self, record: &MafRecord) -> Result<()> {
        let row = record.to_row();
        writeln!(self.inner, "{}", row.join("\t"))?;
        Ok(())
    }

    pub fn flush(&mut self) -> Result<()> {
        self.inner.flush()?;
        Ok(())
    }
}
