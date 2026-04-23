use super::record::{MafRecord, STANDARD_COLUMNS};
use anyhow::Result;
use std::io::Write;

pub struct MafWriter<W: Write> {
    inner: W,
}

impl<W: Write> MafWriter<W> {
    pub fn new(mut inner: W, extra_columns: Vec<String>) -> Result<Self> {
        writeln!(inner, "#version 2.4")?;
        let mut headers: Vec<&str> = STANDARD_COLUMNS.to_vec();
        let extra_refs: Vec<&str> = extra_columns.iter().map(|s| s.as_str()).collect();
        headers.extend_from_slice(&extra_refs);
        writeln!(inner, "{}", headers.join("\t"))?;
        Ok(Self { inner })
    }

    /// Write one MAF record directly to the underlying writer.
    /// Avoids Vec<String> + join by writing each field as raw bytes separated by tabs.
    pub fn write_record(&mut self, record: &MafRecord) -> Result<()> {
        let w = &mut self.inner;
        // Write fields with leading tab separator (first field has no leading tab).
        macro_rules! s {
            ($f:expr) => {{
                w.write_all($f.as_bytes())?;
                w.write_all(b"\t")?;
            }};
        }
        macro_rules! n {
            ($v:expr) => {
                write!(w, "{}\t", $v)?;
            };
        }
        s!(record.hugo_symbol);
        s!(record.entrez_gene_id);
        s!(record.center);
        s!(record.ncbi_build);
        s!(record.chromosome);
        n!(record.start_position);
        n!(record.end_position);
        s!(record.strand);
        s!(record.variant_classification);
        s!(record.variant_type);
        s!(record.reference_allele);
        s!(record.tumor_seq_allele1);
        s!(record.tumor_seq_allele2);
        s!(record.dbsnp_rs);
        s!(record.dbsnp_val_status);
        s!(record.tumor_sample_barcode);
        s!(record.matched_norm_sample_barcode);
        s!(record.match_norm_seq_allele1);
        s!(record.match_norm_seq_allele2);
        s!(record.tumor_validation_allele1);
        s!(record.tumor_validation_allele2);
        s!(record.match_norm_validation_allele1);
        s!(record.match_norm_validation_allele2);
        s!(record.verification_status);
        s!(record.validation_status);
        s!(record.mutation_status);
        s!(record.sequencing_phase);
        s!(record.sequence_source);
        s!(record.validation_method);
        s!(record.score);
        s!(record.bam_file);
        s!(record.sequencer);
        s!(record.tumor_sample_uuid);
        s!(record.matched_norm_sample_uuid);
        s!(record.hgvsc);
        s!(record.hgvsp);
        s!(record.hgvsp_short);
        s!(record.transcript_id);
        s!(record.exon_number);
        s!(record.t_depth);
        s!(record.t_ref_count);
        s!(record.t_alt_count);
        s!(record.n_depth);
        s!(record.n_ref_count);
        s!(record.n_alt_count);
        s!(record.all_effects);
        s!(record.vep_allele);
        s!(record.vep_gene);
        s!(record.vep_feature);
        s!(record.vep_biotype);
        s!(record.vep_canonical);
        s!(record.vep_sift);
        // Last standard field — no trailing tab.
        w.write_all(record.vep_polyphen.as_bytes())?;
        // Extra retained-annotation columns.
        for v in record.extra.values() {
            w.write_all(b"\t")?;
            w.write_all(v.as_bytes())?;
        }
        w.write_all(b"\n")?;
        Ok(())
    }

    pub fn flush(&mut self) -> Result<()> {
        self.inner.flush()?;
        Ok(())
    }
}
