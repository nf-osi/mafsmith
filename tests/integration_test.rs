/// Integration tests for mafsmith: vcf2maf, maf2vcf, vcf2vcf, maf2maf.
///
/// # Running the tests
///
/// Fast tests (no binary needed):
///   cargo test
///
/// Full end-to-end tests against a running mafsmith binary:
///   cargo build && cargo test -- --include-ignored
///
/// Comparison against vcf2maf reference output:
///   scripts/run_comparison.sh
use std::{
    collections::HashMap,
    fs,
    io::{BufRead, BufReader},
    path::{Path, PathBuf},
    process::Command,
};

fn fixtures_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures")
}

fn expected_dir() -> PathBuf {
    fixtures_dir().join("expected")
}

fn mafsmith_bin() -> PathBuf {
    let release = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("target/release/mafsmith");
    let debug = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("target/debug/mafsmith");
    if release.exists() {
        release
    } else {
        debug
    }
}

/// Parse a TSV file (with optional leading `#version` comment) into rows.
fn parse_tsv(path: &Path) -> Vec<HashMap<String, String>> {
    let f = fs::File::open(path).unwrap_or_else(|e| panic!("Cannot open {}: {e}", path.display()));
    let reader = BufReader::new(f);
    let mut headers: Option<Vec<String>> = None;
    let mut rows = Vec::new();

    for line in reader.lines() {
        let line = line.unwrap();
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        if trimmed.starts_with('#') {
            continue; // skip version comment and column-comment lines
        }
        if headers.is_none() {
            headers = Some(trimmed.split('\t').map(|s| s.to_owned()).collect());
            continue;
        }
        let hdrs = headers.as_ref().unwrap();
        let vals: Vec<&str> = trimmed.split('\t').collect();
        let row: HashMap<String, String> = hdrs
            .iter()
            .zip(vals.iter().chain(std::iter::repeat(&"")))
            .map(|(k, v)| (k.clone(), v.to_string()))
            .collect();
        rows.push(row);
    }
    rows
}

/// Columns we compare strictly (depend only on our code, not VEP database versions).
const KEY_COLUMNS: &[&str] = &[
    "Hugo_Symbol",
    "Chromosome",
    "Start_Position",
    "End_Position",
    "Variant_Classification",
    "Variant_Type",
    "Reference_Allele",
    "Tumor_Seq_Allele2",
    "Transcript_ID",
    "Exon_Number",
    "t_depth",
    "t_ref_count",
    "t_alt_count",
    "n_depth",
    "n_ref_count",
    "n_alt_count",
];

// ── Fast fixture-only tests (no binary required) ─────────────────────────────

#[test]
fn expected_tsv_parseable() {
    let rows = parse_tsv(&expected_dir().join("test_b38_key_cols.tsv"));
    assert_eq!(rows.len(), 6, "Expected 6 rows in key-column fixture");

    let r = &rows[0];
    assert_eq!(r["Hugo_Symbol"], "RUNX1");
    assert_eq!(r["Variant_Classification"], "Missense_Mutation");
    assert_eq!(r["Start_Position"], "34792292");
    assert_eq!(r["Reference_Allele"], "A");
    assert_eq!(r["t_alt_count"], "21");

    // Insertion row: position convention and "-" allele
    let r = &rows[1];
    assert_eq!(r["Variant_Type"], "INS");
    assert_eq!(r["Reference_Allele"], "-");
    assert_eq!(r["Start_Position"], "34880555");
    assert_eq!(r["End_Position"], "34880556");

    // Deletion row
    let r = &rows[3];
    assert_eq!(r["Variant_Type"], "DEL");
    assert_eq!(r["Tumor_Seq_Allele2"], "-");
    assert_eq!(r["Start_Position"], "41479183");
    assert_eq!(r["End_Position"], "41479183");

    // Nonsense mutation
    let r = &rows[5];
    assert_eq!(r["Variant_Classification"], "Nonsense_Mutation");
}

#[test]
fn annotated_vcf_fixture_has_csq() {
    let path = fixtures_dir().join("test_b38_annotated.vcf");
    let content =
        fs::read_to_string(&path).unwrap_or_else(|e| panic!("Cannot read annotated VCF: {e}"));
    assert!(
        content.contains("CSQ="),
        "annotated VCF must contain CSQ annotations"
    );
    let data_lines: Vec<&str> = content.lines().filter(|l| !l.starts_with('#')).collect();
    assert_eq!(
        data_lines.len(),
        6,
        "Expected 6 variants in annotated VCF fixture"
    );
}

#[test]
fn test_b38_vcf_fixture_has_25_variants() {
    let path = fixtures_dir().join("test_b38.vcf");
    let content = fs::read_to_string(&path).unwrap_or_else(|e| panic!("Cannot read test VCF: {e}"));
    let count = content.lines().filter(|l| !l.starts_with('#')).count();
    assert_eq!(count, 25);
}

// ── End-to-end tests (require compiled mafsmith binary) ──────────────────────

/// Run mafsmith vcf2maf on the pre-annotated fixture and compare key columns
/// to the expected TSV.
///
/// Requires `cargo build` to have been run first.
/// Enable with: `cargo test -- --include-ignored`
#[test]
#[ignore = "requires compiled mafsmith binary; run `cargo build` first"]
fn vcf2maf_key_columns_match_expected() {
    let bin = mafsmith_bin();
    assert!(
        bin.exists(),
        "mafsmith binary not found at {}; run `cargo build`",
        bin.display()
    );

    let input_vcf = fixtures_dir().join("test_b38_annotated.vcf");
    let expected_tsv = expected_dir().join("test_b38_key_cols.tsv");
    let tmp = tempfile::NamedTempFile::new().unwrap();

    let status = Command::new(&bin)
        .args([
            "vcf2maf",
            "--input-vcf",
            input_vcf.to_str().unwrap(),
            "--output-maf",
            tmp.path().to_str().unwrap(),
            "--genome",
            "grch38",
            "--tumor-id",
            "TUMOR",
            "--normal-id",
            "NORMAL",
            "--skip-annotation",
        ])
        .status()
        .expect("Failed to run mafsmith");

    assert!(status.success(), "mafsmith vcf2maf exited with error");

    let got = parse_tsv(tmp.path());
    let expected = parse_tsv(&expected_tsv);

    assert_eq!(
        got.len(),
        expected.len(),
        "Row count mismatch: got {} rows, expected {}",
        got.len(),
        expected.len()
    );

    for (i, (got_row, exp_row)) in got.iter().zip(expected.iter()).enumerate() {
        for col in KEY_COLUMNS {
            let got_val = got_row.get(*col).map(|s| s.as_str()).unwrap_or("");
            let exp_val = exp_row.get(*col).map(|s| s.as_str()).unwrap_or("");
            assert_eq!(
                got_val, exp_val,
                "Row {i} ({hugo} {chrom}:{pos}): column '{col}' mismatch\n  got:      '{got_val}'\n  expected: '{exp_val}'",
                hugo = exp_row.get("Hugo_Symbol").map(|s| s.as_str()).unwrap_or("?"),
                chrom = exp_row.get("Chromosome").map(|s| s.as_str()).unwrap_or("?"),
                pos = exp_row.get("Start_Position").map(|s| s.as_str()).unwrap_or("?"),
            );
        }
    }
}

/// Regression: GT=0/1 het call with high VAF (≥0.7) must keep Tumor_Seq_Allele1=REF.
///
/// vcf2maf.pl respects an explicit het GT and does not override Allele1 to ALT even
/// when VAF is high (e.g. AD=2,18 → VAF=0.9). mafsmith must match this behaviour.
/// Without the fix, the VAF override would fire and produce Allele1='G' instead of 'A'.
#[test]
#[ignore = "requires compiled mafsmith binary; run `cargo build` first"]
fn het_high_vaf_allele1_stays_ref() {
    let bin = mafsmith_bin();
    assert!(bin.exists(), "mafsmith binary not found; run `cargo build`");

    let input_vcf = fixtures_dir().join("integration_annotated.vcf");
    let tmp = tempfile::NamedTempFile::new().unwrap();

    let status = Command::new(&bin)
        .args([
            "vcf2maf",
            "--input-vcf",
            input_vcf.to_str().unwrap(),
            "--output-maf",
            tmp.path().to_str().unwrap(),
            "--vcf-tumor-id",
            "test_tumor",
            "--tumor-id",
            "test_tumor",
            "--vcf-normal-id",
            "test_normal",
            "--normal-id",
            "test_normal",
            "--genome",
            "grch38",
            "--skip-annotation",
        ])
        .status()
        .expect("Failed to run mafsmith");
    assert!(status.success(), "mafsmith vcf2maf exited with error");

    let rows = parse_tsv(tmp.path());
    let row = rows
        .iter()
        .find(|r| {
            r.get("Start_Position")
                .map(|s| s == "46300000")
                .unwrap_or(false)
        })
        .expect("Regression fixture row 21:46300000 not found in MAF output");

    assert_eq!(
        row.get("Tumor_Seq_Allele1").map(|s| s.as_str()).unwrap_or(""),
        "A",
        "GT=0/1 het with VAF=0.9: Tumor_Seq_Allele1 must stay REF ('A'), not be overridden to ALT ('G')"
    );
    assert_eq!(
        row.get("Tumor_Seq_Allele2")
            .map(|s| s.as_str())
            .unwrap_or(""),
        "G",
        "Tumor_Seq_Allele2 must be ALT ('G')"
    );
}

/// Full end-to-end: run vcf2maf on a real VCF with fastVEP annotation and compare
/// key columns against the vcf2maf reference MAF.
///
/// Requires: mafsmith binary, fastVEP installed, GRCh38 chr21 reference data,
/// and the vcf2maf reference MAF at tests/fixtures/expected/test_b38_vcf2maf.maf.
#[test]
#[ignore = "requires fastVEP, reference data, and vcf2maf reference output"]
fn vcf2maf_vs_vcf2maf_reference() {
    let bin = mafsmith_bin();
    let input_vcf = fixtures_dir().join("test_b38.vcf");
    let reference_maf = expected_dir().join("test_b38_vcf2maf.maf");

    if !reference_maf.exists() {
        eprintln!(
            "Reference MAF not found at {}.\n\
             Run `scripts/generate_expected.sh` to create it using vcf2maf.",
            reference_maf.display()
        );
        return;
    }

    let tmp = tempfile::NamedTempFile::new().unwrap();

    let status = Command::new(&bin)
        .args([
            "vcf2maf",
            "--input-vcf",
            input_vcf.to_str().unwrap(),
            "--output-maf",
            tmp.path().to_str().unwrap(),
            "--genome",
            "grch38",
            "--tumor-id",
            "TUMOR",
            "--normal-id",
            "NORMAL",
        ])
        .status()
        .expect("Failed to run mafsmith");

    assert!(status.success(), "mafsmith vcf2maf exited with error");

    let got = parse_tsv(tmp.path());
    let expected = parse_tsv(&reference_maf);

    // Variant count should match (allow for minor differences in filtering)
    let diff = (got.len() as isize - expected.len() as isize).abs();
    assert!(
        diff <= 2,
        "Row count differs by {diff}: got {}, expected {}",
        got.len(),
        expected.len()
    );

    // Match rows by (Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2)
    let key_of = |row: &HashMap<String, String>| {
        (
            row.get("Chromosome").cloned().unwrap_or_default(),
            row.get("Start_Position").cloned().unwrap_or_default(),
            row.get("Reference_Allele").cloned().unwrap_or_default(),
            row.get("Tumor_Seq_Allele2").cloned().unwrap_or_default(),
        )
    };

    let expected_by_key: HashMap<_, _> = expected.iter().map(|r| (key_of(r), r)).collect();

    let mut mismatches = 0usize;
    for got_row in &got {
        let k = key_of(got_row);
        let Some(exp_row) = expected_by_key.get(&k) else {
            eprintln!("WARNING: variant {:?} not found in reference MAF", k);
            continue;
        };

        for col in KEY_COLUMNS {
            let got_val = got_row.get(*col).map(|s| s.as_str()).unwrap_or("");
            let exp_val = exp_row.get(*col).map(|s| s.as_str()).unwrap_or("");
            if got_val != exp_val {
                eprintln!(
                    "MISMATCH {:?} col '{col}': got '{got_val}' expected '{exp_val}'",
                    k
                );
                mismatches += 1;
            }
        }
    }

    assert_eq!(
        mismatches, 0,
        "{mismatches} column mismatches vs vcf2maf reference"
    );
}

// ── vcf2vcf regression tests ─────────────────────────────────────────────────

/// Regression: vcf2vcf must count only written variants, not total input records.
///
/// Before the fix, `count` incremented for every record read (including filtered ones).
/// vcf2vcf passes non-PASS variants through (matching vcf2vcf.pl); only ref-only
/// records (ALT=".") are dropped.
#[test]
#[ignore = "requires compiled mafsmith binary; run `cargo build` first"]
fn vcf2vcf_filters_and_count_correct() {
    let bin = mafsmith_bin();
    assert!(bin.exists(), "mafsmith binary not found; run `cargo build`");

    let input_vcf = fixtures_dir().join("test_vcf2vcf.vcf");
    let tmp = tempfile::NamedTempFile::with_suffix(".vcf").unwrap();

    let out = Command::new(&bin)
        .args([
            "vcf2vcf",
            "--input-vcf",
            input_vcf.to_str().unwrap(),
            "--output-vcf",
            tmp.path().to_str().unwrap(),
        ])
        .output()
        .expect("Failed to run mafsmith vcf2vcf");

    assert!(
        out.status.success(),
        "mafsmith vcf2vcf exited with error: {}",
        String::from_utf8_lossy(&out.stderr)
    );

    // Log must report 6 written variants (all 7 input records minus the one ref-only ALT=".")
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("wrote 6 variants"),
        "Expected 'wrote 6 variants' in log, got: {stderr}"
    );

    // Count non-header lines in output VCF
    let content = fs::read_to_string(tmp.path()).unwrap();
    let data_lines: Vec<&str> = content.lines().filter(|l| !l.starts_with('#')).collect();
    assert_eq!(
        data_lines.len(),
        6,
        "Output VCF must have exactly 6 data lines, got {}",
        data_lines.len()
    );

    // ref-only variant at chr1:400000 (ALT=".") must be excluded — the only skip criterion
    assert!(
        !content.contains("\t400000\t"),
        "Ref-only variant at pos 400000 must be filtered out"
    );

    // PASS variant at chr1:100000 must be present
    assert!(
        content.contains("\t100000\t"),
        "PASS variant at pos 100000 must be in output"
    );

    // non-PASS (LOWQ) variant at chr1:200000 must be present — vcf2vcf passes non-PASS through
    assert!(
        content.contains("\t200000\t"),
        "LOWQ variant at pos 200000 must be passed through (vcf2vcf does not filter non-PASS)"
    );

    // filter="." variant at chr1:500000 must be present
    assert!(
        content.contains("\t500000\t"),
        "filter='.' variant at pos 500000 must be in output"
    );

    // LOWQ;DP variant at chr1:600000 must be present — vcf2vcf passes non-PASS through
    assert!(
        content.contains("\t600000\t"),
        "LOWQ;DP variant at pos 600000 must be passed through (vcf2vcf does not filter non-PASS)"
    );
}

// ── maf2vcf regression tests ─────────────────────────────────────────────────

/// Regression: maf2vcf must emit per-sample GT values in VCF output.
///
/// Before the fix, each data line ended after the FORMAT field with no sample columns,
/// producing invalid VCF that lacked genotype information.
///
/// Also verifies indel allele conversion (MAF "-" → VCF padding base).
#[test]
#[ignore = "requires compiled mafsmith binary; run `cargo build` first"]
fn maf2vcf_has_sample_columns() {
    let bin = mafsmith_bin();
    assert!(bin.exists(), "mafsmith binary not found; run `cargo build`");

    let input_maf = fixtures_dir().join("test_maf2vcf.maf");
    let tmp = tempfile::NamedTempFile::with_suffix(".vcf").unwrap();

    let status = Command::new(&bin)
        .args([
            "maf2vcf",
            "--input-maf",
            input_maf.to_str().unwrap(),
            "--output-vcf",
            tmp.path().to_str().unwrap(),
            "--genome",
            "grch38",
        ])
        .status()
        .expect("Failed to run mafsmith maf2vcf");

    assert!(status.success(), "mafsmith maf2vcf exited with error");

    let content = fs::read_to_string(tmp.path()).unwrap();
    let header_line = content
        .lines()
        .find(|l| l.starts_with("#CHROM"))
        .expect("VCF must contain a #CHROM header line");

    // Header must declare both sample columns
    assert!(
        header_line.contains("TUMOR_SAMPLE_1"),
        "VCF header must list TUMOR_SAMPLE_1 sample column"
    );
    assert!(
        header_line.contains("NORMAL_SAMPLE_1"),
        "VCF header must list NORMAL_SAMPLE_1 sample column"
    );

    let data_lines: Vec<&str> = content.lines().filter(|l| !l.starts_with('#')).collect();
    assert_eq!(
        data_lines.len(),
        5,
        "Expected 5 variant records, got {}",
        data_lines.len()
    );

    for (i, line) in data_lines.iter().enumerate() {
        let fields: Vec<&str> = line.split('\t').collect();
        // 9 standard VCF fields + 2 sample columns = 11 total
        assert_eq!(
            fields.len(),
            11,
            "Data line {i} must have 11 tab-separated fields (9 VCF + 2 samples), got {}: {line}",
            fields.len()
        );
        // FORMAT column must start with "GT" (may be "GT:AD:DP" when MAF has depth columns)
        assert!(
            fields[8].starts_with("GT"),
            "FORMAT field must start with 'GT' in line {i}, got '{}'",
            fields[8]
        );
        // Sample GT values must be non-empty
        assert!(
            !fields[9].is_empty(),
            "Tumor GT must not be empty in line {i}"
        );
        assert!(
            !fields[10].is_empty(),
            "Normal GT must not be empty in line {i}"
        );
    }

    // hom-alt inference: TP53 row has TSA1=G which != REF (C), so tumor GT must be 1/1
    let tp53_line = data_lines
        .iter()
        .find(|l| l.contains("\t7674220\t"))
        .expect("TP53 row at chr17:7674220 must be present");
    let fields: Vec<&str> = tp53_line.split('\t').collect();
    // Tumor sample is first value of the GT field (split on ':')
    let tumor_gt = fields[9].split(':').next().unwrap_or("");
    assert_eq!(
        tumor_gt, "1/1",
        "TP53 with TSA1!=REF must have tumor GT=1/1"
    );

    // Deletion: BRCA1 AT>- must produce anchor-padded VCF DEL (anchor+AT > anchor)
    let brca1_line = data_lines
        .iter()
        .find(|l| l.contains("chr17") && l.contains("\t43082"))
        .expect("BRCA1 deletion row must be present");
    let fields: Vec<&str> = brca1_line.split('\t').collect();
    let brca1_ref = fields[3];
    let brca1_alt = fields[4];
    assert!(
        brca1_ref.ends_with("AT") && brca1_ref.len() == 3,
        "BRCA1 deletion REF must be anchor+'AT' (3 chars), got '{brca1_ref}'"
    );
    assert_eq!(
        brca1_alt.len(),
        1,
        "BRCA1 deletion ALT must be 1-char anchor base, got '{brca1_alt}'"
    );
    assert_eq!(
        &brca1_ref[..1],
        brca1_alt,
        "BRCA1 deletion REF and ALT must share the same anchor base"
    );

    // Insertion: PIK3CA ->GTC must produce anchor-padded VCF INS (anchor > anchor+GTC)
    let pik3ca_line = data_lines
        .iter()
        .find(|l| l.contains("chr3"))
        .expect("PIK3CA insertion row must be present");
    let fields: Vec<&str> = pik3ca_line.split('\t').collect();
    let pik3ca_ref = fields[3];
    let pik3ca_alt = fields[4];
    assert_eq!(
        pik3ca_ref.len(),
        1,
        "PIK3CA insertion REF must be 1-char anchor base, got '{pik3ca_ref}'"
    );
    assert!(
        pik3ca_alt.ends_with("GTC") && pik3ca_alt.len() == 4,
        "PIK3CA insertion ALT must be anchor+'GTC' (4 chars), got '{pik3ca_alt}'"
    );
    assert_eq!(
        pik3ca_ref,
        &pik3ca_alt[..1],
        "PIK3CA insertion REF and ALT must share the same anchor base"
    );
}

// ── maf2maf regression test ───────────────────────────────────────────────────

/// End-to-end maf2maf round-trip: MAF → VCF → fastVEP annotation → MAF.
///
/// Requires: fastVEP binary, GRCh38 reference FASTA and GFF3 (via `mafsmith fetch`).
/// Run with: cargo test -- --include-ignored maf2maf_roundtrip
#[test]
#[ignore = "requires fastVEP and GRCh38 reference data; run `mafsmith fetch --genome grch38` first"]
fn maf2maf_roundtrip() {
    let bin = mafsmith_bin();
    assert!(bin.exists(), "mafsmith binary not found; run `cargo build`");

    let input_maf = fixtures_dir().join("test_maf2vcf.maf");
    let tmp = tempfile::NamedTempFile::with_suffix(".maf").unwrap();

    let status = Command::new(&bin)
        .args([
            "maf2maf",
            "--input-maf",
            input_maf.to_str().unwrap(),
            "--output-maf",
            tmp.path().to_str().unwrap(),
            "--genome",
            "grch38",
        ])
        .status()
        .expect("Failed to run mafsmith maf2maf");

    assert!(status.success(), "mafsmith maf2maf exited with error");

    let rows = parse_tsv(tmp.path());
    // All 5 variants should survive the round-trip (some may be multi-allelic-split)
    assert!(
        rows.len() >= 5,
        "maf2maf output must have at least 5 rows, got {}",
        rows.len()
    );

    // Every row must have a non-empty Chromosome field
    for row in &rows {
        let chrom = row.get("Chromosome").map(|s| s.as_str()).unwrap_or("");
        assert!(
            !chrom.is_empty(),
            "Chromosome must not be empty after maf2maf round-trip"
        );
    }

    // EGFR SNP must survive
    let egfr = rows
        .iter()
        .find(|r| r.get("Hugo_Symbol").map(|s| s == "EGFR").unwrap_or(false));
    assert!(
        egfr.is_some(),
        "EGFR variant must be present in maf2maf output"
    );
}
