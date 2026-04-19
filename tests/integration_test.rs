/// Integration tests for mafsmith VCF→MAF conversion.
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
    if release.exists() { release } else { debug }
}

/// Parse a TSV file (with optional leading `#version` comment) into rows.
fn parse_tsv(path: &Path) -> Vec<HashMap<String, String>> {
    let f = fs::File::open(path)
        .unwrap_or_else(|e| panic!("Cannot open {}: {e}", path.display()));
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
    let content = fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("Cannot read annotated VCF: {e}"));
    assert!(content.contains("CSQ="), "annotated VCF must contain CSQ annotations");
    let data_lines: Vec<&str> = content
        .lines()
        .filter(|l| !l.starts_with('#'))
        .collect();
    assert_eq!(data_lines.len(), 6, "Expected 6 variants in annotated VCF fixture");
}

#[test]
fn test_b38_vcf_fixture_has_25_variants() {
    let path = fixtures_dir().join("test_b38.vcf");
    let content = fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("Cannot read test VCF: {e}"));
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
    assert!(diff <= 2, "Row count differs by {diff}: got {}, expected {}", got.len(), expected.len());

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

    assert_eq!(mismatches, 0, "{mismatches} column mismatches vs vcf2maf reference");
}
