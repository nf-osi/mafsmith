use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use mafsmith::{
    annotation::{
        consequence::{so_to_variant_classification, variant_type},
        csq::{shorten_hgvsp, CsqFormat},
        transcript::select_transcript,
    },
    vcf::{
        normalization::{maf_positions, normalize},
        VcfReader,
    },
};
use std::{
    io::{BufReader, Cursor},
    path::PathBuf,
};

// ── VCF parsing ───────────────────────────────────────────────────────────────

fn vcf_parse_fixture(c: &mut Criterion) {
    let fixtures = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures");
    let annotated_vcf = fixtures.join("test_b38_annotated.vcf");

    let content = std::fs::read_to_string(&annotated_vcf)
        .expect("annotated VCF fixture not found; run cargo test first");

    let variant_count = content.lines().filter(|l| !l.starts_with('#')).count();

    let mut g = c.benchmark_group("vcf_parse");
    g.throughput(Throughput::Elements(variant_count as u64));

    g.bench_function("parse_annotated_vcf", |b| {
        b.iter(|| {
            let cursor = Cursor::new(content.as_bytes());
            let mut reader = VcfReader::new(BufReader::new(cursor));
            let mut count = 0usize;
            while let Ok(Some(rec)) = reader.next_record() {
                black_box(&rec);
                count += 1;
            }
            count
        })
    });

    g.finish();
}

// ── CSQ parsing ──────────────────────────────────────────────────────────────

const CSQ_HEADER: &str = r#"##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|TSL|MANE">"#;

const CSQ_VALUE: &str = "C|missense_variant|MODERATE|RUNX1|RUNX1|Transcript|ENST00000675419|protein_coding|9/9||c.1286T>G|p.Arg429Gly|1480/5971|1286/1032|429/343|R/G|cGg/cGc|||-1|||HGNC|10471|YES||MANE_Select";

fn csq_parse(c: &mut Criterion) {
    let fmt = CsqFormat::from_header_description(CSQ_HEADER).unwrap();

    let mut g = c.benchmark_group("csq");
    g.throughput(Throughput::Elements(1));

    g.bench_function("parse_single_entry", |b| {
        b.iter(|| {
            let entry = fmt.parse_entry(black_box(CSQ_VALUE), &[]);
            black_box(entry)
        })
    });

    // Simulate a variant with 10 transcript entries (typical in cancer VCFs)
    let multi_csq = std::iter::repeat(CSQ_VALUE)
        .take(10)
        .collect::<Vec<_>>()
        .join(",");
    g.throughput(Throughput::Elements(10));
    g.bench_function("parse_10_entries", |b| {
        b.iter(|| {
            let entries = fmt.parse_all(black_box(&multi_csq), &[]);
            black_box(entries)
        })
    });

    g.finish();
}

// ── Transcript selection ──────────────────────────────────────────────────────

fn transcript_select(c: &mut Criterion) {
    let fmt = CsqFormat::from_header_description(CSQ_HEADER).unwrap();
    let multi_csq = std::iter::repeat(CSQ_VALUE)
        .take(20)
        .collect::<Vec<_>>()
        .join(",");
    let entries = fmt.parse_all(&multi_csq, &[]);

    c.bench_function("transcript_select/20_entries", |b| {
        b.iter(|| {
            let t = select_transcript(black_box(&entries), None);
            black_box(t)
        })
    });
}

// ── Allele normalization ──────────────────────────────────────────────────────

fn normalization(c: &mut Criterion) {
    let cases: &[(&str, &str, &str, u64)] = &[
        ("SNP", "A", "T", 12345),
        ("INS", "A", "ACCTCTT", 34880555),
        ("DEL", "GT", "G", 41479182),
        ("DNP", "AC", "GT", 100000),
    ];

    let mut g = c.benchmark_group("normalization");
    for (name, ref_a, alt_a, pos) in cases {
        g.bench_with_input(BenchmarkId::from_parameter(name), name, |b, _| {
            b.iter(|| {
                let n = normalize(black_box(*pos), black_box(*ref_a), black_box(*alt_a));
                let (s, e) = maf_positions(&n);
                black_box((s, e))
            })
        });
    }
    g.finish();
}

// ── Consequence classification ────────────────────────────────────────────────

fn consequence_classification(c: &mut Criterion) {
    let cases: &[(&str, &[&str], &str, &str)] = &[
        ("missense", &["missense_variant"], "A", "T"),
        ("frameshift_del", &["frameshift_variant"], "AT", "A"),
        (
            "splice_site",
            &["splice_donor_variant", "intron_variant"],
            "A",
            "T",
        ),
        (
            "multi_csq",
            &[
                "missense_variant",
                "splice_region_variant",
                "synonymous_variant",
            ],
            "A",
            "T",
        ),
    ];

    let mut g = c.benchmark_group("consequence");
    for (name, csqs, ref_a, alt_a) in cases {
        g.bench_with_input(BenchmarkId::from_parameter(name), name, |b, _| {
            b.iter(|| {
                let cls = so_to_variant_classification(
                    black_box(csqs),
                    black_box(*ref_a),
                    black_box(*alt_a),
                );
                let vt = variant_type(black_box(*ref_a), black_box(*alt_a));
                black_box((cls, vt))
            })
        });
    }
    g.finish();
}

// ── HGVSp shortening ─────────────────────────────────────────────────────────

fn hgvsp_shorten(c: &mut Criterion) {
    let cases = &[
        ("short", "p.Arg429Gly"),
        ("long", "p.Glu1234_Lys1236delinsTrpAlaVal"),
        ("stop", "p.Gln123Ter"),
    ];

    let mut g = c.benchmark_group("hgvsp");
    for (name, hgvsp) in cases {
        g.bench_with_input(BenchmarkId::from_parameter(name), name, |b, _| {
            b.iter(|| {
                let s = shorten_hgvsp(black_box(hgvsp));
                black_box(s)
            })
        });
    }
    g.finish();
}

// ── End-to-end VCF→MAF conversion (fixture) ──────────────────────────────────

fn end_to_end_fixture(c: &mut Criterion) {
    let fixtures = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures");
    let annotated_vcf = fixtures.join("test_b38_annotated.vcf");
    let content = match std::fs::read_to_string(&annotated_vcf) {
        Ok(s) => s,
        Err(_) => return, // fixture not present; skip
    };
    let variant_count = content.lines().filter(|l| !l.starts_with('#')).count() as u64;

    let fmt = CsqFormat::from_header_description(CSQ_HEADER).unwrap();

    let mut g = c.benchmark_group("end_to_end");
    g.throughput(Throughput::Elements(variant_count));

    g.bench_function("parse_and_classify_fixture", |b| {
        b.iter(|| {
            let cursor = Cursor::new(content.as_bytes());
            let mut reader = VcfReader::new(BufReader::new(cursor));
            let mut results = Vec::new();
            while let Ok(Some(rec)) = reader.next_record() {
                if let Some(csq_val) = rec.csq_value() {
                    let entries = fmt.parse_all(csq_val, &[]);
                    if let Some(t) = select_transcript(&entries, None) {
                        let csqs: Vec<&str> = t.consequences.iter().map(|s| s.as_str()).collect();
                        let n = normalize(rec.pos, &rec.ref_allele, &rec.alt_allele);
                        let cls = so_to_variant_classification(&csqs, &n.ref_allele, &n.alt_allele);
                        results.push(cls);
                    }
                }
            }
            black_box(results)
        })
    });

    g.finish();
}

criterion_group!(
    benches,
    vcf_parse_fixture,
    csq_parse,
    transcript_select,
    normalization,
    consequence_classification,
    hgvsp_shorten,
    end_to_end_fixture,
);
criterion_main!(benches);
