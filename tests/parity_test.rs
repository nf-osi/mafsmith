/// Parity suite: mafsmith vcf2maf vs vcf2maf.pl, cell for cell.
///
/// Each fixture in tests/fixtures/parity/manifest.tsv is converted by mafsmith
/// with `--skip-annotation` and compared against the MAF that vcf2maf.pl
/// `--inhibit-vep` wrote for the same VCF, sample IDs and `--retain-info` /
/// `--retain-fmt` / `--retain-ann` options (tests/fixtures/parity/expected/).
/// The fixtures cover the depth conventions of the callers mafsmith is validated
/// on — Mutect2, Strelka2 SNV/indel, VarScan2, SomaticSniper, Ion Torrent,
/// DeepVariant, Delly SVs — plus a CSQ-annotated VCF for the annotation columns.
///
/// Every column the two MAFs share is compared, apart from the documented
/// exceptions in `KNOWN_DIFFERENCES`. vcf2maf.pl writes ~60 more VEP columns that
/// mafsmith never emits; those are not compared, but every mafsmith column must
/// exist in vcf2maf.pl's header, in the same relative order.
///
/// mafsmith runs with `--strict`, which reproduces vcf2maf.pl for truncated AD
/// arrays; without it mafsmith deliberately reports the counts it can recover.
///
/// After editing a fixture, regenerate its inputs and the references:
///   python3 scripts/build_parity_fixtures.py
///   bash scripts/generate_parity_expected.sh
use std::{
    fs,
    path::{Path, PathBuf},
    process::Command,
};

/// Columns allowed to differ, as (column, reason). Keep this list empty if at all
/// possible: an entry here is a known gap in `--strict`, not a matter of taste.
const KNOWN_DIFFERENCES: &[(&str, &str)] = &[];

fn parity_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/parity")
}

struct Fixture {
    name: String,
    tumor: String,
    normal: Option<String>,
    retain_info: Option<String>,
    retain_fmt: Option<String>,
    retain_ann: Option<String>,
}

fn manifest() -> Vec<Fixture> {
    let text = fs::read_to_string(parity_dir().join("manifest.tsv")).unwrap();
    text.lines()
        .filter(|l| !l.trim().is_empty() && !l.starts_with('#'))
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            assert_eq!(
                f.len(),
                6,
                "manifest line needs 6 tab-separated fields: {l}"
            );
            let opt = |s: &str| (s != "-").then(|| s.to_owned());
            Fixture {
                name: f[0].to_owned(),
                tumor: f[1].to_owned(),
                normal: opt(f[2]),
                retain_info: opt(f[3]),
                retain_fmt: opt(f[4]),
                retain_ann: opt(f[5]),
            }
        })
        .collect()
}

struct Maf {
    header: Vec<String>,
    rows: Vec<Vec<String>>,
}

impl Maf {
    fn read(path: &Path) -> Maf {
        let text = fs::read_to_string(path)
            .unwrap_or_else(|e| panic!("cannot read {}: {e}", path.display()));
        let mut lines = text.lines().filter(|l| !l.starts_with('#'));
        let header: Vec<String> = lines
            .next()
            .unwrap()
            .split('\t')
            .map(String::from)
            .collect();
        let rows = lines
            .map(|l| l.split('\t').map(String::from).collect::<Vec<_>>())
            .collect();
        Maf { header, rows }
    }

    fn col(&self, name: &str) -> usize {
        self.header.iter().position(|h| h == name).unwrap()
    }
}

fn run_mafsmith(fx: &Fixture, out: &Path) {
    let dir = parity_dir();
    let mut args: Vec<String> = vec![
        "vcf2maf".into(),
        "--input-vcf".into(),
        dir.join(format!("{}.vcf", fx.name)).display().to_string(),
        "--output-maf".into(),
        out.display().to_string(),
        "--genome".into(),
        "grch38".into(),
        "--skip-annotation".into(),
        "--strict".into(),
        "--ref-fasta".into(),
        dir.join("reference.fa").display().to_string(),
        "--tumor-id".into(),
        fx.tumor.clone(),
        "--vcf-tumor-id".into(),
        fx.tumor.clone(),
    ];
    if let Some(n) = &fx.normal {
        args.extend([
            "--normal-id".into(),
            n.clone(),
            "--vcf-normal-id".into(),
            n.clone(),
        ]);
    }
    for (flag, val) in [
        ("--retain-info", &fx.retain_info),
        ("--retain-fmt", &fx.retain_fmt),
        ("--retain-ann", &fx.retain_ann),
    ] {
        if let Some(v) = val {
            args.extend([flag.into(), v.clone()]);
        }
    }
    let output = Command::new(env!("CARGO_BIN_EXE_mafsmith"))
        .args(&args)
        .output()
        .expect("failed to run mafsmith");
    assert!(
        output.status.success(),
        "mafsmith failed on {}:\n{}",
        fx.name,
        String::from_utf8_lossy(&output.stderr)
    );
}

/// Compare one fixture; returns a description of every problem found.
fn compare(fx: &Fixture, expected: &Maf, got: &Maf) -> Vec<String> {
    let mut problems = Vec::new();
    let tag = &fx.name;

    // Header: every mafsmith column is a vcf2maf.pl column, in the same order.
    let mut last = None;
    for h in &got.header {
        match expected.header.iter().position(|e| e == h) {
            None => problems.push(format!("{tag}: column {h} is not in vcf2maf.pl's header")),
            Some(i) if last.is_some_and(|l| i < l) => {
                problems.push(format!("{tag}: column {h} is out of vcf2maf.pl's order"))
            }
            Some(i) => last = Some(i),
        }
    }
    // The columns after vcf_pos (retained INFO, SV summary, retained FORMAT) must
    // be exactly vcf2maf.pl's.
    let after_vcf_pos = |m: &Maf| -> Vec<String> { m.header[m.col("vcf_pos") + 1..].to_vec() };
    if after_vcf_pos(expected) != after_vcf_pos(got) {
        problems.push(format!(
            "{tag}: retained columns differ\n  vcf2maf.pl: {:?}\n  mafsmith:   {:?}",
            after_vcf_pos(expected),
            after_vcf_pos(got)
        ));
    }

    if expected.rows.len() != got.rows.len() {
        problems.push(format!(
            "{tag}: {} rows from vcf2maf.pl, {} from mafsmith",
            expected.rows.len(),
            got.rows.len()
        ));
        return problems;
    }

    let shared: Vec<(&String, usize, usize)> = got
        .header
        .iter()
        .enumerate()
        .filter(|(_, h)| !KNOWN_DIFFERENCES.iter().any(|(c, _)| c == h))
        .filter_map(|(gi, h)| {
            expected
                .header
                .iter()
                .position(|e| e == h)
                .map(|ei| (h, ei, gi))
        })
        .collect();
    let (chrom, pos) = (expected.col("Chromosome"), expected.col("vcf_pos"));
    for (r, (e, g)) in expected.rows.iter().zip(&got.rows).enumerate() {
        for (name, ei, gi) in &shared {
            let ev = e.get(*ei).map_or("", String::as_str);
            let gv = g.get(*gi).map_or("", String::as_str);
            if ev != gv {
                problems.push(format!(
                    "{tag} row {r} ({}:{}) {name}: vcf2maf.pl={ev:?} mafsmith={gv:?}",
                    e[chrom], e[pos]
                ));
            }
        }
    }
    problems
}

#[test]
fn vcf2maf_matches_vcf2maf_pl() {
    let fixtures = manifest();
    assert!(!fixtures.is_empty());
    let tmp = tempfile::tempdir().unwrap();

    let mut problems = Vec::new();
    let mut rows = 0usize;
    for fx in &fixtures {
        let out = tmp.path().join(format!("{}.maf", fx.name));
        run_mafsmith(fx, &out);
        let expected = Maf::read(&parity_dir().join(format!("expected/{}.vcf2maf.maf", fx.name)));
        let got = Maf::read(&out);
        rows += got.rows.len();
        problems.extend(compare(fx, &expected, &got));
    }

    assert!(
        problems.is_empty(),
        "{} differences from vcf2maf.pl across {} fixtures:\n{}",
        problems.len(),
        fixtures.len(),
        problems.join("\n")
    );
    eprintln!(
        "{} fixtures, {rows} MAF rows match vcf2maf.pl",
        fixtures.len()
    );
}

/// The suite must actually exercise what it claims to: retained columns from
/// INFO and FORMAT, a column dropped for clashing with vcf2maf.pl's header,
/// populated all_effects, and the SV summary columns with a Fusion pair.
#[test]
fn parity_fixtures_cover_retain_options() {
    let fixtures = manifest();
    assert!(fixtures.iter().any(|f| f.retain_info.is_some()));
    assert!(fixtures.iter().any(|f| f.retain_fmt.is_some()));
    assert!(fixtures.iter().any(|f| f.retain_ann.is_some()));
    let mutect = Maf::read(&parity_dir().join("expected/mutect2_tn.vcf2maf.maf"));
    // AF was requested via --retain-info but is already a vcf2maf.pl column.
    assert_eq!(mutect.header.iter().filter(|h| *h == "AF").count(), 1);
    assert!(mutect.header.iter().any(|h| h == "t_AD"));

    let annotated = Maf::read(&parity_dir().join("expected/annotated.vcf2maf.maf"));
    let effects = annotated.col("all_effects");
    assert!(annotated
        .rows
        .iter()
        .any(|r| r[effects].matches(';').count() > 1));

    let sv = Maf::read(&parity_dir().join("expected/delly_sv.vcf2maf.maf"));
    let fusion = sv.col("Fusion");
    assert!(sv.rows.iter().all(|r| r[fusion].ends_with(" fusion")));
}
