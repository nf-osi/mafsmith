//! Ports of vcf2maf.pl's per-record VCF handling, for output that has to match it
//! exactly: the SV preprocessing pass that rewrites INFO, the `%info` parse, and
//! the tumor/normal FORMAT handling — variant-allele selection, `FixAlleleDepths`
//! and the genotype fix-ups.
//!
//! The FORMAT side is used for the `--retain-fmt` columns, which vcf2maf.pl reads
//! from its rewritten per-sample hashes rather than from the raw VCF, and as the
//! fallback for the depth columns when `depth::extract_depth` finds no allele
//! counts. Perl semantics are kept deliberately — `"0"`/`""` falsiness, string to
//! number coercion, `%.15g` number formatting and `split` dropping trailing empty
//! fields — because each of them shows up in the output.

use crate::vcf::VcfRecord;
use std::collections::{HashMap, HashSet};

/// vcf2maf.pl's default `--min-hom-vaf`: the VAF at which an untyped tumor
/// genotype is called hom-alt.
const MIN_HOM_VAF: f64 = 0.7;

// ── INFO ─────────────────────────────────────────────────────────────────────

/// vcf2maf.pl's SV preprocessing pass rewrites INFO (and ALT) before the main
/// conversion reads them, so retained INFO values must come from the rewritten line.
/// Returns the INFO line to parse and, for split SVs, the `<SVTYPE>` ALT.
pub(crate) fn split_sv_preprocess(
    info: &str,
    ref_allele: &str,
    alts: &[String],
) -> (String, Option<String>) {
    // First pass: plain key→value map, last occurrence wins, flags map to "1".
    let mut svtype: Option<&str> = None;
    for tok in info.split(';') {
        match tok.split_once('=') {
            Some(("SVTYPE", v)) => svtype = Some(v),
            None if tok == "SVTYPE" => svtype = Some("1"),
            _ => {}
        }
    }
    let Some(svtype) = svtype.filter(|v| perl_truthy(v)) else {
        return (info.to_owned(), None);
    };

    let alt = alts.join(",");
    let is_bases = |s: &str, extra: &str| {
        !s.is_empty()
            && s.chars()
                .all(|c| "ACGTNacgtn".contains(c) || extra.contains(c))
    };
    if is_bases(ref_allele, "") && is_bases(&alt, ",") {
        (strip_svtype(info), None)
    } else if matches!(svtype, "BND" | "TRA" | "DEL" | "DUP" | "INV") {
        let line = strip_svtype(info);
        let line = replace_first(&line, "CT=", |rest| {
            let b = rest.as_bytes();
            (b.len() >= 4
                && matches!(b[0], b'3' | b'5')
                && &rest[1..3] == "to"
                && matches!(b[3], b'3' | b'5'))
            .then(|| (format!("Frame={}", &rest[..4]), 4))
        });
        let line = replace_first(&line, "SVMETHOD=", |rest| {
            let n = rest
                .find(|c: char| !(c.is_ascii_alphanumeric() || c == '_' || c == '.'))
                .unwrap_or(rest.len());
            (n > 0).then(|| (format!("Method={}", &rest[..n]), n))
        });
        (line, Some(format!("<{svtype}>")))
    } else {
        (info.to_owned(), None)
    }
}

/// `s/(SVTYPE=\w+;|;SVTYPE=\w+|SVTYPE=\w+)//` — remove the leftmost match.
fn strip_svtype(info: &str) -> String {
    let word_len = |s: &str| {
        s.find(|c: char| !(c.is_ascii_alphanumeric() || c == '_'))
            .unwrap_or(s.len())
    };
    for (i, _) in info.char_indices() {
        let rest = &info[i..];
        if let Some(v) = rest.strip_prefix("SVTYPE=") {
            let n = word_len(v);
            if n > 0 {
                let end = 7 + n;
                // Prefer consuming the trailing ';' (first alternative).
                let end = if rest[end..].starts_with(';') {
                    end + 1
                } else {
                    end
                };
                return format!("{}{}", &info[..i], &info[i + end..]);
            }
        } else if let Some(v) = rest.strip_prefix(";SVTYPE=") {
            let n = word_len(v);
            if n > 0 {
                return format!("{}{}", &info[..i], &info[i + 8 + n..]);
            }
        }
    }
    info.to_owned()
}

/// Replace the leftmost `prefix` occurrence that `tail` accepts. `tail` gets the
/// text after `prefix` and returns the replacement plus how much of it to consume.
fn replace_first(s: &str, prefix: &str, tail: impl Fn(&str) -> Option<(String, usize)>) -> String {
    for (i, _) in s.match_indices(prefix) {
        if let Some((rep, n)) = tail(&s[i + prefix.len()..]) {
            return format!("{}{}{}", &s[..i], rep, &s[i + prefix.len() + n..]);
        }
    }
    s.to_owned()
}

/// Parse INFO into vcf2maf.pl's `%info`: flags become "1", and a repeated key
/// accumulates its distinct values comma-joined.
pub(crate) fn parse_info(info_line: &str) -> HashMap<&str, String> {
    let toks = perl_split(info_line, ';');
    let mut info: HashMap<&str, String> = toks
        .iter()
        .filter(|t| !t.contains('='))
        .map(|t| (*t, "1".to_owned()))
        .collect();
    for tok in toks.iter().filter(|t| t.contains('=')) {
        let (key, val) = tok.split_once('=').unwrap();
        let cur = info.entry(key).or_default();
        if perl_truthy(cur) && perl_split(cur, ',').contains(&val) {
            continue;
        }
        if perl_truthy(cur) {
            cur.push(',');
        }
        cur.push_str(val);
    }
    info
}

// ── FORMAT ───────────────────────────────────────────────────────────────────

pub(crate) type FmtInfo = HashMap<String, String>;

/// One sample's FORMAT data as vcf2maf.pl holds it once it reaches the MAF columns.
#[derive(Debug, Clone)]
pub struct PerlSample {
    /// `%tum_info` / `%nrm_info`: FORMAT key → value, rewritten by allele selection
    /// and `FixAlleleDepths`. Keys absent from the sample are absent here.
    pub fmt: FmtInfo,
    /// `@tum_depths` / `@nrm_depths`: the rewritten AD split on commas.
    depths: Vec<String>,
}

impl PerlSample {
    /// `t_depth` / `n_depth`.
    pub fn depth(&self) -> String {
        self.fmt
            .get("DP")
            .filter(|d| *d != ".")
            .cloned()
            .unwrap_or_default()
    }

    /// `t_ref_count` / `n_ref_count`.
    pub fn ref_count(&self) -> String {
        self.depths.first().cloned().unwrap_or_default()
    }

    /// `t_alt_count` / `n_alt_count`, for vcf2maf.pl's selected ALT.
    pub fn alt_count(&self, var_idx: usize) -> String {
        self.depths.get(var_idx).cloned().unwrap_or_default()
    }
}

/// Both samples of one VCF record, plus the ALT allele vcf2maf.pl selected.
#[derive(Debug, Clone)]
pub struct PerlSamples {
    pub tumor: Option<PerlSample>,
    pub normal: Option<PerlSample>,
    /// VCF allele index of the selected ALT (1 = first ALT).
    pub var_idx: usize,
}

impl PerlSamples {
    /// Run vcf2maf.pl's per-record sample handling on `rec`.
    ///
    /// `tumor_vals` / `normal_vals` are the colon-split sample columns, or `None`
    /// when that sample is not in the VCF.
    pub fn new(rec: &VcfRecord, tumor_vals: Option<&[&str]>, normal_vals: Option<&[&str]>) -> Self {
        let (_, sv_alt) = split_sv_preprocess(&rec.info, &rec.ref_allele, &rec.all_alts);
        let alleles: Vec<String> = std::iter::once(rec.ref_allele.clone())
            .chain(match sv_alt {
                Some(a) => vec![a],
                None => rec.all_alts.clone(),
            })
            .collect();
        sample_infos(&rec.format_keys, tumor_vals, normal_vals, alleles)
    }
}

/// vcf2maf.pl's tumor/normal handling: pick the variant allele from the tumor GT
/// (or its depths), run `FixAlleleDepths` on both samples, and patch up GT.
fn sample_infos(
    format_keys: &[String],
    tumor_vals: Option<&[&str]>,
    normal_vals: Option<&[&str]>,
    mut alleles: Vec<String>,
) -> PerlSamples {
    let parse = |vals: &[&str]| -> FmtInfo {
        let keys = perl_trim(format_keys.iter().map(String::as_str).collect());
        perl_trim(vals.to_vec())
            .into_iter()
            .enumerate()
            .map(|(i, v)| (keys.get(i).copied().unwrap_or("").to_owned(), v.to_owned()))
            .collect()
    };
    let gt_called = |f: &FmtInfo| matches!(f.get("GT"), Some(g) if g != "." && g != "./.");
    let split_ad = |f: &FmtInfo| -> Vec<String> {
        perl_split(&f["AD"], ',')
            .into_iter()
            .map(str::to_owned)
            .collect()
    };

    let mut nrm = normal_vals.map(parse);
    let mut tum = tumor_vals.map(parse);
    let mut var_idx = 1usize;
    let mut tum_depths = Vec::new();

    if let Some(t) = tum.as_mut() {
        if gt_called(t) {
            let tum_gt = perl_split_gt(&t["GT"]);
            // First non-REF tumor allele, preferring one absent from the normal GT.
            let mut pick = tum_gt.iter().find(|g| **g != "0").copied();
            if let Some(n) = nrm.as_ref().filter(|n| gt_called(n)) {
                let nrm_gt: HashSet<&str> = perl_split_gt(&n["GT"]).into_iter().collect();
                pick = tum_gt
                    .iter()
                    .find(|g| **g != "0" && !nrm_gt.contains(**g))
                    .copied();
            }
            let idx = pick
                .filter(|p| !p.is_empty() && p.bytes().all(|b| b.is_ascii_digit()))
                .and_then(|p| p.parse::<usize>().ok())
                .filter(|&i| i < alleles.len());
            match idx {
                Some(i) => var_idx = i,
                None => {
                    var_idx = 1;
                    t.insert("GT".into(), "./.".into());
                }
            }
        }

        fix_allele_depths(&mut alleles, var_idx, t);
        tum_depths = split_ad(t);

        if !gt_called(t) {
            // Untyped tumor: take the deepest ALT, then call het/hom-alt from its VAF.
            for i in 1..tum_depths.len() {
                if perl_truthy(&tum_depths[i])
                    && perl_num(&tum_depths[i])
                        > tum_depths.get(var_idx).map_or(0.0, |d| perl_num(d))
                {
                    var_idx = i;
                }
            }
            let mut gt = "./.";
            if let Some(dp) = t.get("DP").filter(|d| *d != ".") {
                let dp = perl_num(dp);
                if dp != 0.0 {
                    if let Some(d) = tum_depths.get(var_idx) {
                        gt = if perl_num(d) / dp < MIN_HOM_VAF {
                            "0/1"
                        } else {
                            "1/1"
                        };
                    }
                }
            }
            t.insert("GT".into(), gt.into());
        }
    }

    let mut nrm_depths = Vec::new();
    if let Some(n) = nrm.as_mut() {
        fix_allele_depths(&mut alleles, var_idx, n);
        nrm_depths = split_ad(n);
        if !matches!(n.get("GT"), Some(g) if g != ".") {
            n.insert("GT".into(), "./.".into());
        }
    }

    PerlSamples {
        tumor: tum.map(|fmt| PerlSample {
            fmt,
            depths: tum_depths,
        }),
        normal: nrm.map(|fmt| PerlSample {
            fmt,
            depths: nrm_depths,
        }),
        var_idx,
    }
}

/// Port of vcf2maf.pl's `FixAlleleDepths`: derive per-allele depths from whichever
/// caller convention the sample uses, reconcile DP with them, and store the
/// result back as AD (one entry per allele, "." where unknown).
fn fix_allele_depths(alleles: &mut [String], var: usize, f: &mut FmtInfo) {
    let n = alleles.len();
    let has = |f: &FmtInfo, k: &str| f.contains_key(k);
    let num = |f: &FmtInfo, k: &str| f.get(k).map_or(0.0, |v| perl_num(v));
    let first = |f: &FmtInfo, k: &str| -> Option<String> {
        f.get(k)
            .and_then(|v| perl_split(v, ',').first().map(|s| s.to_string()))
    };
    let blank = || vec![Some(String::new()); n];
    // Assigning past the end of a Perl array leaves undefined holes.
    let set = |d: &mut Vec<Option<String>>, i: usize, v: Option<String>| {
        if i >= d.len() {
            d.resize(i + 1, None);
        }
        d[i] = v;
    };
    let digits_or_blank = |v: &str| {
        if !v.is_empty() && v.bytes().all(|b| b.is_ascii_digit()) {
            v.to_owned()
        } else {
            String::new()
        }
    };

    let mut depths: Vec<Option<String>> = Vec::new();
    if let Some(ad) = f.get("AD").filter(|a| *a != ".") {
        depths = perl_split(ad, ',')
            .into_iter()
            .map(|v| Some(digits_or_blank(v)))
            .collect();
    }
    let no_ad = !has(f, "AD");

    if depths.len() == 1 && has(f, "RD") {
        // VarScan: AD is the ALT depth alone, REF depth is in RD.
        depths = blank();
        depths[0] = f.get("RD").cloned();
        set(&mut depths, var, f.get("AD").cloned());
    } else if no_ad && has(f, "BCOUNT") {
        // SomaticSniper: per-base A,C,G,T counts.
        let bcount: Vec<String> = perl_split(&f["BCOUNT"], ',')
            .iter()
            .map(|s| s.to_string())
            .collect();
        depths = alleles
            .iter()
            .map(|a| {
                let idx = match a.as_str() {
                    "A" => Some(0),
                    "C" => Some(1),
                    "G" => Some(2),
                    "T" => Some(3),
                    _ => None,
                };
                Some(idx.and_then(|i| bcount.get(i).cloned()).unwrap_or_default())
            })
            .collect();
    } else if no_ad && ["AU", "CU", "GU", "TU"].iter().all(|k| has(f, k)) {
        // Strelka SNVs: keep tier-1 counts only.
        for b in ["A", "C", "G", "T"] {
            let k = format!("{b}U");
            match first(f, &k) {
                Some(v) => {
                    f.insert(k, v);
                }
                None => {
                    f.remove(&k);
                }
            }
        }
        if n == 2 && alleles[1] == "N" {
            let depth_of = |b: &str| f.get(&format!("{b}U")).map_or(0.0, |v| perl_num(v));
            let mut deepest = ["A", "C", "G", "T"];
            deepest.sort_by(|a, b| depth_of(b).total_cmp(&depth_of(a)));
            alleles[1] = if deepest[0] != alleles[0] {
                deepest[0]
            } else {
                deepest[1]
            }
            .to_owned();
        }
        depths = alleles
            .iter()
            .map(|a| Some(f.get(&format!("{a}U")).cloned().unwrap_or_default()))
            .collect();
    } else if no_ad && has(f, "TIR") && has(f, "TAR") {
        // Strelka indels: tier-1 REF depth in TAR, ALT depth in TIR.
        depths = blank();
        depths[0] = first(f, "TAR");
        set(&mut depths, var, first(f, "TIR"));
    } else if no_ad
        && ["FAZ", "FCZ", "FGZ", "FTZ", "RAZ", "RCZ", "RGZ", "RTZ"]
            .iter()
            .all(|k| has(f, k))
    {
        // CaVEMan: forward + reverse strand counts per base.
        for b in ["A", "C", "G", "T"] {
            let sum = num(f, &format!("F{b}Z")) + num(f, &format!("R{b}Z"));
            f.insert(format!("{b}Z"), perl_fmt_num(sum));
        }
        depths = alleles
            .iter()
            .map(|a| Some(f.get(&format!("{a}Z")).cloned().unwrap_or_default()))
            .collect();
    } else if no_ad && has(f, "AO") && has(f, "RO") {
        // Ion Torrent / FreeBayes without AD: REF in RO, ALTs in AO.
        depths = std::iter::once(f.get("RO").cloned())
            .chain(
                perl_split(&f["AO"], ',')
                    .into_iter()
                    .map(|v| Some(digits_or_blank(v))),
            )
            .collect();
    } else if no_ad && has(f, "RR") && has(f, "RV") {
        // Delly: SV junction read counts.
        depths = blank();
        depths[0] = f.get("RR").cloned();
        set(&mut depths, var, f.get("RV").cloned());
    } else if no_ad && has(f, "DP4") && perl_split(&f["DP4"], ',').len() == 4 {
        let c: Vec<f64> = perl_split(&f["DP4"], ',')
            .iter()
            .map(|v| perl_num(v))
            .collect();
        depths = blank();
        depths[0] = Some(perl_fmt_num(c[0] + c[1]));
        set(&mut depths, var, Some(perl_fmt_num(c[2] + c[3])));
    } else if no_ad && ["PP", "NP", "PR", "NR"].iter().all(|k| has(f, k)) {
        // cgpPindel.
        depths = blank();
        set(
            &mut depths,
            var,
            Some(perl_fmt_num(num(f, "PP") + num(f, "NP"))),
        );
        let dp = perl_fmt_num(num(f, "PR") + num(f, "NR"));
        f.insert("DP".into(), dp);
    } else if no_ad && has(f, "FA") && has(f, "DP") && f["DP"] != "." {
        // ALT allele fraction × total depth.
        depths = blank();
        set(
            &mut depths,
            var,
            Some(format!("{:.0}", num(f, "FA") * num(f, "DP"))),
        );
    } else if no_ad && has(f, "DV") && has(f, "DP") {
        // bcftools/mpileup.
        depths = blank();
        set(&mut depths, var, f.get("DV").cloned());
    } else if has(f, "AD") && depths.len() == 1 {
        // A lone AD value is taken to be the ALT depth.
        let ad = f.get("AD").cloned();
        depths = blank();
        set(&mut depths, var, ad);
    } else if !depths.is_empty() && depths.len() != n {
        depths = blank();
    }

    let d0 = depths.first().cloned().flatten();
    let dv = depths.get(var).cloned().flatten();
    let sum = || {
        perl_fmt_num(
            depths
                .iter()
                .flatten()
                .filter(|d| perl_truthy(d) && *d != ".")
                .map(|d| perl_num(d))
                .sum(),
        )
    };

    // DP must cover the REF and ALT depths; otherwise rebuild it from all depths.
    if let Some(dp) = f.get("DP").filter(|d| *d != ".").map(|d| perl_num(d)) {
        let t0 = d0.as_deref().filter(|d| perl_truthy(d)).map(perl_num);
        let tv = dv.as_deref().filter(|d| perl_truthy(d)).map(perl_num);
        if t0.is_some_and(|x| x > dp)
            || tv.is_some_and(|x| x > dp)
            || matches!((t0, tv), (Some(x), Some(y)) if x + y > dp)
        {
            f.insert("DP".into(), sum());
        }
    }
    if d0.is_some() && dv.is_some() && f.get("DP").is_none_or(|d| d == ".") {
        f.insert("DP".into(), sum());
    }

    let ad = depths
        .iter()
        .map(|d| match d {
            Some(s) if !s.is_empty() => s.as_str(),
            _ => ".",
        })
        .collect::<Vec<_>>()
        .join(",");
    f.insert("AD".into(), ad);
}

// ── Perl semantics ───────────────────────────────────────────────────────────

/// Perl's `split(/sep/, $s)`: trailing empty fields are dropped.
fn perl_split(s: &str, sep: char) -> Vec<&str> {
    perl_trim(s.split(sep).collect())
}

/// `split(/[\/|]/, $gt)`.
fn perl_split_gt(s: &str) -> Vec<&str> {
    perl_trim(s.split(['/', '|']).collect())
}

fn perl_trim(mut v: Vec<&str>) -> Vec<&str> {
    while v.last() == Some(&"") {
        v.pop();
    }
    v
}

/// Perl truthiness of a defined string: false only for "" and "0".
fn perl_truthy(s: &str) -> bool {
    !s.is_empty() && s != "0"
}

/// Perl's string→number conversion: the longest leading decimal number, else 0.
fn perl_num(s: &str) -> f64 {
    let s = s.trim_start();
    let b = s.as_bytes();
    let mut i = 0;
    if i < b.len() && (b[i] == b'+' || b[i] == b'-') {
        i += 1;
    }
    let int_start = i;
    while i < b.len() && b[i].is_ascii_digit() {
        i += 1;
    }
    let mut mantissa_digits = i - int_start;
    if i < b.len() && b[i] == b'.' {
        let frac_start = i + 1;
        let mut j = frac_start;
        while j < b.len() && b[j].is_ascii_digit() {
            j += 1;
        }
        if mantissa_digits > 0 || j > frac_start {
            mantissa_digits += j - frac_start;
            i = j;
        }
    }
    if mantissa_digits == 0 {
        return 0.0;
    }
    if i < b.len() && (b[i] == b'e' || b[i] == b'E') {
        let mut j = i + 1;
        if j < b.len() && (b[j] == b'+' || b[j] == b'-') {
            j += 1;
        }
        let exp_start = j;
        while j < b.len() && b[j].is_ascii_digit() {
            j += 1;
        }
        if j > exp_start {
            i = j;
        }
    }
    s[..i].parse().unwrap_or(0.0)
}

/// Perl's number→string conversion (`%.15g`).
fn perl_fmt_num(x: f64) -> String {
    if x == 0.0 {
        return "0".into();
    }
    if !x.is_finite() {
        return if x.is_nan() {
            "NaN"
        } else if x > 0.0 {
            "Inf"
        } else {
            "-Inf"
        }
        .into();
    }
    let sci = format!("{x:.14e}");
    let (mant, exp) = sci.split_once('e').unwrap();
    let exp: i32 = exp.parse().unwrap();
    let trim = |s: &str| -> String {
        if s.contains('.') {
            s.trim_end_matches('0').trim_end_matches('.').to_owned()
        } else {
            s.to_owned()
        }
    };
    if (-4..15).contains(&exp) {
        trim(&format!("{:.*}", (14 - exp) as usize, x))
    } else {
        let sign = if exp < 0 { '-' } else { '+' };
        format!("{}e{}{:02}", trim(mant), sign, exp.abs())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn s(v: &[&str]) -> Vec<String> {
        v.iter().map(|x| x.to_string()).collect()
    }

    #[test]
    fn info_flags_and_repeats() {
        let info = parse_info("DP=10;SOMATIC;X=1;X=2;X=1;Y=");
        assert_eq!(info["DP"], "10");
        assert_eq!(info["SOMATIC"], "1");
        assert_eq!(info["X"], "1,2");
        assert_eq!(info["Y"], "");
    }

    #[test]
    fn sv_info_rewrite() {
        let (line, alt) = split_sv_preprocess(
            "IMPRECISE;SVTYPE=DEL;END=500;CT=3to5;SVMETHOD=EMBL.DELLYv0.7.8",
            "A",
            &s(&["<DEL>"]),
        );
        assert_eq!(line, "IMPRECISE;END=500;Frame=3to5;Method=EMBL.DELLYv0.7.8");
        assert_eq!(alt.as_deref(), Some("<DEL>"));
        assert_eq!(strip_svtype("SVTYPE=DEL;END=5"), "END=5");
    }

    #[test]
    fn strelka_snv_depths_and_gt() {
        let keys = s(&["DP", "FDP", "AU", "CU", "GU", "TU"]);
        let p = sample_infos(
            &keys,
            Some(&["40", "0", "0,0", "10,12", "0,0", "30,31"]),
            Some(&["35", "0", "0,0", "35,36", "0,0", "0,0"]),
            s(&["C", "T"]),
        );
        let t = p.tumor.unwrap();
        assert_eq!(t.fmt["AD"], "10,30");
        assert_eq!(t.fmt["CU"], "10");
        assert_eq!(t.fmt["GT"], "1/1"); // 30/40 >= 0.7
        assert_eq!(
            (t.depth(), t.ref_count(), t.alt_count(p.var_idx)),
            ("40".into(), "10".into(), "30".into())
        );
        let n = p.normal.unwrap();
        assert_eq!(n.fmt["AD"], "35,0");
        assert_eq!(n.fmt["GT"], "./.");
    }

    #[test]
    fn dp_rebuilt_when_too_small() {
        let keys = s(&["GT", "AD", "DP"]);
        let p = sample_infos(&keys, Some(&["0/1", "8,9", "10"]), None, s(&["A", "G"]));
        assert_eq!(p.tumor.unwrap().fmt["DP"], "17");
    }

    #[test]
    fn no_depth_source_leaves_columns_blank() {
        let keys = s(&["GT", "PR", "SR"]);
        let p = sample_infos(
            &keys,
            Some(&["0/1", "10,8", "5,4"]),
            None,
            s(&["G", "<DEL>"]),
        );
        let t = p.tumor.unwrap();
        assert_eq!(
            (t.depth(), t.ref_count(), t.alt_count(p.var_idx)),
            (String::new(), String::new(), String::new())
        );
    }

    #[test]
    fn perl_numbers() {
        assert_eq!(perl_num("12abc"), 12.0);
        assert_eq!(perl_num("."), 0.0);
        assert_eq!(perl_num("1.5e2x"), 150.0);
        assert_eq!(perl_fmt_num(17.0), "17");
        assert_eq!(perl_fmt_num(0.1 + 0.2), "0.3");
        assert_eq!(perl_fmt_num(1e20), "1e+20");
    }
}
