use crate::{
    annotation::{
        consequence::{so_to_variant_classification, variant_type},
        csq::{shorten_hgvsp, splice_hgvsp_short, CsqEntry, CsqFormat, CsqLightEntry},
        depth::extract_depth,
        transcript::select_transcript_light,
    },
    cli::{Annotator, Vcf2mafArgs},
    maf::{record::MafRecord, writer::MafWriter},
    vcf::{normalization::{maf_positions, normalize}, VcfReader},
};
use anyhow::{bail, Context, Result};
use flate2::read::MultiGzDecoder;
use rayon::prelude::*;
use std::{
    collections::HashSet,
    fmt::Write as FmtWrite,
    fs,
    io::{BufRead, BufReader, BufWriter, Read, Seek, SeekFrom},
    path::{Path, PathBuf},
    process::Command,
};
use tempfile::NamedTempFile;
use tracing::{info, warn};

pub async fn run(args: Vcf2mafArgs) -> Result<()> {
    let data_dir = crate::commands::fetch::data_dir(None)?;
    let genome_str = args.genome.ncbi_build();

    let custom_enst = load_custom_enst(args.custom_enst.as_deref())?;

    // Either run fastVEP → temp file, or read the already-annotated input VCF directly.
    let mut sv_ref_fasta: Option<PathBuf> = args.ref_fasta.clone();
    let (annotated_path, _tmp): (PathBuf, Option<NamedTempFile>) = if args.skip_annotation {
        // Auto-resolve ref-fasta even in skip-annotation mode so we can validate the genome build.
        // If it's not available (user hasn't run `mafsmith fetch`), warn and skip validation.
        let resolved = args.ref_fasta.clone().or_else(|| {
            let p = data_dir.join(genome_str).join("reference.fa");
            if p.exists() { Some(p) } else { None }
        });
        match resolved {
            Some(ref rf) => {
                sv_ref_fasta = Some(rf.clone());
                validate_ref_fasta(rf, &args.input_vcf)?;
            }
            None => warn!(
                "No reference FASTA available — genome build not validated. \
                 Run `mafsmith fetch --genome {}` or pass --ref-fasta.",
                args.genome.ncbi_build().to_lowercase()
            ),
        }
        (args.input_vcf.clone(), None)
    } else {
        match args.annotator {
            Annotator::Fastvep => {
                let ref_fasta = resolve_ref_fasta(&args.ref_fasta, &data_dir, genome_str)?;
                sv_ref_fasta = Some(ref_fasta.clone());
                validate_ref_fasta(&ref_fasta, &args.input_vcf)?;
                let gff3 = resolve_gff3(&args.gff3, &data_dir, genome_str)?;
                let fastvep = resolve_fastvep(&args.fastvep_path, &data_dir)?;
                info!("Annotating variants with fastVEP...");
                let tmp = run_fastvep(&fastvep, &args.input_vcf, &ref_fasta, &gff3)?;
                let path = tmp.path().to_owned();
                (path, Some(tmp))
            }
            Annotator::Vep => {
                // Ref FASTA is optional in VEP mode (VEP manages its own cache).
                // Resolve it for SV base lookup and build validation if available.
                let resolved_fasta = args.ref_fasta.clone().or_else(|| {
                    let p = data_dir.join(genome_str).join("reference.fa");
                    if p.exists() { Some(p) } else { None }
                });
                if let Some(ref rf) = resolved_fasta {
                    sv_ref_fasta = Some(rf.clone());
                    validate_ref_fasta(rf, &args.input_vcf)?;
                } else {
                    warn!(
                        "No reference FASTA available — genome build not validated. \
                         Pass --ref-fasta if needed for SV secondary row base lookup."
                    );
                }
                let vep = resolve_vep(&args.vep_path)?;
                let vep_data = args.vep_data.clone().unwrap_or_else(|| {
                    PathBuf::from(std::env::var("HOME").unwrap_or_default()).join(".vep")
                });
                info!("Annotating variants with Ensembl VEP...");
                let tmp = run_vep(&vep, &args.input_vcf, &vep_data, genome_str, args.vep_forks)?;
                let path = tmp.path().to_owned();
                (path, Some(tmp))
            }
        }
    };

    // Parse and convert
    info!("Converting annotated VCF to MAF...");
    let out_file = fs::File::create(&args.output_maf)
        .with_context(|| format!("Cannot create output MAF: {}", args.output_maf.display()))?;
    let mut writer = MafWriter::new(BufWriter::new(out_file), args.retain_ann.clone())?;

    let in_file = fs::File::open(&annotated_path)
        .with_context(|| format!("Cannot open annotated VCF: {}", annotated_path.display()))?;
    let is_gz = annotated_path.extension().map_or(false, |e| e == "gz");
    let mut vcf = if is_gz {
        let dec = MultiGzDecoder::new(in_file);
        VcfReader::new(BufReader::new(Box::new(dec) as Box<dyn std::io::Read>))
    } else {
        VcfReader::new(BufReader::new(Box::new(in_file) as Box<dyn std::io::Read>))
    };

    // Pull CSQ FORMAT from header before iterating records
    let mut csq_format: Option<CsqFormat> = None;

    // We need to read header first by initializing the reader
    // The reader inits lazily on first next_record call, so prime it:
    let first_rec = vcf.next_record()?;

    // Parse CSQ format from header
    for hline in &vcf.header_lines {
        if hline.contains("ID=CSQ") {
            csq_format = CsqFormat::from_header_description(hline).ok();
            break;
        }
    }
    // When --skip-annotation is set and no CSQ header is present, proceed without
    // annotation: every record will fall through to the Targeted_Region default,
    // matching vcf2maf.pl --inhibit-vep behavior on unannotated VCFs.
    let csq_format = if let Some(fmt) = csq_format {
        fmt
    } else if args.skip_annotation {
        CsqFormat::empty()
    } else {
        anyhow::bail!(
            "fastVEP output did not contain a CSQ INFO header. \
             Was fastVEP run with annotation enabled?"
        )
    };

    let sample_names = vcf.sample_names.clone();
    let tumor_col = resolve_sample_col(
        args.vcf_tumor_id.as_deref(),
        &sample_names,
        "tumor",
        true,
    )?;
    let normal_col = resolve_sample_col(
        args.vcf_normal_id.as_deref(),
        &sample_names,
        "normal",
        false,
    )?;

    let tumor_barcode = args.tumor_id.as_deref().unwrap_or(tumor_col.as_deref().unwrap_or("TUMOR"));
    let normal_barcode = args.normal_id.as_deref().unwrap_or(normal_col.as_deref().unwrap_or("NORMAL"));

    // Pre-compute sample column indices once so each record avoids a linear name scan.
    let tumor_idx = tumor_col.as_deref()
        .and_then(|col| sample_names.iter().position(|s| s == col));
    let normal_idx = normal_col.as_deref()
        .and_then(|col| sample_names.iter().position(|s| s == col));

    // Convert one VCF record to MAF row(s). Returns Vec with optional secondary SV row first.
    let convert_record = |rec: &crate::vcf::VcfRecord| -> Result<Vec<MafRecord>> {
        // Lazy CSQ parsing: parse lightweight entries for all transcripts (no String allocations),
        // then fully parse only the selected transcript.
        let light_entries: Vec<(CsqLightEntry<'_>, &str)> = rec.csq_value()
            .map(|v| csq_format.parse_all_light(v))
            .unwrap_or_default();

        // Unannotated variants:
        // - If the VCF has no CSQ header at all (csq_format is empty), we're in --skip-annotation
        //   mode on an unannotated VCF. Emit every variant as Targeted_Region, matching
        //   vcf2maf.pl --inhibit-vep behavior.
        // - If the VCF was annotated (CSQ header present) but this record has no CSQ entry,
        //   only emit rows for recognized SV types; drop everything else (matching vcf2maf.pl
        //   behavior when VEP was run but left a variant unannotated).
        let default_entry: Option<CsqEntry> = if light_entries.is_empty() {
            if csq_format.fields.is_empty() {
                Some(CsqEntry::default())
            } else {
                match rec.info_field("SVTYPE") {
                    Some("BND") | Some("TRA") | Some("DEL") | Some("DUP") | Some("INV") => {
                        Some(CsqEntry::default())
                    }
                    _ => None,
                }
            }
        } else {
            None
        };

        // Select the best transcript from lightweight entries (cheap sort/scan, no allocation).
        let transcript: CsqEntry = if !light_entries.is_empty() {
            match select_transcript_light(&light_entries, custom_enst.as_ref()) {
                Some(idx) => csq_format.parse_entry(light_entries[idx].1, &args.retain_ann),
                None => match default_entry {
                    Some(e) => e,
                    None => return Ok(vec![]),
                },
            }
        } else {
            match default_entry {
                Some(e) => e,
                None => return Ok(vec![]),
            }
        };

        // Pre-split sample fields once per record to avoid redundant string splits.
        let tumor_vals: Option<Vec<&str>> = tumor_idx
            .and_then(|i| rec.samples_raw.get(i))
            .map(|s| s.split(':').collect());
        let normal_vals: Option<Vec<&str>> = normal_idx
            .and_then(|i| rec.samples_raw.get(i))
            .map(|s| s.split(':').collect());

        // Select the effective tumor ALT allele.
        // - GT has a non-ref allele not seen in normal GT: use it (somatic allele, matches vcf2maf.pl).
        // - GT has only non-ref alleles all seen in normal GT: fall back to first non-ref in GT.
        // - GT is all-ref or missing (hom-ref call in GVCF / no-call): pick the alt with the
        //   highest read depth in AD. vcf2maf.pl uses the same highest-depth fallback, which
        //   is why multi-allelic GVCF sites can differ from simply defaulting to ALT1.
        let nrm_gt_set_for_alt_sel: std::collections::HashSet<usize> = normal_vals.as_deref()
            .and_then(|nrm| {
                let gi = rec.format_keys.iter().position(|k| k == "GT")?;
                let nrm_gt = nrm.get(gi)?;
                if *nrm_gt == "./." || *nrm_gt == "." { return None; }
                Some(nrm_gt.split(|c| c == '/' || c == '|')
                    .filter_map(|s: &str| s.parse().ok())
                    .collect())
            })
            .unwrap_or_default();
        let (effective_alt, effective_alt_vcf_idx) = tumor_vals.as_deref()
            .and_then(|vals| {
                let gt_idx = rec.format_keys.iter().position(|k| k == "GT")?;
                let gt = vals.get(gt_idx)?;
                let idxs: Vec<usize> = gt
                    .split(|c| c == '/' || c == '|')
                    .filter_map(|s| s.parse().ok())
                    .collect();

                // GT has a non-ref allele — pick the somatic one (not seen in normal GT) if possible.
                // vcf2maf.pl: grep {ne "0" and !nrm_gt{$_}} @tum_gt → first allele absent from normal.
                // When normal GT is valid: if no somatic allele found (all shared), vcf2maf.pl resets
                // var_allele_idx=1 and uses depth-based selection — fall through to AD path.
                // When normal GT is no-call/missing: skip the normal-GT filter, take first non-ref.
                let non_ref: Vec<usize> = idxs.iter().filter(|&&i| i > 0).cloned().collect();
                if !non_ref.is_empty() {
                    let somatic_idx = if !nrm_gt_set_for_alt_sel.is_empty() {
                        non_ref.iter().find(|idx| !nrm_gt_set_for_alt_sel.contains(*idx)).copied()
                    } else {
                        // Normal GT is no-call or there is no normal sample: take first non-ref from
                        // tumor GT without filtering (vcf2maf.pl line 640 fallback).
                        non_ref.first().copied()
                    };
                    if let Some(alt_vcf_idx) = somatic_idx {
                        let alt_allele = rec.all_alts.get(alt_vcf_idx - 1)?.clone();
                        return Some((alt_allele, alt_vcf_idx));
                    }
                    // All non-ref GT alleles are in a valid normal GT — no somatic allele found.
                    // vcf2maf.pl: reset var_allele_idx=1, then depth-loop upgrades only when AD
                    // has the full complement of values. Truncated AD → stay at idx=1.
                    let n_alleles = 1 + rec.all_alts.len();
                    let best_idx = rec.format_keys.iter().position(|k| k == "AD")
                        .and_then(|ai| vals.get(ai))
                        .and_then(|ad_str| {
                            let ad_vals: Vec<u32> = ad_str
                                .split(',').filter_map(|v| v.parse().ok()).collect();
                            if ad_vals.len() < n_alleles { return None; } // truncated → keep idx=1
                            // Full AD: find first-maximum depth alt (strict >, same as vcf2maf.pl loop).
                            ad_vals.iter().enumerate().skip(1)
                                .fold(None::<(usize, u32)>, |best, (i, &c)| Some(match best {
                                    None => (i, c),
                                    Some((_, bc)) if c > bc => (i, c),
                                    Some(b) => b,
                                }))
                                .map(|(i, _)| i)
                        })
                        .unwrap_or(1); // default var_allele_idx=1 (vcf2maf.pl fallback)
                    let alt_allele = rec.all_alts.get(best_idx - 1)?.clone();
                    return Some((alt_allele, best_idx));
                }

                // GT is all-ref or missing: pick highest-depth alt from AD.
                let ad_idx = rec.format_keys.iter().position(|k| k == "AD")?;
                let ad_vals: Vec<u32> = vals.get(ad_idx)?
                    .split(',').filter_map(|v| v.parse().ok()).collect();
                // Use first-maximum: when depths tie, prefer the lower allele index (vcf2maf.pl behavior).
                // Rust's max_by_key returns the LAST maximum; fold with strict > keeps the FIRST.
                let best_idx = ad_vals.iter().enumerate().skip(1)
                    .fold(None::<(usize, u32)>, |best, (i, &c)| Some(match best {
                        None => (i, c),
                        Some((_, bc)) if c > bc => (i, c),
                        Some(b) => b,
                    }))
                    .map(|(i, _)| i)?;
                let alt_allele = rec.all_alts.get(best_idx - 1)?.clone();
                Some((alt_allele, best_idx))
            })
            .unwrap_or_else(|| (rec.alt_allele.clone(), 1));

        // Normalize SV ALT alleles to <SVTYPE> form.
        // vcf2maf.pl does `$cols[4] = "<" . $info{SVTYPE} . ">"` for BND/TRA/DEL/DUP/INV,
        // which normalizes <DUP:TANDEM>→<DUP>, BND breakend notation→<BND>, etc.
        //
        // Fast path: SVTYPE only appears in SV records. For annotated non-SV records, INFO starts
        // with "CSQ=" and we can skip the full info_field scan (O(|info|)) in O(7) characters.
        let svtype_tag: Option<&str> = if rec.info.starts_with("SVTYPE=") || rec.info.contains(";SVTYPE=") {
            rec.info_field("SVTYPE")
        } else {
            None
        };
        let is_sv_to_split = matches!(svtype_tag, Some("BND") | Some("TRA") | Some("DEL") | Some("DUP") | Some("INV"));

        // Drop records with unrecognized symbolic ALT alleles (<INS>, <CNV>, etc.) —
        // vcf2maf.pl only processes BND/TRA/DEL/DUP/INV as SV types and silently drops others.
        if effective_alt.starts_with('<') && !is_sv_to_split {
            return Ok(vec![]);
        }

        // Stash raw ALT only for SV records — Manta BND records may omit CHR2/END from INFO,
        // encoding the partner locus in the ALT notation instead.
        let raw_sv_alt = is_sv_to_split.then(|| effective_alt.clone());
        let effective_alt = if is_sv_to_split
            && (effective_alt.starts_with('<') || effective_alt.contains('[') || effective_alt.contains(']'))
        {
            format!("<{}>", svtype_tag.unwrap())
        } else if !effective_alt.starts_with('<')
            && (effective_alt.contains('[') || effective_alt.contains(']'))
        {
            // BND notation without a recognized SVTYPE — fall back to <BND>
            format!("<{}>", svtype_tag.unwrap_or("BND"))
        } else {
            effective_alt
        };

        // Normalize alleles: strip shared prefix/suffix, compute MAF positions.
        let norm = normalize(rec.pos, &rec.ref_allele, &effective_alt);
        let (start_pos, end_pos) = maf_positions(&norm);

        // Detect whether we need to emit a secondary SV breakpoint row.
        // vcf2maf.pl writes a second row at (CHR2, END) for all BND/TRA/DEL/DUP/INV SVs,
        // using samtools faidx to fetch the reference base at the partner location.
        // We replicate this by reading the FASTA index (.fai) directly.
        let sv_secondary: Option<(String, u64, String)> = if is_sv_to_split && effective_alt.starts_with('<') {
            // CHR2 defaults to the same chromosome when absent (vcf2maf.pl behavior for DEL/DUP/INV).
            let chr2 = rec.info_field("CHR2")
                .filter(|s| !s.is_empty())
                .map(|s| s.to_owned())
                .unwrap_or_else(|| rec.chrom.clone());
            let chr2_end = rec.info_field("END")
                .and_then(|s| s.parse::<u64>().ok())
                .map(|end_val| (chr2, end_val))
                .or_else(|| raw_sv_alt.as_deref().and_then(parse_bnd_alt));
            chr2_end.map(|(chr2, end_val)| {
                let ref2 = sv_ref_fasta.as_deref()
                    .and_then(|fa| fasta_fetch_base(fa, &chr2, end_val))
                    .map(|c| c.to_string())
                    .unwrap_or_else(|| rec.ref_allele.clone());
                (chr2, end_val, ref2)
            })
        } else {
            None
        };

        // Allele depths — use effective alt index for correct AD field in multi-allelic records.
        let tumor_depth = tumor_vals.as_deref()
            .map(|vals| extract_depth(rec.format_keys.as_slice(),vals, &rec.ref_allele, &effective_alt, effective_alt_vcf_idx));
        let normal_depth = normal_vals.as_deref()
            .map(|vals| extract_depth(rec.format_keys.as_slice(),vals, &rec.ref_allele, &rec.alt_allele, 1));

        let var_class =
            so_to_variant_classification(&transcript.consequences, &norm.ref_allele, &norm.alt_allele);
        let var_type = variant_type(&norm.ref_allele, &norm.alt_allele);

        // dbSNP: populate from VEP Existing_variation only (matching vcf2maf.pl behavior).
        // vcf2maf.pl does not use the VCF ID column for this field.
        let dbsnp = transcript.existing_variation
            .split(',')
            .filter(|s| s.starts_with("rs") && s[2..].chars().all(|c| c.is_ascii_digit()))
            .collect::<Vec<_>>()
            .join(",");

        // all_effects: semicolon-separated list for every annotated transcript.
        // Use light entries — avoids a full CsqEntry parse for each transcript.
        // hgvsp prefix is stripped inline (rfind(':')) without allocating.
        let mut all_effects = String::new();
        for (light, _) in light_entries.iter().filter(|(l, _)| !l.feature.is_empty()) {
            if !all_effects.is_empty() { all_effects.push(';'); }
            let cls = so_to_variant_classification(&light.consequences, &rec.ref_allele, &rec.alt_allele);
            write!(all_effects, "{},", light.symbol).unwrap();
            for (i, csq) in light.consequences.iter().enumerate() {
                if i > 0 { all_effects.push('&'); }
                all_effects.push_str(csq);
            }
            let hgvsp = light.hgvsp_raw;
            let hgvsp_stripped = hgvsp.rfind(':').map(|p| &hgvsp[p+1..]).unwrap_or(hgvsp);
            write!(all_effects, ",{},{}", cls, hgvsp_stripped).unwrap();
        }

        let hgvsp_short = if transcript.hgvsp.is_empty() && var_class == "Splice_Site" {
            splice_hgvsp_short(&transcript.hgvsc)
        } else {
            shorten_hgvsp(&transcript.hgvsp)
        };

        let fmt_depth = |d: Option<u32>, truncated: bool| -> String {
            if truncated { ".".to_owned() } else { d.map(|v| v.to_string()).unwrap_or_default() }
        };
        let fmt_normal_depth = |d: Option<u32>, truncated: bool| -> String {
            if normal_col.is_none() { String::new() }
            else if truncated { ".".to_owned() }
            else { d.map(|v| v.to_string()).unwrap_or_default() }
        };

        // n_expected_ad: AD needs exactly REF + all ALTs entries for a complete read.
        let n_expected_ad = 1 + rec.all_alts.len();
        // strict mode: when AD has fewer values than expected, suppress depth-based allele
        // calling and report '.' for depth fields — matching vcf2maf.pl behavior.
        let tumor_ad_truncated = args.strict && tumor_vals.as_deref()
            .and_then(|vals| rec.format_keys.iter().position(|k| k == "AD")
                .and_then(|ai| vals.get(ai)))
            .map(|ad_str| ad_str.split(',').count() < n_expected_ad)
            .unwrap_or(false);

        // Tumor_Seq_Allele1: sort GT allele indices and use the smallest, matching vcf2maf.
        // 0/1 or 1|0 → min=0 → ref; 1/1 → min=1 → alt.
        // No GT field (Strelka2 somatic callers): use depth-based VAF to infer het vs hom-alt.
        // GT=1/2 or 0/2 with multi-allelic: the "other" GT allele drives Tumor_Seq_Allele1.
        let tumor_seq_allele1 = tumor_vals.as_deref()
            .and_then(|vals| {
                // When no GT field is present (e.g. Strelka2 somatic), use depth-based VAF to
                // determine het vs hom-alt, matching vcf2maf.pl behavior.
                // VAF >= 0.7 (min_hom_vaf default) → hom-alt → TSA1=alt; else het → TSA1=ref.
                let gt_idx = rec.format_keys.iter().position(|k| k == "GT");
                if gt_idx.is_none() {
                    let ad = extract_depth(rec.format_keys.as_slice(),vals, &rec.ref_allele, &effective_alt, effective_alt_vcf_idx);
                    let hom_alt = match (ad.alt_count, ad.depth()) {
                        (Some(alt), Some(total)) if total > 0 => alt as f64 / total as f64 >= 0.7,
                        _ => false,
                    };
                    return Some(if hom_alt { norm.alt_allele.clone() } else { norm.ref_allele.clone() });
                }

                let gt_idx = gt_idx.unwrap();
                let gt = vals.get(gt_idx)?;
                let mut idxs: Vec<usize> = gt
                    .split(|c| c == '/' || c == '|')
                    .filter_map(|s| s.parse().ok())
                    .collect();

                // Depth-based fallback for missing or hom-ref GT (mirrors vcf2maf.pl):
                // vcf2maf.pl overrides GT=0/0 with depth-based logic — when all GT alleles are
                // reference (or GT is missing), use AD/DP to determine hom-alt vs het.
                // VAF >= 0.7 (min_hom_vaf default) → hom-alt (TSA1=alt); else het (TSA1=ref).
                // Also falls back when every non-ref tumor GT allele is also present in normal GT
                // (no somatic allele identified) — vcf2maf.pl sets GT=./. in this case.
                let all_in_normal_gt = {
                    let non_ref: Vec<usize> = idxs.iter().filter(|&&i| i > 0).cloned().collect();
                    !non_ref.is_empty() && normal_vals.as_deref().map(|nrm| {
                        let nrm_gt_set: std::collections::HashSet<usize> = rec.format_keys.iter()
                            .position(|k| k == "GT")
                            .and_then(|gi| nrm.get(gi))
                            .filter(|s| **s != "./." && **s != ".")
                            .map(|nrm_gt| nrm_gt.split(|c| c == '/' || c == '|')
                                .filter_map(|s: &str| s.parse().ok())
                                .collect::<std::collections::HashSet<usize>>())
                            .unwrap_or_default();
                        !nrm_gt_set.is_empty() && non_ref.iter().all(|idx| nrm_gt_set.contains(idx))
                    }).unwrap_or(false)
                };
                let needs_depth_fallback = idxs.is_empty() || idxs.iter().all(|&i| i == 0) || all_in_normal_gt;
                if needs_depth_fallback {
                    // When GT is a no-call (./.) or hom-ref (0/0), vcf2maf.pl infers
                    // het vs hom-alt from depth: VAF >= min_hom_vaf → TSA1=alt; else REF.
                    // This fires for all VCFs (single-sample and paired alike).
                    // Suppress only when AD is truncated in strict mode.
                    let hom_alt = if tumor_ad_truncated {
                        false
                    } else {
                        (|| -> Option<bool> {
                            let ad_idx = rec.format_keys.iter().position(|k| k == "AD")?;
                            let ad_vals: Vec<u32> = vals.get(ad_idx)?
                                .split(',').filter_map(|v| v.parse().ok()).collect();
                            let var_depth = *ad_vals.get(effective_alt_vcf_idx)?;
                            if var_depth == 0 { return Some(false); }
                            let dp = rec.format_keys.iter().position(|k| k == "DP")
                                .and_then(|i| vals.get(i))
                                .and_then(|v| v.parse::<u32>().ok())
                                .unwrap_or_else(|| ad_vals.iter().sum());
                            if dp == 0 { return None; }
                            Some(var_depth as f64 / dp as f64 >= 0.7)
                        })().unwrap_or(false)
                    };
                    return Some(if hom_alt { norm.alt_allele.clone() } else { norm.ref_allele.clone() });
                }
                idxs.sort();
                let min_idx = idxs[0];
                let _max_idx = idxs[idxs.len() - 1];

                // Find the "other" GT allele — the one that is NOT the effective alt.
                let other_idx = idxs.iter().find(|&&i| i != effective_alt_vcf_idx).copied();

                match other_idx {
                    Some(0) => Some(norm.ref_allele.clone()),
                    Some(i) if i != effective_alt_vcf_idx => {
                        // Other allele is a different alt — strip same prefix as normalization.
                        let other_raw = rec.all_alts.get(i - 1)?;
                        let prefix = rec.ref_allele.as_bytes()
                            .iter().zip(effective_alt.as_bytes().iter())
                            .take_while(|(r, a)| r == a).count();
                        let stripped = &other_raw[prefix.min(other_raw.len())..];
                        Some(if stripped.is_empty() { "-".to_owned() } else { stripped.to_owned() })
                    }
                    _ => {
                        // Both GT alleles are the same (hom-alt: min==max, both == effective_alt_vcf_idx)
                        if min_idx == 0 {
                            Some(norm.ref_allele.clone())
                        } else {
                            Some(norm.alt_allele.clone())
                        }
                    }
                }
            })
            .unwrap_or_else(|| norm.ref_allele.clone());

        // vcf2maf.pl VAF override: when GT is absent, no-call (./.), or hom-ref (0/0), and
        // alt VAF >= min_hom_vaf, the caller under-called; override Tumor_Seq_Allele1 to ALT.
        // vcf2maf.pl respects an explicit het GT (0/1) and does NOT override it even at
        // high VAF — mafsmith matches this behaviour by skipping the override when the tumor
        // GT contains both a ref allele (index 0) and at least one alt allele (index > 0).
        // Only fire when AD has the full complement of allele values (strict count, matching
        // vcf2maf.pl: it skips the override when AD is shorter than expected — e.g. GATK
        // multi-allelic sites where AD is trimmed to only the called allele).
        // Suppress for single-sample VCFs (no normal column): vcf2maf.pl skips the override
        // when there is no matched normal, matching the behavior for absent normal samples.
        let gt_is_explicit_het = tumor_vals.as_deref().and_then(|vals| {
            let gi = rec.format_keys.iter().position(|k| k == "GT")?;
            let gt = vals.get(gi)?;
            let idxs: Vec<usize> = gt.split(|c: char| c == '/' || c == '|')
                .filter_map(|s| s.parse().ok()).collect();
            Some(idxs.iter().any(|&i| i == 0) && idxs.iter().any(|&i| i > 0))
        }).unwrap_or(false);
        let tumor_seq_allele1 = if tumor_seq_allele1 == norm.ref_allele && normal_col.is_some() && !gt_is_explicit_het {
            if let Some(vals) = tumor_vals.as_deref() {
                let ad_full = rec.format_keys.iter().position(|k| k == "AD")
                    .and_then(|i| vals.get(i))
                    .map(|s| s.split(',').count() >= n_expected_ad)
                    .unwrap_or(false);
                if ad_full {
                    let ad = extract_depth(rec.format_keys.as_slice(), vals, &rec.ref_allele, &effective_alt, effective_alt_vcf_idx);
                    match (ad.alt_count, ad.depth()) {
                        (Some(alt), Some(total)) if total > 0 && alt as f64 / total as f64 >= 0.7 => {
                            norm.alt_allele.clone()
                        }
                        _ => tumor_seq_allele1,
                    }
                } else {
                    tumor_seq_allele1
                }
            } else {
                tumor_seq_allele1
            }
        } else {
            tumor_seq_allele1
        };

        let normal_ad_truncated = args.strict && normal_vals.as_deref()
            .and_then(|vals| rec.format_keys.iter().position(|k| k == "AD")
                .and_then(|ai| vals.get(ai)))
            .map(|ad_str| ad_str.split(',').count() < n_expected_ad)
            .unwrap_or(false);

        let mut extra = indexmap::IndexMap::new();
        for field in &args.retain_ann {
            let val = transcript.extra.get(field).cloned().unwrap_or_default();
            extra.insert(field.clone(), val);
        }

        // Match_Norm_Seq_Allele1/2: the actual alleles carried by the normal sample,
        // derived from the normal's GT (like Tumor_Seq_Allele1/2 is derived from tumor GT).
        // vcf2maf.pl maps each normal GT allele index to the normalized allele representation.
        // For "other" alleles (not REF, not effective_alt), apply the SAME prefix offset that was
        // computed for the effective_alt normalization — matching vcf2maf.pl's single-pass strip.
        let effective_prefix = (norm.pos - rec.pos) as usize;
        let (match_norm_seq_allele1, match_norm_seq_allele2) = normal_vals.as_deref()
            .and_then(|vals| {
                let gt_idx = rec.format_keys.iter().position(|k| k == "GT")?;
                let gt = vals.get(gt_idx)?;
                let idxs: Vec<usize> = gt
                    .split(|c| c == '/' || c == '|')
                    .filter_map(|s: &str| s.parse().ok())
                    .collect();
                if idxs.len() < 2 { return None; }
                let resolve = |idx: usize| -> String {
                    if idx == 0 {
                        norm.ref_allele.clone()
                    } else if idx == effective_alt_vcf_idx {
                        norm.alt_allele.clone()
                    } else {
                        rec.all_alts.get(idx - 1)
                            .map(|a| {
                                let stripped = &a[effective_prefix.min(a.len())..];
                                if stripped.is_empty() { "-".to_owned() } else { stripped.to_owned() }
                            })
                            .unwrap_or_else(|| norm.ref_allele.clone())
                    }
                };
                Some((resolve(idxs[0]), resolve(idxs[1])))
            })
            .unwrap_or_else(|| (norm.ref_allele.clone(), norm.ref_allele.clone()));

        let record = MafRecord {
            // vcf2maf.pl: Hugo_Symbol = SYMBOL → Transcript_ID → "Unknown"
            hugo_symbol: if !transcript.symbol.is_empty() {
                transcript.symbol.clone()
            } else if !transcript.feature.is_empty() {
                transcript.feature.clone()
            } else {
                String::from("Unknown")
            },
            entrez_gene_id: String::from("0"),
            center: args.maf_center.clone(),
            ncbi_build: genome_str.to_owned(),
            chromosome: rec.chrom.clone(),
            start_position: start_pos,
            end_position: end_pos,
            strand: String::from("+"),
            variant_classification: var_class.to_owned(),
            variant_type: var_type.to_owned(),
            reference_allele: norm.ref_allele.clone(),
            tumor_seq_allele1: tumor_seq_allele1.clone(),
            tumor_seq_allele2: norm.alt_allele.clone(),
            dbsnp_rs: dbsnp,
            dbsnp_val_status: String::new(),
            tumor_sample_barcode: tumor_barcode.to_owned(),
            matched_norm_sample_barcode: normal_barcode.to_owned(),
            match_norm_seq_allele1,
            match_norm_seq_allele2,
            mutation_status: String::new(),
            hgvsc: transcript.hgvsc.clone(),
            hgvsp: transcript.hgvsp.clone(),
            hgvsp_short,
            transcript_id: transcript.feature.clone(),
            exon_number: transcript.exon.clone(),
            t_depth: fmt_depth(tumor_depth.as_ref().and_then(|d| d.depth()), tumor_ad_truncated),
            t_ref_count: fmt_depth(tumor_depth.as_ref().and_then(|d| d.ref_count), tumor_ad_truncated),
            t_alt_count: fmt_depth(tumor_depth.as_ref().and_then(|d| d.alt_count), tumor_ad_truncated),
            n_depth: fmt_normal_depth(normal_depth.as_ref().and_then(|d| d.depth()), normal_ad_truncated),
            n_ref_count: fmt_normal_depth(normal_depth.as_ref().and_then(|d| d.ref_count), normal_ad_truncated),
            n_alt_count: fmt_normal_depth(normal_depth.as_ref().and_then(|d| d.alt_count), normal_ad_truncated),
            all_effects,
            vep_allele: transcript.allele.clone(),
            vep_gene: transcript.gene.clone(),
            vep_feature: transcript.feature.clone(),
            vep_biotype: transcript.biotype.clone(),
            vep_canonical: if transcript.canonical {
                String::from("YES")
            } else {
                String::new()
            },
            vep_sift: transcript.sift.clone(),
            vep_polyphen: transcript.polyphen.clone(),
            extra,
            ..Default::default()
        };

        // Build output: secondary SV breakpoint row first (matches vcf2maf.pl output order).
        let mut output = Vec::with_capacity(2);
        if let Some((chr2, bnd_pos, ref2)) = &sv_secondary {
            let sec_norm = normalize(*bnd_pos, ref2, &effective_alt);
            let (sec_start, sec_end) = maf_positions(&sec_norm);
            let sec_var_class = so_to_variant_classification(&transcript.consequences, &sec_norm.ref_allele, &sec_norm.alt_allele);
            let sec_var_type = variant_type(&sec_norm.ref_allele, &sec_norm.alt_allele);
            let sec_tsa1 = if tumor_seq_allele1 == norm.ref_allele {
                sec_norm.ref_allele.clone()
            } else {
                sec_norm.alt_allele.clone()
            };
            let mut secondary = record.clone();
            secondary.chromosome = chr2.clone();
            secondary.start_position = sec_start;
            secondary.end_position = sec_end;
            secondary.variant_classification = sec_var_class.to_owned();
            secondary.variant_type = sec_var_type.to_owned();
            secondary.reference_allele = sec_norm.ref_allele.clone();
            secondary.tumor_seq_allele1 = sec_tsa1;
            secondary.tumor_seq_allele2 = sec_norm.alt_allele.clone();
            secondary.match_norm_seq_allele1 = sec_norm.ref_allele.clone();
            secondary.match_norm_seq_allele2 = sec_norm.ref_allele.clone();
            output.push(secondary);
        }
        output.push(record);
        Ok(output)
    };

    // Parallel batch processing: read records in chunks, convert in parallel, write sequentially.
    // Preserves output order while utilizing all CPU cores for the conversion work.
    const BATCH: usize = 50_000;
    let mut batch: Vec<crate::vcf::VcfRecord> = Vec::with_capacity(BATCH);

    let drain = |batch: &[crate::vcf::VcfRecord], writer: &mut MafWriter<BufWriter<fs::File>>| -> Result<()> {
        let results: Vec<Result<Vec<MafRecord>>> = batch
            .par_iter()
            .map(|rec| convert_record(rec))
            .collect();
        for result in results {
            for maf_rec in result? {
                writer.write_record(&maf_rec)?;
            }
        }
        Ok(())
    };

    if let Some(rec) = first_rec {
        batch.push(rec);
    }
    loop {
        match vcf.next_record()? {
            Some(rec) => {
                batch.push(rec);
                if batch.len() >= BATCH {
                    drain(&batch, &mut writer)?;
                    batch.clear();
                }
            }
            None => break,
        }
    }
    if !batch.is_empty() {
        drain(&batch, &mut writer)?;
    }

    writer.flush()?;
    info!("MAF written to {}", args.output_maf.display());
    Ok(())
}

/// Parse CHR2 and END from a BND ALT allele string.
/// Handles all four VCF 4.x BND notations where the partner locus is between
/// matching `]` or `[` bracket characters: `]chr:pos]base`, `base]chr:pos]`,
/// `[chr:pos[base`, `base[chr:pos[`.
fn parse_bnd_alt(alt: &str) -> Option<(String, u64)> {
    let bracket = if alt.contains(']') { ']' } else { '[' };
    let mut positions = alt.match_indices(bracket).map(|(i, _)| i);
    let first = positions.next()?;
    let second = positions.next()?;
    let content = &alt[first + 1..second];
    let colon = content.rfind(':')?;
    let chr2 = content[..colon].to_owned();
    let pos: u64 = content[colon + 1..].parse().ok()?;
    Some((chr2, pos))
}

/// Check that the input VCF's chromosomes appear in the FASTA index (.fai).
/// Samples up to 10 non-symbolic variants; if none of their chromosomes are found
/// in the .fai, the ref-fasta is almost certainly a different genome build.
/// Mirrors vcf2maf.pl's behavior: it dies when samtools fetches 0 flanking sequences.
fn validate_ref_fasta(fasta: &Path, vcf: &Path) -> Result<()> {
    let fai_path = PathBuf::from(format!("{}.fai", fasta.display()));
    let fai_content = fs::read_to_string(&fai_path)
        .with_context(|| format!("Cannot read FASTA index {}", fai_path.display()))?;
    let fai_chroms: HashSet<&str> = fai_content
        .lines()
        .filter_map(|l| l.split('\t').next())
        .collect();

    let f = fs::File::open(vcf)?;
    let is_gz = vcf.extension().map_or(false, |e| e == "gz");
    let lines: Box<dyn Iterator<Item = String>> = if is_gz {
        Box::new(BufReader::new(MultiGzDecoder::new(f)).lines().flatten())
    } else {
        Box::new(BufReader::new(f).lines().flatten())
    };

    let mut sampled: Vec<String> = Vec::new();
    for line in lines {
        if line.starts_with('#') { continue; }
        let mut cols = line.splitn(6, '\t');
        let chrom = match cols.next() { Some(c) => c.to_owned(), None => continue };
        // Skip the ALT check — just sample chromosomes from any data line
        sampled.push(chrom);
        if sampled.len() >= 10 { break; }
    }

    if sampled.is_empty() {
        return Ok(());
    }

    if sampled.iter().any(|c| fai_chroms.contains(c.as_str())) {
        return Ok(());
    }

    bail!(
        "--ref-fasta '{}' does not match the genome build of input VCF '{}'. \
         VCF chromosome '{}' was not found in the FASTA index. \
         Check that the genome build and chromosome naming (chr-prefix) are consistent.",
        fasta.display(), vcf.display(), sampled[0]
    );
}

fn fasta_chrom_has_chr(fasta: &Path) -> bool {
    fs::File::open(fasta)
        .ok()
        .and_then(|f| {
            BufReader::new(f).lines().next().and_then(|l| l.ok())
        })
        .map(|l| l.starts_with(">chr"))
        .unwrap_or(false)
}

fn vcf_chrom_has_chr(vcf: &Path) -> bool {
    let f = match fs::File::open(vcf) {
        Ok(f) => f,
        Err(_) => return false,
    };
    let is_gz = vcf.extension().and_then(|e| e.to_str()) == Some("gz");
    let lines: Box<dyn Iterator<Item = String>> = if is_gz {
        Box::new(BufReader::new(MultiGzDecoder::new(f)).lines().flatten())
    } else {
        Box::new(BufReader::new(f).lines().flatten())
    };
    for line in lines {
        if !line.starts_with('#') {
            return line.starts_with("chr");
        }
    }
    false
}

fn strip_chr_prefix<R: std::io::BufRead, W: std::io::Write>(reader: R, mut writer: W) -> Result<()> {
    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') {
            writeln!(writer, "{}", line)?;
        } else if let Some(rest) = line.strip_prefix("chr") {
            writeln!(writer, "{}", rest)?;
        } else {
            writeln!(writer, "{}", line)?;
        }
    }
    Ok(())
}

fn run_fastvep(
    fastvep: &Path,
    input_vcf: &Path,
    ref_fasta: &Path,
    gff3: &Path,
) -> Result<NamedTempFile> {
    // fastVEP cannot read gzip-compressed GFF3 — decompress to a sibling file so that
    // fastVEP's transcript cache (.fastvep.cache) persists across runs.
    let gff3_owned: PathBuf;
    let gff3_path: &Path = if gff3.extension().and_then(|e| e.to_str()) == Some("gz") {
        let stem = gff3.file_stem().unwrap_or_default();
        let sibling = gff3.with_file_name(stem);
        if !sibling.exists() {
            let src = fs::File::open(gff3)
                .with_context(|| format!("Cannot open GFF3: {}", gff3.display()))?;
            let mut dec = MultiGzDecoder::new(BufReader::new(src));
            let mut out = fs::File::create(&sibling)
                .with_context(|| format!("Cannot write decompressed GFF3 to {}", sibling.display()))?;
            std::io::copy(&mut dec, &mut out).context("Failed to decompress GFF3")?;
        }
        gff3_owned = sibling;
        &gff3_owned
    } else {
        gff3
    };

    // fastVEP cannot read gzip-compressed VCF, and cannot match "chr1" against "1".
    // Pre-process the input VCF into a plain temp file when either condition applies.
    let fasta_has_chr = fasta_chrom_has_chr(ref_fasta);
    let vcf_has_chr = vcf_chrom_has_chr(input_vcf);
    let is_gz = input_vcf.extension().and_then(|e| e.to_str()) == Some("gz");
    let strip_chr = vcf_has_chr && !fasta_has_chr;
    let _stripped_vcf: Option<NamedTempFile>;
    let input_vcf_path: &Path = if strip_chr || is_gz {
        let mut tmp_vcf = NamedTempFile::with_suffix(".vcf")
            .context("Cannot create temp file for VCF preprocessing")?;
        let src = fs::File::open(input_vcf)?;
        let writer = BufWriter::new(tmp_vcf.as_file_mut());
        if strip_chr {
            if is_gz {
                strip_chr_prefix(BufReader::new(MultiGzDecoder::new(src)), writer)?;
            } else {
                strip_chr_prefix(BufReader::new(src), writer)?;
            }
        } else {
            // Decompress only — fastVEP requires plain VCF.
            let mut dec = BufReader::new(MultiGzDecoder::new(src));
            let mut w = writer;
            std::io::copy(&mut dec, &mut w).context("Failed to decompress VCF for fastVEP")?;
        }
        _stripped_vcf = Some(tmp_vcf);
        _stripped_vcf.as_ref().unwrap().path()
    } else {
        _stripped_vcf = None;
        input_vcf
    };

    let tmp = NamedTempFile::new().context("Cannot create temp file for fastVEP output")?;
    let status = Command::new(fastvep)
        .args([
            "annotate",
            "-i",
            input_vcf_path.to_str().unwrap(),
            "-o",
            tmp.path().to_str().unwrap(),
            "--fasta",
            ref_fasta.to_str().unwrap(),
            "--gff3",
            gff3_path.to_str().unwrap(),
            "--hgvs",
            "--output-format",
            "vcf",
        ])
        .status()
        .context("Failed to run fastVEP. Run `mafsmith fetch` to install it.")?;

    if !status.success() {
        bail!("fastVEP exited with status {}", status);
    }
    Ok(tmp)
}

fn resolve_vep(cli: &Option<PathBuf>) -> Result<PathBuf> {
    if let Some(p) = cli {
        return Ok(p.clone());
    }
    let out = Command::new("which").arg("vep").output();
    if let Ok(o) = out {
        if o.status.success() {
            let path = String::from_utf8_lossy(&o.stdout).trim().to_owned();
            return Ok(PathBuf::from(path));
        }
    }
    bail!(
        "Ensembl VEP not found in PATH. Install it (e.g. `conda install -c bioconda ensembl-vep`) \
         or pass --vep-path."
    );
}

fn run_vep(
    vep: &Path,
    input_vcf: &Path,
    vep_data: &Path,
    assembly: &str,
    forks: u32,
) -> Result<NamedTempFile> {
    let tmp = NamedTempFile::new().context("Cannot create temp file for VEP output")?;

    let mut cmd = Command::new(vep);
    cmd.args([
        "--input_file",  input_vcf.to_str().unwrap(),
        "--output_file", tmp.path().to_str().unwrap(),
        "--format",      "vcf",
        "--vcf",
        "--cache",
        "--offline",
        "--dir_cache",   vep_data.to_str().unwrap(),
        "--assembly",    assembly,
        "--hgvs",
        "--hgvsg",
        "--symbol",
        "--biotype",
        "--canonical",
        "--mane",
        "--no_stats",
        "--buffer_size", "5000",
        "--force_overwrite",
    ]);
    if forks > 0 {
        cmd.args(["--fork", &forks.to_string()]);
    }

    let status = cmd.status()
        .context("Failed to run VEP. Is it installed and in PATH?")?;

    if !status.success() {
        bail!("VEP exited with status {}", status);
    }
    Ok(tmp)
}

fn resolve_ref_fasta(
    cli: &Option<PathBuf>,
    data_dir: &Path,
    genome: &str,
) -> Result<PathBuf> {
    if let Some(p) = cli {
        return Ok(p.clone());
    }
    let p = data_dir.join(genome).join("reference.fa");
    if !p.exists() {
        bail!(
            "Reference FASTA not found at {}. Run `mafsmith fetch --genome {}` first.",
            p.display(),
            genome.to_lowercase()
        );
    }
    Ok(p)
}

fn resolve_gff3(cli: &Option<PathBuf>, data_dir: &Path, genome: &str) -> Result<PathBuf> {
    if let Some(p) = cli {
        return Ok(p.clone());
    }
    // Accept both compressed and uncompressed
    for name in &["genes.gff3.gz", "genes.gff3"] {
        let p = data_dir.join(genome).join(name);
        if p.exists() {
            return Ok(p);
        }
    }
    bail!(
        "GFF3 not found under {}. Run `mafsmith fetch --genome {}` first.",
        data_dir.join(genome).display(),
        genome.to_lowercase()
    );
}

fn resolve_fastvep(cli: &Option<PathBuf>, data_dir: &Path) -> Result<PathBuf> {
    if let Some(p) = cli {
        return Ok(p.clone());
    }
    // Check our managed bin dir first
    let managed = data_dir.join("bin").join("fastvep");
    if managed.exists() {
        return Ok(managed);
    }
    // Fall back to $PATH
    if let Ok(path) = which_fastvep() {
        return Ok(path);
    }
    bail!("fastVEP binary not found. Run `mafsmith fetch` to install it.");
}

fn which_fastvep() -> Result<PathBuf> {
    let out = Command::new("which").arg("fastvep").output()?;
    if out.status.success() {
        let path = String::from_utf8_lossy(&out.stdout).trim().to_owned();
        Ok(PathBuf::from(path))
    } else {
        bail!("fastvep not in PATH")
    }
}

fn load_custom_enst(path: Option<&Path>) -> Result<Option<HashSet<String>>> {
    let Some(p) = path else { return Ok(None) };
    let content = fs::read_to_string(p)
        .with_context(|| format!("Cannot read custom ENST file: {}", p.display()))?;
    let set: HashSet<String> = content
        .lines()
        .map(|l| l.trim())
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .map(|l| l.split_whitespace().next().unwrap_or(l).to_owned())
        .collect();
    Ok(Some(set))
}

/// Determine which VCF sample column to use for tumor or normal.
/// Returns None if no samples present (single-sample VCF without tumor/normal designation).
fn resolve_sample_col(
    cli: Option<&str>,
    sample_names: &[String],
    role: &str,
    is_tumor: bool,
) -> Result<Option<String>> {
    if let Some(name) = cli {
        if sample_names.contains(&name.to_owned()) {
            return Ok(Some(name.to_owned()));
        }
        warn!(
            "Specified {} sample '{}' not found in VCF samples {:?}",
            role, name, sample_names
        );
        return Ok(None);
    }
    if sample_names.is_empty() {
        return Ok(None);
    }
    // Heuristic: if two samples, first = tumor, second = normal
    if sample_names.len() >= 2 {
        return Ok(Some(sample_names[if is_tumor { 0 } else { 1 }].clone()));
    }
    if is_tumor {
        Ok(Some(sample_names[0].clone()))
    } else {
        Ok(None)
    }
}

/// Look up a single reference base from an indexed FASTA file.
///
/// Reads the `.fai` index (same path as fasta + ".fai") to locate the sequence,
/// then seeks directly to the byte that holds the requested 1-based position.
/// Tries the chromosome name as given, then with/without "chr" prefix if not found.
fn fasta_fetch_base(fasta_path: &Path, chrom: &str, pos: u64) -> Option<char> {
    let fai_path = PathBuf::from(format!("{}.fai", fasta_path.display()));
    let fai = fs::read_to_string(&fai_path).ok()?;

    let lookup = |name: &str| -> Option<(u64, u64, u64)> {
        for line in fai.lines() {
            let cols: Vec<&str> = line.split('\t').collect();
            if cols.len() >= 5 && cols[0] == name {
                let offset: u64 = cols[2].parse().ok()?;
                let bases_per_line: u64 = cols[3].parse().ok()?;
                let bytes_per_line: u64 = cols[4].parse().ok()?;
                return Some((offset, bases_per_line, bytes_per_line));
            }
        }
        None
    };

    let (offset, bases_per_line, bytes_per_line) = lookup(chrom)
        .or_else(|| {
            // Try toggling the "chr" prefix
            if let Some(stripped) = chrom.strip_prefix("chr") {
                lookup(stripped)
            } else {
                lookup(&format!("chr{}", chrom))
            }
        })?;

    if bases_per_line == 0 || pos == 0 { return None; }
    let line_num = (pos - 1) / bases_per_line;
    let col = (pos - 1) % bases_per_line;
    let byte_pos = offset + line_num * bytes_per_line + col;

    let mut f = fs::File::open(fasta_path).ok()?;
    f.seek(SeekFrom::Start(byte_pos)).ok()?;
    let mut buf = [0u8; 1];
    f.read_exact(&mut buf).ok()?;
    let c = buf[0] as char;
    if c.is_ascii_alphabetic() { Some(c.to_ascii_uppercase()) } else { None }
}
