# rewrites.bio principles analysis

Assessment of how well mafsmith implements the [rewrites.bio](https://rewrites.bio) best practice
principles for AI-assisted bioinformatics rewrites.

---

## Summary

| Section | Principle | Status |
|---------|-----------|--------|
| Philosophy | Credit the Original Authors | Partial |
| Philosophy | Emulate Exactly | Strong |
| Philosophy | Be Transparent About AI | Missing |
| Planning | Think Big | Strong |
| Planning | Work Small | Strong |
| Building | Test and Benchmark with Real Data | Strong |
| Building | Build Only What You Need | Strong |
| Building | Pin Versions and Document | Partial |
| Stewardship | Maintain and Govern | Partial |
| Stewardship | Preserve Compatibility | Strong |
| Stewardship | Release as Open Source | Strong |
| Stewardship | Contribute Upstream Responsibly | Not addressed |

---

## Philosophy

### 1. Credit the Original Authors — Partial

**What's done well:**

The README has a prominent Acknowledgements section that names the original tool
([vcf2maf](https://github.com/mskcc/vcf2maf)), links to the repo, and credits
[@ckandoth](https://github.com/ckandoth) with his ORCID
([0000-0002-1345-3573](https://orcid.org/0000-0002-1345-3573)). The repo description explicitly
frames the project as an adaptation of vcf2maf's "design, field conventions, and edge-case
handling."

**Gaps:**

- No `CITATION.cff` file. Without one, GitHub's "Cite this repository" button doesn't know to
  surface the original paper or guide users to cite Kandoth et al.
- No link to the original vcf2maf publication or preprint in the README.
- No explicit statement asking downstream users to also cite the original tool. The rewrites.bio
  principle specifically calls for this: "When users cite your rewrite, they should know to cite
  the original too."

**Recommendation:** Add a `CITATION.cff` and include a sentence in the README directing users to
cite the original vcf2maf paper alongside mafsmith.

---

### 2. Emulate Exactly — Strong

**What's done well:**

This is the clearest strength of the project. mafsmith achieves 0 conversion-field differences
against `vcf2maf.pl --inhibit-vep` in `--strict` mode across 21 caller types and 6+ reference
datasets (GIAB, SEQC2, ICGC PCAWG, DepMap CCLE, COSMIC). All four subcommands (`vcf2maf`,
`maf2vcf`, `vcf2vcf`, `maf2maf`) are validated. Column ordering, headers, field naming, and
edge-case behavior (multi-allelic, GVCF, SV) are all reproduced.

A "Known intentional differences" section in the README is transparent about the handful of places
where mafsmith deliberately diverges from vcf2maf.pl — and in most cases mafsmith's behaviour is
*more* correct (e.g. SV secondary rows with actual partner chromosome vs. blank fields in
vcf2maf.pl).

**Gaps:**

None substantive. This principle is met thoroughly.

---

### 3. Be Transparent About AI — Missing

**What's done well:**

The rewrites.bio badge has been added to the README.

**Gaps:**

There is no disclosure anywhere in the repository (README, docs, or CLAUDE.md) stating that:

- Claude (or any LLM) was used to write the code
- Which AI tools were used and in what role (code generation, debugging, manuscript drafting, etc.)
- How AI-generated code was validated against the original tool

This is the most significant gap relative to the rewrites.bio principles. The entire rationale for
the badge is to signal transparency, which requires the disclosure to actually be present.

**Recommendation:** Add an "AI assistance" section to the README documenting that mafsmith was
written with Claude (Anthropic), describing the human-in-the-loop validation workflow
(iterative comparison against vcf2maf.pl output on real VCFs, bug-fix-driven fixture expansion),
and noting that all correctness was verified by output comparison against reference tools, not
by trusting AI-generated logic directly.

---

## Planning

### 1. Think Big — Strong

**What's done well:**

mafsmith does not just rewrite vcf2maf — it also integrates fastVEP as the annotation backend, provides a `fetch` subcommand that downloads and installs all reference data automatically, and supports the full four-subcommand surface of the vcf2maf suite. The design explicitly considers the upstream/downstream pipeline context (fastVEP replacing VEP, drop-in compatibility enabling one-line pipeline config changes).

**Gaps:**

None substantive.

---

### 2. Work Small — Strong

**What's done well:**

The CLAUDE.md development log shows a highly iterative process: each bug fix is tied to a specific
VCF record added to the integration fixture, with the root cause, fix location, and validating
dataset all documented. No large speculative features were implemented. The fixture grows only
when a real edge case requires coverage.

**Gaps:**

This process is not surfaced in user-facing documentation. A brief note in the README or
CONTRIBUTING.md about the development methodology would make the validation story clearer for
contributors and users.

---

## Building

### 1. Test and Benchmark with Real Data — Strong

**What's done well:**

Validation spans 21 caller types, 4 subcommands, GRCh37 and GRCh38, paired T/N and single-sample,
germline and somatic, and multiple real published cohorts (GIAB NIST v4.2.1, SEQC2 HCC1395, ICGC
PCAWG, DepMap CCLE WGS, COSMIC). Benchmarks are documented with hardware (AWS c6a.4xlarge, AMD
EPYC 7R13, 16 vCPU, 30 GiB RAM), exact commands, mean ± std timing, and per-sample breakdowns.
Cost and carbon savings are estimated.

**Gaps:**

None substantive.

---

### 2. Build Only What You Need — Strong

**What's done well:**

GRCm39 is implemented but explicitly flagged as "Available; not yet validated against vcf2maf.pl."
Unrecognized symbolic ALTs (`<INS>`, `<CNV>`) are dropped, matching vcf2maf.pl behaviour and
documented in "Known intentional differences." No speculative features were added.

**Gaps:**

The README does not explicitly state what happens when unsupported features are invoked (i.e., that
they fail loudly rather than silently producing wrong output). A brief note would strengthen this.

---

### 3. Pin Versions and Document — Partial

**What's done well:**

The README specifies which tool was validated against (`vcf2maf.pl --inhibit-vep`, VEP 112/115),
which datasets were used (with Synapse IDs, GIAB version, SEQC2 cohort), and which exact commands
were run. The CLAUDE.md development log pins individual VCF files by Synapse ID.

**Gaps:**

- The exact version of `vcf2maf.pl` used for validation is not stated in the README. Equivalence
  claims implicitly depend on a specific version.
- No date of last validation run is recorded anywhere user-facing.
- No stated policy for what happens when vcf2maf.pl releases a new version (re-validation
  trigger, regression testing process).

**Recommendation:** Add a "Validated against" line to the README: `vcf2maf.pl vX.Y.Z, as of
YYYY-MM`. Add a brief policy on re-validation after upstream releases.

---

## Stewardship

### 1. Maintain and Govern — Partial

**What's done well:**

CI is configured (`.github/workflows/ci.yml`) with build, test, clippy, and fmt jobs. Versioned
tags were introduced with the v0.1.0 release. The CLAUDE.md maintains a detailed internal
development log.

**Gaps:**

- No `CONTRIBUTING.md`. Contributors have no guidance on how to submit issues or PRs, what
  validation is required for new caller support, or how edge cases should be documented.
- No `CHANGELOG.md`. There is no user-facing record of what changed between releases.
- No documented issue acknowledgment or triage process.

**Recommendation:** Create a minimal `CONTRIBUTING.md` covering: how to report a new caller
discrepancy, what a fix needs to include (regression test VCF record + CLAUDE.md entry), and how
to run validation locally. Create a `CHANGELOG.md` starting from v0.1.0.

---

### 2. Preserve Compatibility — Strong

**What's done well:**

mafsmith preserves MAF v2.4 column headers and ordering exactly. The `--strict` flag provides
bit-for-bit output compatibility with vcf2maf.pl for downstream scripts and MultiQC modules. The
README explicitly positions the tool as a drop-in replacement, and the "Known intentional
differences" section makes deviations opt-in and documented.

**Gaps:**

None substantive.

---

### 3. Release as Open Source — Strong

**What's done well:**

Apache 2.0 license. Public GitHub repository under `nf-osi`. The permissive license allows HPC
use, forking, and integration into any downstream tool or pipeline.

**Gaps:**

None substantive.

---

### 4. Contribute Upstream Responsibly — Not addressed

**What's done well:**

The CLAUDE.md documents several cases where vcf2maf.pl has known bugs that mafsmith fixes (e.g.
SV secondary rows with blank Chromosome, first-consequence short-circuit for multi-consequence
variants). These are documented but not acted on upstream.

**Gaps:**

There is no documented policy on how discoveries of vcf2maf.pl bugs are handled — whether they are
reported upstream, and if so, with what validation process (manually verified, minimal reproducible
example, no AI-generated test cases). The rewrites.bio principle specifically warns against
automating bug reports and requires manually verifiable examples.

**Recommendation:** Add a brief section to `CONTRIBUTING.md` (or the README) stating: "If
mafsmith reveals what appears to be a bug in vcf2maf.pl, verify it manually with the original
tool on clean data before filing an upstream issue. Create a minimal reproducible VCF example.
Do not use mafsmith output alone as evidence."

---

## Priority action items

1. **Add AI transparency disclosure to README** (principle: Be Transparent About AI) — highest
   priority given the badge is now present. State that Claude was used, describe the validation
   methodology, note that correctness was verified by output comparison, not by trusting generated
   logic.

2. **Add CONTRIBUTING.md** (principle: Maintain and Govern) — cover caller validation requirements,
   fixture contribution process, and upstream bug-reporting policy.

3. **Add CHANGELOG.md** (principle: Maintain and Govern) — start from v0.1.0.

4. **Pin vcf2maf.pl version in README** (principle: Pin Versions and Document) — "Validated
   against vcf2maf.pl vX.Y.Z, last verified YYYY-MM."

5. **Add CITATION.cff and cite-original instruction** (principle: Credit the Original Authors) —
   surface the vcf2maf publication and direct users to cite it alongside mafsmith.

6. **Document upstream contribution policy** (principle: Contribute Upstream Responsibly) — brief
   note on manual verification before filing bugs.
