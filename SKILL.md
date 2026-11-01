---
name: dia-evidence-extraction
description: Use for extracting observed MS2 spectra from DIA raw files with RT-anchored targeted extraction, rendering isoPSV-style evidence pages (ion ladder, observed spectrum, inverted AlphaPeptDeep prediction, ppm panel, RT heatmap), generating Table 1 cross-reference PDFs, and auditing results. Triggers on "extract DIA spectra", "evidence page", "spectrum evidence", "RT heatmap", "Table 1 PDF", "audit evidence", or paths to DIA raw directories with search results.
user-invocable: true
argument-hint: "[search results path] [raw directory] [options]"
---

# DIA spectrum evidence extraction

## When to use

- DIA search results (with peptide sequences, charges, m/z, per-sample RTs) need to be backed by actual observed MS2 spectra from raw files.
- Reviewers request spectrum-level evidence for specific peptides (e.g., pathogen-derived proteins in a clinical study).
- A manuscript table (Table 1) lists peptides that need per-peptide evidence pages.
- Evidence PDFs need to be audited for data accuracy and presentation quality.

## When not to use

- DDA data (use `ms.s1.acquisition-eval` for parameter evaluation instead).
- No per-sample RT information is available (the RT anchor is essential for targeted extraction).
- Quantitative analysis is the goal (use DIA-NN / Spectronaut directly).

## Preconditions & evidence

| Requirement | Evidence source |
|---|---|
| DIA search results table | Spectronaut report or S14-style table with seq/charge/m/z/sample/RT |
| Thermo .raw files | Directory of raw DIA acquisitions |
| isoPSV R library | `risoPSV/` directory with compat.R, ReadMGF_script.R, chimeric_check.R, AA_mass_dataframe.txt |
| AlphaPeptDeep | `peptdeep` executable or Python package |
| Table 1 xlsx (optional) | Manuscript table for cross-reference PDF |

Probe once:
1. Does the search results table contain per-sample RT values? (Column check)
2. Do raw files exist and can RawFileReader open them? (rawdump probe)
3. Does the isoPSV R library load under the system R? (`Rscript smoke.R`)

## Workflow

The pipeline is a chain of numbered scripts. Each takes the previous output as input.

### Phase 1: Target building and extraction
1. **70_build_s14_targets.py** — Parse search results → target table (seq/z/mz/sample/RT/qty).
   Key: resolves modification states from m/z deltas; maps sample names to raw files.
2. **71_extract_targeted.py** — RT-anchored targeted extraction. For each (precursor, sample):
   theoretical b/y → window assignment → RT ± 0.6 min candidate scans → sequence-ladder
   apex scoring → best scan export + MS1 precursor envelope check.
   Env vars: `TGT` (target parquet), `OUTD` (output directory).
3. **73_profiles.py** — Full-RT fragment XIC profiles for heatmaps.
   Three-tier background removal: ion-level 10% apex floor → spectrum-correlation
   (cos × pres) → stored as per-scan profiles.
   Env vars: `TGT`, `ROWSX`, `OUTD`.

### Phase 2: Evidence page rendering
4. **76_build_isopsv_inputs.py** — Build isoPSV R engine inputs (pFind-style table + MGF).
5. **78_s14_export.R** — isoPSV engine renders top block (info lines, ion ladder,
   observed spectrum, inverted prediction, ppm panel). Uses R base graphics only.
6. **79_compose_a4.py** — Compose A4 page: info text (Arial, auto-shrink) + R top block +
   RT heatmap (matplotlib PNG, γ² sharpening + per-row 20% cutoff + 1.5× absolute
   background line + blank-row quantity annotation).

### Phase 3: Supplementary chains
7. **82_extra_targets.py** — For peptides missing from search results: use v9-master
   extraction RTCorrected as anchor.
8. **84_extra_inputs.py** — Build isoPSV inputs for the extra chain.
9. **81_table1_pdf.py** — Table 1 cross-reference PDF (ordered by manuscript table).
10. **86_predict_missing.py** — On-the-fly AlphaPeptDeep for peptides without predictions.

### Phase 4: Audit
11. **85_audit.py** — Three-layer audit:
    A: data reconciliation (S14→targets→extraction→inputs→MGF→Info fields→cosine→heatmap stats)
    B: Table 1 cross-reference (page order, headers, chain provenance)
    C: presentation (key text presence, heatmap non-blank, text overflow)

## Decision rules

| Decision | Rule | Evidence |
|---|---|---|
| RT anchor source | S14 per-sample RT if available; v9-master RTCorrected for missing peptides | Search table vs master table coverage |
| Apex scoring | (bLad+yLad, nmatched, intensity) — sequence ladder first, not raw intensity | DIA mixed spectra contain strong co-isolated interference |
| Background removal | Three tiers: ion-level 10% apex floor → row-level 2× baseline → 1.5× absolute line | Correlation profiles measure background texture for weak peptides |
| Heatmap display | γ² sharpen → 20% row cutoff → 1.5× absolute cutoff; blank rows show S14 quantity | Sparse low-m/z windows have ~956 peaks/scan background |
| Prediction fallback | AlphaPeptDeep on-the-fly for peptides outside the original prediction set | pkl keys derived from Spectronaut report coverage |

## Validation gates

- After 71: `ok_rate == 1.0` and `abs_drt_med < 0.1 min` in summary.json.
- After 73: profile npz row count matches targets × samples.
- After 78: R output page count == input table row count.
- After 79: `total_pages == 1 + n_peptides + n_index_pages`.
- After 85: `errors == 0`.

## Failure modes & recovery

| Failure | Detection | Recovery |
|---|---|---|
| RawFileReader can't open raw | `open_raw` retries 3× then fails | Check .NET Core 3.1 + DLL paths |
| R script fails on Chinese path | `setwd` error | Copy risoPSV/ to ASCII-only path |
| R locale breaks UTF-8 comments | `invalid input found on input connection` | `env -u LC_CTYPE -u LC_ALL -u LANG Rscript` |
| PDF output Permission denied | `pymupdf.FzErrorSystem: cannot remove file` | Close the PDF viewer that holds the lock |
| insert_image silent no-resize | Heatmap same size after rect change | Match figsize aspect to rect aspect exactly |
| Background wash in heatmap | Row coverage > 50% of RT range | Three-tier background removal (already built in) |
| Peptides missing from S14 | 81 outputs `no_spectrum` placeholders | Run 82→71→73→84 for the extra chain |

## Output contract

| Artifact | Format | Schema | Description |
|---|---|---|---|
| `evidence_pdf` | `.pdf` | — | Per-peptide A4 evidence pages (isoPSV layout + RT heatmap) |
| `table1_pdf` | `.pdf` | — | Table 1 cross-reference evidence PDF |
| `audit_report` | `.md` + `.json` | `schemas/output.json` | Three-layer audit results |
| `rows_all.parquet` | `.parquet` | — | Per-(precursor,sample) extraction results |
| `spectra_all.npz` | `.npz` | — | Extracted apex spectra (concatenated arrays) |
| `profiles.npz` | `.npz` | — | Full-RT correlation XIC profiles |
| `pages_manifest.tsv` | `.tsv` | — | Page-to-peptide index |

## Safety & approvals

- **Network**: none (except optional AlphaPeptDeep model download on first run).
- **Filesystem**: read-write (writes to `work/output/` only; does not modify inputs).
- **Secrets**: none required.
- **Destructive**: no. Overwrites previous outputs in the same directory.

## References

- [Configuration reference](docs/configuration.md) — env vars and script parameters
- [Troubleshooting](docs/troubleshooting.md) — 15+ known issues and fixes
- [Architecture](docs/architecture.md) — script chain diagram and data flow
- [isoPSV](https://github.com/Songhao-Jiang/isoPSV) — spectrum visualization engine
- [eval_ms_pipeline](https://github.com/AI4MSOmics/eval_ms_pipeline_repo) — companion DDA parameter evaluation

## Examples

### Example 1: Full chain from Spectronaut report

```bash
# 1. Build targets from search results
python scripts/70_build_s14_targets.py

# 2. RT-anchored targeted extraction
python scripts/71_extract_targeted.py

# 3. Full-RT profiles for heatmap
python scripts/73_profiles.py

# 4. Build isoPSV inputs
python scripts/76_build_isopsv_inputs.py

# 5. R engine renders top blocks
cd risoPSV && env -u LC_CTYPE -u LC_ALL -u LANG \
  OUT_PDF=../out/s14_topblock.pdf Rscript 78_s14_export.R

# 6. Compose final A4 PDFs
python scripts/79_compose_a4.py

# 7. Audit
python scripts/85_audit.py
```

### Example 2: Table 1 cross-reference PDF

```bash
# Extra chain for peptides missing from S14
python scripts/82_extra_targets.py
TGT=work/tmp/s14_targets_extra.parquet OUTD=work/output/extract_extra \
  python scripts/71_extract_targeted.py
TGT=work/tmp/s14_targets_extra.parquet \
  ROWSX=work/output/extract_extra/rows_all.parquet \
  OUTD=work/output/profiles_extra \
  python scripts/73_profiles.py
python scripts/84_extra_inputs.py
# R render extra topblock
cd risoPSV && env -u LC_CTYPE -u LC_ALL -u LANG \
  IN_TSV=in_extra/s14_spectra.tsv \
  OUT_PDF=out/topblock_extra.pdf Rscript 78_s14_export.R
# Table 1 PDF
python scripts/81_table1_pdf.py
```
