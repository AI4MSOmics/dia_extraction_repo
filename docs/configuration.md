# Configuration reference

## Language

**Default: English.** All user-facing text (PDF info lines, heatmap titles, cover/index pages, audit report, R panel labels) is in English by default.

To switch to Chinese, set the `DIA_LANG` environment variable:

```bash
# Chinese output
export DIA_LANG=zh
python scripts/79_compose_a4.py

# R rendering (also picks up DIA_LANG via lang_strings.R)
DIA_LANG=zh env -u LC_CTYPE Rscript risoPSV/78_s14_export.R

# Back to English (default)
unset DIA_LANG
# or
export DIA_LANG=en
```

Translation is handled by `lang.py` (Python) and `risoPSV/lang_strings.R` (R), both reading the same `DIA_LANG` variable. Code comments and script logs remain in English (they don't affect deliverables).

## Environment variables

| Variable | Used by | Default | Description |
|----------|---------|---------|-------------|
| `DIA_LANG` | all | `en` | Output language: `en` (default) or `zh` |
| `TGT` | 71, 73 | `work/tmp/s14_targets.parquet` | Target table path (parquet) |
| `OUTD` | 71, 73 | `work/output/s14_evidence/extract` (71) or `.../profiles` (73) | Output directory |
| `ROWSX` | 73 | `work/output/s14_evidence/extract/rows_all.parquet` | Extraction results path (for profile apex determination) |

### R script (78_s14_export.R)

| Variable | Default | Description |
|----------|---------|-------------|
| `IN_TSV` | `risoPSV/in/s14_spectra.tsv` | Input table path |
| `OUT_PDF` | `risoPSV/out/s14_topblock.pdf` | Output PDF path |
| `SEQ` | (empty = all) | Filter to single sequence |

## Script parameters (hardcoded in script headers, edit to change)

### 71_extract_targeted.py

| Parameter | Default | Description |
|-----------|---------|-------------|
| RT tolerance | 0.6 min (widened to 1.2 if < 2 scans) | Apex search window around S14 RT |
| Fragment ppm | 20 | Fragment m/z match tolerance |
| MS1 ppm | 10 | Precursor envelope match tolerance |
| Min ions for candidate | bLad≥2, yLad≥2, nm≥4 | Sequence-ladder threshold |

### 73_profiles.py

| Parameter | Default | Description |
|-----------|---------|-------------|
| Ion-level background floor | 0.10 | Ion counts only if ≥ 10% of its apex intensity |
| Profile matching ppm | 10 (profiles) / 20 (apex search) | Narrower for XIC to reduce noise |

### 79_compose_a4.py

| Parameter | Default | Description |
|-----------|---------|-------------|
| Heatmap figsize | 8.0 × 5.87 in | Must match rect aspect ratio exactly |
| Heatmap axes | [0.135, 0.108, 0.775, 0.816] | Three-band allocation (title/band/xlabel) |
| γ sharpening | 2 (squared) | Correlation value gamma for contrast |
| Per-row display cutoff | 20% | Bins below 20% of row max are blanked |
| Row quant threshold | 2.0× | Elution window peak must exceed baseline |
| Absolute background line | 1.5× | Bin must exceed row baseline |
| Heatmap rect on A4 | (20, 431, 575, 838) pt | Position within A4 page |
| Top block rect | (2, top_h, 593, 442) pt | R block position |

### 82_extra_targets.py

| Parameter | Default | Description |
|-----------|---------|-------------|
| Precursor selection | tier A > RTQC=0 > ms1/consensus > NMatchedFrags | Best precursor per missing sequence |

## Platform requirements

| Component | Version | Notes |
|-----------|---------|-------|
| Python | ≥3.8 | numpy, pandas, matplotlib, pymupdf, openpyxl, h5py |
| R | ≥4.0 | base graphics only (no ggplot2) |
| .NET Core | 3.1 | For rawdump.exe (ThermoFisher RawFileReader) |
| AlphaPeptDeep | 1.5.0+ | `peptdeep` CLI or Python package |
| ThermoFisher DLLs | latest | RawFileReader.dll + Data.dll |

## Known platform quirks

1. **R + Chinese paths**: `setwd()` fails → copy `risoPSV/` to ASCII-only path.
2. **R + Git Bash locale**: `LC_CTYPE=C.UTF-8` breaks UTF-8 R source → `env -u LC_CTYPE -u LC_ALL -u LANG Rscript`.
3. **PDF viewer lock**: `pymupdf.FzErrorSystem: cannot remove file` → close the PDF viewer before re-saving.
4. **insert_image aspect**: `pymupdf.insert_image` silently preserves aspect ratio → figsize must match rect aspect exactly.
