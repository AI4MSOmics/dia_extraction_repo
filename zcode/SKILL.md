---
name: dia-evidence-extraction
description: Extract observed MS2 spectra from DIA raw files with RT-anchored targeted extraction. Renders isoPSV-style evidence pages (ion ladder, observed spectrum, inverted AlphaPeptDeep prediction, ppm panel, RT heatmap) and Table 1 cross-reference PDFs. Triggers on "extract DIA spectra", "evidence page", "RT heatmap", "Table 1 evidence PDF".
user-invocable: true
argument-hint: "[search results] [raw directory] [output directory]"
---

# DIA Spectrum Evidence Extraction

## Quick start

```bash
# Full chain (7 scripts, ~1 hour for 30 raw files)
python scripts/70_build_s14_targets.py
python scripts/71_extract_targeted.py
python scripts/73_profiles.py
python scripts/76_build_isopsv_inputs.py
cd risoPSV && env -u LC_CTYPE -u LC_ALL -u LANG Rscript 78_s14_export.R && cd ..
python scripts/79_compose_a4.py
python scripts/85_audit.py
```

## What you get

Per-peptide A4 evidence PDF pages:

| Section | Content |
|---------|---------|
| Info lines (3) | Protein, flanked peptide (isoPSV format), charge, mods, RT range |
| Ion ladder | Sequence with b (green, below) / y (red, above) matched sites |
| Observed spectrum | Best-PSM apex scan, grey background + colored matched peaks |
| Inverted prediction | AlphaPeptDeep theoretical spectrum (mirrored), cosine annotation |
| ppm panel | Mass error of matched ions (±20 ppm) |
| RT heatmap | Per-sample fragment XIC, white→red, 1-min bins |

## Script quick reference

| # | Name | Input → Output | Time |
|---|------|---------------|------|
| 70 | build_s14_targets | search xlsx/tsv → targets.parquet | seconds |
| 71 | extract_targeted | targets + raw → rows.parquet + spectra.npz | ~5 min |
| 73 | profiles | targets + raw → profiles.npz | ~35 min |
| 76 | build_isopsv_inputs | rows + pkl → tsv + mgf | seconds |
| 78 | s14_export.R | tsv + mgf → topblock.pdf | ~2 min |
| 79 | compose_a4 | topblock + profiles → evidence.pdf | ~10 min |
| 82 | extra_targets | master → extra targets | seconds |
| 84 | extra_inputs | extra rows → tsv + mgf | seconds |
| 81 | table1_pdf | Table1.xlsx + chains → T1.pdf | ~10 min |
| 85 | audit | PDFs + data → report | ~2 min |
| 86 | predict_missing | tsv → pkl update | ~1 min |

## Key decisions

1. **RT anchor**: Use search result's per-sample RT (not global average)
2. **Apex scoring**: Sequence ladder first (bLad+yLad), not intensity
3. **Background removal**: Three tiers (ion 10%, row 2×, bin 1.5×)
4. **Fallback chains**: v9-master for missing peptides; AlphaPeptDeep on-the-fly for missing predictions
5. **Audit**: Always run 85 after any change

## Common issues

| Problem | Fix |
|---------|-----|
| R fails on Chinese path | Copy risoPSV/ to ASCII path |
| R locale breaks UTF-8 | `env -u LC_CTYPE -u LC_ALL -u LANG Rscript` |
| PDF save Permission denied | Close the PDF viewer |
| Heatmap same size after change | Match figsize aspect to rect aspect exactly |
| Wide RT coverage in heatmap | Three-tier background removal (already built in) |

## Detailed documentation

- [SKILL.md](../SKILL.md) — Full MSHub skill card with decision rules and validation gates
- [docs/architecture.md](../docs/architecture.md) — Script chain diagram
- [docs/configuration.md](../docs/configuration.md) — All parameters
- [docs/troubleshooting.md](../docs/troubleshooting.md) — 15+ issues
