# Architecture

## Script chain

```
┌─────────────────────────────────────────────────────────────────┐
│  Phase 1: Target & Extraction                                    │
├─────────────────────────────────────────────────────────────────┤
│                                                                  │
│  70_build_s14_targets.py                                         │
│    search results xlsx/tsv → s14_targets.parquet                │
│    (seq/z/mz/sample/RT/qty/mod_state/pid)                        │
│         │                                                        │
│         ▼                                                        │
│  71_extract_targeted.py                                          │
│    targets + raw files → rows_all.parquet + spectra_all.npz     │
│    (RT-anchored: theory b/y → window → ±0.6min → ladder apex)   │
│         │                                                        │
│         ▼                                                        │
│  73_profiles.py                                                  │
│    targets + raw files → profiles.npz                            │
│    (full-RT XIC: ion floor 10% → cos×pres per scan)             │
│                                                                  │
├─────────────────────────────────────────────────────────────────┤
│  Phase 2: Evidence Rendering                                     │
├─────────────────────────────────────────────────────────────────┤
│                                                                  │
│  76_build_isopsv_inputs.py                                       │
│    rows + targets + pkl → in/{s14_spectra.tsv, best.mgf}        │
│         │                                                        │
│         ▼                                                        │
│  78_s14_export.R (isoPSV engine)                                │
│    tsv + mgf → out/s14_topblock.pdf                              │
│    (info lines / ion ladder / observed spectrum /               │
│     inverted prediction / ppm panel)                             │
│         │                                                        │
│         ▼                                                        │
│  79_compose_a4.py                                                │
│    topblock + profiles → S14_peptide_evidence_A4_isoPSV.pdf    │
│    (pymupdf: info text + R block + heatmap PNG)                 │
│                                                                  │
├─────────────────────────────────────────────────────────────────┤
│  Phase 3: Supplementary Chains                                   │
├─────────────────────────────────────────────────────────────────┤
│                                                                  │
│  82_extra_targets.py     (peptides missing from S14)            │
│    v9 master → s14_targets_extra.parquet                        │
│         │                                                        │
│         ▼  (TGT=... OUTD=... 71/73 with env overrides)          │
│                                                                  │
│  84_extra_inputs.py      (build inputs for extra chain)         │
│         │                                                        │
│         ▼                                                        │
│  86_predict_missing.py   (on-the-fly AlphaPeptDeep)             │
│         │                                                        │
│         ▼                                                        │
│  81_table1_pdf.py        (Table 1 cross-reference PDF)          │
│    Table1.xlsx + dual chains → Table1_peptide_evidence_A4.pdf  │
│                                                                  │
├─────────────────────────────────────────────────────────────────┤
│  Phase 4: Audit                                                  │
├─────────────────────────────────────────────────────────────────┤
│                                                                  │
│  85_audit.py                                                     │
│    data → presentation → audit_report.md/.json                  │
│    A: data reconciliation (S14→targets→extraction→info)         │
│    B: Table 1 crossref (page order, headers, provenance)        │
│    C: presentation (text, heatmap, overflow)                    │
│                                                                  │
└─────────────────────────────────────────────────────────────────┘
```

## Data flow

```
Search results (xlsx/tsv)          Raw files (.raw)         FASTA
        │                               │                     │
        ▼                               ▼                     ▼
  70: targets.parquet    71: rawdump (C#)          70: flanking
        │                       │                          │
        │         ┌─────────────┤                          │
        │         │             │                          │
        ▼         ▼             ▼                          ▼
  71: extract     73: profiles    76: isoPSV inputs    76: flanks in Info1
        │              │              │
        │              │              ▼
        │              │        78: R topblock.pdf
        │              │              │
        │              ▼              ▼
        │        79: heatmap PNG   (pymupdf compose)
        │              │              │
        └──────────────┴──────────────┘
                              │
                              ▼
                 Evidence PDF (A4, per-peptide)
```

## Heatmap background removal (three tiers)

```
Tier 1 (73, ion level):
  ion i at scan j counts only if int >= 10% × apex_int(i)
  → kills random ±10ppm first-hits in dense windows

Tier 2 (79, row level):
  row's elution-window peak < 2× out-of-window baseline → entire row blank
  → kills rows where fragments are buried in background

Tier 3 (79, bin level):
  bin raw value < 1.5× row baseline → blank
  → kills bins that survive the 20% relative cutoff in flat rows

Result: γ² → 20% row cutoff → 1.5× absolute cutoff
```
