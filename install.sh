#!/bin/bash
# DIA Spectrum Evidence Extraction — one-command installer
set -e
cd "$(dirname "$0")"
echo "=== DIA Evidence Extraction Pipeline ==="
echo "Checking dependencies..."
python -c "import numpy, pandas, matplotlib, pymupdf, openpyxl, h5py" 2>/dev/null \
  && echo "  [OK] Python packages" \
  || { echo "  [!] Install: pip install numpy pandas matplotlib pymupdf openpyxl h5py"; exit 1; }
Rscript -e "cat('R OK\n')" 2>/dev/null \
  && echo "  [OK] R" \
  || echo "  [!] R not found (required for isoPSV rendering)"
[ -f "risoPSV/78_s14_export.R" ] \
  && echo "  [OK] isoPSV R engine" \
  || echo "  [!] risoPSV/ incomplete"
[ -f "tools/rawdump.exe" ] \
  && echo "  [OK] rawdump" \
  || echo "  [!] tools/rawdump.exe missing"
echo "=== Installation complete ==="
echo "Quick start: python scripts/70_build_s14_targets.py"
