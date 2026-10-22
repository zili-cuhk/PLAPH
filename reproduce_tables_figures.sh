#!/bin/bash
# Level A: regenerate every numeric table and every data figure of the article
# from the result files in results/, and check the tables against the published values
# in reference_tables/.
# Requires Python 3 with numpy, pandas and matplotlib. Takes a few minutes.
set -euo pipefail
cd "$(dirname "$0")"
OUT=${1:-output}
mkdir -p "$OUT"
python3 code/python/reproduce_tables.py results/simulation results/tcga results/external/gse72094_result.json \
  reference_tables/main_tables.tex reference_tables/supp_tables.tex "$OUT"
PLAPH_TCGA_DIR=results/tcga python3 code/python/make_figures.py results/simulation "$OUT/figures"
python3 code/python/kmic_geometry.py "$OUT/figures/main/kinked_mic_geometry.pdf"
python3 code/python/pilot_summary.py results/pilot > "$OUT/development_pilot.txt"
echo "Tables: $OUT/tables/   check: $OUT/table_check.txt"
echo "Figures: $OUT/figures/  development pilot numbers: $OUT/development_pilot.txt"
