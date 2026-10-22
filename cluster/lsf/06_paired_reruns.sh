#!/bin/bash
# Paired-design runs (2026-10-22). The first copula_family run gave every
# working family its own data sets; F2-F4, F6, F8 and F10 are now refitted on
# the data sets of F1, F5, F7 and F9 (whose results are unchanged and kept).
# The xi experiment refits the same data sets at every xi.
# Usage: ./06_paired_reruns.sh ORDER CL KAPPA     e.g.  ./06_paired_reruns.sh 3 1.0 0.10
set -euo pipefail
cd "$(dirname "$0")"
source config.sh
ORDER=${1:?spline order}; CL=${2:?FixedCL for xi}; KAPPA=${3:?FixedKappa for xi}
# keep the unpaired fits for the record (a finished chunk would otherwise be skipped)
mkdir -p "$OUT_ROOT/copula_family_unpaired"
for s in 32 33 34 36 38 40; do
  for d in "$OUT_ROOT"/copula_family/V2_copula_family_S${s}_*; do
    [ -d "$d" ] && mv "$d" "$OUT_ROOT/copula_family_unpaired/"
  done
done
if [ -d "$OUT_ROOT/xi_sensitivity" ]; then mv "$OUT_ROOT/xi_sensitivity" "$OUT_ROOT/xi_sensitivity_unpaired_$(date +%Y%m%d%H%M)"; fi
# rows 2-4, 6, 8, 10 of copula_family, 40 chunks each
IDX="41-160,201-240,281-320,361-400" ./submit_array.sh copula_family 400 150 12:00 "" "'SplineOrder',$ORDER"
# 2 sample sizes x 10 xi values x 40 chunks
./submit_array.sh xi_sensitivity 800 150 10:00 "" "'SplineOrder',$ORDER,'FixedCL',$CL,'FixedKappa',$KAPPA"
