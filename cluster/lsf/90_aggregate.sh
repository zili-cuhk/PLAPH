#!/bin/bash
set -euo pipefail
source "$(dirname "$0")/config.sh"
module load ${MATLAB_MODULE}
matlab -singleCompThread -batch "addpath('${CODE_DIR}'); plaph_v2_aggregate('${OUT_ROOT}'); plaph_v2_tcga_summary('${OUT_ROOT}/tcga')"
echo "Summaries in $OUT_ROOT/summary and $OUT_ROOT/tcga/summary"
