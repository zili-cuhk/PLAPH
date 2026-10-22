#!/bin/bash
# Stage 0 (about 30-60 min): 2 replications of scenario 5 with the full EBIC
# grid for PLAPH and PL, data export + grpsurv on the same 2 data sets, and
# one TCGA fit. Inspect the logs and summary before any large submission.
set -euo pipefail
source "$(dirname "$0")/config.sh"
OUT="$OUT_ROOT/smoke"; mkdir -p "$OUT/logs"
bsub <<EOT
#!/bin/bash
#BSUB -J plaph_smoke
#BSUB -q ${QUEUE}
#BSUB -W 6:00
#BSUB -n 1
#BSUB -M 8G
#BSUB -R "rusage[mem=8] span[hosts=1]"
#BSUB -o ${OUT}/logs/smoke_%J.out
#BSUB -e ${OUT}/logs/smoke_%J.err
module load ${MATLAB_MODULE}
matlab -singleCompThread -batch "addpath('${CODE_DIR}'); plaph_v2_smoke('${OUT}', '${WS_DIR}')"
module load ${R_MODULE}
D=\$(ls -d ${OUT}/grpsurv_data/V2_smoke_S*)
Rscript ${CODE_DIR}/grpsurv_comparator.R \$D 1 2 ${OUT}/grpsurv_results/smoke.csv
EOT
