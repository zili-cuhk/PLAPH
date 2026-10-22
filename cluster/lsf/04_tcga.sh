#!/bin/bash
# Stage 4: TCGA-LUAD. Stage A = 20 fixed-tau EBIC grids; Stage B (after A)
# = 4 profile fits + 200 subsamples; Stage C (after A) = additive-threshold path,
# weight exponent and screened-dimension sensitivity at tau=0.
# Optional IDX="8-9" resubmits only those task indices.
# Usage: ./04_tcga.sh A ORDER  |  ./04_tcga.sh B ORDER  |  ./04_tcga.sh C ORDER
set -euo pipefail
source "$(dirname "$0")/config.sh"
STAGE=${1:?A or B}; ORDER=${2:?spline order}
OUT="$OUT_ROOT/tcga"; mkdir -p "$OUT/logs"
if [ "$STAGE" = "A" ]; then N=20; S=1; WALL=20:00; elif [ "$STAGE" = "B" ]; then N=204; S=2; WALL=6:00; else N=9; S=3; WALL=20:00; fi
bsub <<EOT
#!/bin/bash
#BSUB -J plaph_tcga_${STAGE}[${IDX:-1-${N}}]%100
#BSUB -q ${QUEUE}
#BSUB -W ${WALL}
#BSUB -n 1
#BSUB -M 8G
#BSUB -R "rusage[mem=8] span[hosts=1]"
#BSUB -o ${OUT}/logs/tcga_${STAGE}_%J_%I.out
#BSUB -e ${OUT}/logs/tcga_${STAGE}_%J_%I.err
module load ${MATLAB_MODULE}
matlab -singleCompThread -batch "addpath('${CODE_DIR}'); t=plaph_v2_tcga_tasks(${S}); plaph_v2_tcga_task(t(\${LSB_JOBINDEX}), '${WS_DIR}', '${OUT}', 'SplineOrder', ${ORDER})"
EOT
