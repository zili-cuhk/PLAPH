#!/bin/bash
# 2026-11-17. Three steps; run them in order.
#   ./10_numerics_and_new_runs.sh smoke    one job (about 1-3 h): checks that LegacyNumerics=true
#                                          reproduces the stored main D6 replication 1, and runs the
#                                          corrected numerics on the same data. Send the .out file.
#   ./10_numerics_and_new_runs.sh tau      tau_misspec: D5 data fitted at tau 0.2 and 0.8, D6 data
#                                          fitted at tau 0.5; 3 x 200 replications, 120 tasks.
#   ./10_numerics_and_new_runs.sh tcga     TCGA: the eight working copulas whose EBIC minimum was at
#                                          the grid edge c_L = 6, refitted on c_L in {6,7,8,10,12,15,20}.
set -euo pipefail
cd "$(dirname "$0")"
source ./config.sh
MODE=${1:?give smoke, tau or tcga}
case "$MODE" in
smoke)
  bsub <<EOT
#!/bin/bash
#BSUB -J plaph_smoke_numerics
#BSUB -q ${QUEUE}
#BSUB -W 8:00
#BSUB -n 1
#BSUB -M 8G
#BSUB -R "rusage[mem=8] span[hosts=1]"
#BSUB -o ${OUT_ROOT}/logs/smoke_numerics_%J.out
#BSUB -e ${OUT_ROOT}/logs/smoke_numerics_%J.err
module load ${MATLAB_MODULE}
matlab -singleCompThread -batch "addpath('${CODE_DIR}'); plaph_v2_smoke_numerics('${OUT_ROOT}', '${OUT_ROOT}/smoke_numerics')"
EOT
  ;;
tau)
  ./submit_array.sh tau_misspec 120 120 12:00 "" "'SplineOrder',3"
  ;;
tcga)
  OUT="$OUT_ROOT/tcga_gridext"; mkdir -p "$OUT/logs"
  bsub <<EOT
#!/bin/bash
#BSUB -J plaph_tcga_gridext[1-8]
#BSUB -q ${QUEUE}
#BSUB -W 20:00
#BSUB -n 1
#BSUB -M 8G
#BSUB -R "rusage[mem=8] span[hosts=1]"
#BSUB -o ${OUT}/logs/gridext_%J_%I.out
#BSUB -e ${OUT}/logs/gridext_%J_%I.err
module load ${MATLAB_MODULE}
matlab -singleCompThread -batch "addpath('${CODE_DIR}'); t=[\"grid_Clayton_0.30\";\"grid_Frank_0.30\";\"grid_Frank_0.50\";\"grid_Frank_-0.70\";\"grid_Gaussian_0.50\";\"grid_Gaussian_-0.50\";\"grid_Gaussian_-0.70\";\"grid_Gumbel_0.50\"]; plaph_v2_tcga_task(t(\${LSB_JOBINDEX}), '${WS_DIR}', '${OUT}', 'SplineOrder', 3, 'GridCL', [6 7 8 10 12 15 20])"
EOT
  ;;
*) echo "MODE must be smoke, tau or tcga"; exit 1 ;;
esac
