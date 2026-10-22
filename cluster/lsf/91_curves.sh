#!/bin/bash
# Function-estimation summaries (ISB^2, IVar, IMSE; mean curves and bands)
# from the curves stored in the simulation .mat files. Run after the
# simulation experiments have finished; a few minutes.
set -euo pipefail
source "$(dirname "$0")/config.sh"
EXPS=${1:-"'main','ablation_tau0','highcorr','copula_family'"}
bsub <<EOT
#!/bin/bash
#BSUB -J plaph_curves
#BSUB -q ${QUEUE}
#BSUB -W 4:00
#BSUB -n 1
#BSUB -M 16G
#BSUB -R "rusage[mem=16] span[hosts=1]"
#BSUB -o ${OUT_ROOT}/logs/curves_%J.out
#BSUB -e ${OUT_ROOT}/logs/curves_%J.err
module load ${MATLAB_MODULE}
matlab -singleCompThread -batch "addpath('${CODE_DIR}'); plaph_v2_curve_summary('${OUT_ROOT}', 'Experiments', {${EXPS}})"
EOT
