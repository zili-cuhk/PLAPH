#!/bin/bash
# Usage: ./submit_array.sh EXPERIMENT N_TASKS MAX_CONCURRENT WALLTIME [OUT_SUBDIR] [EXTRA_MATLAB_ARGS]
# Optional: IDX="41-160,201-240" ./submit_array.sh ... submits only those task indices.
# Example: ./submit_array.sh main 520 200 10:00
#          ./submit_array.sh spline_pilot 12 12 10:00 pilot_order3 '"SplineOrder",3'
set -euo pipefail
source "$(dirname "$0")/config.sh"
EXP=$1; N=$2; CONC=$3; WALL=$4; SUB=${5:-}; EXTRA=${6:-}
SPEC=${IDX:-1-${N}}
OUT="$OUT_ROOT${SUB:+/$SUB}"
mkdir -p "$OUT/logs"
EXTRA_ARGS=""
if [ -n "$EXTRA" ]; then EXTRA_ARGS=", ${EXTRA//\"/\\\"}"; fi
bsub <<EOT
#!/bin/bash
#BSUB -J plaph_${EXP}${SUB:+_$SUB}[${SPEC}]%${CONC}
#BSUB -q ${QUEUE}
#BSUB -W ${WALL}
#BSUB -n 1
#BSUB -M 8G
#BSUB -R "rusage[mem=8] span[hosts=1]"
#BSUB -o ${OUT}/logs/${EXP}_%J_%I.out
#BSUB -e ${OUT}/logs/${EXP}_%J_%I.err
module load ${MATLAB_MODULE}
matlab -singleCompThread -batch "addpath('${CODE_DIR}'); plaph_v2_run_task('${EXP}', \${LSB_JOBINDEX}, '${OUT}'${EXTRA_ARGS})"
EOT
