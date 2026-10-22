#!/bin/bash
# Stage 3: grpsurv on exported data (after grpsurv_export has finished).
# One task per (scenario, block of 20 data sets): 13 x 10 = 130 tasks.
# First time only: CRAN is unreachable from seadragon; download grpreg_<ver>.tar.gz on another machine, then
#   module load R/4.5.1; R CMD INSTALL -l "$(Rscript -e 'cat(Sys.getenv("R_LIBS_USER"))')" grpreg_<ver>.tar.gz   (3.6.0 used)
set -euo pipefail
source "$(dirname "$0")/config.sh"
DATA="$OUT_ROOT/grpsurv_data"; RES="$OUT_ROOT/grpsurv_results"
mkdir -p "$RES" "$OUT_ROOT/logs"
PATTERN=${1:-V2_grpsurv_export_S*}   # e.g. ./03_grpsurv.sh "V2_grpsurv_export_S21_*"
LIST="$RES/scenario_dirs_$(date +%Y%m%d%H%M%S).txt"   # one list per submission
ls -d "$DATA"/$PATTERN > "$LIST"
NS=$(wc -l < "$LIST"); N=$((NS*10))
bsub <<EOT
#!/bin/bash
#BSUB -J plaph_grpsurv[1-${N}]%100
#BSUB -q ${QUEUE}
#BSUB -W 12:00
#BSUB -n 1
#BSUB -M 8G
#BSUB -R "rusage[mem=8] span[hosts=1]"
#BSUB -o ${OUT_ROOT}/logs/grpsurv_%J_%I.out
#BSUB -e ${OUT_ROOT}/logs/grpsurv_%J_%I.err
module load ${R_MODULE}
I=\${LSB_JOBINDEX}
S=\$(( (I-1)/10 + 1 )); B=\$(( (I-1)%10 ))
DIR=\$(sed -n "\${S}p" ${LIST})
NAME=\$(basename \$DIR)
FIRST=\$(( B*20 + 1 )); LAST=\$(( B*20 + 20 ))
Rscript ${CODE_DIR}/grpsurv_comparator.R \$DIR \$FIRST \$LAST ${RES}/\$NAME/rows_\${FIRST}_\${LAST}.csv
EOT
