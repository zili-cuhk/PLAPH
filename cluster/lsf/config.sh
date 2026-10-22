# PLAPH v2 campaign configuration for an LSF cluster. Edit the paths once.
export CODE_DIR=/path/to/PLAPH_reproducibility/code/matlab
export OUT_ROOT=/path/to/plaph_v2_runs            # all experiment outputs are written here
export WS_DIR=/path/to/tcga_working_set           # TCGA working set (not distributed; see data/README_data.md)
export MATLAB_MODULE=matlab                       # name of the MATLAB environment module
export R_MODULE=R/4.5.1                           # R with grpreg (3.6.0 used) and survival
export QUEUE=medium                               # queue used here required walltime > 3 h, limit 24 h
mkdir -p "$OUT_ROOT/logs"
