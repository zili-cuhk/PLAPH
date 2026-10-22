# Reproducibility materials

**High-Dimensional Partially Linear Additive Cox Models with a Working Copula for Dependent Censoring (PLAPH)**

This archive contains the code, the numerical results and the scripts that produce every
numeric table and every data figure of the main article and the Supplementary Material.
Reproduction is possible at three levels:

| Level | What is reproduced | Needs | Time |
|---|---|---|---|
| **A** | All tables and data figures from the stored results, with an automatic check against the reported values | Python 3 (numpy, pandas, matplotlib) | minutes |
| **B** | The stored simulation results themselves (all experiments, seeds fixed) | MATLAB, R, an LSF cluster or a long local run | about 1,500 core-hours |
| **C** | The TCGA-LUAD analysis and the external evaluation | Level B software, data downloaded from cBioPortal and GEO | hours |

## Level A: tables and figures from the stored results

```bash
./reproduce_tables_figures.sh output
```

This regenerates the 22 numeric tables (written to `output/tables/<label>.tex`) and the nine
data figures (`output/figures/`), recomputes the development-pilot numbers quoted in the
Supplement (`output/development_pilot.txt`), and compares every numeric table cell with the
values reported in the article, stored in `reference_tables/` (the 22 numeric tables as
LaTeX, without captions or notes). The comparison is written to `output/table_check.txt`;
for the submitted version it reports `22 tables, 5468 numeric cells compared, 0 discrepancies`.
The script reads `reference_tables/` only to obtain the row labels of each table; all values
are blanked and recomputed from `results/`.

Software used for Level A: Python 3.13, numpy 2.5, pandas 3.0, matplotlib 3.11 (any recent
versions should work).

## Level B: rerunning the simulations

All estimators are implemented in MATLAB (`code/matlab/`); the published comparator
grpsurv is run in R (`code/R/grpsurv_comparator.R`). The runs reported in the article used
MATLAB R2026a (26.1.0.3251617, Update 2) with the Statistics and Machine Learning, Optimization and Curve Fitting
Toolboxes, and R 4.5.1 with grpreg 3.6.0, on an LSF cluster (one core and 8 GB per task).

**Entry points.**

* `plaph_v2_manifest(experiment)` is the frozen table of designs, seeds and methods for each
  experiment (`main`, `ablation_tau0`, `grpsurv_export`, `highcorr`, `copula_family`,
  `xi_sensitivity`, `spline_pilot`, `smoke`).
  The seed of replication `r` is `seed_root + 1000*seed_scenario + r`; designs that must share
  data (the working families of one generating family; all kink coefficients at one sample
  size; PLAPH and PLAPH-ind) share `seed_scenario`.
* `plaph_v2_run_task(experiment, k, out_root, ...)` runs the `k`-th chunk of five
  replications of an experiment and writes `scenario_001_replicates.csv`,
  `scenario_001_tuning_path.csv` and a `.mat` file with the fitted curves. A chunk whose
  results exist is skipped, so interrupted runs can be resubmitted.
  `plaph_v2_task_count(experiment)` gives the number of chunks.
* `plaph_run_one_scenario_impl.m` is the estimation engine (likelihood, screening, adaptive
  pilot, K-MIC and group hard-thresholding updates, stopping gate, EBIC grid with the tie
  rule, joint refit, model-based inference, and the aligned PL comparator).

**Small local check (about 15 minutes).**

```matlab
addpath('code/matlab')
plaph_v2_smoke('smoke_out')            % 2 replications of design D5 with PLAPH and PL, prints PASS/FAIL lines
```

**Full campaign on an LSF cluster.** Edit `cluster/lsf/config.sh` (paths, modules, queue),
then submit in this order; each script prints the job identifiers.

| Step | Script | Content |
|---|---|---|
| 0 | `00_smoke.sh` | smoke test |
| 1 | `01_spline_pilot.sh`, `05_grid_pilot.sh` | development pilot (basis order, grid confirmation; development seeds) |
| 2 | `02_main_campaign.sh 3` | D1–D13 (PLAPH, PL), PLAPH-ind ablation, data export for grpsurv, high correlation, copula families |
| 3 | `06_paired_reruns.sh 3 1.0 0.10` | paired copula-family refits and the kink-coefficient experiment |
| 4 | `07_highcorr_extras.sh 3` | PLAPH-ind and grpsurv data for the high-correlation design |
| 5 | `03_grpsurv.sh` and `03_grpsurv.sh "V2_grpsurv_export_S21_*"` | grpsurv (R) on the exported data |
| 6 | `90_aggregate.sh`, `91_curves.sh` | summaries (`summary/*.csv`) and function-error summaries from the stored curves |

`02_main_campaign.sh` also submits the copula-family experiment in full; step 3 then refits
the misspecified families on the shared data sets, as reported. After step 6, copy
`OUT_ROOT/{main,ablation_tau0,highcorr,copula_family,xi_sensitivity,grpsurv_results,summary}`
to `results/simulation/` and run Level A.

**Without a cluster.** Every chunk is an ordinary MATLAB call, e.g.
`for k = 1:plaph_v2_task_count("main"), plaph_v2_run_task("main", k, "out"); end`.
A chunk of five replications takes about 20 minutes to 6 hours on one core, depending on the
design. The complete set of simulation runs reported here (including runs later superseded by the
paired copula-family refits) used about 1,500 core-hours.

**Numerical notes.** A fit is *accepted* only if the selection gate and the refit KKT residual
are at most 1e-4; acceptance rates are reported and all summaries use all replications.
Supports are taken from the frozen selection masks before refitting.

## Level C: TCGA-LUAD and external evaluation

The patient-level TCGA working matrices are not distributed. `data/README_data.md` describes
how to obtain the data and the format of the working set expected by
`run_PLAPH_real_data.m`. With the working set in `WS_DIR`:

| Step | Script | Content |
|---|---|---|
| A | `04_tcga.sh A 3` | 20 fixed-association EBIC grids (Frank, Gaussian, Clayton, Gumbel) |
| B | `04_tcga.sh B 3` | 4 profile fits and 200 stratified 80% subsamples (after A) |
| C | `04_tcga.sh C 3` | additive-threshold path, weight exponent and screened-dimension sensitivity (after A) |

The stored TCGA outputs in `results/tcga/` contain gene-level coefficients and fit summaries
only. The external evaluation in GSE72094 is reproduced with

```bash
python3 code/python/extract_gse72094.py GSE72094_series_matrix.txt.gz GPL15048_table.txt gse.csv
python3 code/python/external_gse72094.py gse.csv external.json
Rscript code/R/external_gse72094_check.R gse.csv     # optional cross-check with survival::coxph
```

## Contents

```
code/matlab     estimation engine, experiment manifest and task runner, aggregation
code/R          grpsurv comparator, coxph cross-check of the external evaluation
code/python     table regeneration and check, figures, pilot summary, external evaluation
cluster/lsf     LSF submission scripts (edit config.sh)
results/        simulation, pilot, TCGA and external results used in the article
reference_tables/ the numeric tables of the article, without captions (read by the table check)
data/           data-access notes
MANIFEST.md     table and figure -> script and input map
LICENSE         MIT License
```

Code is released under the MIT License (`LICENSE`); the stored results are provided for verification of the article.

Contact: Zili Liu (ZLiu20@mdanderson.org).
