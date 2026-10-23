# Data access

## Simulations

All simulated data are generated inside `plaph_run_one_scenario_impl.m` from the seeds in
`plaph_v2_manifest.m`; no data files are needed. The grpsurv comparator reads data sets that
the MATLAB engine exports (`grpsurv_export` experiment) and that are therefore identical to
those analyzed by PLAPH.

## TCGA-LUAD

Source: cBioPortal study `luad_tcga_pan_can_atlas_2018` (TCGA PanCancer Atlas), with overall
survival from the TCGA Pan-Cancer Clinical Data Resource. Access is public at
https://www.cbioportal.org/study/summary?id=luad_tcga_pan_can_atlas_2018.

The analysis code (`run_PLAPH_real_data.m`) reads a working set from `WorkingSetDir` with the
files

```
primary_p500_q200_clinical.csv            patient_id, age, sex, stage, os_months, os_event, ...
primary_p500_q200_Z_RNA_p500.csv          500 screened RNA z-scores (patients x features)
primary_p500_q200_Z_RNA_p500_features.csv feature names, in screening order
primary_p500_q200_X_CNA_q200.csv          200 screened continuous log2 CNA values
primary_p500_q200_X_CNA_q200_features.csv feature names, in screening order
```

These patient-level matrices are not distributed with this archive. They were built by the
preprocessing described in Supplementary Section S5.1 in two stages.

1. *Extraction and variance filter.* RNA z-scores and continuous log2 CNA values from the
   cBioPortal download were linked by patient identifier, feature-level missing values were
   replaced by the feature median, features were standardized, and a variance filter kept
   2,000 RNA and 1,000 CNA variables (`Z_RNA_var_top2000_os.csv`, `X_log2CNA_var_top1000_os.csv`,
   `luad_pilot_clinical.csv`). The script for this stage was not retained.
2. *Screen and export* (`code/R/tcga_preprocessing/`, R with `survival` and `splines`).
   `run_luad_pilot_dimension_experiment.R` runs, on the patients with positive follow-up and
   complete age, sex and stage, death- and censoring-margin Cox screens adjusted for age, sex
   and stage: a Wald test for each RNA variable and a likelihood-ratio test of a three-degree-
   of-freedom natural spline for each CNA variable. Each variable is ranked by the larger of its
   two -log10 p-values. `export_luad_pilot_working_sets.R` then writes the leading 500 RNA and
   200 CNA variables as the working set `primary_p500_q200_*` (it also writes two sensitivity
   sets that are not used in the article). Set `PLAPH_LUAD_ROOT` to the folder holding the
   stage-1 files.

`code/python/check_tcga_screen.py` re-runs the stage-2 screens on a working set and compares
them with the screen statistics stored in its feature files. For the working set used in the
article it reproduces all 2,000 stored RNA p-values and coefficients to a relative difference
below 1e-6, all 800 stored CNA p-values to below 1e-11, and the stored ranking exactly
(482 patients, 174 deaths).

The working set, which is derived from publicly available TCGA data, can be provided to the
editor and reviewers on request; all analyses downstream of it (`04_tcga.sh` and the scripts in
`code/`) are included.

## GSE72094

NCBI GEO series GSE72094 (platform GPL15048), public:

```
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE72nnn/GSE72094/matrix/GSE72094_series_matrix.txt.gz
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL15048&targ=self&form=text&view=data
```

`code/python/extract_gse72094.py` extracts the clinical fields and the probe set of CYP17A1;
`code/python/external_gse72094.py` computes the prespecified external evaluation.
