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
preprocessing described in Supplementary Section S5.1 (median imputation and standardization,
variance filter to 2,000 RNA and 1,000 CNA variables, margin-specific Cox screens adjusted
for age, sex and stage). The script that produced the working set from the cBioPortal
download is not included: the original preprocessing script was not retained. The working set
itself, which is derived from publicly available TCGA data, can be provided to the editor and
reviewers on request, and all analyses downstream of it (`04_tcga.sh` and the scripts in
`code/`) are included.

## GSE72094

NCBI GEO series GSE72094 (platform GPL15048), public:

```
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE72nnn/GSE72094/matrix/GSE72094_series_matrix.txt.gz
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL15048&targ=self&form=text&view=data
```

`code/python/extract_gse72094.py` extracts the clinical fields and the probe set of CYP17A1;
`code/python/external_gse72094.py` computes the prespecified external evaluation.
