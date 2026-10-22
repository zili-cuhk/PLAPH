#!/bin/bash
# Stage 2: formal runs with the spline order chosen in stage 1.
# Usage: ./02_main_campaign.sh 4     (or 3)
set -euo pipefail
cd "$(dirname "$0")"
ORDER=${1:?give the spline order chosen by plaph_v2_spline_decision}
X="'SplineOrder',$ORDER"
./submit_array.sh main            520 250 12:00 "" "$X"   # 13 scenarios x 40 chunks
./submit_array.sh ablation_tau0   360 150 12:00 "" "$X"   # 9 tau>0 scenarios x 40 chunks, working tau=0
./submit_array.sh grpsurv_export  520 100  4:00 "" "$X"   # data only, for grpsurv
./submit_array.sh highcorr         40  40 12:00 "" "$X"
./submit_array.sh copula_family   400 150 12:00 "" "$X"
echo "After 'main' finishes and 90_aggregate.sh has run: choose FixedCL/FixedKappa as the interior"
echo "grid point with the largest optimal_share for scenarios 5 and 8 in"
echo "summary/tuning_optimal_share_main.csv, then run ./06_paired_reruns.sh $ORDER <CL> <KAPPA>"
