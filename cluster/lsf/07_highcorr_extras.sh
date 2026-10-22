#!/bin/bash
# PLAPH-ind and grpsurv data for the high-correlation design (S21), appended
# as the last rows of ablation_tau0 and grpsurv_export (2026-10-22).
# Usage: ./07_highcorr_extras.sh ORDER
# When the export has finished:  ./03_grpsurv.sh "V2_grpsurv_export_S21_*"
set -euo pipefail
cd "$(dirname "$0")"
ORDER=${1:?spline order}
IDX="361-400" ./submit_array.sh ablation_tau0  400 40 12:00 "" "'SplineOrder',$ORDER"
IDX="521-560" ./submit_array.sh grpsurv_export 560 40  4:00 "" "'SplineOrder',$ORDER"
