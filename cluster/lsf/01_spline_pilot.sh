#!/bin/bash
# Stage 1: spline-order pilot on development seeds (disjoint from reported
# seeds): scenarios 1, 4, 6, R=20, cubic vs quadratic. Then run
#   matlab -batch "addpath('$CODE_DIR'); plaph_v2_spline_decision('$OUT_ROOT')"
set -euo pipefail
cd "$(dirname "$0")"
./submit_array.sh spline_pilot 12 12 10:00 pilot_order4 "'SplineOrder',4"
./submit_array.sh spline_pilot 12 12 10:00 pilot_order3 "'SplineOrder',3"
