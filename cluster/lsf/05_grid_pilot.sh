#!/bin/bash
# Confirmation pilot (development seeds): quadratic basis, widened EBIC grid
# c_L in {0.5,1,1.5,2,2.5} x kappa in {0.03,0.05,0.075,0.10,0.15}, budget 6000.
set -euo pipefail
cd "$(dirname "$0")"
./submit_array.sh spline_pilot 12 12 4:00 pilot_grid2 "'SplineOrder',3"
