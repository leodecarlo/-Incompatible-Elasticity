#!/bin/bash
export OMP_NUM_THREADS=1
## This assumes that the Firedrake venv is already activated

## Example 11.1 ## TODO repeat now with this larger scale_displacement
python3 hct_traction.py --problem "traction" --clscale 0.04 --n_iterates 20 --scale_displacement 10.0 -options_left 0

## Example 11.2 and its h-refinement. n_iterates needs to be only around 10 for it to converge
python3 hct_traction.py --problem "necking" --clscale 0.018 --n_iterates 50 --scale_displacement 5 -options_left 0
python3 hct_traction.py --problem "necking" --clscale 0.0132 --n_iterates 50 --scale_displacement 5 -options_left 0

## Example 11.3
python3 hct_traction.py --problem "inclusion" --clscale 0.018 --n_iterates 10 -options_left 0
python3 hct_traction.py --problem "inclusion" --clscale 0.1 --n_iterates 10 -options_left 0


