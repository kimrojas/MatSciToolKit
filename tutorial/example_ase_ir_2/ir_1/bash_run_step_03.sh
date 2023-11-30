#!/bin/bash

export ESPRESSO_COMMAND="mpirun -np 4 pw.x"
export ESPRESSO_PSEUDO="$(pwd)"
export ESPRESSO_TMPDIR="tmpdir"

for i in {1..19}
do
  python 03_run_dft_calculation.py --id $i
done