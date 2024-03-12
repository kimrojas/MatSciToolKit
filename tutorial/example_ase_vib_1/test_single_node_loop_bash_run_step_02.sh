#!/bin/bash

export ESPRESSO_COMMAND="mpirun -np 4 pw.x"
export ESPRESSO_PSEUDO="$(pwd)"
export ESPRESSO_TMPDIR="tmpdir"

export OMP_NUM_THREADS=1

for i in {1..25}
do
  python 02_run_dft_calculation.py --id $i
done