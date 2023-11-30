#!/bin/bash

export ESPRESSO_COMMAND="mpirun -np 2 pw.x"
export ESPRESSO_PSEUDO="$(pwd)"
export ESPRESSO_TMPDIR="tmpdir"

python relax.py
