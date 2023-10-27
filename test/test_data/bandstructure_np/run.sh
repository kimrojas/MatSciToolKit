#!/bin/bash

NPROC=1
mpirun -np $NPROC pw.x -in espresso_pwscf.inp > espresso_pwscf.out
mpirun -np $NPROC pw.x -in espresso_pwbands.inp > espresso_pwbands.out
mpirun -np $NPROC bands.x -in espresso_bands.inp > espresso_bands.out