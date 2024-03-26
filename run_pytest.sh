#!/bin/bash

OMP_NUM_THREADS=1 pytest --durations=0  -v --parallel=6
