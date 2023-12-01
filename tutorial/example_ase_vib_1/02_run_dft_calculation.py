from matscitoolkit.ase_vib_tools import DFTrunner
import os

import argparse
argparser = argparse.ArgumentParser()
argparser.add_argument("-i", "--id", type=int)
args = argparser.parse_args()

"""
In this step, we run the DFT calculation. 
This is the most time-consuming step.
However, we can reduce the time by running the calculation in parallel. 
In this case, the parallelization is based on displaced images. 
The parallelization is controlled by the "id" argument which should be a 
mapping of JOB ARRAY ID and displaced image ID.

NOTE: computation of emaxpos can be handled automatically by the code.
"""




input_data = {
    "control": {},
    "system": {
        "ecutwfc": 60,
        "ecutrho": 480,
        "occupations": "smearing",
        "smearing": "gaussian",
        "degauss": 0.01,
        "input_dft": "vdw-df2-b86r",
    },
    "electrons": {
        "electron_maxstep": 300,
        "conv_thr": 1.0e-09,
        "mixing_beta": 0.7,
    },
}

pseudopotentials = {"B": "b_pbe_v1.4.uspp.F.UPF", "H": "h_pbe_v1.4.uspp.F.UPF"}

x = DFTrunner(
    system_id=args.id,
    input_data=input_data,
    pseudopotentials=pseudopotentials,
    kpts=[12, 12, 1],
    dirname="dft_calc/SUFFIX",
    espresso_command=os.environ["ESPRESSO_COMMAND"].split(),
)
x.run()
