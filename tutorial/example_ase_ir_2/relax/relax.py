from ase import Atoms
from ase.optimize import BFGS
from ase.calculators.espresso import Espresso, EspressoProfile
from ase.io import read, write
import os

atoms_obj = read("water.vasp")
input_data = {
    "control": {
        "tprnfor": True,
    },
    "system": {
        "ecutwfc": 60,
        "ecutrho": 480,
        "occupations": "smearing",
        "smearing": "gaussian",
        "degauss": 0.01,
    },
    "electrons": {
        "electron_maxstep": 100,
        "conv_thr": 1.0e-07,
        "mixing_beta": 0.7,
    },
}

atoms_obj.calc = Espresso(
    input_data=input_data,
    pseudopotentials={"O": "o_pbe_v1.2.uspp.F.UPF", "H": "h_pbe_v1.4.uspp.F.UPF"},
    kpts=[1, 1, 1],
    directory="dft_calc",
    profile=EspressoProfile(os.environ["ESPRESSO_COMMAND"].split())
)

import numpy as np

dyn = BFGS(atoms_obj, trajectory='H2O.traj', logfile='H2O.log')
dyn.run(fmax=0.0005)

atoms_obj.center()
write("rlx_H2O.vasp", atoms_obj)
