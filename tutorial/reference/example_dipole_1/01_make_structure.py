from ase.io import read, write
from ase.build import molecule
from ase.calculators.espresso import Espresso, EspressoProfile
from ase.optimize import LBFGS
import numpy as np


atoms = molecule("H2O", cell=10 * np.eye(3))
atoms.center()

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
        "conv_thr": 1.0e-06,
    },
}

pseudopotentials = {"O": "o_pbe_v1.2.uspp.F.UPF", "H": "h_pbe_v1.4.uspp.F.UPF"}

calc = Espresso(
    input_data=input_data,
    kpts=(1, 1, 1),
    pseudopotentials=pseudopotentials,
    profile=EspressoProfile("mpirun pw.x".split()),
    directory="01_dftcalc",    
)

atoms.set_calculator(calc)

LBFGS(atoms, trajectory="01_lbfgs.traj", logfile="01_lbfgs.log").run(fmax=0.01)

atoms.center()
write("h2o.vasp", atoms)
