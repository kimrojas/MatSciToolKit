from ase.io import read, write
from ase.build import molecule
from ase.calculators.espresso import Espresso, EspressoProfile
from ase.optimize import LBFGS
import numpy as np
import matplotlib.pyplot as plt
from ase.visualize.plot import plot_atoms
from ase.units import create_units
units = create_units('2006')

atoms = read("h2o.vasp")


input_data = {
    "control": {
        "tprnfor": True,
        "tefield": True,
        "dipfield": True,
    },
    "system": {
        "ecutwfc": 60,
        "ecutrho": 480,
        "occupations": "smearing",
        "smearing": "gaussian",
        "degauss": 0.01,
        "edir": 3, # Specific for dipole moment parameters
        "eamp": 0.00, # Specific for dipole moment parameters
        "eopreg": 0.0001, # Specific for dipole moment parameters
        "emaxpos": 0.0001 # Specific for dipole moment parameters
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
    directory="02_dftcalc",
)

atoms.set_calculator(calc)

# Get Potential energy
e = atoms.get_potential_energy()
print(f"Potential energy (eV): {e:.6f}")

# Get dipole moment
dm = atoms.get_dipole_moment()
dm = [_ for _ in dm]
print(f"Dipole moment (ASE unit): {dm}", )
dm_debye =  [_ / units['Debye'] for _ in dm]
print(f"Dipole moment (Debye): {dm_debye}")


# Visualization
atoms = read("h2o.vasp")

fig, ax = plt.subplots()
plot_atoms(atoms, ax, radii=0.9, rotation=("-90x,90y,0z"))

ax.arrow(5.25, 7, 0, -3, head_width=0.5, head_length=0.1, fc='black', ec='black', ls="-")
ax.text(5.25, 3, f"Dipole moment (Debye): {dm_debye[2]}", fontsize=12, ha='center', va='center')
ax.set_ylabel("z-axis")
ax.set_xlabel("y-axis")

fig.savefig("02_dipolemoment.png", dpi=300)