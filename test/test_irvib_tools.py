from ase.io import read, write
from ase.build import molecule
from pathlib import Path
from matscitoolkit import ase_vib_tools, ase_ir_tools
import os
import shutil
import pytest
from ase.calculators.espresso import Espresso, EspressoProfile

from temporary_handler import TemporaryEnvironment


@pytest.mark.parametrize("name,expected", [("H2O", 19), ("C2H6", 49)])
@pytest.mark.parametrize("ext", ["vasp", "xyz", "traj"])
@pytest.mark.parametrize("tool", [ase_vib_tools, ase_ir_tools])
def test_generateDisplacedStructures(name, expected, tool, ext):
    tmp = TemporaryEnvironment(name, expected, tool.__name__, ext)

    mol = molecule(name, cell=[5, 5, 5])
    mol.center()
    write(f"{name}.{ext}", mol)

    kwargs = {
        "filepath": f"{name}.{ext}",
        "disp_dir": "displacement_dir",
        "filetype": ext,
    }

    tool.generateDisplacedStructures(**kwargs)

    nfiles = len(list(Path("displacement_dir").glob(f"*.{ext}")))

    # -- Return to main
    tmp.return_to_main()

    # -- Evaluation
    assert nfiles == expected, "Number of files is not correct"

    # -- Clean up
    tmp.clean()


def test_dipolemoment_calculation():
    pseudo_dir = str(Path("test/test_data_irvib/espresso_pseudopotential").absolute())
    tmp = TemporaryEnvironment("dipolemoment")
    
    atoms = molecule("H2O", cell=[10, 10, 10])
    atoms.center()
    input_data = {
        "control": {
            "tefield": True,
            "dipfield": True,
            "pseudo_dir": pseudo_dir,
        },
        "system": {
            "occupations": "smearing",
            "smearing": "fermi-dirac",
            "degauss": 0.02,
            "edir": 3,
            "eamp": 0.00,
            "eopreg": 0.0001,
            "emaxpos": 0.0001,
            "ecutwfc": 30,
            "ecutrho": 240,
        },
    }
    pseudopotential = {"H": "h_pbe_v1.4.uspp.F.UPF", "O": "o_pbe_v1.2.uspp.F.UPF"}
    profile = EspressoProfile("mpirun -np 1 pw.x".split())
    atoms.calc = Espresso(
        input_data=input_data,
        pseudopotentials=pseudopotential,
        kpts=(1, 1, 1),
        profile=profile,
    )

    atoms.get_potential_energy()
    dipol_arr = atoms.get_dipole_moment().tolist()
    # [-0.0, -0.0, -0.36960743277046937]
    expected_dipol_arr = [0, 0, -0.36991972]
    with open("dipolemoment.txt", "w") as f:
        f.write(str(dipol_arr))
    assert dipol_arr == pytest.approx(expected_dipol_arr, abs=0.02)

    tmp.return_to_main()
    tmp.clean()
