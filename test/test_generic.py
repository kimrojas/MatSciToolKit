import pytest
from matscitoolkit.analysis_workflow.ABC_workflow import WorkflowBaseClass
from matscitoolkit.analysis_workflow.espresso import VibrationalAnalysisWorkflow
from ase.io import read, write
from ase.build import molecule
from temporary_handler import TemporaryEnvironment
from pathlib import Path

# def test_WorkflowBaseClass():
#     tmp = TemporaryEnvironment("generic")

#     water = molecule("H2O", cell=[1, 1, 1])
#     write("water.traj", water)

#     for i in range(1, 19+1):
#         a = WorkflowBaseClass(filepath="water.traj", jobnumber=i, cache="cache", debug=True)
#         a.get_displaced_structure()
#         a.close_logger()


#     tmp.return_to_main()
pseudo_dir = str(Path("test/test_data_irvib/espresso_pseudopotential").absolute())
espresso_command = f"mpirun -np 1 pw.x".split()
input_data_template = {
    # "tprnfor": False,
    # "control": {"tprnfor": False},
    "pseudo_dir": pseudo_dir,
    "occupations": "smearing",
    "smearing": "fermi-dirac",
    "degauss": 0.02,
    "ecutwfc": 20,
    "ecutrho": 160,
}

pseudopotentials = {
    "H2O": {"H": "h_pbe_v1.4.uspp.F.UPF", "O": "o_pbe_v1.2.uspp.F.UPF"},
}


def test_WorkflowBaseClass():
    tmp = TemporaryEnvironment("generic")

    water = molecule("H2O", cell=[5, 5, 5])
    write("water.traj", water)


    for i in range(1, 19 + 1):
        a = VibrationalAnalysisWorkflow(filepath="water.traj", jobnumber=i, cache="cache", debug=True)
        a.get_displaced_structure()

        input_data = input_data_template

        a.run(
            input_data=input_data,
            pseudopotentials=pseudopotentials["H2O"],
            kpts=(1, 1, 1),
            profile=espresso_command,
        )
        
        a.close_logger()

    tmp.return_to_main()
