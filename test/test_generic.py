import pytest
from matscitoolkit.analysis_workflow.ABC_workflow import WorkflowBaseClass
from matscitoolkit.analysis_workflow.espresso import VibrationalAnalysisWorkflow
from ase.io import read, write
from ase.build import molecule
from temporary_handler import TemporaryEnvironment

# def test_WorkflowBaseClass():
#     tmp = TemporaryEnvironment("generic")
    
#     water = molecule("H2O", cell=[1, 1, 1])
#     write("water.traj", water)
    
#     for i in range(1, 19+1):
#         a = WorkflowBaseClass(filepath="water.traj", jobnumber=i, cache="cache", debug=True)
#         a.get_displaced_structure()
#         a.close_logger()
    

#     tmp.return_to_main()


def test_WorkflowBaseClass():
    tmp = TemporaryEnvironment("generic")
    
    water = molecule("H2O", cell=[1, 1, 1])
    water.set_masses([10, 10, 160])
    write("water.traj", water)
    
    for i in range(1, 19+1):
        a = VibrationalAnalysisWorkflow(filepath="water.traj", jobnumber=i, cache="cache", debug=True)
        a.get_displaced_structure()
        a.close_logger()
    

    tmp.return_to_main()
