import pytest
from matscitoolkit.analysis_workflow.ABC_workflow import WorkflowBaseClass
from ase.io import read, write
from ase.build import molecule
from temporary_handler import TemporaryEnvironment

def test_WorkflowBaseClass():
    tmp = TemporaryEnvironment("generic")
    
    water = molecule("H2O", cell=[1, 1, 1])
    write("water.traj", water)
    
    for i in range(1, 19):
        a = WorkflowBaseClass(filepath="water.traj", jobnumber=i, cache="cache", debug=True)
        a.close_logger()
    

    tmp.return_to_main()
