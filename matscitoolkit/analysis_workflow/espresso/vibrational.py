from matscitoolkit.analysis_workflow.ABC_workflow import WorkflowBaseClass
from matscitoolkit.utils.ensure_key import ensure_key
from ase.vibrations import Vibrations, Infrared
from ase.calculators.espresso import Espresso, EspressoProfile
import json


class VibrationalAnalysisWorkflow(WorkflowBaseClass):

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def run(self, input_data, pseudopotentials, kpts, profile=["mpirun", "pw.x"], directory="dft",):        
        # Work directory
        directory = self.cache / directory / self.job["fullname"]
        self.goto_workdir(directory)
        
        self.log.info(f"Input data detected: \n{json.dumps(input_data, indent=4)}")
        
        
        input_data = ensure_key(input_data, "tprnfor", True, logger=self.log)
        self.log.info(f"Input data after ensure_key: \n{json.dumps(input_data, indent=4)}")

        # Define calculator
        calculator = Espresso(
            input_data=input_data,
            pseudopotentials=pseudopotentials,
            kpts=kpts,
            profile=EspressoProfile(profile),
        )

        # Run the calculation
        super().run(calculator=calculator, directory=directory)      

        # Return to main directory
        self.goto_maindir()

    def clean(self, directory=None):
        pass
