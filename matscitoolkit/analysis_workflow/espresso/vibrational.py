from matscitoolkit.analysis_workflow.ABC_workflow import WorkflowBaseClass
from ase.vibrations import Vibrations, Infrared
from ase.calculators.espresso import Espresso, EspressoProfile


class VibrationalAnalysisWorkflow(WorkflowBaseClass):

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def run(self, input_data, pseudopotentials, kpts, profile=["mpirun", "pw.x"], directory="dft",):        
        # Work directory
        directory = self.cache / directory / self.job["fullname"]
        self.goto_workdir(directory)
        
        self.log.info(f"Input data: {input_data}")
        input_data = self.ensure_key(input_data, "tprnfor", True)
        self.log.info(f"Input data: {input_data}")

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
