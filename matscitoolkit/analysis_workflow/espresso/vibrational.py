from matscitoolkit.analysis_workflow.ABC_workflow import WorkflowBaseClass


class VibrationalAnalysisWorkflow(WorkflowBaseClass):
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def run(self):
        pass
    
    def clean(self, directory=None):
        pass