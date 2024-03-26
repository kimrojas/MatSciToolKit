from abc import ABC, abstractmethod
from ase.io import read, write
from pathlib import Path
from matscitoolkit.analysis_workflow.logger import logger


def test_logger():
    log = logger()
    # Example logging messages
    log.debug("This is a debug message")
    log.info("This is an info message")
    log.warning("This is a warning message")
    log.error("This is an error message")
    log.critical("This is a critical message")
    log.assert_(True, "This is an assertion message")
    # log.assert_(False, "This is an assertion message")


class WorkflowBaseClass(ABC):

    def __init__(self, filepath=None, jobnumber=None, cache="cache", debug=True):
        # Reference structure
        self.filepath = Path(filepath)
        self.filename = self.filepath.name
        self.filestem = self.filepath.stem
        self.filetype = self.filepath.suffix

        # Initialize logger
        nfiles = 1 + (len(read(self.filepath)) * 6)
        dim = len(str(nfiles))
        self.jobnumber = int(jobnumber)
        self.log = logger(logfile=f"{self.filestem}_{self.jobnumber:0{dim}d}.log", debug=debug)
        self.log.info(f"Subclass name: {self.__class__.__name__}s")

        # Initialize cache directory
        self.cache = Path(cache)
        self.cache.mkdir(exist_ok=True)
        self.log.info(f"Cache directory created: {str(self.cache)}")

        # Input checks
        self.log.info(f"Reference structure: '{self.filename}'")
        self.log.info(f"Expected number of displaced structures: {nfiles}")
        self.log.info(f"Job number: {self.jobnumber}")
        self.log.assert_(1 <= self.jobnumber <= nfiles, f"Job number {self.jobnumber} is out of range [1-{nfiles}]")
    
    def close_logger(self):
        self.log.close()

