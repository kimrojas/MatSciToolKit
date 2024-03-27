from abc import ABC, abstractmethod
from ase.io import read, write
from ase.vibrations import Vibrations, Infrared
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
        self.jobnumber = int(jobnumber)
        self.log = logger(logfile=f"{self.filestem}_{self.jobnumber}.log", debug=debug)
        self.log.info(f"Subclass name: {self.__class__.__name__}s")

        # Initialize cache directory
        self.cache = Path(cache)
        self.cache.mkdir(exist_ok=True)
        self.log.info(f"Cache directory created: {str(self.cache)}")

    def get_displaced_structure(self, generatefile=True, directory=None, methodkwargs={}):
        """Produces the displaced structure for a given job number"""

        # Add default values for methodkwargs
        default_methodkwargs = {"indices": None, "delta": 0.01, "nfree": 2, "directions": None}
        for k, v in default_methodkwargs.items():
            methodkwargs.setdefault(k, v)

        # Compute all displaced structures
        self.log.info(f"Reference structure: '{self.filename}'")
        displaced_structures = dict(enumerate(Infrared(read(self.filepath), **methodkwargs).iterdisplace(), start=1))
        self.nfiles = len(displaced_structures)
        self.dim = len(str(self.nfiles))
        self.log.info(f"Expected number of displaced structures: {self.nfiles}")

        # Select/Map the job number to the displaced structure
        self.log.info(f"Job number: {self.jobnumber}")
        self.log.assert_(
            1 <= self.jobnumber <= self.nfiles, f"Job number {self.jobnumber} is out of range [1-{self.nfiles}]"
        )
        disp, atm = displaced_structures[self.jobnumber]

        # Organize job information
        self.log.info(f"Job name: {disp.name}")
        self.job = {"number": self.jobnumber, "name": disp.name, "structure": atm}

        # Print/Output structure file for displaced structure
        if directory is None:
            directory = self.cache / "displaced_structures"
        else:
            directory = self.cache / directory

        if generatefile:
            dispfile = directory / f"{self.jobnumber:0{self.dim}d}.{self.job['name']}{self.filetype}"
            directory.mkdir(exist_ok=True)  # Initialize directory
            write(dispfile, self.job["structure"])  # Write structure file
            self.log.info(f"Displaced structure saved in: {str(dispfile)}")

    @abstractmethod
    def run(self):
        """RUN DFT CALCULATION on self.job['structure']"""
        pass

    @abstractmethod
    def clean(self, directory=None):
        """CLEAN TEMPORARY DIRECTORY"""
        pass

    def close_logger(self):
        self.log.close()
