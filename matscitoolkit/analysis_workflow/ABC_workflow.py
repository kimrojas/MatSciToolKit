from abc import ABC, abstractmethod
from ase.io import read, write
from ase.vibrations import Vibrations, Infrared
from pathlib import Path
from matscitoolkit.analysis_workflow.logger import logger
from matscitoolkit.utils.ensure_key import ensure_key
import os
from copy import copy, deepcopy


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

        # Main directory
        self.main_path = Path.cwd()

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
        self.job = {
            "number": self.jobnumber,
            "name": disp.name,
            "structure": atm,
            "fullname": f"{self.jobnumber:0{self.dim}d}.{disp.name}",
        }

        # Print/Output structure file for displaced structure
        if directory is None:
            directory = self.cache / "displaced_structures"
        else:
            directory = self.cache / directory

        if generatefile:
            dispfile = directory / f"{self.job['fullname']}{self.filetype}"
            directory.mkdir(exist_ok=True)  # Initialize directory
            write(dispfile, self.job["structure"])  # Write structure file
            self.log.info(f"Displaced structure saved in: {str(dispfile)}")

    def irun(self, calculator, directory, tag=""):
        """RUN DFT CALCULATION on self.job['structure']"""
        self.log.info(f"Running '{self.job['fullname']}' {tag} in {directory}/")
        
        # Attach calculator
        self.job["structure"].calc = calculator
        
        # Start energy calculation
        try:
            self.job["structure"].get_potential_energy()
        except Exception as e:
            self.log.error(f"Error: {e}")
            self.goto_maindir()
            raise Exception(f"Error: {e}")
        else:
            self.log.info(f"Energy calculation successful")

        # Start energy calculation
        try:
            self.job["structure"].get_forces()
        except Exception as e:
            self.log.error(f"Error: {e}")
            self.goto_maindir()
            raise Exception(f"Error: {e}")
        else:
            self.log.info(f"Force calculation successful")

    @abstractmethod
    def clean(self, directory=None):
        """CLEAN TEMPORARY DIRECTORY"""
        pass
    
    # def collect(self):
        

    def goto_workdir(self, directory):
        directory = Path(directory)
        directory.mkdir(exist_ok=True, parents=True)
        os.chdir(directory)

    def goto_maindir(self):
        os.chdir(self.main_path)

    def close_logger(self):
        self.log.close()
