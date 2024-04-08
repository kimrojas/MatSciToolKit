from abc import ABC, abstractmethod
from ase.io import read, write
from ase.vibrations import Vibrations, Infrared
from pathlib import Path
from matscitoolkit.analysis_workflow.logger import logger
from matscitoolkit.utils.ensure_key import ensure_key
import os
from copy import copy, deepcopy
from matscitoolkit.utils.serial import serialize, deserialize


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

    def __init__(self, filepath=None, jobnumber=None, cache="cache", debug=True, logfile=None):
        # Reference structure
        self.filepath = Path(filepath)
        self.filename = self.filepath.name
        self.filestem = self.filepath.stem
        self.filetype = self.filepath.suffix

        # Initialize logger
        self.jobnumber = int(jobnumber)
        if logfile is None:
            logfile = f"{self.filestem}_{self.jobnumber}.log"
        self.log = logger(logfile=logfile, debug=debug)
        self.log.info(f"Subclass name: {self.__class__.__name__}s")

        # Initialize cache directory
        self.cache = Path(cache)
        self.cache.mkdir(exist_ok=True)
        self.log.info(f"Cache directory created: {str(self.cache)}")

        # Main directory
        self.main_path = Path.cwd()

    def get_displaced_structure(self, generatefile=True, directory=None, method="vib", methodkwargs={}):
        """Produces the displaced structure for a given job number"""
        self.method = method
        if method == "vib":
            self.methodclass = Vibrations
            default_methodkwargs = {"indices": None, "delta": 0.01, "nfree": 2}
        if method == "ir":
            self.methodclass = Infrared
            default_methodkwargs = {"indices": None, "delta": 0.01, "nfree": 2, "directions": None}

        # Add default values for methodkwargs
        for k, v in default_methodkwargs.items():
            methodkwargs.setdefault(k, v)

        # Compute all displaced structures
        self.log.info(f"Reference structure: '{self.filename}'")
        displaced_structures = dict(
            enumerate(self.methodclass(read(self.filepath), **methodkwargs).iterdisplace(), start=1)
        )
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

    def collect(self, method="vib", cache_source="dft", outputpattern="espresso.pwo"):
        self.log.info(":: COLLECT module ::")
        source = self.cache / cache_source
        destination = method
        self.log.info(f"Collecting files from {source} to {destination}")
        self.log.info(f"Destination directory: {destination}")
        self.log.info(f"Source directory: {source}")

        # Create destination directory
        filelist = sorted(list(Path(source).rglob(outputpattern)))
        filelisttxt = "\n" + "\n".join([str(f) for f in filelist])
        parentslist = [str(f.parent).split("/") for f in filelist]
        keylist = [p[2].split(".")[-1] for p in parentslist]
        keylisttxt = "\n" + "\n".join(keylist)
        Path(destination).mkdir(exist_ok=True)
        
        # Check file depth
        filedepth = len(parentslist[0])
        # filedepth = 3 means, dipole are in just 1 file
        # filedepth = 4 means, dipole separated by dimensions in multiple files.
        
        self.log.info(f"Files to be collected: {filelisttxt}")
        self.log.info(f"Key list: {keylisttxt}")
        
        # Read data (force/dipole) and serialize
        if filedepth == 3:
            for f, key in zip(filelist, keylist):
                atoms_obj = read(f)
                
                serialize_input = {}
                serialize_input["force_arr"] = atoms_obj.get_forces()
                serialize_input["output_file"] = Path(destination) / f"cache.{key}.json"
                
                if method == "ir":
                    serialize_input["dipole_arr"] = atoms_obj.get_dipole()

                serialize(**serialize_input)
        
        if filedepth == 4:
            pass
        
        
        
        

    def goto_workdir(self, directory):
        directory = Path(directory)
        directory.mkdir(exist_ok=True, parents=True)
        os.chdir(directory)

    def goto_maindir(self):
        os.chdir(self.main_path)

    def close_logger(self):
        self.log.close()
