import ase
from ase.io import read, write
from ase.vibrations import Vibrations, Infrared
from ase.calculators.espresso import Espresso, EspressoProfile

import os
import json
import pickle
import numpy as np
from datetime import datetime
from pprint import pprint
from glob import glob
from tqdm import tqdm

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.style as mplstyle
from matplotlib.animation import FuncAnimation
from matplotlib.ticker import MultipleLocator, AutoMinorLocator

from matscitoolkit.simple_logger import print_log
from matscitoolkit.utils.parallel import process_map
from matscitoolkit.utils.plot import plot_function, plot_spectra


def generateDisplacedStructures(filepath: str, disp_dir: str = "displacement_dir", filetype: str = "traj"):
    """
    Generate displaced structures from a given structure file.

    Parameters
    ----------
    filepath : str
        Path to the structure file.
    disp_dir : str, optional
        Path to the directory where the displaced structures will be saved, by default "displacement_dir".
    """

    print_log("============================================================")
    print_log(f"Generating displaced structures \t {str(datetime.now())}")
    print_log("============================================================")

    # Initialize
    atoms_obj = read(filepath)
    ir_obj = Vibrations(atoms_obj)
    os.makedirs(disp_dir, exist_ok=True)

    # Generate displaced structures
    images = list(ir_obj.iterdisplace())
    print_log(f"\tTotal items: {len(images)}")
    width = len(str(len(images)))

    for i, image in enumerate(images, start=1):
        disp, atm = image
        data = {"name": disp.name, "mode": disp.a, "direction": disp.i, "sign": disp.sign, "ndisp": disp.ndisp}

        filename = f"{i:0{width}d}.{data['name']}.{filetype}"
        write(os.path.join(disp_dir, filename), atm)
        print_log(f"\t{i:>{width}d} {data['name']:10s} :: {data}")

    print_log("\n\n")


class DFTrunner:
    """
    Run DFT calculation using ASE Espresso calculator.

    Parameters
    ----------
    system_id : int
        ID of the system to be calculated.
    input_data : dict
        Input data for Espresso calculator.
    pseudopotentials : dict
        Pseudopotentials for Espresso calculator.
    kpts : list
        K-points for Espresso calculator.
    displacement_dir : str, optional
        Path to the directory where the displaced structures are saved, by default "displacement_dir".
    dirname : str, optional
        Path to the directory where the DFT calculation will be saved, by default "dft_calc/SUFFIX".
    espresso_command : list, optional
        Command to run Espresso, by default ["mpirun", "pw.x"].
    field_directions : list, optional
        List of field directions, by default [1, 2, 3].
    emaxpos : str, optional
        EMAXPOS for Espresso calculator, by default "auto".
    eamp : float, optional
        EAMP for Espresso calculator, by default 0.005.
    eopreg : float, optional
        EOPREG for Espresso calculator, by default 0.01.
    calc_check : str, optional
        Path to the directory where the calculation check files will be saved, by default "calc_check".
    debug : bool, optional
        Print debug messages, by default False.
    """

    def __init__(
        self,
        system_id: int,
        input_data: dict,
        pseudopotentials: dict,
        kpts: list,
        displacement_dir="displacement_dir",
        dirname="dft_calc/SUFFIX",
        espresso_command=["mpirun", "pw.x"],
        calc_check="calc_check",
        debug=False,
        filetype="traj",
    ):
        self.system_id = system_id
        self.input_data = input_data
        self.pseudopotentials = pseudopotentials
        self.kpts = kpts
        self.displacement_dir = displacement_dir
        self.dirname = dirname
        self.espresso_command = espresso_command
        self.calc_check = calc_check
        self.debug = debug
        self.filetype = filetype

        # Initialize system data
        filelist = glob(os.path.join(self.displacement_dir, f"*.{self.filetype}"))
        filedict = {int(os.path.basename(f).split(".")[0]): f for f in filelist}
        filepath = filedict[self.system_id]
        self.filename = os.path.basename(filepath)
        self.system_data = read(filepath)
        base_filepath = filedict[1]

    def run(self):
        """Run DFT calculation."""

        if "ESPRESSO_TMPDIR" not in os.environ:
            os.environ["ESPRESSO_TMPDIR"] = f"tmpdir"

        os.makedirs(self.calc_check, exist_ok=True)
        _fn = self.filename.removesuffix(f".{self.filetype}")
        calc_check_file = os.path.join(self.calc_check, f"{_fn}.log")
        fw = open(calc_check_file, "w")

        fw.write("============================================================\n")
        fw.write(f"Running DFT calculation :  \t {str(datetime.now())}\n")
        fw.write("============================================================\n")

        i_espresso_args = self.make_espresso_args()
        try:
            calc = Espresso(**i_espresso_args)
            # self.system_data.set_calculator(calc)
            self.system_data.calc = calc
            self.system_data.get_potential_energy()
            energy = self.system_data.get_potential_energy()
            force = self.system_data.get_forces()
            fw.write(f"results ::  success\n")
        except Exception as e:
            fw.write(f"results ::  failed\n")
            fw.write(f"error :: {e}\n")

        self.remove_tmpdir()

    def make_input_data(self):
        """Make input data for Espresso calculator."""
        input_data_template = {"control": {}, "system": {}, "electrons": {}}
        input_data_template.update(self.input_data)

        input_data_template["control"].update({"tprnfor": True})

        return input_data_template

    def make_espresso_args(self):
        """Make Espresso arguments."""
        espresso_args = {
            "input_data": self.make_input_data(),
            "pseudopotentials": self.pseudopotentials,
            "kpts": self.kpts,
            "directory": self.dirname,
            "profile": EspressoProfile(argv=self.espresso_command),
        }

        # Change directory name
        FileID = self.filename.removesuffix(f".{self.filetype}")
        self.directory_loc = espresso_args["directory"].replace("SUFFIX", f"/{FileID}")
        espresso_args["directory"] = self.directory_loc

        return espresso_args

    def check_dict(self, d):
        pprint(d, indent=2)

    def remove_tmpdir(self):
        # check if ESPRESSO_TMPDIR is set
        os.system(f"rm -r {self.directory_loc}/{os.environ['ESPRESSO_TMPDIR']}")


class Collector:
    """
    Collect Forces from DFT calculation.

    Parameters
    ----------
    dft_calc_dir : str, optional
        Path to the directory where the DFT calculation is saved, by default "dft_calc".
    """

    def __init__(self, dft_calc_dir="dft_calc"):
        self.dft_calc_dir = dft_calc_dir
        self.dir_list = sorted(glob(f"{dft_calc_dir}/*"))

        pprint(self.dir_list)

        self.dir_dict = {}
        for i in self.dir_list:
            dispname = os.path.basename(i).split(".")[-1]
            self.dir_dict[dispname] = glob(os.path.join(i, "*pwo"))[0]

        pprint(self.dir_dict)

        # Initialize data dictionary
        self.data_dict = {}

        # Initialize directory
        os.makedirs("vib", exist_ok=True)

    def collect(self):
        """Collect forces from DFT calculation."""
        for key, val in self.dir_dict.items():
            full_dispname = f"cache.{key}.json"

            atoms_obj = read(val)
            forces = atoms_obj.get_forces()

            output_file = f"vib/{full_dispname}"
            self.serialize(forces, output_file)

        # # SAVE ALL COLLECTED DATA
        # with open("collected_data.json", "w") as f:
        #     json.dump(self.data_dict, f, indent=4, default=self.convert_np_to_list)

    def convert_np_to_list(self, obj):
        """Convert numpy array to list."""
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return obj

    def custom_deserializer(self, dct):
        """Custom deserializer for JSON."""
        if "__ndarray__" in dct:
            shape, dtype, data = dct["__ndarray__"]
            return np.array(data, dtype=dtype).reshape(shape)
        return dct

    def deserialize(self, input_file):
        """Deserialize JSON file."""
        with open(input_file, "r") as f:
            return json.load(f, object_hook=self.custom_deserializer)

    def custom_serializer(self, obj):
        """Custom serializer for JSON."""
        if isinstance(obj, np.ndarray):
            return {"__ndarray__": [obj.shape, obj.dtype.name, obj.flatten().tolist()]}
        raise TypeError(f"Type {type(obj)} not serializable")

    def serialize(self, force_arr, output_file):
        """Serialize data to JSON."""
        data = {"forces": force_arr}
        with open(output_file, "w") as f:
            json.dump(data, f, default=self.custom_serializer)
        return json.dumps(data, default=self.custom_serializer)


class PostProcess:
    """
    Post-process DFT calculation using ASE Vibrations and Infrared.

    Parameters
    ----------
    structure_file : str, optional
        Path to the structure file, by default None.
    nproc : int, optional
        Number of processors to be used, by default 1.
    plot_title : str, optional
        Title of the spectra plot, by default "Infrared Plot".
    debug : bool, optional
        Print debug messages, by default True.
    use_cache : bool, optional
        Use cached data, by default False.
    """

    def __init__(
        self, structure_file=None, nproc=1, plot_title="Vibration Plot", debug=True, use_cache=False, filetype="traj"
    ):
        # Initialize all parameters

        self.MODES_DIR = "vib_modes"
        self.HIGH_MODES_DIR = "vib_high_modes"
        self.SUMMARY_FILE = "report_summary.dat"
        self.SPECTRA_FILE = "report_spectra.dat"
        self.SPECTRA_PLOT_FILE = "report_spectra.png"
        self.SPECTRA_PLOT_TITLE = plot_title
        self.CACHE_FILE = "report_vibration.cache.pkl"
        self.PARALLEL_NPROC = nproc
        self.debug = debug
        self.filetype = filetype

        if self.debug:
            print_log("============================================================")
            print_log(f"Post-processing data \t {str(datetime.now())}")
            print_log("============================================================")
            print_log("\n----- GENERAL SETTINGS...")
            for var_name, var_value in self.__dict__.items():
                print_log(f"\t{var_name}: {var_value}")

        # Initialize placeholders
        self.vib = None
        self.atoms_obj = None

        # Cache mode
        if use_cache:
            with open(self.CACHE_FILE, "rb") as f:
                self.vib = pickle.load(f)
        else:
            if structure_file is None:
                self.atoms_obj = read(glob(f"displacement_dir/*.{self.filetype}")[0])
            else:
                self.atoms_obj = read(structure_file)
            self.vib = Vibrations(self.atoms_obj, name="vib")
            self.vib.get_vibrations()
            with open(self.CACHE_FILE, "wb") as f:
                pickle.dump(self.vib, f)

    def generate_summary(self):
        # STEP 1: Generate VIB summary
        if self.debug:
            print_log("\n----- GENERATING SUMMARY...")
            print_log(f"\t Summary file: {self.SUMMARY_FILE}")

        with open(self.SUMMARY_FILE, "w") as f:
            self.vib.summary(log=f)

    def generate_spectra(self, start=200, end=5000, width=30, method="Lorentzian"):
        """
        Generate IR spectra.

        Parameters
        ----------
        start : int, optional
            Starting frequency, by default 200.
        end : int, optional
            Ending frequency, by default 5000.
        """

        # STEP 2: Generate IR spectra
        if self.debug:
            print_log("\n----- GENERATING SPECTRA...")
            print_log(f"\t Spectra file: {self.SPECTRA_FILE}")
            print_log(f"\t Width: {width}")
            print_log(f"\t Method: {method}")

        self.vib.write_dos(
            out=self.SPECTRA_FILE,
            start=start,
            end=end,
            width=width,
            type=method,
            # intensity_unit="(D/A)2/amu",
        )

        self.SMEARING_WIDTH = width

    def load_spectra(self):
        if self.debug:
            print_log("\n----- LOADING SPECTRA...")
        self.spectra = np.genfromtxt(self.SPECTRA_FILE, names=["wavenumber", "intensity"])

        return self.spectra

    def load_summary(self):
        # STEP 3: Summary: load, clean, filter
        ## 3.1 Load summary
        modes = np.genfromtxt(
            self.SUMMARY_FILE,
            names=["mode", "frequency_mev", "frequency_cm"],
            skip_header=3,
            skip_footer=2,
            dtype=[int, float, float],
        )
        self.modes = modes

        if self.debug:
            print_log("\n----- LOADING SUMMARY...")
            print_log(f"\t Summary dataset shape:  {modes.shape}")
            print_log(f"\t Summary dataset columns:  {modes.dtype.names}")

        pprint(self.modes)

        ## 3.2 Clean
        ## remove rows with nan values
        nan_mask = np.isnan(modes["frequency_mev"]) | np.isnan(modes["frequency_cm"])
        self.clean_modes = modes[~nan_mask]

        pprint(self.clean_modes)

        if self.debug:
            print_log(f"\t Clean dataset shape:  {self.clean_modes.shape}")
            print_log(f"\t Clean dataset columns:  {self.clean_modes.dtype.names}")

    def generate_mode_traj(self):
        if self.debug:
            print_log("\n----- GENERATING MODE TRAJECTORY FILES...")
        os.makedirs(self.MODES_DIR, exist_ok=True)

        def write_traj(mode):
            self.vib.write_mode(n=mode)

        parent = os.getcwd()
        os.chdir(self.MODES_DIR)
        results = process_map(
            write_traj,
            self.modes["mode"],
            desc="generating mode trajectory files",
            max_workers=self.PARALLEL_NPROC,
        )
        os.chdir(parent)

    def generate_mode_gif(self, mode_type="all"):
        """
        Generate mode GIF files.

        Parameters
        ----------
        mode_type : str, optional
            Type of mode, by default "all".
        """
        if self.debug:
            print_log("\n----- GENERATING MODE GIF FILES...")

        if mode_type == "all":
            mode_type_dir = self.MODES_DIR

        def write_gif(traj):
            atoms_obj = read(traj, index=":")
            mode = traj.split("/")[-1].removesuffix(".traj")
            plot_function(os.path.join(mode_type_dir, mode + ".gif"), atoms_obj, title=f"IR mode: {mode}")

        trajfiles = sorted(glob(os.path.join(mode_type_dir, "*.traj")))

        # for t in tqdm(trajfiles):
        # write_gif(t)
        plot_results = process_map(
            write_gif,
            trajfiles,
            max_workers=self.PARALLEL_NPROC,
            desc="generating mode GIF files",
        )

    def generate_spectra_plot(self, figsize=(4.4, 3.52), fontsize=9, dpi=150):
        """
        Generate spectra plot.

        Parameters
        ----------
        figsize : tuple, optional
            Figure size, by default (4.4, 3.52).
        fontsize : int, optional
            Font size, by default 9.
        dpi : int, optional
            DPI, by default 150.
        """
        if self.debug:
            print_log("\n----- GENERATING SPECTRA PLOT... ")
            print_log("\t Plotting spectra...")

        plot_spectra(
            figsize=(4.4, 3.52),
            fontsize=9,
            dpi=150,
            spectra_xy=[self.spectra["wavenumber"], self.spectra["intensity"]],
            spectra_label=f"Vibrations DOS spectra \nwith smearing width: {self.SMEARING_WIDTH}",
            exact_modes_x=self.clean_modes["frequency_cm"],
            exact_modes_y=1,
            exact_modes_label="exact modes",
            exact_modes_color="g",
            xlabel="Wavenumber, cm$^{-1}$",
            ylabel="Vibrational DOS, arb. units.",
            xrange=[0, 4000],
            title=self.SPECTRA_PLOT_TITLE,
            plotfile=self.SPECTRA_PLOT_FILE,
        )

        if self.debug:
            print_log("\t Saving spectra plot...")

        # fig.savefig(self.SPECTRA_PLOT_FILE, dpi=dpi)
