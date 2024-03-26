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
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.style as mplstyle
from matplotlib.animation import FuncAnimation
from matplotlib.ticker import MultipleLocator, AutoMinorLocator

matplotlib.use("Agg")
mplstyle.use("fast")

from matscitoolkit.simple_logger import print_log
from matscitoolkit.utils.parallel import process_map
from matscitoolkit.utils.plot import plot_function, plot_spectra


def generateDisplacedStructures(
    filepath: str,
    disp_dir: str = "displacement_dir",
    filetype: str = "traj",
):
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
    ir_obj = Infrared(atoms_obj)
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


def computeEMAXPOS(filepath, field_direction=[1, 2, 3], debug=False):
    """
    Compute EMAXPOS from a given structure file.

    Parameters
    ----------
    filepath : str
        Path to the structure file.
    field_direction : list, optional
        List of field directions, by default [1, 2, 3].
    debug : bool, optional
        Print debug messages, by default False.
    """

    if debug:
        print_log("============================================================")
        print_log(f"Computing EMAXPOS \t {str(datetime.now())}")
        print_log("============================================================")

    atoms_obj = read(filepath)
    pos = atoms_obj.get_positions()
    x, y, z = pos.T

    cell = atoms_obj.get_cell()
    cx, cy, cz = cell.diagonal()

    x0, x1 = min(x) / cx, max(x) / cx
    y0, y1 = min(y) / cy, max(y) / cy
    z0, z1 = min(z) / cz, max(z) / cz

    emaxpos_x = (x0 + x1 + 1) / 2
    emaxpos_y = (y0 + y1 + 1) / 2
    emaxpos_z = (z0 + z1 + 1) / 2

    emaxpos = {1: emaxpos_x, 2: emaxpos_y, 3: emaxpos_z}

    if debug:
        # print results
        print_log(f"\tatomic space - x: {x0:.3f} {x1:.3f}")
        print_log(f"\tatomic space - y: {y0:.3f} {y1:.3f}")
        print_log(f"\tatomic space - z: {z0:.3f} {z1:.3f}")

        print_log(f"\tEMAXPOS (along directions): X={emaxpos_x:.3f} | Y={emaxpos_y:.3f} | Z={emaxpos_z:.3f}")
        print_log(f"\t**note: EMAXPOS should only be used in directions with vacuum")

    return [emaxpos[i] for i in field_direction]


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
        field_directions=[1, 2, 3],
        emaxpos="auto",
        eamp=0.0,
        eopreg=0.0001,
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
        self.field_directions = field_directions
        self.emaxpos = emaxpos
        self.eamp = eamp
        self.eopreg = eopreg
        self.calc_check = calc_check
        self.debug = debug
        self.filetype = filetype

        # Initialize system data
        filelist = glob(os.path.join(self.displacement_dir, f"*.{self.filetype}"))
        filedict = {int(os.path.basename(f).split(".")[0]): f for f in filelist}
        self.filepath = Path(filedict[self.system_id])
        self.filename = self.filepath.name
        self.system_data = read(self.filepath)
        base_filepath = filedict[1]

        # Initialize emaxpos
        if self.emaxpos == "auto":
            self.emaxpos = computeEMAXPOS(base_filepath, self.field_directions)
        else:
            assert len(self.emaxpos) == len(
                self.field_directions
            ), "emaxpos and field_directions must have the same length"


    def run(self):
        """Run DFT calculation."""

        if 'ESPRESSO_TMPDIR' not in os.environ:
            os.environ['ESPRESSO_TMPDIR'] = f"tmpdir"

        os.makedirs(self.calc_check, exist_ok=True)
        # calc_check_file = os.path.join(self.calc_check, f"{self.filename.removesuffix('.vasp')}.log")
        calc_check_file = os.path.join(self.calc_check, f"{self.filepath.stem}.log")
        fw = open(calc_check_file, "w")

        for i_edir, i_emaxpos in zip(self.field_directions, self.emaxpos):
            i_espresso_args = self.make_espresso_args(i_edir, i_emaxpos)
            if self.debug:
                self.check_dict(i_espresso_args)

            fw.write("============================================================\n")
            fw.write(f"Running DFT calculation : edir = {i_edir} \t {str(datetime.now())}\n")
            fw.write("============================================================\n")

            try:
                calc = Espresso(**i_espresso_args)
                # self.system_data.set_calculator(calc)
                self.system_data.calc = calc
                self.system_data.get_potential_energy()
                energy = self.system_data.get_potential_energy()
                force = self.system_data.get_forces()
                dipole = self.system_data.get_dipole_moment()
                fw.write(f"results :: edir {i_edir} success\n")
            except Exception as e:
                fw.write(f"results :: edir {i_edir} failed\n")
                fw.write(f"error :: {e}\n")

            self.remove_tmpdir()

    def make_input_data(self):
        """Make input data for Espresso calculator."""
        input_data_template = {"control": {}, "system": {}, "electrons": {}}
        input_data_template.update(self.input_data)

        input_data_template["control"].update({"tprnfor": True, "tefield": True, "dipfield": True})
        input_data_template["system"].update({"eamp": self.eamp, "eopreg": self.eopreg})

        return input_data_template

    def make_espresso_args(self, i_edir, i_emaxpos):
        """Make Espresso arguments."""
        espresso_args = {
            "input_data": self.make_input_data(),
            "pseudopotentials": self.pseudopotentials,
            "kpts": self.kpts,
            "directory": self.dirname,
            "profile": EspressoProfile(argv=self.espresso_command),
        }

        # Add emaxpos and edir
        espresso_args["input_data"]["system"].update({"edir": i_edir, "emaxpos": i_emaxpos})

        # Change directory name
        # FileID = f"{self.filename.removesuffix('.vasp')}"
        FileID = f"{self.filepath.stem}"
        self.directory_loc = espresso_args["directory"].replace("SUFFIX", f"/edir_{i_edir}/{FileID}")
        espresso_args["directory"] = self.directory_loc

        return espresso_args

    def check_dict(self, d):
        pprint(d, indent=2)

    def remove_tmpdir(self):
        # check if ESPRESSO_TMPDIR is set
        os.system(f"rm -r {self.directory_loc}/{os.environ['ESPRESSO_TMPDIR']}")



def check_dict(d):
    pprint(d, indent=2)


class DipoleCollector:
    """
    Collect dipole moment from DFT calculation.
    
    Parameters
    ----------
    dft_calc_dir : str, optional
        Path to the directory where the DFT calculation is saved, by default "dft_calc".
    """
    
    
    def __init__(self, dft_calc_dir="dft_calc"):
        self.dft_calc_dir = dft_calc_dir
        self.edir_list = sorted(glob(os.path.join(dft_calc_dir, "edir_*")))

        self.edir_dict = {}
        for i in self.edir_list:
            self.edir_dict[os.path.basename(i)] = sorted(glob(os.path.join(i, "*/*pwo")))

        # Initialize data dictionary
        self.data_dict = {}

        # Initialize directory
        os.makedirs("ir", exist_ok=True)

    def collect(self):
        """Collect dipole moment from DFT calculation."""
        def init_entry(full_dispname):
            if full_dispname not in self.data_dict:
                self.data_dict[full_dispname] = {
                    "all_forces": [],
                    "all_dipoles": [],
                    "net_force": [],
                    "net_dipole": [],
                }

        for key, val in self.edir_dict.items():
            for ival in tqdm(val, desc=f"Collecting from {key}"):
                filepath = ival
                _, edir, full_dispname, _ = filepath.split("/")
                dispname = full_dispname.split(".")[-1]
                # print(edir, full_dispname)

                init_entry(full_dispname)

                D = self.data_dict[full_dispname]
                # File path
                D[f"{edir}_file"] = filepath
                # Read force and dipole moment
                atoms_obj = read(filepath)
                D[f"{edir}_force"] = atoms_obj.get_forces()
                D[f"{edir}_dipole"] = atoms_obj.get_dipole_moment()

                D["all_forces"].append(D[f"{edir}_force"])
                D["all_dipoles"].append(D[f"{edir}_dipole"])

                # Compute net force and dipole moment
                D["net_force"] = np.average(D["all_forces"], axis=0)
                D["net_dipole"] = np.sum(D["all_dipoles"], axis=0)

                # JSON OUTPUT
                D["json_output_file"] = f"ir/cache.{dispname}.json"

                # WRITE SERIALIZED DATA
                self.serialize(D["net_force"], D["net_dipole"], D["json_output_file"])

        # SAVE ALL COLLECTED DATA
        with open("collected_data.json", "w") as f:
            json.dump(self.data_dict, f, indent=4, default=self.convert_np_to_list)

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

    def serialize(self, force_arr, dipole_arr, output_file):
        """Serialize data to JSON."""
        data = {"forces": force_arr, "dipole": dipole_arr}
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
    
    def __init__(self, structure_file=None, suffix="", infrared_options={}, nproc=1, plot_title="Infrared Plot", debug=True, use_cache=False):
        # Initialize all parameters

        self.INFRARED_OPTIONS = infrared_options
        if suffix != "":
            suffix = "_" + suffix

        self.MODES_DIR = "ir_modes{}".format(suffix)
        self.HIGH_MODES_DIR = "ir_high_modes{}".format(suffix)
        self.SUMMARY_FILE = "report_summary{}.dat".format(suffix)
        self.SPECTRA_FILE = "report_spectra{}.dat".format(suffix)
        self.SPECTRA_PLOT_FILE = "report_spectra{}.png".format(suffix)
        self.SPECTRA_PLOT_TITLE = plot_title
        self.CACHE_FILE = "report_infrared{}.cache.pkl".format(suffix)
        self.PARALLEL_NPROC = nproc
        self.debug = debug

        if self.debug:
            print_log("============================================================")
            print_log(f"Post-processing data \t {str(datetime.now())}")
            print_log("============================================================")
            print_log("\n----- GENERAL SETTINGS...")
            for var_name, var_value in self.__dict__.items():
                print_log(f"\t{var_name}: {var_value}")

        # Initialize placeholders
        self.ir = None
        self.atoms_obj = None
        
        # Infrared class options
        if "name" not in self.INFRARED_OPTIONS:
            self.INFRARED_OPTIONS["name"] = "ir"

        # Cache mode
        if use_cache:
            with open(self.CACHE_FILE, "rb") as f:
                self.ir = pickle.load(f)
        else:
            if structure_file is None:
                ref_file = Path("displacement_dir")
                self.atoms_obj = read(next(ref_file.glob("*.eq.*")))
                # self.atoms_obj = read(glob("displacement_dir/*.vasp")[0])
            else:
                self.atoms_obj = read(structure_file)
            
            
            self.ir = Infrared(self.atoms_obj, **self.INFRARED_OPTIONS)
            self.ir.get_vibrations()
            with open(self.CACHE_FILE, "wb") as f:
                pickle.dump(self.ir, f)

    def generate_summary(self):
        # STEP 1: Generate IR summary
        if self.debug:
            print_log("\n----- GENERATING SUMMARY...")
            print_log("\t Summary file: {self.SUMMARY_FILE}")

        with open(self.SUMMARY_FILE, "w") as f:
            self.ir.summary(log=f)

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

        self.ir.write_spectra(
            out=self.SPECTRA_FILE,
            start=start,
            end=end,
            width=width,
            type=method,
            intensity_unit="(D/A)2/amu",
        )

        self.SMEARING_WIDTH = width

    def load_summary(self, ntop=4):
        # STEP 3: Summary: load, clean, filter
        ## 3.1 Load summary
        modes = np.genfromtxt(
            self.SUMMARY_FILE,
            names=["mode", "frequency_mev", "frequency_cm", "intensity"],
            skip_header=4,
            skip_footer=4,
            dtype=[int, float, float, float],
        )
        self.modes = modes

        if self.debug:
            print_log("\n----- LOADING SUMMARY...")
            print_log(f"\t Summary dataset shape:  {modes.shape}")
            print_log(f"\t Summary dataset columns:  {modes.dtype.names}")

        ## 3.2 Clean
        ## remove rows with nan values
        nan_mask = np.isnan(modes["frequency_mev"]) | np.isnan(modes["frequency_cm"])
        self.clean_modes = modes[~nan_mask]

        ## 3.3 Filter
        ## filter by intensity
        a = self.clean_modes["intensity"]
        sorted_indices = np.argsort(a)
        top_index = sorted_indices[-ntop:][::-1]
        self.top_modes = self.clean_modes[top_index]

        if self.debug:
            print_log(f"\t Clean dataset shape:  {self.clean_modes.shape}")
            print_log(f"\t Clean dataset columns:  {self.clean_modes.dtype.names}")
            print_log(f"\t Top modes:")
            for _m, _mev, _cm, _int in self.top_modes:
                print_log(f"\t\t {_m:>5d}  {_mev:>.4}  {_cm:>.4f}  {_int:<.8f}")

        return self.clean_modes, self.top_modes

    def load_spectra(self):
        if self.debug:
            print_log("\n----- LOADING SPECTRA...")
        self.spectra = np.genfromtxt(self.SPECTRA_FILE, names=["wavenumber", "abs_ir_intensity", "absorbance_scaled"])
        return self.spectra

    def generate_mode_traj(self):
        if self.debug:
            print_log("\n----- GENERATING MODE TRAJECTORY FILES...")
        os.makedirs(self.MODES_DIR, exist_ok=True)
        os.makedirs(self.HIGH_MODES_DIR, exist_ok=True)

        def write_traj(mode):
            self.ir.write_mode(n=mode)

        parent = os.getcwd()
        os.chdir(self.MODES_DIR)
        results = process_map(
            write_traj,
            self.modes["mode"],
            desc="generating mode trajectory files",
            max_workers=self.PARALLEL_NPROC,
        )
        os.chdir(parent)

        parent = os.getcwd()
        os.chdir(self.HIGH_MODES_DIR)
        results = process_map(
            write_traj,
            self.top_modes["mode"],
            desc="generating high mode trajectory files",
            max_workers=self.PARALLEL_NPROC,
        )
        os.chdir(parent)

    def generate_mode_gif(self, mode_type="top"):
        """
        Generate mode GIF files.
        
        Parameters
        ----------
        mode_type : str, optional
            Type of mode, by default "top".        
        """
        if self.debug:
            print_log("\n----- GENERATING MODE GIF FILES...")

        if mode_type == "top":
            mode_type_dir = self.HIGH_MODES_DIR
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
            write_gif, trajfiles, max_workers=self.PARALLEL_NPROC, desc="generating mode GIF files"
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
            spectra_xy=[self.spectra["wavenumber"], self.spectra["abs_ir_intensity"]],
            spectra_label=f"Infrared spectra \nwith smearing width: {self.SMEARING_WIDTH}",
            exact_modes_x=self.clean_modes["frequency_cm"],
            exact_modes_y=self.clean_modes["intensity"],
            exact_modes_label="exact modes",
            exact_modes_color="r",
            xlabel="Frequency, cm$^{-1}$",
            ylabel="Intensity, (D/Å)$^2$ amu$^{-1}$",
            xrange=[0, 4000],
            title=self.SPECTRA_PLOT_TITLE,
            plotfile=self.SPECTRA_PLOT_FILE,
        )

        if self.debug:
            print_log("\t Saving spectra plot...")
            
        # if self.debug:
        #     print_log("\n----- GENERATING SPECTRA PLOT... ")
        #     print_log("\t Plotting spectra...")

        # fig, ax = plt.subplots(figsize=figsize)
        # ax.fill_between(
        #     self.spectra["wavenumber"],
        #     self.spectra["abs_ir_intensity"],
        #     label=f"Infrared spectra \nwith smearing width: {self.SMEARING_WIDTH}",
        #     color="k",
        #     alpha=0.8,
        # )
        # ax.vlines(
        #     x=self.clean_modes["frequency_cm"],
        #     ymin=0,
        #     ymax=self.clean_modes["intensity"],
        #     label="exact modes",
        #     color="r",
        #     alpha=1,
        #     linewidth=1,
        # )

        # ax.set_xlabel("Frequency, cm$^{-1}$", fontsize=fontsize)
        # ax.set_ylabel("Intensity, (D/Å)$^2$ amu$^{-1}$", fontsize=fontsize)
        # ax.tick_params(axis="x", labelsize=fontsize)
        # ax.tick_params(axis="y", labelsize=fontsize)

        # low, high = ax.get_ybound()
        # freq_limit = [0, 4000]
        # ax.set_ylim([low, high * 1.2])
        # ax.set_xlim(freq_limit)

        # ax.set_title(self.SPECTRA_PLOT_TITLE)
        # ax.xaxis.set_minor_locator(MultipleLocator(100))

        # fig.tight_layout()

        # if self.debug:
        #     print_log("\t Saving spectra plot...")

        # fig.savefig(self.SPECTRA_PLOT_FILE, dpi=dpi)
