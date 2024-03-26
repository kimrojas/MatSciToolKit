from ase.io import read, write
from ase.build import molecule
from pathlib import Path
from matscitoolkit import ase_vib_tools, ase_ir_tools
import os
import shutil
import pytest
from pathos import multiprocessing
from ase.calculators.espresso import Espresso, EspressoProfile


from temporary_handler import TemporaryEnvironment

# GLOBAL VARIABLES
pseudo_dir = str(Path("test/test_data_irvib/espresso_pseudopotential").absolute())
espresso_command = f"mpirun -np 1 pw.x".split()
profile = EspressoProfile(argv=espresso_command)


input_data_template = {
    "tprnfor": True,
    "pseudo_dir": pseudo_dir,
    "occupations": "smearing",
    "smearing": "fermi-dirac",
    "degauss": 0.02,
    "ecutwfc": 20,
    "ecutrho": 160,
}

pseudopotentials = {
    "H2O": {"H": "h_pbe_v1.4.uspp.F.UPF", "O": "o_pbe_v1.2.uspp.F.UPF"},
}


@pytest.fixture
def parallel(request):
    return request.config.getoption("parallel")


@pytest.mark.parametrize("name,expected", [("H2O", 19), ("C2H6", 49)])
@pytest.mark.parametrize("ext", ["vasp", "xyz", "traj"])
@pytest.mark.parametrize("tool", [ase_vib_tools, ase_ir_tools])
def test_generateDisplacedStructures(name, expected, tool, ext):
    tmp = TemporaryEnvironment(name, expected, tool.__name__, ext)

    mol = molecule(name, cell=[5, 5, 5])
    mol.center()
    write(f"{name}.{ext}", mol)

    kwargs = {
        "filepath": f"{name}.{ext}",
        "disp_dir": "displacement_dir",
        "filetype": ext,
    }

    tool.generateDisplacedStructures(**kwargs)

    nfiles = len(list(Path("displacement_dir").glob(f"*.{ext}")))

    # -- Return to main
    tmp.return_to_main()

    # -- Evaluation
    assert nfiles == expected, "Number of files is not correct"

    # -- Clean up
    tmp.clean()


def test_dipolemoment_calculation():

    tmp = TemporaryEnvironment("dipolemoment")

    atoms = molecule("H2O", cell=[5, 5, 5])
    atoms.center()
    input_data = input_data_template.copy()
    input_data.update(
        {
            "tefield": True,
            "dipfield": True,
            "edir": 3,
            "eamp": 0.00,
            "eopreg": 0.0001,
            "emaxpos": 0.0001,
        }
    )

    atoms.calc = Espresso(
        input_data=input_data,
        pseudopotentials=pseudopotentials["H2O"],
        kpts=(1, 1, 1),
        profile=profile,
    )

    atoms.get_potential_energy()
    dipol_arr = atoms.get_dipole_moment().tolist()
    # [-0.0, -0.0, -0.3299047698800686]
    expected_dipol_arr = [0, 0, -0.329]
    with open("dipolemoment.txt", "w") as f:
        f.write(str(dipol_arr))
    assert dipol_arr == pytest.approx(expected_dipol_arr, abs=0.1)

    tmp.return_to_main()
    tmp.clean()


class TestVibToolsWorkflow:
    tmp = TemporaryEnvironment("vibtools_workflow")

    atoms = molecule("H2O", cell=[5, 5, 5])
    atoms.center()
    filepath = "water.traj"
    write(filepath, atoms)

    nfiles = len(atoms) * 6 + 1

    tmp.return_to_main()

    def test_generate_displaced_structures(self):
        self.tmp.goto_tmp()

        structure_filepath = self.filepath
        ase_vib_tools.generateDisplacedStructures(structure_filepath)

        displacement_path = Path("displacement_dir")
        assert displacement_path.is_dir(), "[DisplacedGenerator] `displacement_dir` directory was not created"

        ndisplaced = len(list(displacement_path.glob("*.traj")))
        assert ndisplaced == 19, "[DisplacedGenerator] Number of displaced structures is not correct"

        self.tmp.return_to_main()

    def test_DFT_calculation(self, parallel):
        self.tmp.goto_tmp()

        input_data = input_data_template.copy()

        # -- Run DFT calculation (in parallel)
        def run_dft(args_id):
            dft = ase_vib_tools.DFTrunner(
                system_id=args_id,  # JOB array id -> displaced image id
                input_data=input_data,  # No need for specific parameter for dipole moment
                pseudopotentials=pseudopotentials["H2O"],
                kpts=[1, 1, 1],
                dirname="dft_calc/SUFFIX",
                espresso_command=espresso_command,
            )
            dft.run()

        with multiprocessing.Pool(processes=parallel) as pool:
            pool.map(run_dft, range(1, self.nfiles + 1))

        # -- [DFTrunner] evaluation
        dft_path = Path("dft_calc")
        assert dft_path.is_dir(), "[DFTrunner] `dft_calc` directory was not created"

        ndft = len(list(dft_path.glob("*/")))
        assert ndft == self.nfiles, "[DFTrunner] Number of DFT calculation is not correct"

        self.tmp.return_to_main()

    def test_collector(self):
        self.tmp.goto_tmp()

        # Step 3: Collector
        collector = ase_vib_tools.Collector()
        collector.collect()

        # -- [Collector] evaluation
        vib_path = Path("vib")
        assert vib_path.is_dir(), "[Collector] `vib` directory was not created"

        ncache = len(list(vib_path.glob("*.json")))
        assert ncache == self.nfiles, "[Collector] Number of cache files is not correct"

        self.tmp.return_to_main()

    def test_vibrational_analysis(self, parallel):
        self.tmp.goto_tmp()

        # Step 4: Vibrational analysis
        pp = ase_vib_tools.PostProcess(nproc=parallel)
        pp.generate_summary()
        pp.generate_spectra()
        pp.load_summary()
        pp.load_spectra()
        pp.generate_mode_traj()
        pp.generate_mode_gif()
        pp.generate_spectra_plot()

        self.tmp.return_to_main()

    # def test_clean(self):
    #     self.tmp.clean()


class TestIRToolsWorkflow:
    tmp = TemporaryEnvironment("ir_workflow")

    atoms = molecule("H2O", cell=[5, 5, 5])
    atoms.center()
    filepath = "water.traj"
    write(filepath, atoms)

    nfiles = len(atoms) * 6 + 1

    tmp.return_to_main()

    def test_generate_displaced_structures(self):
        self.tmp.goto_tmp()

        structure_filepath = self.filepath
        ase_ir_tools.generateDisplacedStructures(structure_filepath)

        displacement_path = Path("displacement_dir")
        assert displacement_path.is_dir(), "[DisplacedGenerator] `displacement_dir` directory was not created"

        ndisplaced = len(list(displacement_path.glob("*.traj")))
        assert ndisplaced == 19, "[DisplacedGenerator] Number of displaced structures is not correct"

        self.tmp.return_to_main()

    def test_computeEMAXPOS(self):
        self.tmp.goto_tmp()

        structure_filepath = self.filepath
        output = ase_ir_tools.computeEMAXPOS(structure_filepath, debug=True)

        assert 3 == len(output), "[computeEMAXPOS] Number of elements in the output is not correct"
        assert isinstance(output[0], float), "[computeEMAXPOS] First element is not float"
        assert isinstance(output[1], float), "[computeEMAXPOS] Second element is not float"
        assert isinstance(output[2], float), "[computeEMAXPOS] Third element is not float"
        
        self.tmp.return_to_main()

    def test_DFT_calculation(self, parallel):
        self.tmp.goto_tmp()

        input_data = input_data_template.copy()

        field_kwargs = {
            "field_directions": [1, 2, 3],
            "emaxpos": "auto",
            "eamp": 0.0,
            "eopreg": 0.0001,
        }

        # -- Run DFT calculation (in parallel)
        def run_dft(args_id):
            dft = ase_ir_tools.DFTrunner(
                system_id=args_id,  # JOB array id -> displaced image id
                input_data=input_data,  # No need for specific parameter for dipole moment
                pseudopotentials=pseudopotentials["H2O"],
                kpts=[1, 1, 1],
                dirname="dft_calc/SUFFIX",
                espresso_command=espresso_command,
                **field_kwargs,
            )
            dft.run()

        with multiprocessing.Pool(processes=parallel) as pool:
            pool.map(run_dft, range(1, self.nfiles + 1))

        # -- [DFTrunner] evaluation
        dft_path = Path("dft_calc")
        assert dft_path.is_dir(), "[DFTrunner] `dft_calc` directory was not created"

        ndft = len(list(dft_path.rglob("*.pwo")))
        assert ndft == self.nfiles*3, f"[DFTrunner] Number of DFT calculation is not correct {ndft} {self.nfiles*3}"

        self.tmp.return_to_main()

    def test_collector(self):
        self.tmp.goto_tmp()

        # Step 3: Collector
        collector = ase_ir_tools.DipoleCollector()
        collector.collect()

        # -- [Collector] evaluation
        ir_path = Path("ir")
        assert ir_path.is_dir(), "[Collector] `ir` directory was not created"

        ncache = len(list(ir_path.glob("*.json")))
        assert ncache == self.nfiles, "[Collector] Number of cache files is not correct"

        self.tmp.return_to_main()

    def test_infrared_analysis(self, parallel):
        self.tmp.goto_tmp()

        # Step 4: Vibrational analysis
        pp = ase_ir_tools.PostProcess(nproc=parallel)
        pp.generate_summary()
        pp.generate_spectra()
        pp.load_summary()
        pp.load_spectra()
        pp.generate_mode_traj()
        pp.generate_mode_gif()
        pp.generate_spectra_plot()

        self.tmp.return_to_main()

    # def test_clean(self):
    #     self.tmp.clean()
