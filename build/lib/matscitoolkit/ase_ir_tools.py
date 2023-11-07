from ase.io import read, write
from ase.vibrations import Vibrations, Infrared
from ase.calculators.espresso import Espresso, EspressoProfile

import os
import json
from datetime import datetime
from matscitoolkit.simple_logger import print_log


def generateDisplacedStructures(
    filepath: str,
    calcname: str = "ir",
    disp_dir: str = "displacement_dir",
    logfile: str = "MatciSciToolkit.log",
):
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

        filename = f"{i:0{width}d}.{data['name']}.vasp"
        write(os.path.join(disp_dir, filename), atm)
        print_log(f"\t{i:>{width}d} {data['name']:10s} :: {data}")

    print_log("\n\n")


def computeEMAXPOS(filepath, field_direction=[1, 2, 3]):
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

    # print results
    print_log(f"\tatomic space - x: {x0:.3f} {x1:.3f}")
    print_log(f"\tatomic space - y: {y0:.3f} {y1:.3f}")
    print_log(f"\tatomic space - z: {z0:.3f} {z1:.3f}")

    emaxpos_x = (x0 + x1 + 1) / 2
    emaxpos_y = (y0 + y1 + 1) / 2
    emaxpos_z = (z0 + z1 + 1) / 2

    emaxpos = {1: emaxpos_x, 2: emaxpos_y, 3: emaxpos_z}

    print_log(f"\tEMAXPOS (along directions): X={emaxpos_x:.3f} | Y={emaxpos_y:.3f} | Z={emaxpos_z:.3f}")
    print_log(f"\t**note: EMAXPOS should only be used in directions with vacuum")

    return [emaxpos[i] for i in field_direction]


class DFTrunner:
    def __init__(
        self,
        system_id: str,
        input_data: dict,
        pseudopotentials: dict,
        kpts: list,
        dirname="dft_calc_SUFFIX",
        espresso_command=["mpirun", "pw.x"],
        field_directions=[1, 2, 3],
        emaxpos="auto",
        eamp=0.005,
        eopreg=0.01,
    ):
        self.system_id = system_id
        self.input_data = input_data
        self.pseudopotentials = pseudopotentials
        self.kpts = kpts
        self.dirname = dirname
        self.espresso_command
        self.field_directions = field_directions
        self.emaxpos = emaxpos
        self.eamp = eamp
        self.eopreg = eopreg
       
        # Espresso args
        espresso_args = self.make_espresso_args()
        self.check_dict(espresso_args)

    def make_input_data(self):
        input_data_addon = {
            "control": {"tprnfor": True, "tefield": True, "dipfield": True},
            "system": {"eamp": self.eamp, "eopreg": self.eopreg},
        }
        
        input_data_template = self.input_data.copy()
        for key, val in input_data_addon.items():
            if key not in input_data_template.keys():
                input_data_template[key] = val
            else:
                input_data_template[key].update(val)
        
        return input_data_template

    def make_espresso_args(self):
        espresso_args = {
            "input_data": self.make_input_data(),
            "pseudopotentials": self.pseudopotentials,
            "kpts": self.kpts,
            "directory": self.dirname,
            "profile": EspressoProfile(argv=self.espresso_command),
        }
        
        return espresso_args
    
    def check_dict(self, d):
        s = json.dumps(d, indent=2)
        print(s)
        








        # EMAXPOS initialization
        # if emaxpos is 'auto':
        #     emaxpos = computeEMAXPOS()
        # assert len(emaxpos) == len(field_directions), "emaxpos and field_directions must have the same length"

        # Update input data
        # input_data.update(
        #     {
        #         "control": {
        #             "tprnfor": True,
        #             "tefield": True,
        #             "dipfield": True,
        #         },
        #         "system": {
        #             "eamp": eamp,
        #             "eopreg": eopreg,
        #         },
        #     }
        # )

        # self.input_data = input_data
        # s = json.dumps(input_data, indent=2)
        # print(s)
        # "edir": 2,
        # "emaxpos": EMAXPOS,
