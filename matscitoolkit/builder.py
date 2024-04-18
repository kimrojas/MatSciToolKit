from ase.io import read, write
from ase import Atoms
from pathlib import Path
import numpy as np
from ase.visualize import view
from ase.data import covalent_radii


class AddAdsorbateMethod:
    def __init__(self):
        self.substrate = None
        self.adsorbate = None

        # -- Adsorbate properties --
        self.adsorbate_origin = None
        self.adsorbate_rotation = None

        # -- Substrate properties --
        self.substrate_reference = None

        # -- Method properties --
        self.height = None

        # -- Output --
        self.structure_complex = None

    def set_adsorbate(self, adsorbate, adsorbate_origin="COG", adsorbate_rotation=[]):
        """
        Sets the adsorbate for the builder.

        Parameters:
            adsorbate (str or Path or Atoms): The adsorbate to set. It can be a path to a file, a string representing a file path, or an Atoms object.
            adsorbate_origin (str or int or list, optional): The origin of the adsorbate. It can be one of the following:
                - "COM" (default): Center of mass of the adsorbate.
                - "COG": Center of geometry of the adsorbate.
                - int: User-defined atom index representing the origin.
                - list: List of user-defined atom indices representing the center of all indices.
            adsorbate_rotation (list or tuple, optional): The rotation of the adsorbate. Each rotation is specified as a tuple with the format (angle, axis vector).
                - if the list is empty, no rotation is applied.
                - if the list is not empty, each tuple should have two elements: angle and axis vector.
                - Example: [[90, [1, 0, 0]], [180, [0, 1, 0]]] will rotate the adsorbate 90 degrees around the x-axis and THEN 180 degrees around the y-axis.
                - Example: [[90, 'x'], [180, 'y']] will rotate the adsorbate 90 degrees around the x-axis and THEN 180 degrees around the y-axis.

        Raises:
            ValueError: If the adsorbate is not an Atoms object or a valid file path.

        Returns:
            None
        """

        # -- Set adsorbate --
        if isinstance(adsorbate, str) or isinstance(adsorbate, Path):
            self.adsorbate = read(adsorbate)
        elif isinstance(adsorbate, Atoms):
            self.adsorbate = adsorbate
        else:
            raise ValueError("Adsorbate should be an Atoms object or a path to a file.")

        # -- Set adsorbate origin --
        if adsorbate_origin == "COM":  # default to center of mass
            self.adsorbate_origin = self.adsorbate.get_center_of_mass()
        if adsorbate_origin == "COG":  # center of geometry
            self.adsorbate_origin = self.adsorbate.get_positions().mean(axis=0)
        if isinstance(adsorbate_origin, int):  # user defined atom index
            self.adsorbate_origin = self.adsorbate.get_positions()[adsorbate_origin]
        if isinstance(adsorbate_origin, list):  # center of all user defined atom index
            stack = np.array([self.adsorbate.get_positions()[i] for i in adsorbate_origin])
            self.adsorbate_origin = stack.mean(axis=0)

        self.adsorbate.translate(-self.adsorbate_origin)

        # -- Set adsorbate rotation --
        if isinstance(adsorbate_rotation, list) or isinstance(adsorbate_rotation, tuple):
            for irot in adsorbate_rotation:
                assert len(irot) == 2, "Rotation should be a list of tuples with the format (angle, axis vector)"
                angle, vector_axis = irot
                self.adsorbate.rotate(a=angle, v=vector_axis)

    def set_substrate(self, substrate, substrate_reference):
        """
        Set the substrate and substrate reference for the builder.

        Parameters:
            substrate (str or Path or Atoms): The substrate to set. It can be either a path to a file, a string representing the path, or an Atoms object.
            substrate_reference (int or list): The reference position(s) on the substrate. It can be either an integer representing the index of a single position or a list of integers representing multiple positions.

        Raises:
            ValueError: If the substrate is not an Atoms object or a valid path to a file.

        """
        # -- Set substrate --
        if isinstance(substrate, str) or isinstance(substrate, Path):
            self.substrate = read(substrate)
        elif isinstance(substrate, Atoms):
            self.substrate = substrate
        else:
            raise ValueError("Substrate should be an Atoms object or a path to a file.")

        # -- Set substrate reference --
        if isinstance(substrate_reference, int):
            self.substrate_reference = self.substrate.get_positions()[substrate_reference]
        if isinstance(substrate_reference, list):
            stack = np.array([self.substrate.get_positions()[i] for i in substrate_reference])
            self.substrate_reference = stack.mean(axis=0)

    def set_height(self, height):
        """
        Set the height of the object.

        Parameters:
        height (float or int): The height value to be set.

        Raises:
        ValueError: If the height is not a float or an integer.
        """
        if isinstance(height, float) or isinstance(height, int):
            self.height = float(height)
        else:
            raise ValueError("Height should be a float or an integer.")

    def build(self):
        _adsorbate = self.adsorbate.copy()
        _substrate = self.substrate.copy()
        
        # -- Align origin to reference --
        _adsorbate.translate(self.substrate_reference)

        # -- Adjust height --
        _adsorbate.translate([0, 0, self.height])

        # -- Combine substrate and adsorbate --
        self.substrate_complex = _substrate + _adsorbate
        
        # -- Check for ionic overlap --
        warn = self.warn_overlap(_adsorbate, _substrate)
        return warn

    def get_structure(self):
        return self.substrate_complex.copy()
    
    def warn_overlap(self, ads, sub, scale=0.8):
        ads_pos, ads_sym, ads_num = ads.get_positions(), ads.get_chemical_symbols(), ads.get_atomic_numbers()
        sub_pos, sub_sym, sub_num = sub.get_positions(), sub.get_chemical_symbols(), sub.get_atomic_numbers()
        
        warn_status = False
        
        for iads, (iads_pos, iads_sym, iads_num) in enumerate(zip(ads_pos, ads_sym, ads_num)):
            for isub, (isub_pos, isub_sym, isub_num) in enumerate(zip(sub_pos, sub_sym, sub_num)):
                if np.linalg.norm(iads_pos - isub_pos) < scale * (covalent_radii[iads_num] + covalent_radii[isub_num]):
                    print(f"  Warning: Adsorbate atom [{iads_sym:>2}_{iads:<3d}] and Substrate atom [{isub_sym:>2}_{isub:<3d}] are too close.")
                    warn_status = True
                    
        return warn_status

def add_adsorbate(
    adsorbate,
    substrate,
    height,
    substrate_reference,
    adsorbate_origin="COG",
    adsorbate_rotation=[],
):
    obj = AddAdsorbateMethod()
    obj.set_adsorbate(adsorbate, adsorbate_origin, adsorbate_rotation)
    obj.set_substrate(substrate, substrate_reference)
    
    if isinstance(height, float) or isinstance(height, int):
        height = [float(height)]
    
    images = []
    for h in height:
        print(f"Building structure at height: {h}")
        obj.set_height(h)
        obj.build()
        structure = obj.get_structure()
        images.append(structure)
    
    if len(images) == 1:
        return images[0]
    else:
        return images
