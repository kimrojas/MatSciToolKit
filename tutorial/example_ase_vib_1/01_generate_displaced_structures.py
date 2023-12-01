from matscitoolkit.ase_vib_tools import generateDisplacedStructures

# In this step, we generate all the displaced structures.

structure_filepath = "hb_unitcell.poscar"
generateDisplacedStructures(structure_filepath)