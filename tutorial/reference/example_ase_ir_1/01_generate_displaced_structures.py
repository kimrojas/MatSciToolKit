from matscitoolkit.ase_ir_tools import generateDisplacedStructures
from matscitoolkit.ase_ir_tools import computeEMAXPOS

# In this step, we generate all the displaced structures.
structure_filepath = "hb_unitcell.poscar"
generateDisplacedStructures(structure_filepath)


structure_filepath = "hb_unitcell.poscar"
computeEMAXPOS(structure_filepath, debug=True)