from matscitoolkit.ase_vib_tools import generateDisplacedStructures
from matscitoolkit.ase_ir_tools import computeEMAXPOS

# In this step, we generate all the displaced structures.
structure_filepath = "hb_unitcell.poscar"
generateDisplacedStructures(structure_filepath)


# This next step is not necessary as it is automatically handled by the dft interface
# but for clarity, lets print this out.
structure_filepath = "hb_unitcell.poscar"
computeEMAXPOS(structure_filepath, debug=True)
