from matscitoolkit.ase_ir_tools import generateDisplacedStructures

# In this step, we generate all the displaced structures.

structure_filepath = "rlx_H2O.vasp"
generateDisplacedStructures(structure_filepath)