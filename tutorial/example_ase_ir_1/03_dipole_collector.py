from matscitoolkit.ase_ir_tools import DipoleCollector

"""
This script collects the dipole moment from the DFT calculation.
It also summarizes these into the cache file that ASE can read. 

As you can see, everything is automatically handled by the code.
"""


x = DipoleCollector()
x.collect()