from matscitoolkit.ase_ir_tools import PostProcess

"""
In this part, we simply use ASE to read the cache file and make a a vibrations data

Also, we can already plot it. 
"""

x = PostProcess(nproc=4)
x.generate_summary()
x.generate_spectra()
x.load_summary()
x.load_spectra()
x.generate_mode_traj()
x.generate_mode_gif('all')
x.generate_spectra_plot()

