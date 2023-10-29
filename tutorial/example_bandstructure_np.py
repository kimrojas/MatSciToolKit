from matscitoolkit import gnu_parser
from matscitoolkit import qe_tools

import matplotlib.pyplot as plt

gnufile = "../test/test_data/bandstructure_np/si_bands.dat.gnu"
bandsout = "../test/test_data/bandstructure_np/espresso_bands.out"
pwscfout = "../test/test_data/bandstructure_np/espresso_pwscf.out"

# --- High symmetry points (user-defined)
labels = ["L", "G", "X", "U", "G"]

# --- Fermi energy
efermi = qe_tools.get_fermi(pwscfout)

# --- GNU data
gnudata = gnu_parser.read_gnu_datafile(gnufile)
print("GNU data shape: ", gnudata.shape)

# --- QE data
xcoord = qe_tools.get_high_sympoints(bandsout)

# --- Plot band structure
fig, ax = plt.subplots()
## plot bands
for band in range(gnudata.shape[0]):
    ax.plot(gnudata[band, :, 0], gnudata[band, :, 1], "k-", linewidth=1, alpha=0.7)
## plot high symmetry points
for x in xcoord:
    ax.axvline(x, color="k", linestyle="--", linewidth=0.5)
## plot labels using xticks
ax.set_xticks(xcoord, labels=labels)

ax.set_ylabel("Energy (eV)")

plt.show()
