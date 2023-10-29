from matscitoolkit import gnu_parser
from matscitoolkit import qe_tools
from matscitoolkit import plot_tools

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines

gnufile_up = "../test/test_data/bandstructure_sp/fe_bands_up.dat.gnu"
gnufile_dn = "../test/test_data/bandstructure_sp/fe_bands_dn.dat.gnu"
bandsout = "../test/test_data/bandstructure_sp/espresso_bands_up.out"
pwscfout = "../test/test_data/bandstructure_sp/espresso_pwscf.out"

# --- High symmetry points (user-defined)
labels = ["G", "H", "N", "G", "P", "H", "P", "N"]

# --- Fermi energy
efermi = qe_tools.get_fermi(pwscfout)

# --- GNU data
gnudata_up = gnu_parser.read_gnu_datafile(gnufile_up)
print("GNU data up shape: ", gnudata_up.shape)
gnudata_dn = gnu_parser.read_gnu_datafile(gnufile_dn)
print("GNU data dn shape: ", gnudata_dn.shape)

# --- QE data
xcoord = qe_tools.get_high_sympoints(bandsout)

# --- Plot band structure
fig, (ax1, ax2) = plt.subplots(1, 2, sharex=True, sharey=True)

## plot bands
for i, (ax, gnudata, color, label) in enumerate(
    zip([ax1, ax2], [gnudata_up, gnudata_dn], ["k", "r"], ["Spin-up", "Spin-dn"])
):
    for band in range(gnudata.shape[0]):
        ax.plot(
            gnudata[band, :, 0],
            gnudata[band, :, 1],
            color=color,
            linewidth=1,
            alpha=0.7,
            label=label,
        )

for ax in (ax1, ax2):

    ## plot high symmetry points
    for x in xcoord:
        ax.axvline(x, color="k", linestyle="--", linewidth=0.5)

    ## plot labels using xticks
    processed_locations, processed_labels = plot_tools.tick_combiner(xcoord, labels)
    ax1.set_xticks(processed_locations, labels=processed_labels)

    ## plot fermi energy
    ax.axhline(efermi, color="b", linestyle="--", linewidth=0.5, alpha=0.5)

ax1.set_ylabel("Energy (eV)")
ax1.set_ylim(efermi - 10, efermi + 10)

ax1.set_title("Spin-up")
ax2.set_title("Spin-dn")

plt.show()
