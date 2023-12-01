import itertools
from typing import Any, Dict, List, Optional, Tuple, Union

import ase.data
import ase.data.colors
import numpy as np
from ase.atoms import Atoms
from ase.cell import Cell
from matplotlib.axes import Axes
from matplotlib.collections import PatchCollection, PolyCollection
from matplotlib.patches import Circle
from numpy.typing import NDArray

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.style as mplstyle
from matplotlib.animation import FuncAnimation
from matplotlib.ticker import MultipleLocator, AutoMinorLocator

matplotlib.use("Agg")
mplstyle.use("fast")

"""
SOURCE FROM: https://gitlab.com/agox/agox/-/blob/stable/agox/utils/plot.py
sourced becaused agox needs to be built with more packages, making it large to install.
"""


def plot_atoms(
    ax: Axes,
    atoms: Atoms,
    *,
    plane: str = "xy+",
    radius_factor: float = 1.0,
    repeat: int = 0,
    patch_kwargs: Optional[Dict[str, Any]] = None,
    repeat_patch_kwargs: Optional[Dict[str, Any]] = None,
) -> List[PatchCollection]:
    """Plot a 2D projection of an `Atoms` object into an `Axes` object.

    Parameters
    ----------
    ax : Axes
        Axes to plot the atoms in.
    atoms : Atoms
        Atoms to plot.
    plane : str, optional
        Plane to project the atoms on. See `plane_to_indices` for information
        on valid plane strings, by default 'xy+'.
    radius_factor : float, optional
        Factor to multiply the atomic covalent radii with, by default 1.0.
    repeat : int, optional
        Number of times to repeat the Atoms object in dimensions with periodic
        boundary conditions, by default 0. Repeated copies are not included in
        autoscaling. Note that if none of the two projection axes have periodic
        boundary conditions, this parameter is effectively ignored.
    patch_kwargs : Optional[Dict[str, Any]], optional
        Keyword arguments forwarded to the `Circle` constructor for each atom,
        by default None.
    repeat_patch_kwargs : Optional[Dict[str, Any]], optional
        Keyword arguments additionally forwarded to the `Circle` constructor
        for each repeated atom, by default None.

    Returns
    -------
    List[PatchCollection]
        List of PatchCollections that make up the Atoms object. If the Atoms
        object is not repeated, this list contains one element. If it is
        repeated, this list contains two elements: a patch collection
        containing the non-repeated patches and a patch collection containing
        the repeated patches.
    """

    if patch_kwargs is None:
        patch_kwargs = {}
    if repeat_patch_kwargs is None:
        repeat_patch_kwargs = {}

    d0, d1, d2, kv = plane_to_indices(plane)

    radii = radius_factor * ase.data.covalent_radii[atoms.numbers]
    colors = ase.data.colors.jmol_colors[atoms.numbers]

    if repeat > 0:
        atoms_ = atoms.copy()
        pbc = atoms.get_pbc()

        repeat_range_d0 = range(-repeat, repeat + 1) if pbc[d0] else [0]
        repeat_range_d1 = range(-repeat, repeat + 1) if pbc[d1] else [0]

        for x, y in itertools.product(repeat_range_d0, repeat_range_d1):
            if x == y == 0:
                continue
            offset_uc = np.zeros(3)
            offset_uc[[d0, d1]] = [x, y]
            offset = atoms.cell.cartesian_positions(offset_uc)
            repeat_atoms = atoms.copy()
            repeat_atoms.translate(offset)
            atoms_ += repeat_atoms
    else:
        atoms_ = atoms

    n_atoms = len(atoms)
    patches: List[Circle] = []
    repeat_patches: List[Circle] = []

    for i in np.argsort(kv * atoms_.positions[:, d2]):
        patch_kwargs_ = {**patch_kwargs} if i < n_atoms else {**patch_kwargs, **repeat_patch_kwargs}

        position = atoms_[i].position
        radius = radii[i % n_atoms]
        color = colors[i % n_atoms]

        circle = Circle(position[[d0, d1]], radius=radius, ec="k", fc=color, **patch_kwargs_)

        if i < n_atoms:
            patches.append(circle)
        else:
            repeat_patches.append(circle)

    pcs: List[PatchCollection] = []

    p = PatchCollection(patches, match_original=True)
    ax.add_collection(p)
    pcs.append(p)

    if len(repeat_patches) > 0:
        repeat_p = PatchCollection(repeat_patches, match_original=True)
        ax.add_collection(repeat_p, autolim=False)
        pcs.append(repeat_p)

    ax.set_aspect("equal")
    ax.tick_params(
        bottom=False,
        top=False,
        left=False,
        right=False,
        labelbottom=False,
        labeltop=False,
        labelleft=False,
        labelright=False,
    )

    return pcs


def plot_cell(
    ax: Axes,
    cell: Union[Cell, NDArray],
    *,
    plane: str = "xy+",
    offset: Optional[NDArray] = None,
    collection_kwargs: Dict[str, Any] = None,
):
    """Plot a 2D projection of a `Cell` object into an `Axes` object.

    Parameters
    ----------
    ax : Axes
        Axes to plot the cell in.
    cell : Union[Cell, NDArray]
        Cell to plot.
    plane : str, optional
        Plane to project the cell on. See `plane_to_indices` for information on
        valid plane strings, by default 'xy+'.
    offset : Optional[NDArray], optional
        Cartesian coordinates to apply as offset to the base of the cell, by
        default None (no offset).
    collection_kwargs : Dict[str, Any], optional
        Keyword arguments additionally forwarded to the `PolyCollection`
        constructor for the cell, by default None.

    Raises
    ------
    ValueError
        Raised if `offset` or `cell` has an incorrect shape.
    """

    default_collection_kwargs = dict(edgecolors="k", facecolors="none", linestyles="dotted")

    if offset is None:
        offset = np.zeros(3)
    else:
        offset = np.asarray(offset)

    if offset.shape != (3,):
        raise ValueError(f"`offset` must have shape (3,), input has {offset.shape}")

    if isinstance(cell, Cell):
        cell = cell.complete()

    if cell.shape != (3, 3):
        raise ValueError(f"`cell` must have shape (3, 3), input has {cell.shape}")

    if collection_kwargs is None:
        collection_kwargs = {}

    collection_kwargs = {**default_collection_kwargs, **collection_kwargs}

    d0, d1, _, _ = plane_to_indices(plane)

    verts_uc = np.array(
        [
            [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]],  # xy, z = 0
            [[0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]],  # xy, z = 1
            [[0, 0, 0], [1, 0, 0], [1, 0, 1], [0, 0, 1]],  # xz, y = 0
            [[0, 1, 0], [1, 1, 0], [1, 1, 1], [0, 1, 1]],  # xz, y = 1
            [[0, 0, 0], [0, 1, 0], [0, 1, 1], [0, 0, 1]],  # yz, x = 0
            [[1, 0, 0], [1, 1, 0], [1, 1, 1], [1, 0, 1]],  # yz, x = 1
        ]
    )
    verts_3d = verts_uc @ cell + offset
    verts_2d = verts_3d[:, :, [d0, d1]]

    p = PolyCollection(verts_2d, **collection_kwargs)
    ax.add_collection(p)


def plane_to_indices(plane: str) -> Tuple[int, int, int, int]:
    """Convert a string representation of a plane to a 4-tuple of indices
    usable to perform operations using this plane as a projection.

    Parameters
    ----------
    plane : str
        String representation of the plane, in the format /[xyz]{2}[+-]?/. The
        first two characters define the plane's Cartesian axes (and should not
        be equal), and the third optional character defines whether the
        perpendicular axis is positive (+) or negative (-); the default is
        positive.

        Examples of string representations: 'xy', 'xy+', 'xy-', 'xz', 'zx+'.

    Returns
    -------
    Tuple[int, int, int, int]
        A tuple (d0, d1, d2, kv) describing this plane:
        - d0: index of the first plane axis (0, 1, or 2);
        - d1: index of the second plane axis (0, 1, or 2);
        - d2: index of the axis perpendicular to the plane (0, 1, or 2);
        - kv: sign of the axis perpendicular to the plane (-1 or 1).

    Raises
    ------
    ValueError
        Raised if the plane specification is invalid.
    """
    if not isinstance(plane, str):
        plane = str(plane)

    if len(plane) < 2:
        raise ValueError(f"Invalid plane specification {plane}")

    I = np.eye(3)

    d0 = ord(plane[0]) - ord("x")  # 'x' -> 0, 'y' -> 1, 'z' -> 2
    d1 = ord(plane[1]) - ord("x")  # 'x' -> 0, 'y' -> 1, 'z' -> 2

    if len(plane) > 2:
        s = -(ord(plane[2]) - ord(","))  # '+' -> 1, '-' -> -1
    else:
        s = 1

    if d0 not in [0, 1, 2] or d1 not in [0, 1, 2] or d0 == d1 or s not in [-1, 1]:
        raise ValueError(f"Invalid plane specification {plane}")

    d2 = ({0, 1, 2} - {d0} - {d1}).pop()

    i, j = I[d0, :], I[d1, :]
    k = np.cross(i, j)
    kv = s * int(np.sum(k))  # -1 or 1

    return (d0, d1, d2, kv)


def plot_function(savename, images, title):
    """
    Plot function.

    Parameters
    ----------
    savename : str
        Path to the file where the plot will be saved.
    images : list
        List of images.
    title : str
        Title of the plot.
    """
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12, 5))

    fig.suptitle(title)

    init_image = images[0]

    low_x, low_y, low_z = np.inf, np.inf, np.inf
    high_x, high_y, high_z = -np.inf, -np.inf, -np.inf
    for image in images:
        image.get_positions()
        low_x = min(low_x, image.get_positions()[:, 0].min())
        low_y = min(low_y, image.get_positions()[:, 1].min())
        low_z = min(low_z, image.get_positions()[:, 2].min())
        high_x = max(high_x, image.get_positions()[:, 0].max())
        high_y = max(high_y, image.get_positions()[:, 1].max())
        high_z = max(high_z, image.get_positions()[:, 2].max())

    init_cell = init_image.get_cell()
    cell_x_max = max(init_cell[:, 0].sum(), init_cell[:, 0].max())
    cell_y_max = max(init_cell[:, 1].sum(), init_cell[:, 1].max())
    cell_z_max = max(init_cell[:, 2].sum(), init_cell[:, 2].max())

    maxhigh_x = max(high_x, cell_x_max)
    maxhigh_y = max(high_y, cell_y_max)
    maxhigh_z = max(high_z, cell_z_max)

    minlow_x = min(low_x, 0)
    minlow_y = min(low_y, 0)
    minlow_z = min(low_z, 0)

    max_rad = ase.data.covalent_radii[init_image.numbers].max() * 1.2

    def xy_set_lim():
        ax1.set_xlim(minlow_x - max_rad, maxhigh_x + max_rad)
        ax1.set_ylim(minlow_y - max_rad, maxhigh_y + max_rad)

    def xy_update(image):
        ax1.cla()
        # ax1.set_aspect("equal")
        xy_set_lim()
        plot_atoms(ax=ax1, atoms=image)
        plot_cell(ax=ax1, cell=image.get_cell())
        ax1.set_ylabel("$y$-axis")
        ax1.set_xlabel("$x$-axis")

    def yz_set_lim():
        ax2.set_xlim(low_y - max_rad, high_y + max_rad)
        ax2.set_ylim(low_z - max_rad, high_z + max_rad)

    def yz_update(image):
        ax2.cla()
        # ax2.set_aspect("equal")
        yz_set_lim()
        plot_atoms(ax=ax2, atoms=image, plane="yz")
        plot_cell(ax=ax2, cell=image.get_cell(), plane="yz")
        ax2.set_ylabel("$z$-axis")
        ax2.set_xlabel("$y$-axis")

    def zx_set_lim():
        ax3.set_xlim(low_z - max_rad, high_z + max_rad)
        ax3.set_ylim(low_x - max_rad, high_x + max_rad)

    def zx_update(image):
        ax3.cla()
        # ax3.set_aspect("equal")
        zx_set_lim()
        plot_atoms(ax=ax3, atoms=image, plane="zx")
        plot_cell(ax=ax3, cell=image.get_cell(), plane="zx")
        ax3.set_ylabel("$x$-axis")
        ax3.set_xlabel("$z$-axis")

    def update(image):
        xy_update(image)
        yz_update(image)
        zx_update(image)

    ani = FuncAnimation(fig, update, frames=images[::2], interval=1000 / 30)

    ani.save(savename, dpi=150, writer="ffmpeg")

    plt.close(fig)


def plot_spectra(
    figsize,
    fontsize,
    dpi,
    spectra_xy,
    spectra_label,
    exact_modes_x,
    exact_modes_y,
    exact_modes_label,
    exact_modes_color,
    xlabel,
    ylabel,
    xrange,
    title,
    plotfile,
):
    fig, ax = plt.subplots(figsize=figsize)
    ax.fill_between(
        spectra_xy[0],
        spectra_xy[1],
        label=spectra_label,
        color="k",
        alpha=0.8,
    )
    ax.vlines(
        x=exact_modes_x,
        ymin=0,
        ymax=exact_modes_y,
        label=exact_modes_label,
        color=exact_modes_color,
        alpha=1,
        linewidth=1,
    )

    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.set_ylabel(ylabel, fontsize=fontsize)
    ax.tick_params(axis="x", labelsize=fontsize)
    ax.tick_params(axis="y", labelsize=fontsize)

    low, high = ax.get_ybound()
    freq_limit = xrange
    ax.set_ylim([low, high * 1.2])
    ax.set_xlim(freq_limit)

    ax.set_title(title)
    ax.xaxis.set_minor_locator(MultipleLocator(100))

    fig.tight_layout()

    fig.savefig(plotfile, dpi=dpi)
