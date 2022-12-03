# %%
import numpy as np
import pandas as pd
from pymatgen.core import Lattice, Structure
from pymatviz import plot_structure_2d

from matbench_discovery.plots import plt
from matbench_discovery.structure import perturb_structure

__author__ = "Janosh Riebesell"
__date__ = "2022-12-02"


# %%
ax = pd.Series(np.random.weibull(1.5, 100000)).hist(bins=100)
title = "Distribution of perturbation magnitudes"
ax.set(xlabel="magnitude of perturbation", ylabel="count", title=title)


# %%
struct = Structure(
    lattice=Lattice.cubic(5),
    species=("Fe", "O"),
    coords=((0, 0, 0), (0.5, 0.5, 0.5)),
)

ax = plot_structure_2d(struct)
ax.set(title=f"Original structure: {struct.formula}")
ax.set_aspect("equal")


# %%
fig, axs = plt.subplots(3, 4, figsize=(12, 10))
for idx, ax in enumerate(axs.flat, 1):
    plot_structure_2d(perturb_structure(struct), ax=ax)
    ax.set(title=f"perturbation {idx}")
