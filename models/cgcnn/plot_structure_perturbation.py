# %%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pymatviz as pmv
from pymatgen.core import Lattice, Structure

from matbench_discovery.structure import perturb_structure

__author__ = "Janosh Riebesell"
__date__ = "2022-12-02"

rng = np.random.default_rng(0)


# %%
ax = pd.Series(rng.weibull(1.5, 100_000)).hist(bins=100)
title = "Distribution of perturbation magnitudes"
ax.set(xlabel="magnitude of perturbation", ylabel="count", title=title)


# %%
struct = Structure(
    lattice=Lattice.cubic(5),
    species=("Fe", "O"),
    coords=((0, 0, 0), (0.5, 0.5, 0.5)),
)

ax = pmv.structure_2d(struct)
ax.set(title=f"Original structure: {struct.formula}")
ax.set_aspect("equal")


# %%
fig, axs = plt.subplots(3, 4, figsize=(12, 10))
for idx, ax in enumerate(axs.flat, start=1):
    pmv.structure_2d(perturb_structure(struct), ax=ax)
    ax.set(title=f"perturbation {idx}")
