# %%
import numpy as np
import pandas as pd
import pymatviz as pmv
from pymatgen.core import Lattice, Structure

from matbench_discovery.structure import perturb_structure

rng = np.random.default_rng(seed=0)
pmv.set_plotly_template("pymatviz_dark")


# %%
fig = pd.Series(rng.weibull(1.5, 100_000)).hist(bins=100, backend="plotly")
title = "Distribution of perturbation magnitudes"
fig.layout.update(xaxis_title="Perturbation Magnitude", title=title)
fig.show()


# %%
struct = Structure(
    lattice=Lattice.cubic(5),
    species=("Fe", "Fe"),
    coords=((0, 0, 0), (0.5, 0.5, 0.5)),
)

fig = pmv.structure_2d_plotly(struct)
fig.layout.update(title=f"Original structure: {struct.formula}")
fig.show()


# %%
pmv.structure_2d_plotly(
    [perturb_structure(struct) for _ in range(12)],
    subplot_title=lambda _struct, idx: f"perturbation {idx}",
)
