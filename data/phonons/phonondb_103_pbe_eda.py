"""Exploratory data analysis of the 103 structures in Togo's phononDB PBE dataset.

Used for the anharmonic phonon analysis, specifically the thermal conductivity kappa.

See https://github.com/atztogo/phonondb/blob/bba206/README.md#url-links-to-phono3py-
finite-displacement-method-inputs-of-103-compounds-on-mdr-at-nims-pbe for details.
"""

# %%
from collections import defaultdict

import ase.io
import moyopy
import moyopy.interface
import pymatviz as pmv

from matbench_discovery.data import DataFiles

__date__ = "2025-01-14"


# %%
atoms_list = ase.io.read(DataFiles.phonondb_pbe_structures.path, index=":")


# %% visually inspect first 12 structures
fig = pmv.structure_3d_plotly(atoms_list[:12], n_cols=3, scale=0.5)
fig.show()


# %%
elem_counts: dict[str, int] = defaultdict(int)
for atoms in atoms_list:
    for symb in atoms.symbols:
        elem_counts[symb] += 1


# %%
fig = pmv.ptable_heatmap_plotly(elem_counts, fmt=".0f")
fig.show()


# %% plot spacegroup distribution
spg_nums: dict[str, int] = {}
for atoms in atoms_list:
    moyo_cell = moyopy.interface.MoyoAdapter.from_atoms(atoms).data
    moyo_data = moyopy.MoyoDataset(moyo_cell)
    spg_nums[atoms.info["material_id"]] = moyo_data.number

fig = pmv.spacegroup_sunburst(spg_nums.values(), show_counts="value+percent")
fig.show()
