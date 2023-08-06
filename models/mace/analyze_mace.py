"""Investigate MACE energy underpredictions."""


# %%
import os

import pandas as pd
from pymatviz import density_scatter, ptable_heatmap_plotly, spacegroup_sunburst

from matbench_discovery import plots as plots
from matbench_discovery.data import df_wbm
from matbench_discovery.preds import PRED_FILES

__author__ = "Janosh Riebesell"
__date__ = "2023-07-23"

module_dir = os.path.dirname(__file__)
id_col = "material_id"
target_col = "e_form_per_atom_mp2020_corrected"
pred_col = "e_form_per_atom_mace"


# %%
df_mace = pd.read_csv(PRED_FILES.MACE).set_index(id_col)
df_mace[list(df_wbm)] = df_wbm

wyckoff_col, spg_col = "wyckoff_spglib", "spacegroup"
df_mace[spg_col] = df_wbm[wyckoff_col].str.split("_").str[2].astype(int)


# %%
density_scatter(df=df_mace, x=target_col, y=pred_col)


# %%
df_bad = df_mace.query(f"{target_col} - {pred_col} > 2")

ax = density_scatter(df=df_bad, x=target_col, y=pred_col)
ax.set(title=f"{len(df_bad):,} MACE severe energy underpredictions")


# %%
fig = ptable_heatmap_plotly(df_bad.formula)
title = f"Elements in {len(df_bad):,} MACE severe energy underpredictions"
fig.layout.title.update(text=title, x=0.4, y=0.95)
fig.show()


# %%
fig = spacegroup_sunburst(df_bad[spg_col], title="MACE spacegroups")
title = f"Spacegroup sunburst of {len(df_bad):,} MACE severe energy underpredictions"
fig.layout.title.update(text=title, x=0.5)
fig.show()


"""
Space groups of MACE underpredictions look unremarkable but unusually heavy in Silicon,
Lanthanide and heavy alkali metals.
"""
