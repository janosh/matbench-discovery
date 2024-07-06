"""Investigate MACE energy underpredictions."""

# %%
import os

import pandas as pd
from pymatviz import density_scatter, ptable_heatmap_plotly, spacegroup_sunburst
from pymatviz.enums import Key
from pymatviz.io import save_fig

from matbench_discovery import plots as plots
from matbench_discovery.data import df_wbm
from matbench_discovery.enums import MbdKey
from matbench_discovery.preds import PRED_FILES

__author__ = "Janosh Riebesell"
__date__ = "2023-07-23"

module_dir = os.path.dirname(__file__)
pred_col = "e_form_per_atom_mace"


# %%
df_mace = pd.read_csv(PRED_FILES.MACE).set_index(Key.mat_id)
df_mace[list(df_wbm)] = df_wbm

df_mace[Key.spg_num] = df_wbm[MbdKey.init_wyckoff].str.split("_").str[2].astype(int)


# %%
ax = density_scatter(df=df_mace, x=MbdKey.e_form, y=pred_col)
ax.set(title=f"{len(df_mace):,} MACE severe energy underpredictions")
save_fig(ax, "mace-hull-dist-scatter.pdf")


# %%
df_low = df_mace.query(f"{MbdKey.e_form} - {pred_col} > 2")

ax = density_scatter(df=df_low, x=MbdKey.e_form, y=pred_col)
ax.set(title=f"{len(df_low):,} MACE severe energy underpredictions")
save_fig(ax, "mace-too-low-hull-dist-scatter.pdf")


# %%
fig = ptable_heatmap_plotly(df_low[Key.formula])
title = f"Elements in {len(df_low):,} MACE severe energy underpredictions"
fig.layout.title.update(text=title, x=0.4, y=0.95)
fig.show()

save_fig(fig, "mace-too-low-elements-heatmap.pdf")


# %%
fig = spacegroup_sunburst(df_low[Key.spg_num], title="MACE spacegroups")
title = f"Spacegroup sunburst of {len(df_low):,} MACE severe energy underpredictions"
fig.layout.title.update(text=title, x=0.5)
fig.show()

save_fig(fig, "mace-too-low-spacegroup-sunburst.pdf")


"""
Space groups of MACE underpredictions look unremarkable but unusually heavy in Silicon,
Lanthanide and heavy alkali metals.
"""
