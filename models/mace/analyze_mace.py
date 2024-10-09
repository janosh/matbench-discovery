"""Investigate MACE energy underpredictions."""

# %%
import os

import pandas as pd
import pymatviz as pmv
from pymatviz.enums import Key

from matbench_discovery import SITE_FIGS
from matbench_discovery import plots as plots
from matbench_discovery.data import Model, df_wbm
from matbench_discovery.enums import MbdKey

__author__ = "Janosh Riebesell"
__date__ = "2023-07-23"

module_dir = os.path.dirname(__file__)
e_form_mace_col = "e_form_per_atom_mace"


# %%
df_mace = pd.read_csv(Model.mace.path).set_index(Key.mat_id)
df_mace[list(df_wbm)] = df_wbm

df_mace[Key.spg_num] = df_wbm[MbdKey.init_wyckoff].str.split("_").str[2].astype(int)


# %%
ax = pmv.density_scatter(df=df_mace, x=MbdKey.e_form_dft, y=e_form_mace_col)
ax.set(title=f"{len(df_mace):,} MACE severe energy underpredictions")
pmv.save_fig(ax, "mace-hull-dist-scatter.pdf")


# %%
df_low = df_mace.query(f"{MbdKey.e_form_dft} - {e_form_mace_col} > 2")

ax = pmv.density_scatter(df=df_low, x=MbdKey.e_form_dft, y=e_form_mace_col)
ax.set(title=f"{len(df_low):,} MACE severe energy underpredictions")
pmv.save_fig(ax, "mace-too-low-hull-dist-scatter.pdf")


# %%
fig = pmv.ptable_heatmap_plotly(df_low[Key.formula])
title = f"Elements in {len(df_low):,} MACE severe energy underpredictions"
fig.layout.title.update(text=title, x=0.4, y=0.95)
fig.show()

pmv.save_fig(fig, "mace-too-low-elements-heatmap.pdf")


# %%
fig = pmv.spacegroup_sunburst(df_low[Key.spg_num], title="MACE spacegroups")
title = f"Spacegroup sunburst of {len(df_low):,} MACE severe energy underpredictions"
fig.layout.title.update(text=title, x=0.5)
fig.show()

pmv.save_fig(fig, "mace-too-low-spacegroup-sunburst.pdf")


"""
Space groups of MACE underpredictions look unremarkable but unusually heavy in Silicon,
Lanthanide and heavy alkali metals.
"""

# %%
bad_mask = (df_mace[e_form_mace_col] - df_mace[MbdKey.e_form_dft]) < -5
print(f"{sum(bad_mask)=}")

fig = pmv.density_scatter_plotly(
    df_mace[~bad_mask],
    x=MbdKey.e_form_dft,
    y=e_form_mace_col,
    log_density=(log := True),
)
fig.layout.yaxis.title = MbdKey.e_form_dft.replace("DFT", "MACE")
fig.show()
pmv.save_fig(fig, f"{SITE_FIGS}/mace-wbm-IS2RE-e-form-parity.svelte")


# %%
fig = pmv.density_scatter_plotly(
    df_mace[~bad_mask], x="uncorrected_energy", y="mace_energy", log_density=log
)
fig.layout.yaxis.title = "MACE energy"

fig.show()
pmv.save_fig(fig, f"{SITE_FIGS}/mace-wbm-IS2RE-raw-energy-parity.svelte")
