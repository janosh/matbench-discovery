"""
Investigate NequIP energy underpredictions, templated from the
MACE analysis script from Janosh.
"""

# uses matbench-discovery matbench-discovery commit ID 012ccfe,
# k_srme commit ID 0269a946, pymatviz v0.15.1


# %%
import os

import pandas as pd
import pymatviz as pmv
from pymatviz.enums import Key

from matbench_discovery import SITE_FIGS
from matbench_discovery import plots as plots
from matbench_discovery.data import df_wbm
from matbench_discovery.enums import MbdKey

__author__ = "Janosh Riebesell, Seán Kavanagh"

module_dir = os.path.dirname(__file__)
e_form_nequip_col = "e_form_per_atom_nequip"

filtered_csv_path = "./results/nequip-filtered_preds.csv.gz"
if not os.path.exists(filtered_csv_path):
    filtered_csv_path = "./nequip-filtered_preds.csv.gz"

df_nequip = pd.read_csv(filtered_csv_path)
df_nequip = df_nequip.set_index("material_id")
df_wbm["e_form_per_atom_nequip_error"] = (
    df_wbm["e_form_per_atom_mp2020_corrected"] - df_nequip["e_form_per_atom_nequip"]
)

df_nequip[list(df_wbm)] = df_wbm

df_nequip[Key.spg_num] = (
    df_wbm[MbdKey.init_wyckoff_spglib].str.split("_").str[2].astype(int)
)


# %%
fig = pmv.density_scatter_plotly(df=df_nequip, x=MbdKey.e_form_dft, y=e_form_nequip_col)
fig.layout.title = f"{len(df_nequip):,} Nequipsevere energy underpredictions"
pmv.save_fig(fig, "nequip-hull-dist-scatter.png")


# %%
df_low = df_nequip.query(f"{MbdKey.e_form_dft} - {e_form_nequip_col} > 2")

fig = pmv.density_scatter_plotly(df=df_low, x=MbdKey.e_form_dft, y=e_form_nequip_col)
fig.layout.title = f"{len(df_low):,} Nequip severe energy underpredictions"
pmv.save_fig(fig, "nequip-too-low-hull-dist-scatter.png")


# %%
fig = pmv.ptable_heatmap_plotly(df_low[Key.formula])
title = f"Elements in {len(df_low):,} Nequip severe energy underpredictions"
fig.layout.title.update(text=title, x=0.4, y=0.95)
fig.show()

pmv.save_fig(fig, "nequip-too-low-elements-heatmap.png")


# %%
fig = pmv.spacegroup_sunburst(df_low[Key.spg_num], title="Nequip spacegroups")
title = f"Spacegroup sunburst of {len(df_low):,} Nequip severe energy underpredictions"
fig.layout.title.update(text=title, x=0.5)
fig.show()

pmv.save_fig(fig, "nequip-too-low-spacegroup-sunburst.png")


# %%
bad_mask = (df_nequip[e_form_nequip_col] - df_nequip[MbdKey.e_form_dft]) < -5
print(f"{sum(bad_mask)=}")

fig = pmv.density_scatter_plotly(
    df_nequip[~bad_mask],
    x=MbdKey.e_form_dft,
    y=e_form_nequip_col,
    log_density=(log := True),
)
fig.layout.yaxis.title = MbdKey.e_form_dft.replace("DFT", "Nequip")
fig.show()
pmv.save_fig(fig, f"{SITE_FIGS}/nequip-wbm-IS2RE-e-form-parity.png")


# %%
print(df_nequip.columns)
fig = pmv.density_scatter_plotly(
    df_nequip[~bad_mask], x="uncorrected_energy", y="nequip_energy", log_density=log
)
fig.layout.yaxis.title = "Nequip energy"

fig.show()
pmv.save_fig(fig, f"{SITE_FIGS}/nequip-wbm-IS2RE-raw-energy-parity.png")
