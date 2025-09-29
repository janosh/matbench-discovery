"""
Investigate Allegro energy underpredictions, templated from the
MACE analysis script from Janosh.
"""

# uses matbench-discovery matbench-discovery commit ID 012ccfe,
# k_srme commit ID 0269a946, pymatviz v0.15.1

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
e_form_allegro_col = "e_form_per_atom_allegro"

filtered_csv_path = "./results/allegro-filtered_preds.csv.gz"
if not os.path.isfile(filtered_csv_path):
    filtered_csv_path = "./allegro-filtered_preds.csv.gz"

df_allegro = pd.read_csv(filtered_csv_path)
df_allegro = df_allegro.set_index("material_id")
df_wbm["e_form_per_atom_allegro_error"] = (
    df_wbm["e_form_per_atom_mp2020_corrected"] - df_allegro["e_form_per_atom_allegro"]
)

df_allegro[list(df_wbm)] = df_wbm

df_allegro[Key.spg_num] = (
    df_wbm[MbdKey.init_wyckoff_spglib].str.split("_").str[2].astype(int)
)


fig = pmv.density_scatter(df=df_allegro, x=MbdKey.e_form_dft, y=e_form_allegro_col)
fig.layout.title = f"{len(df_allegro):,} Allegro severe energy underpredictions"
pmv.save_fig(fig, "allegro-hull-dist-scatter.png")


df_low = df_allegro.query(f"{MbdKey.e_form_dft} - {e_form_allegro_col} > 2")

fig = pmv.density_scatter(df=df_low, x=MbdKey.e_form_dft, y=e_form_allegro_col)
df = fig.layout.title = f"{len(df_low):,} Allegro severe energy underpredictions"
pmv.save_fig(fig, "allegro-too-low-hull-dist-scatter.png")


fig = pmv.ptable_heatmap_plotly(df_low[Key.formula])
title = f"Elements in {len(df_low):,} Allegro severe energy underpredictions"
fig.layout.title.update(text=title, x=0.4, y=0.95)
fig.show()

pmv.save_fig(fig, "allegro-too-low-elements-heatmap.png")


fig = pmv.spacegroup_sunburst(df_low[Key.spg_num], title="Allegro spacegroups")
title = f"Spacegroup sunburst of {len(df_low):,} Allegro severe energy underpredictions"
fig.layout.title.update(text=title, x=0.5)
fig.show()

pmv.save_fig(fig, "allegro-too-low-spacegroup-sunburst.png")


bad_mask = (df_allegro[e_form_allegro_col] - df_allegro[MbdKey.e_form_dft]) < -5
print(f"{sum(bad_mask)=}")

fig = pmv.density_scatter(
    df=df_allegro[~bad_mask],
    x=MbdKey.e_form_dft,
    y=e_form_allegro_col,
    log_density=(log := True),
)
fig.layout.yaxis.title = MbdKey.e_form_dft.replace("DFT", "Allegro")
fig.show()
pmv.save_fig(fig, f"{SITE_FIGS}/allegro-wbm-IS2RE-e-form-parity.png")


print(df_allegro.columns)
fig = pmv.density_scatter(
    df=df_allegro[~bad_mask],
    x="uncorrected_energy",
    y="allegro_energy",
    log_density=log,
)
fig.layout.yaxis.title = "Allegro energy"

fig.show()
pmv.save_fig(fig, f"{SITE_FIGS}/allegro-wbm-IS2RE-raw-energy-parity.png")
