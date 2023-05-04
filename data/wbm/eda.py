# %%
import os

import numpy as np
import pandas as pd
import plotly.express as px
from pymatgen.core import Composition
from pymatviz import (
    count_elements,
    ptable_heatmap_plotly,
    spacegroup_sunburst,
)
from pymatviz.utils import save_fig

from matbench_discovery import FIGS, ROOT, today
from matbench_discovery import plots as plots
from matbench_discovery.data import DATA_FILES, df_wbm
from matbench_discovery.energy import mp_elem_reference_entries
from matbench_discovery.preds import df_each_err, each_true_col

__author__ = "Janosh Riebesell"
__date__ = "2023-03-30"

"""
WBM exploratory data analysis (EDA).
Start with comparing MP and WBM elemental prevalence.
"""

module_dir = os.path.dirname(__file__)
about_data_page = f"{ROOT}/site/src/routes/about-the-data"


# %% load MP training set
df_mp = pd.read_csv(DATA_FILES.mp_energies, na_filter=False)


# %%
for count_mode in ["occurrence", "composition"]:
    wbm_elem_counts = count_elements(df_wbm.formula, count_mode=count_mode).astype(int)

    wbm_elem_counts.to_json(f"{about_data_page}/wbm-element-counts-{count_mode}.json")
    mp_elem_counts = count_elements(df_mp.formula_pretty, count_mode=count_mode).astype(
        int
    )
    mp_elem_counts.to_json(f"{about_data_page}/mp-element-counts-{count_mode}.json")


# %% export element counts by WBM step to JSON
df_wbm["step"] = df_wbm.index.str.split("-").str[1].astype(int)
assert df_wbm.step.between(1, 5).all()
for batch in range(1, 6):
    count_elements(df_wbm[df_wbm.step == batch].formula).to_json(
        f"{about_data_page}/wbm-element-counts-{batch=}.json"
    )

# export element counts by arity (how many elements in the formula)
comp_col = "composition"
df_wbm[comp_col] = df_wbm.formula.map(Composition)

for arity, df_mp in df_wbm.groupby(df_wbm[comp_col].map(len)):
    count_elements(df_mp.formula).to_json(
        f"{about_data_page}/wbm-element-counts-{arity=}.json"
    )


# %%
wbm_fig = ptable_heatmap_plotly(
    wbm_elem_counts.drop("Xe"),
    log=True,
    colorscale="RdBu",
    hover_props=dict(atomic_number="atomic number"),
    hover_data=wbm_elem_counts,
)

title = "WBM Elements"
wbm_fig.update_layout(
    title=dict(text=title, x=0.35, y=0.9, font_size=20),
    xaxis=dict(fixedrange=True),
    yaxis=dict(fixedrange=True),
    paper_bgcolor="rgba(0,0,0,0)",
)
wbm_fig.show()


# %%
wbm_fig.write_image(f"{module_dir}/figs/wbm-elements.svg", width=1000, height=500)
# save_fig(wbm_fig, f"{FIGS}/wbm-elements.svelte")


# %%
mp_fig = ptable_heatmap_plotly(
    mp_elem_counts[mp_elem_counts > 1],
    log=True,
    colorscale="RdBu",
    hover_props=dict(atomic_number="atomic number"),
    hover_data=mp_elem_counts,
)

title = "MP Elements"
mp_fig.update_layout(
    title=dict(text=title, x=0.35, y=0.9, font_size=20),
    xaxis=dict(fixedrange=True),
    yaxis=dict(fixedrange=True),
    paper_bgcolor="rgba(0,0,0,0)",
)
mp_fig.show()


# %%
mp_fig.write_image(f"{module_dir}/figs/{today}-mp-elements.svg", width=1000, height=500)
# save_fig(mp_fig, f"{FIGS}/mp-elements.svelte")


# %% histogram of energy above MP convex hull for WBM
col = "e_above_hull_mp2020_corrected_ppd_mp"
# col = "e_form_per_atom_mp2020_corrected"
mean, std = df_wbm[col].mean(), df_wbm[col].std()

range_x = (mean - 2 * std, mean + 2 * std)
counts, bins = np.histogram(df_wbm[col], bins=150, range=range_x)
x_label = "WBM energy above MP convex hull (eV/atom)"
df_hist = pd.DataFrame([counts, bins], index=["count", x_label]).T

fig = df_hist.plot.area(x=x_label, y="count", backend="plotly", range_x=range_x)

if col.startswith("e_above_hull"):
    n_stable = sum(df_wbm[col] <= 0)
    n_unstable = sum(df_wbm[col] > 0)
    assert n_stable + n_unstable == len(df_wbm.dropna())

    dummy_mae = (df_wbm[col] - df_wbm[col].mean()).abs().mean()

    title = (
        f"n={len(df_wbm.dropna()):,} with {n_stable:,} stable + {n_unstable:,} "
        f"unstable, dummy MAE={dummy_mae:.2f}"
    )
    fig.update_layout(title=dict(text=title, x=0.5, y=0.95))

fig.update_layout(showlegend=False)

for x_pos, label in (
    (mean, f"{mean = :.2f}"),
    (mean - std, f"{mean - std = :.2f}"),
    (mean + std, f"{mean + std = :.2f}"),
):
    anno = dict(text=label, yshift=-10, xshift=-5, xanchor="right")
    fig.add_vline(x=x_pos, line=dict(width=1, dash="dash"), annotation=anno)

fig.show()

save_fig(fig, f"{FIGS}/wbm-each-hist.svelte")
save_fig(fig, "./figs/wbm-each-hist.svg", width=1000, height=500)


# %%
e_col, n_atoms_col = "Energy (eV/atom)", "Number of Atoms"
elem_num_col = "Atomic number"
mp_ref_data = [
    {
        "Element": key,
        e_col: entry.energy_per_atom,
        elem_num_col: entry.composition.elements[0].number,
        n_atoms_col: entry.composition.num_atoms,
        "Name": entry.composition.elements[0].long_name,
        "Material ID": entry.entry_id.replace("-GGA", ""),
    }
    for key, entry in mp_elem_reference_entries.items()
]
df_ref = pd.DataFrame(mp_ref_data).sort_values(elem_num_col)


# %% plot MP elemental reference energies vs atomic number
# marker size = number of atoms in reference structure
fig = df_ref.round(2).plot.scatter(
    x=elem_num_col, y=e_col, backend="plotly", hover_data=list(df_ref), size=n_atoms_col
)
fig.update_traces(mode="markers+lines")
fig.layout.margin = dict(l=0, r=0, t=0, b=0)

# add text annotations showing element symbols
for symbol, e_per_atom, num, *_ in df_ref.itertuples(index=False):
    fig.add_annotation(x=num, y=e_per_atom, text=symbol, showarrow=False, font_size=10)

fig.show()

save_fig(fig, f"{FIGS}/mp-elemental-ref-energies.svelte")


# %% plot 2d and 3d t-SNE projections of one-hot encoded element vectors summed by
# weight in each WBM composition. TLDR: no obvious structure in the data
# was hoping to find certain clusters to have higher or lower errors after seeing
# many models struggle on the halogens in per-element error periodic table heatmaps
# https://janosh.github.io/matbench-discovery/models
df_2d_tsne = pd.read_csv(f"{module_dir}/tsne/one-hot-112-composition-2d.csv.gz")
df_2d_tsne = df_2d_tsne.set_index("material_id")

df_3d_tsne = pd.read_csv(f"{module_dir}/tsne/one-hot-112-composition-3d.csv.gz")
model = "Wrenformer"
df_3d_tsne = pd.read_csv(
    f"{module_dir}/tsne/one-hot-112-composition+{model}-each-err-3d-metric=eucl.csv.gz"
)
df_3d_tsne = df_3d_tsne.set_index("material_id")

df_wbm[list(df_2d_tsne)] = df_2d_tsne
df_wbm[list(df_3d_tsne)] = df_3d_tsne
df_wbm[list(df_each_err.add_suffix(" abs EACH error"))] = df_each_err.abs()


# %%
color_col = f"{model} abs EACH error"
clr_range_max = df_wbm[color_col].mean() + df_wbm[color_col].std()


# %%
fig = px.scatter(
    df_wbm,
    x="2d t-SNE 1",
    y="2d t-SNE 2",
    color=color_col,
    hover_name="material_id",
    hover_data=("formula", each_true_col),
    range_color=(0, clr_range_max),
)
fig.show()


# %%
fig = px.scatter_3d(
    df_wbm,
    x="3d t-SNE 1",
    y="3d t-SNE 2",
    z="3d t-SNE 3",
    color=color_col,
    custom_data=["material_id", "formula", each_true_col, color_col],
    range_color=(0, clr_range_max),
)
fig.data[0].hovertemplate = (
    "<b>material_id: %{customdata[0]}</b><br><br>"
    "t-SNE: (%{x:.2f}, %{y:.2f}, %{z:.2f})<br>"
    "Formula: %{customdata[1]}<br>"
    "E<sub>above hull</sub>: %{customdata[2]:.2f}<br>"
    f"{color_col}: %{{customdata[3]:.2f}}<br>"
)
fig.show()


# %%
wyk_col, spg_col = "wyckoff_spglib", "spacegroup"
df_wbm[spg_col] = df_wbm[wyk_col].str.split("_").str[2].astype(int)
df_mp[spg_col] = df_mp[wyk_col].str.split("_").str[2].astype(int)


# %%
fig = spacegroup_sunburst(df_wbm[spg_col], width=350, height=350, show_counts="percent")
fig.layout.title.update(text="WBM Spacegroup Sunburst", x=0.5, font_size=14)
fig.show()
save_fig(fig, f"{FIGS}/spacegroup-sunburst-wbm.svelte")


# %%
fig = spacegroup_sunburst(df_mp[spg_col], width=350, height=350, show_counts="percent")
fig.layout.title.update(text="MP Spacegroup Sunburst", x=0.5, font_size=14)
fig.show()
save_fig(fig, f"{FIGS}/spacegroup-sunburst-mp.svelte")
# would be good to have consistent order of crystal systems between sunbursts but not
# controllable yet
# https://github.com/plotly/plotly.py/issues/4115
# https://github.com/plotly/plotly.js/issues/5341
# https://github.com/plotly/plotly.js/issues/4728
