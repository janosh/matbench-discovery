# %%
import os

import numpy as np
import pandas as pd
from pymatgen.core import Composition
from pymatviz import count_elements, ptable_heatmap_plotly
from pymatviz.utils import save_fig

from matbench_discovery import FIGS, ROOT, today
from matbench_discovery.data import df_wbm
from matbench_discovery.energy import mp_elem_reference_entries
from matbench_discovery.plots import pio

"""
Compare MP and WBM elemental prevalence. Starting with WBM, MP below.
"""

module_dir = os.path.dirname(__file__)
print(f"{pio.templates.default=}")
about_data_page = f"{ROOT}/site/src/routes/about-the-test-set"


# %%
wbm_elem_counts = count_elements(df_wbm.formula).astype(int)

# wbm_elem_counts.to_json(f"{about_data_page}/wbm-element-counts.json")

# export element counts by WBM step to JSON
df_wbm["step"] = df_wbm.index.str.split("-").str[1].astype(int)
assert df_wbm.step.between(1, 5).all()
for batch in range(1, 6):
    count_elements(df_wbm[df_wbm.step == batch].formula).to_json(
        f"{about_data_page}/wbm-element-counts-{batch=}.json"
    )

# export element counts by arity (how many elements in the formula)
comp_col = "composition"
df_wbm[comp_col] = df_wbm.formula.map(Composition)

for arity, df in df_wbm.groupby(df_wbm[comp_col].map(len)):
    count_elements(df.formula).to_json(
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


# %% load MP training set
df = pd.read_json(f"{module_dir}/../mp/2022-08-13-mp-energies.json.gz")
mp_elem_counts = count_elements(df.formula_pretty).astype(int)

# mp_elem_counts.to_json(f"{about_data_page}/mp-element-counts.json")


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
mp_ref_data = [
    {
        "Element": key,
        e_col: ref.energy_per_atom,
        n_atoms_col: ref.composition.num_atoms,
        "Name": ref.composition.elements[0].long_name,
        "Number": ref.composition.elements[0].number,
    }
    for key, ref in mp_elem_reference_entries.items()
]
df_ref = pd.DataFrame(mp_ref_data).sort_values("Number")


# %% plot MP elemental reference energies vs atomic number
# marker size = number of atoms in reference structure
fig = df_ref.round(2).plot.scatter(
    x="Number", y=e_col, backend="plotly", hover_data=list(df_ref), size=n_atoms_col
)
fig.update_traces(mode="markers+lines")
fig.layout.margin = dict(l=0, r=0, t=0, b=0)

# add text annotations showing element symbols
for symbol, e_per_atom, *_, num in df_ref.itertuples(index=False):
    fig.add_annotation(x=num, y=e_per_atom, text=symbol, showarrow=False, font_size=10)

fig.show()

save_fig(fig, f"{FIGS}/mp-elemental-ref-energies.svelte")
