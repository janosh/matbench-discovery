# %%
import os

import pandas as pd
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


# %%
wbm_elem_counts = count_elements(df_wbm.formula).astype(int)

out_elem_counts = f"{ROOT}/site/src/routes/about-the-test-set/wbm-element-counts.json"
# wbm_elem_counts.to_json(out_elem_counts)


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
save_fig(wbm_fig, f"{FIGS}/{today}-wbm-elements.svelte")


# %% load MP training set
df = pd.read_json(f"{module_dir}/../mp/2022-08-13-mp-energies.json.gz")
mp_elem_counts = count_elements(df.formula_pretty).astype(int)

# mp_elem_counts.to_json(
#     f"{ROOT}/site/src/routes/about-the-test-set/{today}-mp-element-counts.json"
# )
mp_elem_counts.describe()


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
# save_fig(mp_fig, f"{FIGS}/{today}-mp-elements.svelte")


# %% histogram of energy above MP convex hull for WBM
col = "e_above_hull_mp2020_corrected_ppd_mp"
# col = "e_form_per_atom_mp2020_corrected"
mean, std = df_wbm[col].mean(), df_wbm[col].std()

fig = df_wbm[col].hist(
    bins=100,
    backend="plotly",
    range_x=[mean - 2 * std, mean + 2 * std],
    template="plotly_dark",
)

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

fig.update_layout(showlegend=False, paper_bgcolor="rgba(0,0,0,0)")
fig.update_xaxes(title="WBM energy above MP convex hull (eV/atom)")

for x_pos, label in zip(
    [mean, mean + std, mean - std],
    [f"{mean = :.2f}", f"{mean + std = :.2f}", f"{mean - std = :.2f}"],
):
    anno = dict(text=label, yshift=-10, xshift=5)
    fig.add_vline(x=x_pos, line=dict(width=1, dash="dash"), annotation=anno)

fig.show()


# subsample x
for trace in fig.data:
    trace.x = trace.x[::8]

save_fig(fig, f"{FIGS}/{today}-wbm-each-hist.svelte")
save_fig(fig, f"./figs/{today}-wbm-each-hist.svg", width=1000, height=500)


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

save_fig(fig, f"{FIGS}/{today}-mp-elemental-ref-energies.svelte")
