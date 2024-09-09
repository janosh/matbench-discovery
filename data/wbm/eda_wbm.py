"""WBM exploratory data analysis (EDA).
Start with comparing MP and WBM elemental prevalence.
"""

# %%
import os

import numpy as np
import pandas as pd
import plotly.express as px
import pymatviz as pmv
from matplotlib.colors import SymLogNorm
from pymatgen.core import Composition, Structure
from pymatviz.enums import Key
from pymatviz.utils import PLOTLY, si_fmt_int

from matbench_discovery import PDF_FIGS, ROOT, SITE_FIGS, STABILITY_THRESHOLD
from matbench_discovery import plots as plots
from matbench_discovery.data import DataFiles, df_wbm
from matbench_discovery.energy import mp_elem_ref_entries
from matbench_discovery.enums import MbdKey
from matbench_discovery.preds import df_each_err

__author__ = "Janosh Riebesell"
__date__ = "2023-03-30"

module_dir = os.path.dirname(__file__)
data_page = f"{ROOT}/site/src/routes/data"


# %% load MP training set
df_mp = pd.read_csv(DataFiles.mp_energies.path, na_filter=False)
df_mp = df_mp.set_index(Key.mat_id)


# %%
wbm_occu_counts = pmv.count_elements(df_wbm[Key.formula], count_mode="occurrence")
wbm_occu_counts = wbm_occu_counts.dropna().astype(int)
wbm_comp_counts = pmv.count_elements(df_wbm[Key.formula], count_mode="composition")
wbm_comp_counts = wbm_comp_counts.dropna().astype(int)

mp_occu_counts = pmv.count_elements(df_mp[Key.formula], count_mode="occurrence")
mp_occu_counts = mp_occu_counts.dropna().astype(int)

mp_comp_counts = pmv.count_elements(df_mp[Key.formula], count_mode="composition")
mp_comp_counts = mp_comp_counts.dropna().astype(int)

all_counts = (
    ("wbm", "occurrence", wbm_occu_counts),
    ("wbm", "composition", wbm_comp_counts),
    ("mp", "occurrence", mp_occu_counts),
    ("mp", "composition", mp_comp_counts),
)


# %% print prevalence of stable structures in full WBM and uniq-prototypes only
print(f"{STABILITY_THRESHOLD=}")
for df, label in (
    (df_wbm, "full WBM"),
    (df_wbm.query(Key.uniq_proto), "WBM unique prototypes"),
):
    n_stable = sum(df[MbdKey.each_true] <= STABILITY_THRESHOLD)
    stable_rate = n_stable / len(df)
    print(f"{label}: {stable_rate=:.1%} ({n_stable:,} out of {len(df):,})")

# on 2024-04-15: STABILITY_THRESHOLD=0
# full WBM: stable_rate=16.7% (42,825 out of 256,963)
# WBM unique prototypes: stable_rate=15.3% (32,942 out of 215,488)


# %%
for dataset, count_mode, elem_counts in all_counts:
    filename = f"{dataset}-element-counts-by-{count_mode}"
    elem_counts.to_json(f"{data_page}/{filename}.json")

    title = f"Number of {dataset.upper()} structures containing each element"
    fig = pmv.ptable_heatmap_plotly(elem_counts, font_size=10, fmt=",.0f")
    fig.layout.title.update(text=title, x=0.4, y=0.9)
    fig.show()

    # saving matplotlib heatmap to PDF mostly for historical reasons, could also use
    # pmv.ptable_heatmap_plotly
    ax_elem_counts = pmv.ptable_heatmap(
        elem_counts,
        cbar_title=f"{dataset.upper()} Element Count",
        log=(log := SymLogNorm(linthresh=100)),
    )
    if log:
        filename += "-symlog" if isinstance(log, SymLogNorm) else "-log"
    pmv.save_fig(ax_elem_counts, f"{PDF_FIGS}/{filename}.pdf")


# %% ratio of WBM to MP counts
normalized = True
ax_ptable = pmv.ptable_heatmap_ratio(
    wbm_occu_counts / (len(df_wbm) if normalized else 1),
    mp_occu_counts / (len(df_mp) if normalized else 1),
    zero_color="#efefef",
    exclude_elements="Xe Th Pa U Np Pu".split(),
)
img_name = "wbm-mp-ratio-element-counts-by-occurrence"
if normalized:
    img_name += "-normalized"
pmv.save_fig(ax_ptable, f"{PDF_FIGS}/{img_name}.pdf")


# %% export element counts by WBM step to JSON
df_wbm["step"] = df_wbm.index.str.split("-").str[1].astype(int)
assert df_wbm.step.between(1, 5).all()
for batch in range(1, 6):
    pmv.count_elements(df_wbm[df_wbm.step == batch][Key.formula]).to_json(
        f"{data_page}/wbm-element-counts-{batch=}.json"
    )

# export element counts by arity (how many elements in the formula)
df_wbm[Key.composition] = df_wbm[Key.formula].map(Composition)

for arity, df_mp in df_wbm.groupby(df_wbm[Key.composition].map(len)):
    pmv.count_elements(df_mp[Key.formula]).to_json(
        f"{data_page}/wbm-element-counts-{arity=}.json"
    )


# %%
for dataset, count_mode, elem_counts in all_counts:
    fig = pmv.ptable_heatmap_plotly(
        elem_counts.drop("Xe")[elem_counts > 1],
        font_size=11,
        color_bar=dict(title=dict(text=f"WBM {count_mode} counts", font_size=24)),
        log=True,
        hover_props=dict(atomic_number="atomic number"),
        hover_data=wbm_occu_counts,
    )

    fig.layout.margin = dict(l=0, r=0, b=0, t=0)
    fig.show()
    svg_path = f"{module_dir}/figs/wbm-elements.svg"
    # pmv.save_fig(fig, svg_path, width=1000, height=500)
    pmv.save_fig(fig, f"{PDF_FIGS}/{dataset}-element-{count_mode}-counts.pdf")


# %% histogram of energy distance to MP convex hull for WBM
e_col = MbdKey.each_true
# e_col = MbdKey.e_form_raw
# e_col = MbdKey.e_form
mean, std = df_wbm[e_col].mean(), df_wbm[e_col].std()

range_x = (mean - 2 * std, mean + 2 * std)
counts, bins = np.histogram(df_wbm[e_col], bins=150, range=range_x)
bins = bins[1:]  # remove left-most bin edge
left_counts = counts[bins < 0]
right_counts = counts[bins >= 0]

assert e_col.startswith(("e_form_per_atom", "e_above_hull"))
x_label = "energy above MP convex hull" if "above" in e_col else "formation energy"
y_label = "Number of Structures"
fig = px.bar(
    x=bins[bins < 0],
    y=left_counts,
    labels={"x": f"WBM {x_label} (eV/atom)", "y": y_label},
)
fig.add_bar(x=bins[bins >= 0], y=right_counts)
fig.update_traces(width=(bins[1] - bins[0]))  # make bars touch

if e_col.startswith("e_above_hull"):
    n_stable = sum(df_wbm[e_col] <= STABILITY_THRESHOLD)
    n_unstable = sum(df_wbm[e_col] > STABILITY_THRESHOLD)
    assert n_stable + n_unstable == len(df_wbm.dropna())

    dummy_mae = (df_wbm[e_col] - df_wbm[e_col].mean()).abs().mean()

    title = (
        f"{si_fmt_int(len(df_wbm.dropna()))} structures with {si_fmt_int(n_stable)} "
        f"stable + {si_fmt_int(n_unstable)} unstable (stable rate="
        f"{n_stable / len(df_wbm):.1%})"
    )
    fig.layout.title = dict(text=title, x=0.5, font_size=16, y=0.95)

    # add red/blue annotations to left and right of mean saying stable/unstable
    for idx, (label, x_pos) in enumerate(
        (("stable", mean - std), ("unstable", mean + std))
    ):
        fig.add_annotation(
            x=x_pos,
            y=0.5,
            text=label,
            showarrow=False,
            font_size=18,
            font_color=px.colors.qualitative.Plotly[idx],
            yref="paper",
            xanchor="right",
            xshift=-40,
        )


fig.layout.margin = dict(l=0, r=0, b=0, t=40)
fig.update_layout(showlegend=False)

for x_pos, label in (
    (mean, f"{mean = :.2f}"),
    (mean - std, f"{mean - std = :.2f}"),
    (mean + std, f"{mean + std = :.2f}"),
):
    anno = dict(text=label, yshift=-10, xshift=-5, xanchor="right")
    line_width = 3 if x_pos == mean else 2
    fig.add_vline(x=x_pos, line=dict(width=line_width, dash="dash"), annotation=anno)

fig.show()
suffix = {
    MbdKey.each_true: "hull-dist",
    MbdKey.e_form_dft: "e-form",
    MbdKey.e_form_raw: "e-form-uncorrected",
}[e_col]
img_name = f"hist-wbm-{suffix}"
pmv.save_fig(fig, f"{SITE_FIGS}/{img_name}.svelte")
# pmv.save_fig(fig, f"./figs/{img_name}.svg", width=800, height=500)
pmv.save_fig(fig, f"{PDF_FIGS}/{img_name}.pdf", width=600, height=300)


# %%
e_col, n_atoms_col = "Energy (eV/atom)", "Number of Atoms"
atom_num_col = "Atomic number"
mp_ref_data = [
    {
        "Element": key,
        e_col: entry.energy_per_atom,
        atom_num_col: entry.composition.elements[0].number,
        n_atoms_col: entry.composition.num_atoms,
        "Name": entry.composition.elements[0].long_name,
        "Material ID": entry.entry_id.replace("-GGA", ""),
    }
    for key, entry in mp_elem_ref_entries.items()
]
df_ref = pd.DataFrame(mp_ref_data).sort_values(atom_num_col)


# %% plot MP elemental reference energies vs atomic number
# marker size = number of atoms in reference structure
fig = df_ref.round(2).plot.scatter(
    x=atom_num_col, y=e_col, backend=PLOTLY, hover_data=list(df_ref), size=n_atoms_col
)
fig.update_traces(mode="markers+lines")
fig.layout.margin = dict(l=0, r=0, t=0, b=0)

# add text annotations showing element symbols
for symbol, e_per_atom, num, *_ in df_ref.itertuples(index=False):
    fig.add_annotation(x=num, y=e_per_atom, text=symbol, showarrow=False, font_size=10)

fig.show()

pmv.save_fig(fig, f"{SITE_FIGS}/mp-elemental-ref-energies.svelte")
pmv.save_fig(fig, f"{PDF_FIGS}/mp-elemental-ref-energies.pdf")


# %% plot 2d and 3d t-SNE projections of one-hot encoded element vectors summed by
# weight in each WBM composition. TLDR: no obvious structure in the data
# was hoping to find certain clusters to have higher or lower errors after seeing
# many models struggle on the halogens in per-element error periodic table heatmaps
# https://janosh.github.io/matbench-discovery/models
df_2d_tsne = pd.read_csv(f"{module_dir}/tsne/one-hot-112-composition-2d.csv.gz")
df_2d_tsne = df_2d_tsne.set_index(Key.mat_id)

df_3d_tsne = pd.read_csv(f"{module_dir}/tsne/one-hot-112-composition-3d.csv.gz")

df_wbm[list(df_2d_tsne)] = df_2d_tsne
df_wbm[list(df_3d_tsne)] = df_3d_tsne
df_wbm[list(df_each_err.add_suffix(" abs EACH error"))] = df_each_err.abs()


# %%
fig = px.scatter(
    df_wbm,
    x="2d t-SNE 1",
    y="2d t-SNE 2",
    color="step",
    hover_name=Key.mat_id,
    hover_data=(Key.formula, MbdKey.each_true),
)
fig.show()


# %%
fig = px.scatter_3d(
    df_wbm,
    x="3d t-SNE 1",
    y="3d t-SNE 2",
    z="3d t-SNE 3",
    color="step",
    custom_data=[Key.mat_id, Key.formula, MbdKey.each_true],
)
fig.data[0].hovertemplate = (
    "<b>material_id: %{customdata[0]}</b><br><br>"
    "t-SNE: (%{x:.2f}, %{y:.2f}, %{z:.2f})<br>"
    "Formula: %{customdata[1]}<br>"
    "E<sub>above hull</sub>: %{customdata[2]:.2f}<br>"
    "WBM step: %{customdata[3]:.2f}<br>"
)
fig.show()


# %%
df_wbm[Key.spg_num] = df_wbm[MbdKey.init_wyckoff].str.split("_").str[2].astype(int)
df_mp[Key.spg_num] = df_mp[f"{Key.wyckoff}_spglib"].str.split("_").str[2].astype(int)


# %%
fig = pmv.spacegroup_sunburst(
    df_wbm[Key.spg_num], width=350, height=350, show_counts="percent"
)
fig.layout.title.update(text="WBM Spacegroup Sunburst", x=0.5, font_size=14)
fig.layout.margin = dict(l=0, r=0, t=30, b=0)
fig.show()
pmv.save_fig(fig, f"{SITE_FIGS}/spacegroup-sunburst-wbm.svelte")
pmv.save_fig(fig, f"{PDF_FIGS}/spacegroup-sunburst-wbm.pdf")


# %%
fig = pmv.spacegroup_sunburst(
    df_mp[Key.spg_num], width=350, height=350, show_counts="percent"
)
fig.layout.title.update(text="MP Spacegroup Sunburst", x=0.5, font_size=14)
fig.layout.margin = dict(l=0, r=0, t=30, b=0)
fig.show()
pmv.save_fig(fig, f"{SITE_FIGS}/spacegroup-sunburst-mp.svelte")
pmv.save_fig(fig, f"{PDF_FIGS}/spacegroup-sunburst-mp.pdf")
# would be good to have consistent order of crystal systems between sunbursts but not
# controllable yet
# https://github.com/plotly/plotly.py/issues/4115
# https://github.com/plotly/plotly.js/issues/5341
# https://github.com/plotly/plotly.js/issues/4728


# %% compute compositional arity histograms
arity_col = "arity"
df_wbm[arity_col] = df_wbm[Key.formula].map(Composition).map(len)
df_mp[arity_col] = df_mp[Key.formula].map(Composition).map(len)

mp_arity_counts = df_mp[arity_col].value_counts().sort_index() / len(df_mp)
wbm_arity_counts = df_wbm[arity_col].value_counts().sort_index() / len(df_wbm)

df_arity = pd.DataFrame({"MP": mp_arity_counts, "WBM": wbm_arity_counts}).query(
    "0 < index < 7"
)

fig = px.bar(df_arity, barmode="group")
fig.layout.legend.update(x=1, y=1, xanchor="right", yanchor="top", title=None)
fig.layout.margin = dict(l=0, r=0, b=0, t=0)
fig.layout.yaxis.title = "Fraction of Structures in Dataset"
fig.layout.xaxis.title = "Number of Elements in Formula"

fig.show()
img_name = "mp-vs-wbm-arity-hist"
pmv.save_fig(fig, f"{SITE_FIGS}/{img_name}.svelte")
pmv.save_fig(fig, f"{PDF_FIGS}/{img_name}.pdf", width=450, height=280)


# %% find large structures that changed symmetry during relaxation
df_sym_change = (
    df_wbm.query(f"{MbdKey.init_wyckoff} != {Key.wyckoff}_spglib")
    .filter(regex="wyckoff|sites")
    .nlargest(10, Key.n_sites)
)


# %%
df_wbm_structs = pd.read_json(DataFiles.wbm_cses_plus_init_structs.path)
df_wbm_structs = df_wbm_structs.set_index(Key.mat_id)


# %%
for wbm_id in df_sym_change.index:
    init_struct = Structure.from_dict(df_wbm_structs.loc[wbm_id][Key.init_struct])
    final_struct = Structure.from_dict(df_wbm_structs.loc[wbm_id][Key.cse]["structure"])
    init_struct.properties[Key.mat_id] = f"{wbm_id}-init"
    final_struct.properties[Key.mat_id] = f"{wbm_id}-final"

    pmv.structure_2d([init_struct, final_struct])


# %% export initial and final structures with symmetry change to CIF
wbm_id = df_sym_change.index[0]

struct = Structure.from_dict(df_wbm_structs.loc[wbm_id][Key.cse]["structure"])
struct.to(f"{module_dir}/{wbm_id}.cif")
struct.to(f"{module_dir}/{wbm_id}.json")

struct = Structure.from_dict(df_wbm_structs.loc[wbm_id][Key.init_struct])
struct.to(f"{module_dir}/{wbm_id}-init.cif")
struct.to(f"{module_dir}/{wbm_id}-init.json")
