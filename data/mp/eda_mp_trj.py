"""MPtrj exploratory data analysis (EDA)."""

# %%
import os
from collections import defaultdict
from typing import Any

import ase
import ase.io.extxyz
import matplotlib.colorbar
import matplotlib.colors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.express as px
import pymatviz as pmv
from matplotlib.colors import SymLogNorm
from pymatgen.core import Composition
from pymatgen.core.tensors import Tensor
from pymatviz.enums import Key
from tqdm import tqdm

from matbench_discovery import MP_DIR, PDF_FIGS, ROOT, SITE_FIGS
from matbench_discovery.data import DataFiles, ase_atoms_from_zip, df_wbm
from matbench_discovery.energy import get_e_form_per_atom
from matbench_discovery.enums import MbdKey

__author__ = "Janosh Riebesell"
__date__ = "2023-11-22"

data_page = f"{ROOT}/site/src/routes/data"


# %% load MP element counts by occurrence to compute ratio with MPtrj
mp_occu_counts = pd.read_json(
    f"{data_page}/mp-element-counts-by-occurrence.json", typ="series"
)
df_mp = pd.read_csv(DataFiles.mp_energies.path, na_filter=False)
df_mp = df_mp.set_index(Key.mat_id)
assert sum(df_mp[Key.formula].isna() | (df_mp[Key.formula] == "")) == 0


# %% --- load preprocessed MPtrj summary data if available ---
mp_trj_summary_path = f"{MP_DIR}/2022-09-16-mp-trj-summary.json.bz2"
if os.path.isfile(mp_trj_summary_path):
    df_mp_trj = pd.read_json(mp_trj_summary_path)
    df_mp_trj.index.name = "frame_id"
else:
    print("MPtrj summary data not found, run cell below to generate")


# %% extract extXYZ files from zipped directory without unpacking the whole archive
# takes ~8 mins on M2 Max
# takes ~5 mins on M3 Max
atoms_list = ase_atoms_from_zip(DataFiles.mp_trj_extxyz.path)

mp_trj_atoms: dict[str, list[ase.Atoms]] = defaultdict(list)
for atoms in atoms_list:
    mp_id = atoms.info.get(Key.mat_id, "no-id")
    assert mp_id.startswith(("mp-", "mvc-"))
    mp_trj_atoms[mp_id].append(atoms)

del atoms_list  # free up memory

assert len(mp_trj_atoms) == 145_923  # number of unique MP IDs


# %%
def info_dict_to_id(info: dict[str, int | str]) -> str:
    """Construct a unique frame ID from the atoms info dict."""
    return f"{info[Key.task_id]}-{info['calc_id']}-{info['ionic_step']}"


df_mp_trj = pd.DataFrame(
    {
        info_dict_to_id(atoms.info): atoms.info
        | {
            Key.forces: atoms.get_forces(),
            Key.stress: atoms.get_stress(),
            Key.magmoms: atoms.get_magnetic_moments()
            if "magmoms" in atoms.calc.results
            else None,
            Key.formula: str(atoms.symbols),
            Key.atom_nums: atoms.symbols,
        }
        for atoms_list in tqdm(mp_trj_atoms.values(), total=len(mp_trj_atoms))
        for atoms in atoms_list
    }
).T.convert_dtypes()  # convert object columns to float/int where possible
df_mp_trj.index.name = "frame_id"
assert len(df_mp_trj) == 1_580_395  # number of total frames
if Key.formula not in df_mp_trj:
    raise KeyError(f"{Key.formula!s} not in {df_mp_trj.columns=}")

# this is the unrelaxed (but MP2020 corrected) formation energy per atom of the actual
# relaxation step
df_mp_trj[MbdKey.e_form_dft] = [
    get_e_form_per_atom(
        {"composition": row[Key.formula], "energy": row["mp2020_corrected_energy"]}
    )
    for _idx, row in tqdm(
        df_mp_trj.iterrows(), total=len(df_mp_trj), desc="Compute formation energies"
    )
]
df_mp_trj[Key.stress_trace] = [
    np.trace(Tensor.from_voigt(stress)) / 3 for stress in tqdm(df_mp_trj[Key.stress])
]


# %%
df_mp_trj.to_json(mp_trj_summary_path)


# %%
def tile_count_anno(hist_vals: list[Any]) -> dict[str, Any]:
    """Annotate each periodic table tile with the number of values in its histogram."""
    face_color = cmap(norm(np.sum(len(hist_vals)))) if hist_vals else "none"
    bbox = dict(facecolor=face_color, alpha=0.4, pad=2, edgecolor="none")
    return dict(text=pmv.si_fmt(len(hist_vals), ".0f"), bbox=bbox)


# %% plot per-element magmom histograms
ptable_magmom_hist_path = f"{MP_DIR}/2022-09-16-mp-trj-elem-magmoms.json.bz2"
srs_mp_trj_elem_magmoms = locals().get("srs_mp_trj_elem_magmoms")

if os.path.isfile(ptable_magmom_hist_path):
    srs_mp_trj_elem_magmoms = pd.read_json(ptable_magmom_hist_path, typ="series")
if srs_mp_trj_elem_magmoms is None:
    # project magmoms onto symbols in dict
    df_mp_trj_elem_magmom = pd.DataFrame(
        [
            dict(zip(elems, magmoms, strict=False))
            for elems, magmoms in df_mp_trj[[Key.atom_nums, Key.magmoms]]
            .dropna()
            .itertuples(index=False)
        ]
    )

    srs_mp_trj_elem_magmoms = {
        col: list(df_mp_trj_elem_magmom[col].dropna()) for col in df_mp_trj_elem_magmom
    }
    pd.Series(srs_mp_trj_elem_magmoms).to_json(ptable_magmom_hist_path)

cmap = plt.get_cmap(color_map := "viridis")
norm = matplotlib.colors.LogNorm(vmin=1, vmax=150_000)

fig_ptable_magmoms = pmv.ptable_hists(
    srs_mp_trj_elem_magmoms,
    symbol_pos=(0.2, 0.8),
    log=True,
    cbar_title="Magmoms ($μ_B$)",
    cbar_title_kwargs=dict(fontsize=16),
    cbar_coords=(0.18, 0.85, 0.42, 0.02),
    # annotate each element with its number of magmoms in MPtrj
    # anno_kwds=tile_count_anno,
    colormap=color_map,
)

cbar_ax = fig_ptable_magmoms.figure.add_axes([0.26, 0.78, 0.25, 0.015])
cbar = matplotlib.colorbar.ColorbarBase(
    cbar_ax, cmap=cmap, norm=norm, orientation="horizontal"
)
pmv.save_fig(fig_ptable_magmoms, f"{PDF_FIGS}/mp-trj-magmoms-ptable-hists.pdf")


# %% plot per-element force histograms
ptable_force_hist_path = f"{MP_DIR}/2022-09-16-mp-trj-elem-forces.json.bz2"
srs_mp_trj_elem_forces = locals().get("srs_mp_trj_elem_forces")

if os.path.isfile(ptable_force_hist_path):
    srs_mp_trj_elem_forces = pd.read_json(ptable_force_hist_path, typ="series")
if srs_mp_trj_elem_forces is None:
    df_mp_trj_elem_forces = pd.DataFrame(
        [
            dict(zip(elems, np.abs(forces).mean(axis=1), strict=False))
            for elems, forces in df_mp_trj[[Key.atom_nums, Key.forces]].itertuples(
                index=False
            )
        ]
    )
    mp_trj_elem_forces = {
        col: list(df_mp_trj_elem_forces[col].dropna()) for col in df_mp_trj_elem_forces
    }
    srs_mp_trj_elem_forces = pd.Series(mp_trj_elem_forces)
    srs_mp_trj_elem_forces.to_json(ptable_force_hist_path)

cmap = plt.get_cmap(color_map := "viridis")
norm = matplotlib.colors.LogNorm(vmin=1, vmax=1_000_000)

max_force = 10  # eV/Å
fig_ptable_forces = pmv.ptable_hists(
    srs_mp_trj_elem_forces.copy().map(lambda x: [val for val in x if val < max_force]),
    symbol_pos=(0.3, 0.8),
    log=True,
    cbar_title="1/3 Σ|Forces| (eV/Å)",
    cbar_title_kwargs=dict(fontsize=16),
    cbar_coords=(0.18, 0.85, 0.42, 0.02),
    x_range=(0, max_force),
    # anno_kwds=tile_count_anno,
    colormap=color_map,
)

cbar_ax = fig_ptable_forces.figure.add_axes([0.26, 0.78, 0.25, 0.015])
cbar = matplotlib.colorbar.ColorbarBase(
    cbar_ax, cmap=cmap, norm=norm, orientation="horizontal"
)

pmv.save_fig(fig_ptable_forces, f"{PDF_FIGS}/mp-trj-forces-ptable-hists.pdf")


# %% plot histogram of number of sites per element
# TODO fix weirdness with 6x10^e0 y axis label on Cl tile
ptable_n_sites_hist_path = f"{MP_DIR}/2022-09-16-mp-trj-elem-n-sites.json.bz2"
srs_mp_trj_elem_n_sites = locals().get("srs_mp_trj_elem_n_sites")

if os.path.isfile(ptable_n_sites_hist_path):
    srs_mp_trj_elem_n_sites = pd.read_json(ptable_n_sites_hist_path, typ="series")
elif srs_mp_trj_elem_n_sites is None:
    # construct a series of lists of site numbers per element (i.e. how often each
    # element appears in a structure with n sites)
    # create all df cols as int dtype
    df_mp_trj_elem_n_sites = pd.DataFrame(
        [
            dict.fromkeys(set(site_nums), len(site_nums))
            for site_nums in df_mp_trj[Key.atom_nums]
        ]
    )
    mp_trj_elem_n_sites = {
        col: list(df_mp_trj_elem_n_sites[col].dropna().astype(int))
        for col in df_mp_trj_elem_n_sites
    }
    srs_mp_trj_elem_n_sites = pd.Series(mp_trj_elem_n_sites).sort_index()
    srs_mp_trj_elem_n_sites.to_json(ptable_n_sites_hist_path)


cmap = plt.get_cmap("Blues")
cbar_ticks = (100, 1_000, 10_000, 100_000, 1_000_000)
norm = matplotlib.colors.LogNorm(vmin=min(cbar_ticks), vmax=max(cbar_ticks))

fig_ptable_sites = pmv.ptable_hists(
    srs_mp_trj_elem_n_sites,
    symbol_pos=(0.8, 0.9),
    log=True,
    cbar_title="Number of Sites",
    cbar_title_kwargs=dict(fontsize=16),
    cbar_coords=(0.18, 0.85, 0.42, 0.02),
    # anno_kwds=lambda hist_vals: dict(
    #     text=pmv.si_fmt(len(hist_vals), ".0f"),
    #     xy=(0.8, 0.6),
    #     bbox=dict(pad=2, edgecolor="none", facecolor="none"),
    # ),
    x_range=(1, 300),
    # hist_kwds=lambda hist_vals: dict(
    #     color=cmap(norm(len(hist_vals))), edgecolor="none"
    # ),
)

# turn off y axis for helium (why is it even there?)
fig_ptable_sites.axes[17].get_yaxis().set_visible(b=False)

cbar_ax = fig_ptable_sites.figure.add_axes([0.23, 0.73, 0.31, 0.025])
cbar = matplotlib.colorbar.ColorbarBase(
    cbar_ax,
    cmap=cmap,
    norm=norm,
    orientation="horizontal",
    ticks=cbar_ticks,
)
cbar.set_label("Number of atoms in MPtrj structures", fontsize=16)
cbar.ax.xaxis.set_label_position("top")

pmv.save_fig(fig_ptable_sites, f"{PDF_FIGS}/mp-trj-n-sites-ptable-hists.pdf")


# %%
elem_counts: dict[str, dict[str, int]] = {}
for count_mode in ("composition", "occurrence"):
    trj_elem_counts = pmv.count_elements(df_mp_trj[Key.formula], count_mode=count_mode)
    elem_counts[count_mode] = trj_elem_counts
    filename = f"mp-trj-element-counts-by-{count_mode}"
    trj_elem_counts.to_json(f"{data_page}/{filename}.json")


# %% TODO https://github.com/janosh/pymatviz/issues/188 font sizes and box sizes
count_mode = "composition"
trj_elem_counts = pd.read_json(
    f"{data_page}/mp-trj-element-counts-by-{count_mode}.json", typ="series"
)

excl_elems = "He Ne Ar Kr Xe".split() if (excl_noble := False) else ()

ax_ptable = pmv.ptable_heatmap(  # matplotlib version looks better for SI
    trj_elem_counts,
    # zero_color="#efefef",
    log=(log := SymLogNorm(linthresh=10_000)),
    exclude_elements=excl_elems,  # drop noble gases
    # cbar_range=None if excl_noble else (10_000, None),
    show_values=(show_vals := True),
    # label_font_size=17 if show_vals else 25,
    # value_font_size=14,
    cbar_title="MPtrj Element Counts",
)

img_name = f"mp-trj-element-counts-by-{count_mode}"
if log:
    img_name += "-symlog" if isinstance(log, SymLogNorm) else "-log"
if excl_noble:
    img_name += "-excl-noble"
pmv.save_fig(ax_ptable, f"{PDF_FIGS}/{img_name}.pdf")


# %%
normalized = True
ax_ptable = pmv.ptable_heatmap_ratio(
    trj_elem_counts / (len(df_mp_trj) if normalized else 1),
    mp_occu_counts / (len(df_mp) if normalized else 1),
    zero_color="#efefef",
    fmt=".2f",
    not_in_denominator=None,
    not_in_numerator=None,
    not_in_either=None,
)

img_name = "mp-trj-mp-ratio-element-counts-by-occurrence"
if normalized:
    img_name += "-normalized"
pmv.save_fig(ax_ptable, f"{PDF_FIGS}/{img_name}.pdf")


# %% plot formation energy per atom distribution
# pdf_kwds defined to use the same figure size for all plots
fig = pmv.histogram(df_mp_trj[MbdKey.e_form_dft], bins=300, opacity=1)
if log := False:
    fig.update_yaxes(type="log")
fig.layout.xaxis.title = "E<sub>form</sub> (eV/atom)"
count_col = "Number of Structures"
fig.layout.yaxis.title = count_col
fig.show()

pdf_kwds = dict(width=500, height=300)
# pmv.save_fig(
#     fig, f"{PDF_FIGS}/mp-trj-e-form-hist{'-log' if log else ''}.pdf", **pdf_kwds
# )
# pmv.save_fig(fig, f"{SITE_FIGS}/mp-trj-e-form-hist.svelte")


# %% plot forces distribution
fig = pmv.histogram(df_mp_trj[Key.forces].explode().explode().abs(), bins=300)
fig.layout.xaxis.title = "|Forces| (eV/Å)"
fig.layout.yaxis.title = count_col
fig.update_yaxes(type="log")
fig.show()

# pmv.save_fig(fig, f"{PDF_FIGS}/mp-trj-forces-hist.pdf", **pdf_kwds)
# pmv.save_fig(fig, f"{SITE_FIGS}/mp-trj-forces-hist.svelte")


# %% plot hydrostatic stress distribution
fig = pmv.histogram(df_mp_trj[Key.stress_trace], bins=300)
fig.layout.xaxis.title = "1/3 Tr(σ) (eV/Å³)"  # noqa: RUF001
fig.layout.yaxis.title = count_col
fig.update_yaxes(type="log")
fig.show()

# pmv.save_fig(fig, f"{PDF_FIGS}/mp-trj-stresses-hist.pdf", **pdf_kwds)
# pmv.save_fig(fig, f"{SITE_FIGS}/mp-trj-stresses-hist.svelte")


# %% plot magmoms distribution
fig = pmv.histogram(df_mp_trj[Key.magmoms].dropna().explode(), bins=300)
fig.layout.xaxis.title = "Magmoms (μB)"
fig.layout.yaxis.title = count_col
fig.update_yaxes(type="log")
fig.show()

# pmv.save_fig(fig, f"{PDF_FIGS}/mp-trj-magmoms-hist.pdf", **pdf_kwds)
# pmv.save_fig(fig, f"{SITE_FIGS}/mp-trj-magmoms-hist.svelte")


# %%
for df in (df_mp_trj, df_mp, df_wbm):
    if Key.arity not in df:
        df[Key.arity] = df[Key.formula].map(Composition).map(len)


# %%
df_arity = pd.DataFrame(
    {
        key: df[Key.arity].value_counts().sort_index() / len(df)
        for key, df in (("MP", df_mp), ("MPtrj", df_mp_trj), ("WBM", df_wbm))
    }
)
df_arity = df_arity.query("0 < index < 7")

fig = px.bar(df_arity, barmode="group")
fig.update_traces(marker_line_width=0)
fig.layout.legend.update(x=1, y=1, xanchor="right", yanchor="top", title=None)
fig.layout.margin = dict(l=0, r=0, b=0, t=0)
fig.layout.yaxis.title = "Fraction of Structures in Dataset"
fig.layout.xaxis.title = "Number of Elements in Formula"

fig.show()
img_name = "mp-vs-mp-trj-vs-wbm-arity-hist"
pmv.save_fig(fig, f"{SITE_FIGS}/{img_name}.svelte")
pmv.save_fig(fig, f"{PDF_FIGS}/{img_name}.pdf", width=450, height=280)


# %% calc n_sites from per-site atomic numbers
df_mp_trj[Key.n_sites] = df_mp_trj[Key.atom_nums].map(len)
n_sites_hist, n_sites_bins = np.histogram(
    df_mp_trj[Key.n_sites], bins=range(1, df_mp_trj[Key.n_sites].max() + 1)
)

n_struct_col = "Number of Structures"
df_n_sites = pd.DataFrame({Key.n_sites: n_sites_bins[:-1], n_struct_col: n_sites_hist})
log_y = False


# %% plot n_sites distribution
fig = px.bar(df_n_sites, x=Key.n_sites, y=n_struct_col, log_y=log_y, range_x=(1, 200))
# add inset plot with log scale
fig.add_bar(
    x=df_n_sites[Key.n_sites],
    y=df_n_sites[n_struct_col],
    showlegend=False,
    xaxis="x2",
    yaxis="y2",
    marker=dict(color=fig.data[0].marker.color),  # same color as main plot
)

bin_width = n_sites_bins[1] - n_sites_bins[0]
fig.update_traces(width=bin_width, marker_line_width=0)
# add cumulative distribution as 2nd y axis
fig.add_scatter(
    x=df_n_sites[Key.n_sites],
    y=df_n_sites[n_struct_col].cumsum() / df_n_sites[n_struct_col].sum(),
    mode="lines",
    name="Cumulative",
    xaxis="x3",
    yaxis="y3",
    hovertemplate="x: %{x}<br>y: %{y:.1%}",
)
# add inset title 'log-scaled to show tail'
inset_domain = [0.4, 1]
fig.layout.xaxis2 = dict(domain=inset_domain, anchor="y2")
fig.layout.yaxis2 = dict(
    domain=inset_domain,
    anchor="x2",
    type="log",
    title="log-scaled to show tail",
    title_standoff=0,
)

fig.layout.yaxis3 = dict(  # move y3 axis to right side of y2
    overlaying="y2", side="right", tickformat=".0%"
)
fig.layout.xaxis3 = dict(overlaying="x2", visible=False)

# project line from 90% cumulative to x axis
x_90 = df_n_sites[Key.n_sites][
    (df_n_sites[n_struct_col].cumsum() / df_n_sites[n_struct_col].sum()) < 0.9
].iloc[-1]
for x0, y0, x1, y1 in (
    (x_90, 0, x_90, 0.9),
    (x_90, 0.9, df_n_sites[Key.n_sites].max(), 0.9),
):
    fig.add_shape(
        type="line",
        **dict(x0=x0, y0=y0, x1=x1, y1=y1),
        line=dict(width=1, dash="dot"),
        xref="x3",
        yref="y3",
    )
fig.layout.yaxis3.update(showgrid=False, rangemode="tozero")

fig.layout.margin = dict(l=5, r=5, b=5, t=5)
fig.layout.legend.update(x=0.96, y=0.18, xanchor="right", bgcolor="rgba(0,0,0,0)")
fig.show()

img_name = "mp-trj-n-sites-hist"
if log_y:
    img_name += "-log"
pmv.save_fig(fig, f"{SITE_FIGS}/{img_name}.svelte")
# pmv.save_fig(fig, f"{PDF_FIGS}/{img_name}.pdf", width=450, height=300)
