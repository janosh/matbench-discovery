"""MPtrj exploratory data analysis (EDA)."""


# %%
import io
import os
from typing import Any
from zipfile import ZipFile

import ase
import ase.io.extxyz
import matplotlib.colors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.express as px
from matplotlib.colors import SymLogNorm
from pymatgen.core import Composition, Element
from pymatviz import count_elements, ptable_heatmap, ptable_heatmap_ratio, ptable_hists
from pymatviz.io import save_fig
from pymatviz.utils import si_fmt
from tqdm import tqdm

from matbench_discovery import (
    DATA_DIR,
    PDF_FIGS,
    ROOT,
    SITE_FIGS,
    formula_col,
    id_col,
    n_sites_col,
    stress_col,
    stress_trace_col,
)
from matbench_discovery.data import DATA_FILES, df_wbm

__author__ = "Janosh Riebesell"
__date__ = "2023-11-22"

data_page = f"{ROOT}/site/src/routes/data"
e_form_per_atom_col = "ef_per_atom"
magmoms_col = "magmoms"
forces_col = "forces"
site_nums_col = "site_nums"


# %% load MP element counts by occurrence to compute ratio with MPtrj
mp_occu_counts = pd.read_json(
    f"{data_page}/mp-element-counts-by-occurrence.json", typ="series"
)
df_mp = pd.read_csv(DATA_FILES.mp_energies, na_filter=False).set_index(id_col)


# %% --- load preprocessed MPtrj summary data if available ---
mp_trj_summary_path = f"{DATA_DIR}/mp/mp-trj-2022-09-summary.json.bz2"
if os.path.isfile(mp_trj_summary_path):
    df_mp_trj = pd.read_json(mp_trj_summary_path)
    df_mp_trj.index.name = "frame_id"
else:
    print("MPtrj summary data not found, run cell below to generate")


# %% downloaded mptrj-gga-ggapu.tar.gz from https://drive.google.com/drive/folders/1JQ-ry1RHvNliVg1Ut5OuyUxne51RHiT_
# and extracted the mptrj-gga-ggapu directory (6.2 GB) to data/mp using macOS Finder
# then zipped it to mp-trj-extxyz.zip (also using Finder, 1.6 GB)
zip_path = f"{DATA_DIR}/mp/2023-11-22-mp-trj-extxyz-by-yuan.zip"
mp_trj_atoms: dict[str, list[ase.Atoms]] = {}

# extract extXYZ files from zipped directory without unpacking the whole archive
# takes ~8 min on M2 Max
for name in tqdm((zip_file := ZipFile(zip_path)).namelist()):
    if name.startswith("mptrj-gga-ggapu/mp-"):
        mp_id = name.split("/")[1].split(".")[0]
        assert mp_id.startswith("mp-")
        assert mp_id not in mp_trj_atoms

        with zip_file.open(name) as file:
            # wrap byte stream with TextIOWrapper to use as file
            text_file = io.TextIOWrapper(file, encoding="utf-8")
            atoms_list = list(ase.io.extxyz.read_xyz(text_file, index=slice(None)))
        mp_trj_atoms[mp_id] = atoms_list


assert len(mp_trj_atoms) == 145_919  # number of unique MP IDs


# %%
info_to_id = lambda info: f"{info['task_id']}-{info['calc_id']}-{info['ionic_step']}"

df_mp_trj = pd.DataFrame(
    {
        info_to_id(atoms.info): atoms.info
        | {key: atoms.arrays.get(key) for key in ("forces", "magmoms")}
        | {"formula": str(atoms.symbols), site_nums_col: atoms.symbols}
        for atoms_list in tqdm(mp_trj_atoms.values(), total=len(mp_trj_atoms))
        for atoms in atoms_list
    }
).T.convert_dtypes()  # convert object columns to float/int where possible
df_mp_trj.index.name = "frame_id"
assert len(df_mp_trj) == 1_580_312  # number of total frames
assert formula_col in df_mp_trj

# this is the unrelaxed (but MP2020 corrected) formation energy per atom of the actual
# relaxation step
df_mp_trj[stress_trace_col] = [
    np.trace(stress) / 3 for stress in tqdm(df_mp_trj[stress_col])
]


# %%
df_mp_trj.to_json(mp_trj_summary_path)


# %%
def tile_count_anno(hist_vals: list[Any]) -> dict[str, Any]:
    """Annotate each periodic table tile with the number of values in its histogram."""
    facecolor = cmap(norm(np.sum(len(hist_vals)))) if hist_vals else "none"
    bbox = dict(facecolor=facecolor, alpha=0.4, pad=2, edgecolor="none")
    return dict(text=si_fmt(len(hist_vals), ".0f"), bbox=bbox)


# %% plot per-element magmom histograms
ptable_magmom_hist_path = f"{DATA_DIR}/mp/mp-trj-2022-09-elem-magmoms.json.bz2"

if os.path.isfile(ptable_magmom_hist_path):
    srs_mp_trj_elem_magmoms = pd.read_json(ptable_magmom_hist_path, typ="series")
elif "srs_mp_trj_elem_magmoms" not in locals():
    # project magmoms onto symbols in dict
    df_mp_trj_elem_magmom = pd.DataFrame(
        [
            dict(zip(elems, magmoms))
            for elems, magmoms in df_mp_trj.set_index(site_nums_col)[magmoms_col]
            .dropna()
            .items()
        ]
    )

    srs_mp_trj_elem_magmoms = {
        col: list(df_mp_trj_elem_magmom[col].dropna()) for col in df_mp_trj_elem_magmom
    }
    pd.Series(srs_mp_trj_elem_magmoms).to_json(ptable_magmom_hist_path)

cmap = plt.get_cmap(color_map := "viridis")
norm = matplotlib.colors.LogNorm(vmin=1, vmax=150_000)

fig_ptable_magmoms = ptable_hists(
    srs_mp_trj_elem_magmoms,
    symbol_pos=(0.2, 0.8),
    log=True,
    cbar_title="Magmoms ($μ_B$)",
    cbar_title_kwds=dict(fontsize=16),
    cbar_coords=(0.18, 0.85, 0.42, 0.02),
    # annotate each element with its number of magmoms in MPtrj
    anno_kwds=tile_count_anno,
    colormap=color_map,
)

cbar_ax = fig_ptable_magmoms.figure.add_axes([0.26, 0.78, 0.25, 0.015])
cbar = matplotlib.colorbar.ColorbarBase(
    cbar_ax, cmap=cmap, norm=norm, orientation="horizontal"
)
save_fig(fig_ptable_magmoms, f"{PDF_FIGS}/mp-trj-magmoms-ptable-hists.pdf")


# %% plot per-element force histograms
ptable_force_hist_path = f"{DATA_DIR}/mp/mp-trj-2022-09-elem-forces.json.bz2"

if os.path.isfile(ptable_force_hist_path):
    srs_mp_trj_elem_forces = pd.read_json(ptable_force_hist_path, typ="series")
elif "srs_mp_trj_elem_forces" not in locals():
    df_mp_trj_elem_forces = pd.DataFrame(
        [
            dict(zip(elems, np.abs(forces).mean(axis=1)))
            for elems, forces in df_mp_trj.set_index(site_nums_col)[forces_col].items()
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
fig_ptable_forces = ptable_hists(
    srs_mp_trj_elem_forces.copy().map(lambda x: [val for val in x if val < max_force]),
    symbol_pos=(0.3, 0.8),
    log=True,
    cbar_title="1/3 Σ|Forces| (eV/Å)",
    cbar_title_kwds=dict(fontsize=16),
    cbar_coords=(0.18, 0.85, 0.42, 0.02),
    x_range=(0, max_force),
    anno_kwds=tile_count_anno,
    colormap=color_map,
)

cbar_ax = fig_ptable_forces.figure.add_axes([0.26, 0.78, 0.25, 0.015])
cbar = matplotlib.colorbar.ColorbarBase(
    cbar_ax, cmap=cmap, norm=norm, orientation="horizontal"
)

save_fig(fig_ptable_forces, f"{PDF_FIGS}/mp-trj-forces-ptable-hists.pdf")


# %% plot histogram of number of sites per element
ptable_n_sites_hist_path = f"{DATA_DIR}/mp/mp-trj-2022-09-elem-n-sites.json.bz2"

if os.path.isfile(ptable_n_sites_hist_path):
    srs_mp_trj_elem_n_sites = pd.read_json(ptable_n_sites_hist_path, typ="series")
elif "mp_trj_elem_n_sites" not in locals():
    # construct a series of lists of site numbers per element (i.e. how often each
    # element appears in a structure with n sites)
    # create all df cols as int dtype
    df_mp_trj_elem_n_sites = pd.DataFrame(
        [
            dict.fromkeys(set(site_nums), len(site_nums))
            for site_nums in df_mp_trj[site_nums_col]
        ]
    ).astype(int)
    mp_trj_elem_n_sites = {
        col: list(df_mp_trj_elem_n_sites[col].dropna())
        for col in df_mp_trj_elem_n_sites
    }
    srs_mp_trj_elem_n_sites = pd.Series(mp_trj_elem_n_sites).sort_index()

    srs_mp_trj_elem_n_sites.index = srs_mp_trj_elem_n_sites.index.map(
        Element.from_Z
    ).map(str)
    srs_mp_trj_elem_n_sites.to_json(ptable_n_sites_hist_path)


cmap = plt.get_cmap("Blues")
cbar_ticks = (100, 1_000, 10_000, 100_000, 1_000_000)
norm = matplotlib.colors.LogNorm(vmin=min(cbar_ticks), vmax=max(cbar_ticks))

fig_ptable_sites = ptable_hists(
    srs_mp_trj_elem_n_sites,
    symbol_pos=(0.8, 0.9),
    log=True,
    cbar_title="Number of Sites",
    cbar_title_kwds=dict(fontsize=16),
    cbar_coords=(0.18, 0.85, 0.42, 0.02),
    anno_kwds=lambda hist_vals: dict(
        text=si_fmt(len(hist_vals), ".0f"),
        xy=(0.8, 0.6),
        bbox=dict(pad=2, edgecolor="none", facecolor="none"),
    ),
    x_range=(1, 300),
    hist_kwds=lambda hist_vals: dict(
        color=cmap(norm(len(hist_vals))), edgecolor="none"
    ),
)

# turn off y axis for helium (why is it even there?)
fig_ptable_sites.axes[17].get_yaxis().set_visible(False)

cbar_ax = fig_ptable_sites.figure.add_axes([0.23, 0.8, 0.31, 0.025])
cbar = matplotlib.colorbar.ColorbarBase(
    cbar_ax,
    cmap=cmap,
    norm=norm,
    orientation="horizontal",
    ticks=cbar_ticks,
)
cbar.set_label("Number of atoms in MPtrj structures", fontsize=16)
cbar.ax.xaxis.set_label_position("top")

save_fig(fig_ptable_sites, f"{PDF_FIGS}/mp-trj-n-sites-ptable-hists.pdf")


# %%
elem_counts: dict[str, dict[str, int]] = {}
for count_mode in ("composition", "occurrence"):
    trj_elem_counts = count_elements(
        df_mp_trj[formula_col], count_mode=count_mode
    ).astype(int)
    elem_counts[count_mode] = trj_elem_counts
    filename = f"mp-trj-element-counts-by-{count_mode}"
    trj_elem_counts.to_json(f"{data_page}/{filename}.json")


# %%
count_mode = "occurrence"
trj_elem_counts = pd.read_json(
    f"{data_page}/mp-trj-element-counts-by-{count_mode}.json", typ="series"
)

excl_elems = "He Ne Ar Kr Xe".split() if (excl_noble := False) else ()

ax_ptable = ptable_heatmap(  # matplotlib version looks better for SI
    trj_elem_counts,
    zero_color="#efefef",
    log=(log := SymLogNorm(linthresh=10_000)),
    exclude_elements=excl_elems,  # drop noble gases
    # cbar_range=None if excl_noble else (10_000, None),
    label_font_size=17,
    value_font_size=14,
    cbar_title="MPtrj Element Counts",
)

img_name = f"mp-trj-element-counts-by-{count_mode}"
if log:
    img_name += "-symlog" if isinstance(log, SymLogNorm) else "-log"
if excl_noble:
    img_name += "-excl-noble"
save_fig(ax_ptable, f"{PDF_FIGS}/{img_name}.pdf")


# %%
normalized = True
ax_ptable = ptable_heatmap_ratio(
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
save_fig(ax_ptable, f"{PDF_FIGS}/{img_name}.pdf")


# %% plot formation energy per atom distribution
count_col = "Number of Structures"
axes_kwds = dict(linewidth=1, ticks="outside")
pdf_kwds = dict(width=500, height=300)

x_col, y_col = "E<sub>form</sub> (eV/atom)", count_col

if "df_e_form" not in locals():  # only compute once for speed
    e_form_hist = np.histogram(df_mp_trj[e_form_per_atom_col], bins=300)
    df_e_form = pd.DataFrame(e_form_hist, index=[y_col, x_col]).T.round(3)

fig = px.bar(df_e_form, x=x_col, y=count_col, log_y=True)

bin_width = df_e_form[x_col].diff().iloc[-1] * 1.2
fig.update_traces(width=bin_width, marker_line_width=0)
fig.layout.xaxis.update(**axes_kwds)
fig.layout.yaxis.update(**axes_kwds)
fig.layout.margin = dict(l=5, r=5, b=5, t=5)
fig.show()
save_fig(fig, f"{PDF_FIGS}/mp-trj-e-form-hist.pdf", **pdf_kwds)
save_fig(fig, f"{SITE_FIGS}/mp-trj-e-form-hist.svelte")


# %% plot forces distribution
# use numpy to pre-compute histogram
x_col, y_col = "|Forces| (eV/Å)", count_col

if "df_forces" not in locals():  # only compute once for speed
    forces_hist = np.histogram(
        df_mp_trj[forces_col].explode().explode().abs(), bins=300
    )
    df_forces = pd.DataFrame(forces_hist, index=[y_col, x_col]).T.round(3)

fig = px.bar(df_forces, x=x_col, y=count_col, log_y=True)

bin_width = df_forces[x_col].diff().iloc[-1] * 1.2
fig.update_traces(width=bin_width, marker_line_width=0)
fig.layout.xaxis.update(**axes_kwds)
fig.layout.yaxis.update(**axes_kwds)
fig.layout.margin = dict(l=5, r=5, b=5, t=5)
fig.show()
save_fig(fig, f"{PDF_FIGS}/mp-trj-forces-hist.pdf", **pdf_kwds)
save_fig(fig, f"{SITE_FIGS}/mp-trj-forces-hist.svelte")


# %% plot hydrostatic stress distribution
x_col, y_col = "1/3 Tr(σ) (eV/Å³)", count_col  # noqa: RUF001

if "df_stresses" not in locals():  # only compute once for speed
    stresses_hist = np.histogram(df_mp_trj[stress_trace_col], bins=300)
    df_stresses = pd.DataFrame(stresses_hist, index=[y_col, x_col]).T.round(3)

fig = px.bar(df_stresses, x=x_col, y=y_col, log_y=True)

bin_width = (df_stresses[x_col].diff().mean()) * 1.2
fig.update_traces(width=bin_width, marker_line_width=0)
fig.layout.xaxis.update(**axes_kwds)
fig.layout.yaxis.update(**axes_kwds)
fig.layout.margin = dict(l=5, r=5, b=5, t=5)
fig.show()

save_fig(fig, f"{PDF_FIGS}/mp-trj-stresses-hist.pdf", **pdf_kwds)
save_fig(fig, f"{SITE_FIGS}/mp-trj-stresses-hist.svelte")


# %% plot magmoms distribution
x_col, y_col = "Magmoms (μ<sub>B</sub>)", count_col

if "df_magmoms" not in locals():  # only compute once for speed
    magmoms_hist = np.histogram(df_mp_trj[magmoms_col].dropna().explode(), bins=300)
    df_magmoms = pd.DataFrame(magmoms_hist, index=[y_col, x_col]).T.round(3)

fig = px.bar(df_magmoms, x=x_col, y=y_col, log_y=True)

bin_width = df_magmoms[x_col].diff().iloc[-1] * 1.2
fig.update_traces(width=bin_width, marker_line_width=0)
fig.layout.xaxis.update(**axes_kwds)
fig.layout.yaxis.update(**axes_kwds)
fig.layout.margin = dict(l=5, r=5, b=5, t=5)
fig.show()
save_fig(fig, f"{PDF_FIGS}/mp-trj-magmoms-hist.pdf", **pdf_kwds)
save_fig(fig, f"{SITE_FIGS}/mp-trj-magmoms-hist.svelte")


# %%
arity_col = "arity"
for df in (df_mp_trj, df_mp, df_wbm):
    if arity_col not in df:
        df[arity_col] = df[formula_col].map(Composition).map(len)


# %%
df_arity = pd.DataFrame(
    {
        key: df[arity_col].value_counts().sort_index() / len(df)
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
save_fig(fig, f"{SITE_FIGS}/{img_name}.svelte")
save_fig(fig, f"{PDF_FIGS}/{img_name}.pdf", width=450, height=280)


# %% calc n_sites from per-site atomic numbers
df_mp_trj[n_sites_col] = df_mp_trj[site_nums_col].map(len)
n_sites_hist, n_sites_bins = np.histogram(
    df_mp_trj[n_sites_col], bins=range(1, df_mp_trj[n_sites_col].max() + 1)
)

n_struct_col = "Number of Structures"
df_n_sites = pd.DataFrame({n_sites_col: n_sites_bins[:-1], n_struct_col: n_sites_hist})
log_y = False


# %% plot n_sites distribution
fig = px.bar(df_n_sites, x=n_sites_col, y=n_struct_col, log_y=log_y, range_x=(1, 200))
# add inset plot with log scale
fig.add_bar(
    x=df_n_sites[n_sites_col],
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
    x=df_n_sites[n_sites_col],
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

fig.layout.margin = dict(l=5, r=5, b=5, t=5)
fig.layout.legend.update(x=0.96, y=0.25, xanchor="right")
fig.show()

img_name = "mp-trj-n-sites-hist"
if log_y:
    img_name += "-log"
save_fig(fig, f"{SITE_FIGS}/{img_name}.svelte")
# save_fig(fig, f"{PDF_FIGS}/{img_name}.pdf", width=450, height=300)
