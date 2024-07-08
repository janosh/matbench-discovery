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
from pymatviz import (
    count_elements,
    plot_histogram,
    ptable_heatmap,
    ptable_heatmap_ratio,
    ptable_hists,
)
from pymatviz.enums import Key
from pymatviz.io import save_fig
from pymatviz.utils import si_fmt
from tqdm import tqdm

from matbench_discovery import MP_DIR, PDF_FIGS, ROOT, SITE_FIGS
from matbench_discovery.data import DATA_FILES, df_wbm
from matbench_discovery.enums import MbdKey

__author__ = "Janosh Riebesell"
__date__ = "2023-11-22"

data_page = f"{ROOT}/site/src/routes/data"


# %% load MP element counts by occurrence to compute ratio with MPtrj
mp_occu_counts = pd.read_json(
    f"{data_page}/mp-element-counts-by-occurrence.json", typ="series"
)
df_mp = pd.read_csv(DATA_FILES.mp_energies, na_filter=False).set_index(Key.mat_id)


# %% --- load preprocessed MPtrj summary data if available ---
mp_trj_summary_path = f"{MP_DIR}/mp-trj-2022-09-summary.json.bz2"
if os.path.isfile(mp_trj_summary_path):
    df_mp_trj = pd.read_json(mp_trj_summary_path)
    df_mp_trj.index.name = "frame_id"
else:
    print("MPtrj summary data not found, run cell below to generate")


# %% downloaded mptrj-gga-ggapu.tar.gz from https://drive.google.com/drive/folders/1JQ-ry1RHvNliVg1Ut5OuyUxne51RHiT_
# and extracted the mptrj-gga-ggapu directory (6.2 GB) to data/mp using macOS Finder
# then zipped it to mp-trj-extxyz.zip (also using Finder, 1.6 GB)
zip_path = f"{MP_DIR}/2023-11-22-mp-trj-extxyz-by-yuan.zip"
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
def info_dict_to_id(info: dict[str, int | str]) -> str:
    """Construct a unique frame ID from the atoms info dict."""
    return f"{info[Key.task_id]}-{info['calc_id']}-{info['ionic_step']}"


df_mp_trj = pd.DataFrame(
    {
        info_dict_to_id(atoms.info): atoms.info
        | {key: atoms.arrays.get(key) for key in ("forces", "magmoms")}
        | {"formula": str(atoms.symbols), Key.atom_nums: atoms.symbols}
        for atoms_list in tqdm(mp_trj_atoms.values(), total=len(mp_trj_atoms))
        for atoms in atoms_list
    }
).T.convert_dtypes()  # convert object columns to float/int where possible
df_mp_trj.index.name = "frame_id"
assert len(df_mp_trj) == 1_580_312  # number of total frames
if Key.formula not in df_mp_trj:
    raise KeyError(f"{Key.formula!s} not in {df_mp_trj.columns=}")

# this is the unrelaxed (but MP2020 corrected) formation energy per atom of the actual
# relaxation step
df_mp_trj = df_mp_trj.rename(columns={"ef_per_atom": MbdKey.e_form_dft})
df_mp_trj[Key.stress_trace] = [
    np.trace(stress) / 3 for stress in tqdm(df_mp_trj[Key.stress])
]


# %%
df_mp_trj.to_json(mp_trj_summary_path)


# %%
def tile_count_anno(hist_vals: list[Any]) -> dict[str, Any]:
    """Annotate each periodic table tile with the number of values in its histogram."""
    face_color = cmap(norm(np.sum(len(hist_vals)))) if hist_vals else "none"
    bbox = dict(facecolor=face_color, alpha=0.4, pad=2, edgecolor="none")
    return dict(text=si_fmt(len(hist_vals), ".0f"), bbox=bbox)


# %% plot per-element magmom histograms
ptable_magmom_hist_path = f"{MP_DIR}/mp-trj-2022-09-elem-magmoms.json.bz2"
srs_mp_trj_elem_magmoms = locals().get("srs_mp_trj_elem_magmoms")

if os.path.isfile(ptable_magmom_hist_path):
    srs_mp_trj_elem_magmoms = pd.read_json(ptable_magmom_hist_path, typ="series")
if srs_mp_trj_elem_magmoms is None:
    # project magmoms onto symbols in dict
    df_mp_trj_elem_magmom = pd.DataFrame(
        [
            dict(zip(elems, magmoms, strict=False))
            for elems, magmoms in df_mp_trj.set_index(Key.atom_nums)[Key.magmoms]
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
ptable_force_hist_path = f"{MP_DIR}/mp-trj-2022-09-elem-forces.json.bz2"
srs_mp_trj_elem_forces = locals().get("srs_mp_trj_elem_forces")

if os.path.isfile(ptable_force_hist_path):
    srs_mp_trj_elem_forces = pd.read_json(ptable_force_hist_path, typ="series")
if srs_mp_trj_elem_forces is None:
    df_mp_trj_elem_forces = pd.DataFrame(
        [
            dict(zip(elems, np.abs(forces).mean(axis=1), strict=False))
            for elems, forces in df_mp_trj.set_index(Key.atom_nums)[Key.forces].items()
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
ptable_n_sites_hist_path = f"{MP_DIR}/mp-trj-2022-09-elem-n-sites.json.bz2"
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
fig_ptable_sites.axes[17].get_yaxis().set_visible(b=False)

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
        df_mp_trj[Key.formula], count_mode=count_mode
    ).astype(int)
    elem_counts[count_mode] = trj_elem_counts
    filename = f"mp-trj-element-counts-by-{count_mode}"
    trj_elem_counts.to_json(f"{data_page}/{filename}.json")


# %%
count_mode = "composition"
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
    show_values=(show_vals := True),
    label_font_size=17 if show_vals else 25,
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
# pdf_kwds defined to use the same figure size for all plots
fig = plot_histogram(df_mp_trj[MbdKey.e_form_dft], bins=300)
# fig.update_yaxes(type="log")
fig.layout.xaxis.title = "E<sub>form</sub> (eV/atom)"
count_col = "Number of Structures"
fig.layout.yaxis.title = count_col
fig.show()

pdf_kwds = dict(width=500, height=300)
# save_fig(fig, f"{PDF_FIGS}/mp-trj-e-form-hist.pdf", **pdf_kwds)
# save_fig(fig, f"{SITE_FIGS}/mp-trj-e-form-hist.svelte")


# %% plot forces distribution
fig = plot_histogram(df_mp_trj[Key.forces].explode().explode().abs(), bins=300)
fig.layout.xaxis.title = "|Forces| (eV/Å)"
fig.layout.yaxis.title = count_col
fig.update_yaxes(type="log")
fig.show()

# save_fig(fig, f"{PDF_FIGS}/mp-trj-forces-hist.pdf", **pdf_kwds)
# save_fig(fig, f"{SITE_FIGS}/mp-trj-forces-hist.svelte")


# %% plot hydrostatic stress distribution
fig = plot_histogram(df_mp_trj[Key.stress_trace], bins=300)
fig.layout.xaxis.title = "1/3 Tr(σ) (eV/Å³)"  # noqa: RUF001
fig.layout.yaxis.title = count_col
fig.update_yaxes(type="log")
fig.show()

# save_fig(fig, f"{PDF_FIGS}/mp-trj-stresses-hist.pdf", **pdf_kwds)
# save_fig(fig, f"{SITE_FIGS}/mp-trj-stresses-hist.svelte")


# %% plot magmoms distribution
fig = plot_histogram(df_mp_trj[Key.magmoms].dropna().explode(), bins=300)
fig.layout.xaxis.title = "Magmoms (μB)"
fig.layout.yaxis.title = count_col
fig.update_yaxes(type="log")
fig.show()

# save_fig(fig, f"{PDF_FIGS}/mp-trj-magmoms-hist.pdf", **pdf_kwds)
# save_fig(fig, f"{SITE_FIGS}/mp-trj-magmoms-hist.svelte")


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
save_fig(fig, f"{SITE_FIGS}/{img_name}.svelte")
save_fig(fig, f"{PDF_FIGS}/{img_name}.pdf", width=450, height=280)


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
save_fig(fig, f"{SITE_FIGS}/{img_name}.svelte")
# save_fig(fig, f"{PDF_FIGS}/{img_name}.pdf", width=450, height=300)
