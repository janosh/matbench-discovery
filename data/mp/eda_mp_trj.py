"""MPtrj exploratory data analysis (EDA)."""


# %%
import io
from zipfile import ZipFile

import ase
import ase.io.extxyz
import numpy as np
import pandas as pd
import plotly.express as px
from pymatviz import count_elements, ptable_heatmap, ptable_heatmap_ratio
from pymatviz.io import save_fig
from pymatviz.utils import si_fmt
from tqdm import tqdm

from matbench_discovery import (
    DATA_DIR,
    PDF_FIGS,
    ROOT,
    SITE_FIGS,
    formula_col,
    stress_col,
    stress_trace_col,
)
from matbench_discovery.data import DATA_FILES
from matbench_discovery.plots import quantity_labels

__author__ = "Janosh Riebesell"
__date__ = "2023-11-22"

data_page = f"{ROOT}/site/src/routes/data"
e_form_per_atom_col = "ef_per_atom"
magmoms_col = "magmoms"
forces_col = "forces"


# %% downloaded mptrj-gga-ggapu.tar.gz from https://drive.google.com/drive/folders/1JQ-ry1RHvNliVg1Ut5OuyUxne51RHiT_
# and extracted the mptrj-gga-ggapu directory (6.2 GB) to data/mp using macOS Finder
# then zipped it to mp-trj-extxyz.zip (also using Finder, 1.6 GB)
zip_path = f"{DATA_DIR}/mp/mp-trj-extxyz-by-yuan.zip"
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


assert len(mp_trj_atoms) == 145_919


# %%
df_mp_trj = pd.DataFrame(
    {
        f"{atm.info['task_id']}-{atm.info['calc_id']}-{atm.info['ionic_step']}": {
            "formula": str(atm.symbols)
        }
        | {key: atm.arrays.get(key) for key in ("forces", "magmoms")}
        | atm.info
        for atoms_list in mp_trj_atoms.values()
        for atm in atoms_list
    }
).T.convert_dtypes()  # convert object columns to float/int where possible
df_mp_trj.index.name = "frame_id"
assert len(df_mp_trj) == 1_580_312
assert formula_col in df_mp_trj

# this is the unrelaxed (but MP2020 corrected) formation energy per atom of the actual
# relaxation step
df_mp_trj[stress_trace_col] = [
    np.trace(stress) / 3 for stress in tqdm(df_mp_trj[stress_col])
]


# %%
df_mp_trj.to_json(f"{DATA_DIR}/mp/mp-trj-2022-09-summary.json.bz2")


# %% load MPtrj summary data
df_mp_trj = pd.read_json(f"{DATA_DIR}/mp/mp-trj-2022-09-summary.json.bz2")
df_mp_trj.index.name = "frame_id"


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
if "trj_elem_counts" not in locals():
    trj_elem_counts = pd.read_json(
        f"{data_page}/mp-trj-element-counts-by-occurrence.json", typ="series"
    )

excl_elems = "He Ne Ar Kr Xe".split() if (excl_noble := True) else ()

ax_ptable = ptable_heatmap(  # matplotlib version looks better for SI
    trj_elem_counts,
    fmt=lambda x, _: si_fmt(x, ".1f"),
    cbar_fmt=lambda x, _: si_fmt(x, ".0f"),
    zero_color="#efefef",
    log=(log := True),
    # drop noble gases
    exclude_elements=excl_elems,
)

img_name = f"mp-trj-element-counts-by-occurrence{'-log' if log else ''}"
if excl_noble:
    img_name += "-excl-noble"
save_fig(ax_ptable, f"{PDF_FIGS}/{img_name}.pdf")


# %% load MP element counts by occurrence to compute ratio with MPtrj
mp_occu_counts = pd.read_json(
    f"{data_page}/mp-element-counts-by-occurrence.json", typ="series"
)
df_mp = pd.read_csv(DATA_FILES.mp_energies, na_filter=False)


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

img_name = (
    f"mp-trj-mp-ratio-element-counts-by-occurrence-{'normalized' if normalized else ''}"
)
save_fig(ax_ptable, f"{PDF_FIGS}/{img_name}.pdf")


# %% plot formation energy per atom distribution
count_col = "Number of structures"
pdf_kwds = dict(width=500, height=300)
if "e_form_per_atom_hist" not in locals():  # only compute once for speed
    e_form_per_atom_hist, e_form_per_atom_bins = np.histogram(
        df_mp_trj[e_form_per_atom_col], bins=300
    )
fig = px.bar(
    x=e_form_per_atom_bins[:-1],
    y=e_form_per_atom_hist,
    log_y=True,
)
bin_width = (e_form_per_atom_bins[1] - e_form_per_atom_bins[0]) * 1.2
fig.update_traces(width=bin_width, marker_line_width=0)
fig.layout.xaxis.update(linecolor="lightgray", title="E<sub>form</sub> (eV/atom)")
fig.layout.yaxis.update(linecolor="lightgray", title=count_col)
fig.layout.margin = dict(l=5, r=5, b=5, t=5)
fig.show()
save_fig(fig, f"{PDF_FIGS}/mp-trj-e-form-hist.pdf", **pdf_kwds)
save_fig(fig, f"{SITE_FIGS}/mp-trj-e-form-hist.svelte")


# %% plot forces distribution
# use numpy to pre-compute histogram
if "forces_hist" not in locals():  # only compute once for speed
    forces_hist, forces_bins = np.histogram(
        df_mp_trj[forces_col].abs().explode().explode(), bins=300
    )
fig = px.bar(x=forces_bins[:-1], y=forces_hist, log_y=True)
bin_width = (forces_bins[1] - forces_bins[0]) * 1.2
fig.update_traces(width=bin_width, marker_line_width=0)
fig.layout.xaxis.update(linecolor="lightgray", title="|Forces| (eV/Å)")
fig.layout.yaxis.update(linecolor="lightgray", title=count_col)
fig.layout.margin = dict(l=5, r=5, b=5, t=5)
fig.show()
save_fig(fig, f"{PDF_FIGS}/mp-trj-forces-hist.pdf", **pdf_kwds)
save_fig(fig, f"{SITE_FIGS}/mp-trj-forces-hist.svelte")


# %% plot hydrostatic stress distribution
if "stress_hist" not in locals():  # only compute once for speed
    stress_hist, stress_bins = np.histogram(df_mp_trj[stress_trace_col], bins=300)
fig = px.bar(x=stress_bins[:-1], y=stress_hist, log_y=True)
bin_width = (stress_bins[1] - stress_bins[0]) * 1.2
fig.update_traces(width=bin_width, marker_line_width=0)
fig.layout.xaxis.update(linecolor="lightgray", title=quantity_labels[stress_trace_col])
fig.layout.yaxis.update(linecolor="lightgray", title=count_col)
fig.layout.margin = dict(l=5, r=5, b=5, t=5)
fig.show()
save_fig(fig, f"{PDF_FIGS}/mp-trj-stresses-hist.pdf", **pdf_kwds)
save_fig(fig, f"{SITE_FIGS}/mp-trj-stresses-hist.svelte")


# %% plot magmoms distribution
if "magmoms_hist" not in locals():  # only compute once for speed
    magmoms_hist, magmoms_bins = np.histogram(
        df_mp_trj[magmoms_col].dropna().explode(), bins=300
    )

fig = px.bar(x=magmoms_bins[:-1], y=magmoms_hist, log_y=True)
bin_width = (magmoms_bins[1] - magmoms_bins[0]) * 1.2
fig.update_traces(width=bin_width, marker_line_width=0)
fig.layout.xaxis.update(linecolor="lightgray", title="Magmoms (μ<sub>B</sub>)")
fig.layout.yaxis.update(linecolor="lightgray", title=count_col)
fig.layout.margin = dict(l=5, r=5, b=5, t=5)
fig.show()
save_fig(fig, f"{PDF_FIGS}/mp-trj-magmoms-hist.pdf", **pdf_kwds)
save_fig(fig, f"{SITE_FIGS}/mp-trj-magmoms-hist.svelte")
