"""Analyze structures and composition with largest mean error across all models.
Maybe there's some chemistry/region of materials space that all models struggle with?
Might point to deficiencies in the data or models architecture.
"""

# %%
import itertools

import matplotlib.pyplot as plt
import pandas as pd
import plotly.express as px
import plotly.graph_objs as go
import pymatviz as pmv
from pymatgen.core import Composition, Structure
from pymatviz.enums import Key
from pymatviz.utils import PLOTLY
from tqdm import tqdm

from matbench_discovery import PDF_FIGS, SITE_FIGS, WBM_DIR
from matbench_discovery.data import DataFiles, df_wbm
from matbench_discovery.enums import MbdKey
from matbench_discovery.metrics import classify_stable
from matbench_discovery.preds import df_each_err, df_each_pred, df_metrics, df_preds

__author__ = "Janosh Riebesell"
__date__ = "2023-02-15"

models = list(df_each_pred)
fp_diff_col = "site_stats_fingerprint_init_final_norm_diff"


# %%
df_cse = pd.read_json(DataFiles.wbm_cses_plus_init_structs.path).set_index(Key.mat_id)


# %% plot the highest and lowest error structures before and after relaxation
n_rows, n_cols = 5, 4
for good_or_bad, init_or_final in itertools.product(
    ("best", "worst"), ("initial", "final")
):
    fig, axs = plt.subplots(n_rows, n_cols, figsize=(4 * n_cols, 3 * n_rows))
    n_structs = len(axs.flat)
    struct_col = {"initial": Key.init_struct, "final": Key.cse}[init_or_final]

    errors = {
        "best": df_each_err[MbdKey.each_err_models].nsmallest(n_structs),
        "worst": df_each_err[MbdKey.each_err_models].nlargest(n_structs),
    }[good_or_bad]
    title = (
        f"{good_or_bad.title()} {len(errors)} {init_or_final} structures (across "
        f"{len(list(df_each_pred))} models)\nErrors in (ev/atom)"
    )
    fig.suptitle(title, fontsize=20, fontweight="bold", y=1.05)

    for idx, (mat_id, error) in enumerate(errors.items(), start=1):
        struct = df_cse[struct_col].loc[mat_id]
        if "structure" in struct:
            struct = struct["structure"]
        struct = Structure.from_dict(struct)
        ax = pmv.structure_2d(struct, ax=axs.flat[idx - 1])
        _, spg_num = struct.get_space_group_info()
        formula = struct.composition.reduced_formula
        ax_title = f"{idx}. {formula} (spg={spg_num})\n{mat_id} {error=:.2f}"
        ax.set_title(ax_title, fontweight="bold")
    out_path = f"{PDF_FIGS}/{good_or_bad}-{len(errors)}-structures-{init_or_final}.webp"
    # fig.savefig(out_path, dpi=300)


# %%
n_structs = 1000
fig = go.Figure()
for idx, model in enumerate((MbdKey.each_err_models, *df_metrics)):
    large_errors = df_each_err[model].abs().nlargest(n_structs)
    small_errors = df_each_err[model].abs().nsmallest(n_structs)
    for label, errors in (("min", large_errors), ("max", small_errors)):
        fig.add_histogram(
            x=df_wbm.loc[errors.index][fp_diff_col].values,
            name=f"{model} err<sub>{label}</sub>",
            visible="legendonly" if idx else True,
            legendgroup=model,
            hovertemplate=("SSFP diff: %{x:.2f}<br>Count: %{y}"),
        )

title = (
    f"Norm-diff between initial/final SiteStatsFingerprint<br>"
    f"of the {n_structs} highest and lowest error structures for each model"
)
fig.layout.title.update(text=title, xanchor="center", x=0.5)
fig.layout.legend.update(
    title="",
    yanchor="top",
    y=0.98,
    xanchor="right",
    x=0.98,
    font_size=12,
    orientation="h",
)
fig.layout.xaxis.title = "|SSFP<sub>initial</sub> - SSFP<sub>final</sub>|"
fig.layout.yaxis.title = "Count"

fig.show()

# pmv.save_fig(fig, f"{FIGS}/hist-largest-each-errors-fp-diff-models.svelte")


# %%
n_structs = 100
fig = go.Figure()
for idx, model in enumerate(df_metrics):
    errors = df_each_err[model].abs().nlargest(n_structs)
    model_mae = errors.mean().round(3)
    fig.add_scatter(
        x=df_wbm.loc[errors.index][fp_diff_col].values,
        y=errors.values,
        mode="markers",
        name=f"{model} · MAE={model_mae:.2f}",
        visible="legendonly" if idx else True,
        hovertemplate=(
            "material ID: %{customdata[0]}<br>"
            "formula: %{customdata[1]}<br>"
            "FP norm diff: %{x}<br>"
            "error: %{y} eV/atom"
        ),
        customdata=df_wbm.loc[errors.index][[Key.mat_id, Key.formula]].values,
        legendrank=model_mae,
    )

title = (
    f"Norm-diff between initial/final SiteStatsFingerprint<br>"
    f"of the {n_structs} highest-error structures for each model"
)
fig.layout.title.update(text=title, xanchor="center", x=0.5)
fig.layout.legend.update(
    title="", yanchor="top", y=0.98, xanchor="right", x=0.98, font_size=14
)
fig.layout.xaxis.title = "|SSFP<sub>initial</sub> - SSFP<sub>final</sub>|"
fig.layout.yaxis.title = "Absolute error (eV/atom)"

fig.show()

# pmv.save_fig(fig, f"{FIGS}/scatter-largest-each-errors-fp-diff-models.svelte")


# %%
df_mp = pd.read_csv(DataFiles.mp_energies.path, na_filter=False).set_index(Key.mat_id)
train_count_col = "MP Occurrences"
df_elem_counts = pmv.count_elements(
    df_mp[Key.formula], count_mode="occurrence"
).to_frame(name=train_count_col)
n_examp_for_rarest_elem_col = "Examples for rarest element in structure"
df_wbm[n_examp_for_rarest_elem_col] = [
    df_elem_counts[train_count_col].loc[list(map(str, Composition(formula)))].min()
    for formula in tqdm(df_wbm[Key.formula])
]
df_preds[n_examp_for_rarest_elem_col] = df_wbm[n_examp_for_rarest_elem_col]


# %% find materials that were misclassified by all models
for model in df_each_pred:
    true_pos, false_neg, false_pos, true_neg = classify_stable(
        df_each_pred[model], df_preds[MbdKey.each_true]
    )
    df_preds[f"{model} true pos"] = true_pos
    df_preds[f"{model} false neg"] = false_neg
    df_preds[f"{model} false pos"] = false_pos
    df_preds[f"{model} true neg"] = true_neg


df_preds["All models true pos"] = df_preds.filter(like=" true pos").all(axis=1)
df_preds["All models false neg"] = df_preds.filter(like=" false neg").all(axis=1)
df_preds["All models false pos"] = df_preds.filter(like=" false pos").all(axis=1)
df_preds["All models true neg"] = df_preds.filter(like=" true neg").all(axis=1)

df_preds.filter(like="All models ").sum()


# %%
normalized = True
elem_counts: dict[str, pd.Series] = {}
for col in ("All models false neg", "All models false pos"):
    elem_counts[col] = elem_counts.get(
        col, pmv.count_elements(df_preds[df_preds[col]][Key.formula])
    )
    fig = pmv.ptable_heatmap_plotly(
        elem_counts[col] / df_elem_counts[train_count_col]
        if normalized
        else elem_counts[col],
        color_bar=dict(title=col),
        fmt=".3f",
        cscale_range=[0, 0.1],
    )
    fig.show()

# TODO plot these for each model individually


# %% map abs EACH model errors onto elements in structure weighted by composition
# fraction and average over all test set structures
df_comp = pd.json_normalize(
    [Composition(comp).as_dict() for comp in tqdm(df_wbm[Key.formula])]
).set_index(df_wbm.index)

# bar plot showing number of structures in MP containing each element
(len(df_comp) - df_comp.isna().sum()).sort_values().plot.bar(backend=PLOTLY)

# df_comp = df_comp.dropna(axis=1, thresh=100)  # remove Xe with only 1 entry


# %% TODO investigate if structures with largest mean error across all models error can
# be attributed to DFT gone wrong. would be cool if models can be run across large
# databases as correctness checkers
df_each_err.abs().mean().sort_values()
df_each_err.abs().mean(axis=1).nlargest(25)


# %%
n_points = 1000
df_largest_fp_diff = df_wbm[fp_diff_col].nlargest(n_points)

fig = go.Figure()
colors = px.colors.qualitative.Plotly

for idx, model in enumerate(df_metrics):
    color = colors[idx]
    model_mae = df_each_err[model].loc[df_largest_fp_diff.index].abs().mean()

    visible = "legendonly" if idx else True
    fig.add_scatter(
        x=df_largest_fp_diff.values,
        y=df_each_err[model].loc[df_largest_fp_diff.index].abs(),
        mode="markers",
        name=f"{model} · MAE={model_mae:.2f}",
        visible=visible,
        hovertemplate=(
            "ID: %{customdata[0]}<br>"
            "formula: %{customdata[1]}<br>"
            "FP diff: %{x}<br>"
            "error: %{y}<extra></extra>"
        ),
        customdata=df_preds[[Key.mat_id, Key.formula]]
        .loc[df_largest_fp_diff.index]
        .values,
        legendgroup=model,
        marker=dict(color=color),
        legendrank=model_mae,
    )
    # add dashed mean line for each model that toggles with the scatter plot
    # fig.add_hline(
    #     y=model_mae,
    #     line=dict(dash="dash"),
    #     annotation=dict(text=f"{model} mean", x=0, xanchor="left", font_size=10),
    # )
    # get color from scatter plot
    fig.add_scatter(
        x=[df_largest_fp_diff.min(), df_largest_fp_diff.max()],
        y=[model_mae, model_mae],
        line=dict(dash="dash", width=2, color=color),
        legendgroup=model,
        showlegend=False,
        visible=visible,
    )


title = (
    f"Absolute errors in model-predicted E<sub>above hull</sub> for {n_points} "
    "structures<br>with largest norm-diff of initial/final SiteStatsFingerprint in WBM "
    "test set"
)
fig.layout.title.update(text=title, xanchor="center", x=0.5)
fig.layout.legend.update(
    title="",
    yanchor="top",
    y=0.98,
    xanchor="right",
    x=0.98,
    font_size=14,
    tracegroupgap=0,
)
fig.layout.xaxis.title = "|SSFP<sub>initial</sub> - SSFP<sub>final</sub>|"
fig.layout.yaxis.title = "|E<sub>above hull</sub> error| (eV/atom)"
fig.layout.margin = dict(t=40, b=0, l=0, r=0)
fig.show()


# %%
pmv.save_fig(fig, f"{SITE_FIGS}/largest-fp-diff-each-error-models.svelte")
pmv.save_fig(fig, f"{PDF_FIGS}/large-fp-diff-vs-each-error.pdf")


# %%
tsne_cols = ["t-SNE 1", "t-SNE 2"]
tsne_comp_2d_path = f"{WBM_DIR}/tsne/onehot-112-composition-2d.csv"
df_wbm[tsne_cols] = pd.read_csv(tsne_comp_2d_path, index_col=0)

df_wbm["wbm_step"] = df_wbm.index.str.split("-").str[1]


# %% t-SNE 2D plot of composition with discrete color based on band gap > 1 eV
fig = px.scatter(
    df_wbm,
    x=tsne_cols[0],
    y=tsne_cols[1],
    color=(df_wbm[Key.bandgap_pbe] > 1).map({True: "band gap", False: "no gap"}),
    hover_name=Key.mat_id,
    hover_data=(Key.formula,),
)
fig.show()

pmv.save_fig(fig, f"{PDF_FIGS}/tsne-2d-composition-by-wbm-step-bandgap.png", scale=3)


# %% violin plot of EACH error for largest norm-diff FP structures for each model
y_label = "E<sub>above hull</sub> error (eV/atom)"
df_melt = (
    df_each_err.loc[df_largest_fp_diff.index]
    .abs()
    .melt(var_name="Model", value_name=y_label)
)
fig = px.violin(df_melt, x="Model", y=y_label, color="Model")
fig.layout.update(showlegend=False)
title = (
    f"Absolute errors in predicted E<sub>above hull</sub> for {len(df_largest_fp_diff)}"
    " structures with largest FP norm-diff before/after relaxation"
)
fig.layout.title.update(text=title, x=0.5, xanchor="center", y=0.95)
fig.show()


# %% violin plot of norm-diff FP in structures with largest EACH error for each model
y_label = "E<sub>above hull</sub> error (eV/atom)"
n_structs = 1000

for label, which in (("min", "nlargest"), ("max", "nsmallest")):
    fig = go.Figure()
    for model in df_metrics:
        errors = getattr(df_each_err[model].abs(), which)(n_structs)
        fig.add_violin(
            x=df_wbm.loc[errors.index][fp_diff_col].values,
            name=f"{model} err<sub>{label}</sub>",
            legendgroup=model,
            hovertemplate=("SSFP diff: %{x:.2f}<br>Count: %{y}"),
            spanmode="hard",
        )
    fig.layout.update(showlegend=False)
    fig.layout.xaxis.title = "SSFP norm-diff before/after relaxation"
    fig.show()


# %% violin plot of EACH error for largest norm-diff FP structures for each model
y_label = "E<sub>above hull</sub> error (eV/atom)"
df_melt = (
    df_each_err.loc[df_largest_fp_diff.index]
    .abs()
    .melt(var_name="Model", value_name=y_label)
)
fig = px.violin(df_melt, x="Model", y=y_label, color="Model")
fig.layout.update(showlegend=False)
title = (
    f"Absolute errors in predicted E<sub>above hull</sub> for {len(df_largest_fp_diff)}"
    " structures with largest FP norm-diff before/after relaxation"
)
fig.layout.title.update(text=title, x=0.5, xanchor="center", y=0.95)
fig.show()
