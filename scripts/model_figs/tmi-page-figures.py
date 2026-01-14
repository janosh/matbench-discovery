"""Generate figures for the TMI (Too Much Information) pages.

Analyzes structures and compositions with largest mean error across all models.
Maybe there's some chemistry/region of materials space that all models struggle with?
Might point to deficiencies in the data or models architecture.
"""

# %%
import math

import pandas as pd
import plotly.express as px
import plotly.graph_objs as go
import pymatviz as pmv
from pymatgen.core import Composition, Element
from pymatviz.enums import Key
from pymatviz.process_data import bin_df_cols
from pymatviz.utils import df_ptable
from tqdm.auto import tqdm

from matbench_discovery import ROOT, SITE_FIGS
from matbench_discovery.cli import cli_args
from matbench_discovery.data import df_wbm, load_df_wbm_with_preds
from matbench_discovery.enums import MbdKey, TestSubset

__author__ = "Janosh Riebesell"
__date__ = "2023-02-15"

test_set_std_col = "Test set standard deviation"
elem_col, size_col = "Element", "group"
fp_diff_col = "site_stats_fingerprint_init_final_norm_diff"


# %%
test_subset = globals().get("test_subset", TestSubset.uniq_protos)
models_to_plot = cli_args.models

# Load predictions for specified models
df_preds = load_df_wbm_with_preds(
    models=models_to_plot, subset=cli_args.test_subset
).round(3)

# Calculate errors for each model
df_each_err = pd.DataFrame(index=df_preds.index)
for model in models_to_plot:
    df_each_err[model.label] = df_preds[model.label] - df_preds[MbdKey.e_form_dft]

# project average model error onto periodic table
df_comp = pd.DataFrame(
    Composition(comp).as_dict() for comp in df_preds[Key.formula]
).set_index(df_preds.index)


# %% compute number of samples per element in training set
# we count raw element occurrence, i.e. not weighted by composition, based on the
# hypothesis that models don't learn more about iron and oxygen from Fe2O3 than from FeO
counts_path = f"{ROOT}/site/src/routes/data/mp-element-counts-by-occurrence.json"
df_elem_err = pd.read_json(counts_path, typ="series")
train_count_col = "MP Occurrences"
df_elem_err = df_elem_err.reset_index(name=train_count_col).set_index("index")
df_elem_err.index.name = "symbol"

# compute std dev of DFT hull dist for each element in test set
df_elem_err[test_set_std_col] = (
    df_comp.where(pd.isna, 1) * df_preds[MbdKey.each_true].to_numpy()[:, None]
).std()


# %% plot number of structures containing each element in MP and WBM
for label, srs in (
    ("MP", df_elem_err[train_count_col]),
    ("WBM", df_comp.where(pd.isna, 1).sum()),
):
    title = f"Number of {label} structures containing each element"
    srs = srs.sort_values().copy()
    srs.index = [f"{len(srs) - idx} {el}" for idx, el in enumerate(srs.index)]
    fig = srs.plot.bar(backend="plotly", title=title)
    fig.layout.update(showlegend=False)
    fig.show()


# %% plot structure counts for each element in MP and WBM in a grouped bar chart
df_struct_counts_base = pd.DataFrame(index=df_elem_err.index)
df_struct_counts_base["MP"] = df_elem_err[train_count_col]
df_struct_counts_base["WBM"] = df_comp.where(pd.isna, 1).sum()
min_count = 10  # only show elements with at least 10 structures
df_struct_counts_base = df_struct_counts_base[
    df_struct_counts_base.sum(axis=1) > min_count
]

for normalized in (False, True):
    df_struct_counts = df_struct_counts_base.copy()
    if normalized:
        df_struct_counts["MP"] /= len(df_preds) / 100
        df_struct_counts["WBM"] /= len(df_preds) / 100
    y_col = "percent" if normalized else "count"

    df_melt = df_struct_counts.reset_index().melt(
        var_name=(clr_col := "Dataset"),
        value_name=y_col,
        id_vars=(symbol_col := "symbol"),
    )
    fig = df_melt.sort_values([y_col, symbol_col]).plot.bar(
        x=symbol_col,
        y=y_col,
        backend="plotly",
        title="Number of structures containing each element",
        color=clr_col,
        barmode="group",
    )

    fig.layout.update(bargap=0.1)
    fig.layout.legend.update(x=0.02, y=0.98, font_size=16)
    fig.show()
    pmv.save_fig(fig, f"{SITE_FIGS}/tmi/bar-element-counts-mp+wbm-{normalized=}.svelte")


# %%
# plot per-element std dev of DFT hull dist
fig = pmv.ptable_heatmap_plotly(
    df_elem_err[test_set_std_col], fmt=".2f", colorscale="Inferno"
)
fig.show()


# %%
# scatter plot error by element against prevalence in training set
# for checking correlation and R2 of elemental prevalence in MP training data vs.
# model error

# compute mean absolute error per element for each model
# weight by element presence (1 if element in structure, 0 otherwise)
df_elem_present = df_comp.where(pd.isna, 1).fillna(0)
for model in models_to_plot:
    model_err = df_each_err[model.label].abs()
    # for each element, compute mean error across structures containing it
    weighted_err = df_elem_present.multiply(model_err, axis=0)
    df_elem_err[model.label] = weighted_err.sum() / df_elem_present.sum()

df_elem_err[elem_col] = [Element(el).long_name for el in df_elem_err.index]

df_melt = df_elem_err.melt(
    id_vars=[train_count_col, test_set_std_col, elem_col],
    value_name=(val_col := "Error"),
    var_name=(clr_col := "Model"),
    ignore_index=False,
)

df_melt[size_col] = df_ptable[size_col].fillna(0)
fig = df_melt.plot.scatter(
    x=train_count_col,
    y=val_col,
    color=clr_col,
    backend="plotly",
    hover_name=elem_col,
    title="Per-element error vs element occurrence in MP training set",
    hover_data={val_col: ":.2f", train_count_col: ":,.0f"},
)
for trace in fig.data:
    trace.visible = "legendonly"
fig.update_traces(textposition="top center")  # place text above scatter points
fig.layout.title.update(xanchor="center", x=0.5)
fig.layout.legend.update(x=1, y=1, xanchor="right", yanchor="top", title="")
fig.show()

pmv.save_fig(fig, f"{SITE_FIGS}/tmi/element-prevalence-vs-error.svelte")


# %%
# plot EACH errors against least prevalent element in structure (by occurrence in
# MP training set). this seems to correlate more with model error
n_of_rarest_elem_col = "Examples for rarest element in structure"
df_preds[n_of_rarest_elem_col] = [
    df_elem_err[train_count_col].loc[list(map(str, Composition(formula)))].min()
    for formula in tqdm(df_preds[Key.formula])
]

df_melt = (
    df_each_err.abs()
    .reset_index()
    .melt(var_name="Model", value_name=Key.each_pred, id_vars=Key.mat_id)
    .set_index(Key.mat_id)
)
df_melt[n_of_rarest_elem_col] = df_preds[n_of_rarest_elem_col]

df_bin = bin_df_cols(
    df_melt, [n_of_rarest_elem_col, Key.each_pred], group_by_cols=["Model"]
)
df_bin = df_bin.reset_index().set_index(Key.mat_id)
df_bin[Key.formula] = df_preds[Key.formula]


# %%
n_cols = 3
n_rows = math.ceil(len(models_to_plot) / n_cols)
facet_col = "Model"

fig = px.scatter(
    df_bin.reset_index(),
    x=n_of_rarest_elem_col,
    y=Key.each_pred,
    color=facet_col,
    facet_col=facet_col,
    facet_col_wrap=n_cols,
    facet_col_spacing=0.02,
    facet_row_spacing=0.02,
    category_orders={facet_col: [m.label for m in models_to_plot]},
    hover_data={Key.mat_id: True, Key.formula: True, facet_col: False},
    width=300 * n_cols,
    height=180 * n_rows,
)

fig.update_traces(marker=dict(size=3))
fig.layout.paper_bgcolor = "rgba(0,0,0,0)"
fig.layout.showlegend = False

# remove "Model=" prefix from subplot titles
for anno in fig.layout.annotations:
    anno.text = anno.text.split("=")[1]

# shared axis titles
axis_titles = dict(xref="paper", yref="paper", showarrow=False, font_size=14)
fig.add_annotation(  # x-axis title
    text="MP occurrence count of least prevalent element in structure",
    x=0.5,
    y=0,
    yshift=-50,
    **axis_titles,
)
fig.add_annotation(  # y-axis title
    text="Absolute error in E<sub>above hull</sub>",
    x=0,
    xshift=-60,
    y=0.5,
    textangle=-90,
    **axis_titles,
)

# standardize margins
portrait = n_rows > n_cols
fig.layout.margin.update(l=60, r=10, t=0 if portrait else 10, b=60 if portrait else 10)

axes_kwargs = dict(matches=None, title_text="", showgrid=True)
fig.update_xaxes(**axes_kwargs)
fig.update_yaxes(**axes_kwargs)

fig.show()
pmv.save_fig(
    fig, f"{SITE_FIGS}/tmi/each-error-vs-least-prevalent-element-in-struct.svelte"
)


# %% --- Fingerprint-based analysis ---
# Analyze correlation between relaxation change (measured by SiteStatsFingerprint diff)
# and model errors

model_labels = [m.label for m in models_to_plot]


# %% histogram of FP diff for structures with largest/smallest errors
n_structs = 1000
fig = go.Figure()
for idx, model in enumerate(model_labels):
    large_errors = df_each_err[model].abs().nlargest(n_structs)
    small_errors = df_each_err[model].abs().nsmallest(n_structs)
    for label, errors in (("min", small_errors), ("max", large_errors)):
        fig.add_histogram(
            x=df_wbm.loc[errors.index][fp_diff_col].values,
            name=f"{model} err<sub>{label}</sub>",
            visible="legendonly" if idx else True,
            legendgroup=model,
            hovertemplate="SSFP diff: %{x:.2f}<br>Count: %{y}",
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

pmv.save_fig(fig, f"{SITE_FIGS}/tmi/hist-largest-each-errors-fp-diff-models.svelte")


# %% scatter plot:
# FP diff vs error for highest-error structures
n_structs = 100
fig = go.Figure()
for idx, model in enumerate(model_labels):
    errors = df_each_err[model].abs().nlargest(n_structs)
    model_mae = errors.mean()
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

pmv.save_fig(fig, f"{SITE_FIGS}/tmi/scatter-largest-each-errors-fp-diff-models.svelte")


# %% scatter plot: errors for structures with largest FP diff (most relaxation change)
n_points = 1000
# filter to only materials in the predictions subset
df_wbm_subset = df_wbm.loc[df_wbm.index.intersection(df_each_err.index)]
df_largest_fp_diff = df_wbm_subset[fp_diff_col].nlargest(n_points)

fig = go.Figure()
colors = px.colors.qualitative.Plotly

for idx, model in enumerate(model_labels):
    color = colors[idx % len(colors)]
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
    # add dashed mean line for each model
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
pmv.save_fig(fig, f"{SITE_FIGS}/tmi/scatter-largest-fp-diff-each-error-models.svelte")
