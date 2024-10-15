"""Analyze structures and composition with largest mean error across all models.
Maybe there's some chemistry/region of materials space that all models struggle with?
Might point to deficiencies in the data or models architecture.
"""

# %%
import pandas as pd
import plotly.express as px
import pymatviz as pmv
from pymatgen.core import Composition, Element
from pymatviz.enums import Key
from pymatviz.utils import PLOTLY, bin_df_cols, df_ptable
from tqdm import tqdm

from matbench_discovery import PDF_FIGS, ROOT, SITE_FIGS
from matbench_discovery.data import df_wbm
from matbench_discovery.enums import MbdKey, TestSubset
from matbench_discovery.preds import (
    Model,
    df_each_err,
    df_metrics,
    df_metrics_uniq_protos,
    df_preds,
)

__author__ = "Janosh Riebesell"
__date__ = "2023-02-15"


# %%
test_subset = globals().get("test_subset", TestSubset.uniq_protos)

if test_subset == TestSubset.uniq_protos:
    df_preds = df_preds.query(Key.uniq_proto)
    df_each_err = df_each_err.loc[df_preds.index]
    df_metrics = df_metrics_uniq_protos
    df_wbm = df_wbm.loc[df_preds.index]


for df in (df_each_err, df_preds):
    df[MbdKey.each_err_models] = df_each_err.abs().mean(axis=1)


# %% project average model error onto periodic table
df_comp = pd.DataFrame(
    Composition(comp).as_dict() for comp in df_wbm[Key.formula]
).set_index(df_wbm.index)


# %% compute number of samples per element in training set
# we count raw element occurrence, i.e. not weighted by composition, based on the
# hypothesis that models don't learn more about iron and oxygen from Fe2O3 than from FeO
counts_path = f"{ROOT}/site/src/routes/data/mp-element-counts-by-occurrence.json"
df_elem_err = pd.read_json(counts_path, typ="series")
train_count_col = "MP Occurrences"
df_elem_err = df_elem_err.reset_index(name=train_count_col).set_index("index")
df_elem_err.index.name = "symbol"


# %% plot number of structures containing each element in MP and WBM
for label, srs in (
    ("MP", df_elem_err[train_count_col]),
    ("WBM", df_comp.where(pd.isna, 1).sum()),
):
    title = f"Number of {label} structures containing each element"
    srs = srs.sort_values().copy()
    srs.index = [f"{len(srs) - idx} {el}" for idx, el in enumerate(srs.index)]
    fig = srs.plot.bar(backend=PLOTLY, title=title)
    fig.layout.update(showlegend=False)
    fig.show()


# %% plot structure counts for each element in MP and WBM in a grouped bar chart
df_struct_counts = pd.DataFrame(index=df_elem_err.index)
df_struct_counts["MP"] = df_elem_err[train_count_col]
df_struct_counts["WBM"] = df_comp.where(pd.isna, 1).sum()
min_count = 10  # only show elements with at least 10 structures
df_struct_counts = df_struct_counts[df_struct_counts.sum(axis=1) > min_count]
normalized = False
if normalized:
    df_struct_counts["MP"] /= len(df_preds) / 100
    df_struct_counts["WBM"] /= len(df_wbm) / 100
y_col = "percent" if normalized else "count"

df_melt = df_struct_counts.reset_index().melt(
    var_name=(clr_col := "Dataset"), value_name=y_col, id_vars=(symbol_col := "symbol")
)
fig = df_melt.sort_values([y_col, symbol_col]).plot.bar(
    x=symbol_col,
    y=y_col,
    backend=PLOTLY,
    title="Number of structures containing each element",
    color=clr_col,
    barmode="group",
)

fig.layout.update(bargap=0.1)
fig.layout.legend.update(x=0.02, y=0.98, font_size=16)
fig.show()
pmv.save_fig(fig, f"{SITE_FIGS}/bar-element-counts-mp+wbm-{normalized=}.svelte")


# %% compute std dev of DFT hull dist for each element in test set
test_set_std_col = "Test set standard deviation"
df_elem_err[test_set_std_col] = (
    df_comp.where(pd.isna, 1) * df_wbm[MbdKey.each_true].to_numpy()[:, None]
).std()


# %% plot per-element std dev of DFT hull dist
fig = pmv.ptable_heatmap_plotly(
    df_elem_err[test_set_std_col], fmt=".2f", colorscale="Inferno"
)
fig.show()


# %%
normalized = True
cs_range = (0, 0.5)  # same range for all plots
# cs_range = (None, None)  # different range for each plot
for model in (*df_metrics, MbdKey.each_err_models):
    df_elem_err[model] = (df_comp * df_each_err[model].abs().to_numpy()[:, None]).mean()
    # don't change series values in place, would change the df
    per_elem_err = df_elem_err[model].copy(deep=True)
    per_elem_err.name = f"{model} (eV/atom)"
    if normalized:
        per_elem_err /= df_elem_err[test_set_std_col]
        per_elem_err.name = f"{model} (normalized by test set std)"
    fig = pmv.ptable_heatmap_plotly(
        per_elem_err, fmt=".2f", colorscale="Inferno", cscale_range=cs_range
    )
    fig.show()


# %%
expected_cols = {
    *"ALIGNN, BOWSR, CGCNN, CGCNN+P, CHGNet, M3GNet, MEGNet, "
    f"{train_count_col}, {MbdKey.each_err_models}, {test_set_std_col}, Voronoi RF, "
    "Wrenformer".split(", ")
}
if missing_cols := expected_cols - {*df_elem_err}:
    raise ValueError(f"{missing_cols=} not in {df_elem_err.columns=}")
if any(df_elem_err.isna().sum() > 35):
    raise ValueError("Too many NaNs in df_elem_err")
df_elem_err.round(4).to_json(f"{SITE_FIGS}/per-element-each-errors.json")


# %% scatter plot error by element against prevalence in training set
# for checking correlation and R2 of elemental prevalence in MP training data vs.
# model error
elem_col = "Element"
df_elem_err[elem_col] = [Element(el).long_name for el in df_elem_err.index]

df_melt = df_elem_err.melt(
    id_vars=[train_count_col, test_set_std_col, elem_col],
    value_name=(val_col := "Error"),
    var_name=(clr_col := "Model"),
    ignore_index=False,
)
size_col = "group"
df_melt[size_col] = df_ptable[size_col].fillna(0)
fig = df_melt.plot.scatter(
    x=train_count_col,
    y=val_col,
    color=clr_col,
    backend=PLOTLY,
    # size=size_col,
    hover_name=elem_col,
    # text=df_melt.index.where(
    #     (df_melt[val_col] > 0.04) | (df_melt[train_count_col] > 6_000)
    # ),
    title="Per-element error vs element occurrence in MP training set",
    hover_data={val_col: ":.2f", train_count_col: ":,.0f"},
)
for trace in fig.data:
    if trace.name in ("CHGNet", "Voronoi RF", MbdKey.each_err_models):
        continue
    trace.visible = "legendonly"
fig.update_traces(textposition="top center")  # place text above scatter points
fig.layout.title.update(xanchor="center", x=0.5)
fig.layout.legend.update(x=1, y=1, xanchor="right", yanchor="top", title="")
fig.show()

# pmv.save_fig(fig, f"{FIGS}/element-prevalence-vs-error.svelte")
pmv.save_fig(fig, f"{PDF_FIGS}/element-prevalence-vs-error.pdf")


# %% plot EACH errors against least prevalent element in structure (by occurrence in
# MP training set). this seems to correlate more with model error
n_examp_for_rarest_elem_col = "Examples for rarest element in structure"
df_wbm[n_examp_for_rarest_elem_col] = [
    df_elem_err[train_count_col].loc[list(map(str, Composition(formula)))].min()
    for formula in tqdm(df_wbm[Key.formula])
]


# %%
df_melt = (
    df_each_err.abs()
    .reset_index()
    .melt(var_name="Model", value_name=Key.each_pred, id_vars=Key.mat_id)
    .set_index(Key.mat_id)
)
df_melt[n_examp_for_rarest_elem_col] = df_wbm[n_examp_for_rarest_elem_col]

df_bin = bin_df_cols(
    df_melt, [n_examp_for_rarest_elem_col, Key.each_pred], group_by_cols=["Model"]
)
df_bin = df_bin.reset_index().set_index(Key.mat_id)
df_bin[Key.formula] = df_wbm[Key.formula]


# %%
fig = px.scatter(
    df_bin.reset_index(),
    x=n_examp_for_rarest_elem_col,
    y=Key.each_pred,
    color="Model",
    facet_col="Model",
    facet_col_wrap=3,
    hover_data={Key.mat_id: True, Key.formula: True, "Model": False},
    title="Absolute errors in model-predicted E<sub>above hull</sub> vs. occurrence "
    "count in MP training set<br>of least prevalent element in structure",
)
fig.layout.update(showlegend=False)
fig.layout.title.update(x=0.5, xanchor="center", y=0.95)
fig.layout.margin.update(t=100)
# remove axis labels
fig.update_xaxes(title="")
fig.update_yaxes(title="")
for anno in fig.layout.annotations:
    anno.text = anno.text.split("=")[1]

fig.add_annotation(
    text="MP occurrence count of least prevalent element in structure",
    x=0.5,
    y=-0.18,
    xref="paper",
    yref="paper",
    showarrow=False,
)
fig.add_annotation(
    text="Absolute error in E<sub>above hull</sub>",
    x=-0.07,
    y=0.5,
    xref="paper",
    yref="paper",
    showarrow=False,
    textangle=-90,
)

fig.show()
pmv.save_fig(fig, f"{SITE_FIGS}/each-error-vs-least-prevalent-element-in-struct.svelte")


# %% plot histogram of model errors for each element
model = Model.mace.label
fig_ptable_each_errors = pmv.ptable_hists(
    df_comp * (df_each_err[model].to_numpy()[:, None]),
    log=True,
    cbar_title=f"{model} convex hull distance errors (eV/atom)",
    cbar_title_kwargs=dict(fontsize=16),
    x_range=(-0.5, 0.5),
    symbol_pos=(0.1, 0.8),
    colormap="viridis",
)


# %%
img_name = f"ptable-each-error-hists-{model.lower().replace(' ', '-')}"
pmv.save_fig(fig_ptable_each_errors, f"{PDF_FIGS}/{img_name}.pdf")
