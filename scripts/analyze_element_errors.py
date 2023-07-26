"""Analyze structures and composition with largest mean error across all models.
Maybe there's some chemistry/region of materials space that all models struggle with?
Might point to deficiencies in the data or models architecture.
"""


# %%
import pandas as pd
import plotly.express as px
from pymatgen.core import Composition, Element
from pymatviz import ptable_heatmap_plotly
from pymatviz.utils import bin_df_cols, df_ptable, save_fig
from tqdm import tqdm

from matbench_discovery import FIGS, PDF_FIGS, ROOT
from matbench_discovery.data import df_wbm
from matbench_discovery.preds import (
    df_each_err,
    df_preds,
    each_pred_col,
    each_true_col,
    model_mean_err_col,
)

__author__ = "Janosh Riebesell"
__date__ = "2023-02-15"

df_each_err[model_mean_err_col] = df_preds[model_mean_err_col] = df_each_err.abs().mean(
    axis=1
)


# %% map average model error onto elements
frac_comp_col = "fractional composition"
df_wbm[frac_comp_col] = [
    Composition(comp).fractional_composition for comp in tqdm(df_wbm.formula)
]

df_frac_comp = pd.DataFrame(comp.as_dict() for comp in df_wbm[frac_comp_col]).set_index(
    df_wbm.index
)
assert all(
    df_frac_comp.sum(axis=1).round(6) == 1
), "composition fractions don't sum to 1"

# df_frac_comp = df_frac_comp.dropna(axis=1, thresh=100)  # remove Xe with only 1 entry


# %% compute number of samples per element in training set
# counting element occurrences not weighted by composition, assuming model don't learn
# much more about iron and oxygen from Fe2O3 than from FeO
df_elem_err = pd.read_json(
    f"{ROOT}/site/src/routes/about-the-data/mp-element-counts-occurrence.json",
    typ="series",
)
train_count_col = "MP Occurrences"
df_elem_err = df_elem_err.reset_index(name=train_count_col).set_index("index")


# %%
fig = ptable_heatmap_plotly(df_elem_err[train_count_col], font_size=10)
title = "Number of MP structures containing each element"
fig.layout.title.update(text=title, x=0.4, y=0.9)
fig.show()


# %%
test_set_std_col = "Test set standard deviation"
df_elem_err[test_set_std_col] = (
    df_frac_comp.where(pd.isna, 1) * df_wbm[each_true_col].to_numpy()[:, None]
).std()


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
    backend="plotly",
    # size=size_col,
    hover_name=elem_col,
    # text=df_melt.index.where(
    #     (df_melt[val_col] > 0.04) | (df_melt[train_count_col] > 6_000)
    # ),
    title="Per-element error vs element-occurrence in MP training set",
    hover_data={val_col: ":.2f", train_count_col: ":,.0f"},
)
for trace in fig.data:
    if trace.name in ("CHGNet", "Voronoi RF", model_mean_err_col):
        continue
    trace.visible = "legendonly"
fig.update_traces(textposition="top center")  # place text above scatter points
fig.layout.title.update(xanchor="center", x=0.5)
fig.layout.legend.update(x=1, y=1, xanchor="right", yanchor="top", title="")
fig.show()

# save_fig(fig, f"{FIGS}/element-prevalence-vs-error.svelte")
save_fig(fig, f"{PDF_FIGS}/element-prevalence-vs-error.pdf")


# %% plot EACH errors against least prevalent element in structure (by occurrence in
# MP training set). this seems to correlate more with model error
n_examp_for_rarest_elem_col = "Examples for rarest element in structure"
df_wbm[n_examp_for_rarest_elem_col] = [
    df_elem_err[train_count_col].loc[list(map(str, Composition(formula)))].min()
    for formula in tqdm(df_wbm.formula)
]


# %%
df_melt = (
    df_each_err.abs()
    .reset_index()
    .melt(var_name="Model", value_name=each_pred_col, id_vars="material_id")
    .set_index("material_id")
)
df_melt[n_examp_for_rarest_elem_col] = df_wbm[n_examp_for_rarest_elem_col]

df_bin = bin_df_cols(df_melt, [n_examp_for_rarest_elem_col, each_pred_col], ["Model"])
df_bin = df_bin.reset_index().set_index("material_id")
df_bin["formula"] = df_wbm.formula


# %%
fig = px.scatter(
    df_bin.reset_index(),
    x=n_examp_for_rarest_elem_col,
    y=each_pred_col,
    color="Model",
    facet_col="Model",
    facet_col_wrap=3,
    hover_data=dict(material_id=True, formula=True, Model=False),
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
save_fig(fig, f"{FIGS}/each-error-vs-least-prevalent-element-in-struct.svelte")
