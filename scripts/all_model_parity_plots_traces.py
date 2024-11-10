"""parity plot of actual vs predicted e_above_hull and e_form_per_atom for all
models. First 2 plots put all models in single figure with selectable traces.
Last plot is split into 2x3 subplots, one for each model.
"""

# %%
from typing import Literal, get_args

import plotly.express as px
import pymatviz as pmv
from pymatviz.enums import Key
from pymatviz.utils import bin_df_cols

from matbench_discovery import SITE_FIGS
from matbench_discovery.enums import MbdKey, TestSubset
from matbench_discovery.preds.discovery import (
    df_metrics,
    df_metrics_uniq_protos,
    df_preds,
)

__author__ = "Janosh Riebesell"
__date__ = "2022-11-28"

legend = dict(x=1, y=0, xanchor="right", yanchor="bottom", title=None)

# toggle between formation energy and energy above convex hull
EnergyType = Literal["e-form", "each"]
use_e_form, use_each = get_args(EnergyType)
which_energy: EnergyType = globals().get("which_energy", use_each)
if which_energy == use_each:
    e_pred_col = Key.each_pred
    e_true_col = MbdKey.each_true
elif which_energy == use_e_form:
    e_true_col = MbdKey.e_form_dft
    e_pred_col = Key.e_form_pred
else:
    raise ValueError(f"Unexpected {which_energy=}")


test_subset = globals().get("test_subset", TestSubset.uniq_protos)

if test_subset == TestSubset.uniq_protos:
    df_preds = df_preds.query(Key.uniq_proto)
    df_metrics = df_metrics_uniq_protos


# %%
facet_col = "Model"
hover_cols = (MbdKey.each_true, Key.formula)
models = list(df_metrics.T.MAE.nsmallest(6).index)  # top 6 models by MAE
models = list(df_metrics)  # all models

df_melt = df_preds.melt(
    id_vars=(df_preds.index.name, MbdKey.e_form_dft, *hover_cols),
    var_name=facet_col,
    value_vars=models,
    value_name=Key.e_form_pred,
)

df_melt[Key.each_pred] = (
    df_melt[MbdKey.each_true] + df_melt[Key.e_form_pred] - df_melt[MbdKey.e_form_dft]
)

df_bin = bin_df_cols(
    df_melt,
    bin_by_cols=[e_true_col, e_pred_col],
    group_by_cols=[facet_col],
    n_bins=300,
    bin_counts_col=(bin_cnt_col := "bin counts"),
)
df_bin = df_bin.reset_index()

# sort legend and facet plots by MAE
legend_order = list(df_metrics.T.MAE.sort_values().index)


# determine each point's classification to color them by
# now unused, can be used to color points by TP/FP/TN/FN
# true_pos, false_neg, false_pos, true_neg = classify_stable(
#     df_bin[e_true_col], df_bin[e_pred_col]
# )
# clf_col = "classified"
# df_bin[clf_col] = np.array(clf_labels)[
#     true_pos * 0 + false_neg * 1 + false_pos * 2 + true_neg * 3
# ]


# %% parity plot of actual vs predicted e_form_per_atom
fig = px.scatter(
    df_bin,
    x=MbdKey.e_form_dft,
    y=Key.e_form_pred,
    color=facet_col,
    hover_data=hover_cols,
    hover_name=df_preds.index.name,
    opacity=0.7,
    category_orders={facet_col: legend_order},
)

for trace in fig.data:
    # initially hide all traces, let users select which models to compare
    trace.visible = "legendonly"
    model = trace.name
    if model not in df_preds:
        print(f"Unexpected {model=}, not in {models=}")
        continue
    MAE, R2 = df_metrics[model][["MAE", "R2"]]
    trace.name = f"{model} · {MAE=:.2f} · R<sup>2</sup>={R2:.2f}"

fig.layout.legend.update(legend)
pmv.powerups.add_identity_line(fig)
fig.show()

img_name = f"{SITE_FIGS}/e-form-parity-models"
# pmv.save_fig(fig, f"{img_path}.svelte")


# %% parity plot of actual vs predicted e_above_hull
fig = px.scatter(
    df_bin,
    x=e_true_col,
    y=e_pred_col,
    color=facet_col,
    hover_data=hover_cols,
    hover_name=df_preds.index.name,
    opacity=0.7,
    category_orders={facet_col: legend_order},
)

for trace in fig.data:
    trace.visible = "legendonly"
    model = trace.name
    if model not in df_preds:
        print(f"Unexpected {model=}, not in {models=}")
        continue
    MAE, R2 = df_metrics[model][["MAE", "R2"]]
    trace.name = f"{model} · {MAE=:.2f} · R<sup>2</sup>={R2:.2f}"

fig.layout.legend.update(legend)
pmv.powerups.add_identity_line(fig)
fig.show()

img_name = f"{SITE_FIGS}/e-above-hull-parity-models"
# pmv.save_fig(fig, f"{img_path}.svelte")
