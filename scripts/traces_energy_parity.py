"""Parity plot of actual vs predicted e_above_hull and e_form_per_atom for all
models. Unlike energy_parity_tiles.py, this script puts all models in a
single figure with hidable traces.
"""

# %%
import plotly.express as px
import pymatviz as pmv
from pymatviz.enums import Key
from pymatviz.utils import bin_df_cols

from matbench_discovery import SITE_FIGS
from matbench_discovery.cli import cli_args
from matbench_discovery.data import load_df_wbm_with_preds
from matbench_discovery.enums import MbdKey, TestSubset
from matbench_discovery.metrics.discovery import dfs_metrics

__author__ = "Janosh Riebesell"
__date__ = "2022-11-28"

legend = dict(x=1, y=0, xanchor="right", yanchor="bottom", title=None)

# toggle between formation energy and energy above convex hull
which_energy = cli_args.energy_type  # 'each' or 'e_form'

if which_energy == Key.each:
    e_pred_col = Key.each_pred
    e_true_col = MbdKey.each_true
elif which_energy == Key.e_form:
    e_true_col = MbdKey.e_form_dft
    e_pred_col = Key.e_form_pred
else:
    raise ValueError(f"Unexpected {which_energy=}")

test_subset = globals().get("test_subset", TestSubset.uniq_protos)

# Get list of models from command line args (after energy type), fall back to all models
models_to_plot = cli_args.models

# Load predictions for specified models
df_preds = load_df_wbm_with_preds(
    models=models_to_plot, subset=cli_args.test_subset
).round(3)


# %%
facet_col = "Model"
hover_cols = (MbdKey.each_true, Key.formula)

df_melt = df_preds.melt(
    id_vars=(df_preds.index.name, MbdKey.e_form_dft, *hover_cols),
    var_name=facet_col,
    value_vars=[model.label for model in models_to_plot],
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
legend_order = list(dfs_metrics[test_subset].T.MAE.sort_values().index)


# %% parity plot of actual vs predicted e_form_per_atom
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
    # initially hide all traces, let users select which models to compare
    trace.visible = "legendonly"
    model = trace.name
    if model not in df_preds:
        print(f"Unexpected {model=}, not in {[m.label for m in models_to_plot]}")
        continue
    MAE, R2 = dfs_metrics[test_subset][model][["MAE", "R2"]]
    trace.name = f"{model} · {MAE=:.2f} · R<sup>2</sup>={R2:.2f}"

fig.layout.legend.update(legend)
pmv.powerups.add_identity_line(fig)
fig.show()

img_name = f"{which_energy}-parity-models"
# pmv.save_fig(fig, f"{SITE_FIGS}/{img_name}.svelte")


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
        print(f"Unexpected {model=}, not in {[m.label for m in models_to_plot]}")
        continue
    MAE, R2 = dfs_metrics[test_subset][model][["MAE", "R2"]]
    trace.name = f"{model} · {MAE=:.2f} · R<sup>2</sup>={R2:.2f}"

fig.layout.legend.update(legend)
pmv.powerups.add_identity_line(fig)
fig.show()

img_path = f"{SITE_FIGS}/e-above-hull-parity-models"
# pmv.save_fig(fig, f"{img_path}.svelte")
