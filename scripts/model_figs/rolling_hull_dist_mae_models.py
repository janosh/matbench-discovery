"""Plot rolling MAE as a function of hull distance for all models."""

# %%
import numpy as np
import pymatviz as pmv

from matbench_discovery import PDF_FIGS, SITE_FIGS
from matbench_discovery.cli import cli_args
from matbench_discovery.enums import MbdKey, TestSubset
from matbench_discovery.metrics.discovery import dfs_metrics
from matbench_discovery.plots import rolling_mae_vs_hull_dist
from matbench_discovery.preds.discovery import df_each_pred, df_preds

__author__ = "Rhys Goodall, Janosh Riebesell"
__date__ = "2022-06-18"

df_err, df_std = None, None  # variables to cache rolling MAE and std


# %%
test_subset = globals().get("test_subset", TestSubset.uniq_protos)

if test_subset == TestSubset.uniq_protos:
    df_preds = df_preds.query(MbdKey.uniq_proto)
    df_each_pred = df_each_pred.loc[df_preds.index]

show_non_compliant = globals().get("show_non_compliant", False)
models_to_plot = [
    model.label
    for model in cli_args.models
    if model.is_complete and (show_non_compliant or model.is_compliant)
]
mae_vals = dfs_metrics[test_subset].loc["MAE", models_to_plot]
model_ranking = mae_vals.sort_values().index[::-1]

fig, df_err, df_std = rolling_mae_vs_hull_dist(
    e_above_hull_true=df_preds[MbdKey.each_true],
    e_above_hull_preds=df_each_pred[models_to_plot][model_ranking],
    with_sem=False,
    df_rolling_err=df_err,
    df_err_std=df_std,
    show_dummy_mae=False,
    legend_loc="default",
    y_lim=(0, 0.1),
)


# Show only the top N models by default
show_n_best_models = 6
for trace in fig.data:
    trace.visible = (
        True if trace.name in model_ranking[-show_n_best_models:] else "legendonly"
    )

# add negligible noise to prevent strange binning artifacts in the marginal plot
small_noise = np.random.default_rng(seed=0).random(len(df_preds)) * 1e-12
counts, bins = np.histogram(
    df_preds[MbdKey.each_true] + small_noise,
    bins=200,  # match the histogram clf plots.
    range=(-0.7, 0.7),
)
fig.add_scatter(
    x=bins, y=counts, name="Density", fill="tozeroy", showlegend=False, yaxis="y2"
)
fig.data[-1].marker.color = "rgba(0, 150, 200, 1)"

# update layout to include marginal plot
fig.layout.update(
    yaxis1=dict(domain=[0, 0.75]),  # main yaxis
    yaxis2=dict(  # marginal yaxis
        domain=[0.8, 1],
        tickformat="s",
        tickvals=[*range(0, 100_000, 2_000)],
        range=[0, 8_000],
    ),
)
fig.show()


# %%
img_suffix = "" if show_non_compliant else "-only-compliant"
img_name = f"rolling-mae-vs-hull-dist-models{img_suffix}"
pmv.save_fig(fig, f"{SITE_FIGS}/{img_name}.svelte")
pmv.save_fig(fig, f"{PDF_FIGS}/{img_name}.pdf")
