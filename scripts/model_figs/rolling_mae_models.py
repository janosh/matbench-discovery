"""Plot rolling MAE as a function of hull distance for all models."""

# %%
import numpy as np
import pymatviz as pmv
from pymatviz.enums import Key

from matbench_discovery import PDF_FIGS, SITE_FIGS
from matbench_discovery.enums import MbdKey, TestSubset
from matbench_discovery.models import MODEL_METADATA, model_is_compliant
from matbench_discovery.plots import rolling_mae_vs_hull_dist
from matbench_discovery.preds.discovery import df_each_pred, df_preds, models

__author__ = "Rhys Goodall, Janosh Riebesell"
__date__ = "2022-06-18"

df_err, df_std = None, None  # variables to cache rolling MAE and std


# %%
test_subset = globals().get("test_subset", TestSubset.uniq_protos)

if test_subset == TestSubset.uniq_protos:
    df_preds = df_preds.query(Key.uniq_proto)
    df_each_pred = df_each_pred.loc[df_preds.index]

show_non_compliant = globals().get("show_non_compliant", False)
models_to_plot = [
    model
    for model in models
    if show_non_compliant or model_is_compliant(MODEL_METADATA[model])
]

fig, df_err, df_std = rolling_mae_vs_hull_dist(
    e_above_hull_true=df_preds[MbdKey.each_true],
    e_above_hull_preds=df_each_pred[models_to_plot],
    with_sem=False,
    df_rolling_err=df_err,
    df_err_std=df_std,
    show_dummy_mae=False,
    width=500,
    height=500,
    legend_loc="below",
)

show_n_best_models = None
for trace in fig.data:
    model = trace.name.split(" MAE=")[0]
    if show_n_best_models and model in models_to_plot[show_n_best_models:]:
        trace.visible = "legendonly"  # show only top models by default

rng = np.random.default_rng()
small_noise = rng.random(len(df_preds)) * 1e-12
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
        domain=[0.8, 1], tickformat="s", tickvals=[*range(0, 100_000, 2000)]
    ),
)
fig.show()


# %%
img_suffix = "" if show_non_compliant else "-only-compliant"
img_name = f"rolling-mae-vs-hull-dist-models{img_suffix}"
pmv.save_fig(fig, f"{SITE_FIGS}/{img_name}.svelte")
pmv.save_fig(fig, f"{PDF_FIGS}/{img_name}.pdf")
