"""Histogram of the energy difference (either according to DFT ground truth [default] or
model predicted energy) to the convex hull for materials in the WBM data set. The
histogram stacks true/false positives/negatives with different colors.
"""

# %%
import math
from typing import Final

import pymatviz as pmv
from pymatviz.enums import Key

from matbench_discovery import PDF_FIGS, SITE_FIGS
from matbench_discovery.cli import cli_args
from matbench_discovery.data import load_df_wbm_with_preds
from matbench_discovery.enums import MbdKey, TestSubset
from matbench_discovery.metrics.discovery import dfs_metrics
from matbench_discovery.plots import hist_classified_stable_vs_hull_dist

__author__ = "Janosh Riebesell"
__date__ = "2022-12-01"

show_non_compliant = cli_args.show_non_compliant
models_to_plot = [
    model
    for model in cli_args.models
    if model.is_complete and (show_non_compliant or model.is_compliant)
]
test_subset: TestSubset = globals().get("test_subset", TestSubset.uniq_protos)
models_to_plot = sorted(  # sort models by F1
    models_to_plot,
    key=lambda model: -dfs_metrics[test_subset][model.label][Key.f1.symbol],
)
df_preds = load_df_wbm_with_preds(models=models_to_plot, subset=test_subset)


# %%
n_cols = 3
if cli_args.use_full_rows:
    # drop last models that don't fit in last row
    n_rows = len(models_to_plot) // n_cols
    models_to_plot = models_to_plot[: n_rows * n_cols]
else:
    n_rows = math.ceil(len(models_to_plot) / n_cols)

hover_cols = (df_preds.index.name, MbdKey.e_form_dft, MbdKey.each_true, Key.formula)
facet_col = "Model"

# 'true' or 'pred': whether to put DFT or model-predicted hull distances on the x-axis
which_energy: Final = "pred"
hist_clf_kwargs = dict(
    facet_col=facet_col,
    facet_col_wrap=n_cols,
    category_orders={facet_col: [m.label for m in models_to_plot]},
    facet_col_spacing=0.04,
    facet_row_spacing=0.04,
)

df_melt = df_preds.melt(
    id_vars=hover_cols,
    value_vars=[model.label for model in models_to_plot],
    var_name=facet_col,
    value_name=Key.e_form_pred,
)

df_melt[Key.each_pred] = (
    df_melt[MbdKey.each_true] + df_melt[Key.e_form_pred] - df_melt[MbdKey.e_form_dft]
)

fig = hist_classified_stable_vs_hull_dist(
    df=df_melt,
    each_true_col=MbdKey.each_true,
    each_pred_col=Key.each_pred,
    which_energy=which_energy,
    rolling_acc=None,
    **hist_clf_kwargs,
)

metrics_in_plot_titles = True
for anno in fig.layout.annotations:
    model_name = anno.text = anno.text.split("=", 1).pop()
    if (
        model_name not in [m.label for m in models_to_plot]
        or not metrics_in_plot_titles
    ):
        continue
    F1, FPR, FNR, DAF = (
        dfs_metrics[test_subset][model_name][x] for x in ("F1", "FPR", "FNR", "DAF")
    )
    # anno.text = f"{model_name} · {F1=:.2f} · {FPR=:.2f} · {FNR=:.2f} · {DAF=:.2f}"
    anno.text = f"{model_name} · {F1=:.2f}"

# set the figure size based on the number of rows and columns
fig.layout.height = 230 * n_rows
fig.layout.width = 280 * n_cols
fig.layout.paper_bgcolor = "rgba(0,0,0,0)"

# set the shared y and x axis ranges
fig.update_yaxes(range=[0, 9_000], title_text=None, matches=None, tickformat="s")
fig.update_xaxes(range=[-0.4, 0.4], title_text=None, matches=None)

axis_titles = dict(xref="paper", yref="paper", showarrow=False, font_size=16)
fig.add_annotation(  # x-axis title
    x=0.5,
    y=0,
    yshift=-50,
    text=MbdKey.each_true.label,
    borderpad=5,
    **axis_titles,
)
fig.add_annotation(  # y-axis title
    x=0,
    xshift=-70,
    y=0.5,
    text=Key.each_pred.label,
    textangle=-90,
    borderpad=5,
    **axis_titles,
)

# place the legend above the subplots
fig.layout.legend.update(
    y=1.03, xanchor="center", x=0.5, bgcolor="rgba(0,0,0,0)", orientation="h"
)
# standardize the margins and template
portrait = n_rows > n_cols
fig.layout.margin.update(l=60, r=10, t=0 if portrait else 10, b=60 if portrait else 10)

fig.show()


# %%
img_suffix = "" if show_non_compliant else "-only-compliant"
img_name = f"hist-clf-{which_energy}-hull-dist-models-{n_rows}x{n_cols}{img_suffix}"
pmv.save_fig(fig, f"{SITE_FIGS}/{img_name}.svelte")
pmv.save_fig(fig, f"{PDF_FIGS}/{img_name}.pdf")
