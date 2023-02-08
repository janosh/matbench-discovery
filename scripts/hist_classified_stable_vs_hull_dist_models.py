# %%
from typing import Final

from pymatviz.utils import save_fig

from matbench_discovery import ROOT, STATIC, today
from matbench_discovery.plots import hist_classified_stable_vs_hull_dist, plt
from matbench_discovery.preds import (
    df_metrics,
    df_wbm,
    e_form_col,
    each_true_col,
    models,
)

__author__ = "Janosh Riebesell"
__date__ = "2022-12-01"

"""
Histogram of the energy difference (either according to DFT ground truth [default] or
model predicted energy) to the convex hull for materials in the WBM data set. The
histogram stacks true/false positives/negatives with different colors.
"""


# %%
hover_cols = (df_wbm.index.name, e_form_col, each_true_col, "formula")
e_form_preds = "e_form_per_atom_pred"
each_pred_col = "e_above_hull_pred"
facet_col = "Model"

df_melt = df_wbm.melt(
    id_vars=hover_cols,
    # value_vars=models,
    value_vars=list(df_metrics),
    var_name=facet_col,
    value_name=e_form_preds,
)

df_melt[each_pred_col] = (
    df_melt[each_true_col] + df_melt[e_form_preds] - df_melt[e_form_col]
)


# %%
backend: Final = "plotly"
rows, cols = len(models) // 2, 2
which_energy: Final = "true"
kwds = (
    dict(facet_col=facet_col, facet_col_wrap=cols)
    if backend == "plotly"
    else dict(by=facet_col, figsize=(20, 20), layout=(rows, cols), bins=500)
)

fig = hist_classified_stable_vs_hull_dist(
    df=df_melt,
    each_true_col=each_true_col,
    each_pred_col=each_pred_col,
    which_energy=which_energy,
    backend=backend,
    rolling_acc=None,
    stability_threshold=None,
    **kwds,  # type: ignore[arg-type]
)


# TODO add line showing the true histogram of the hull distance distribution on each subplot
show_metrics = False
if backend == "matplotlib":
    fig = plt.gcf()
    fig.suptitle(f"{today} {which_energy=}", y=1.04, fontsize=18, fontweight="bold")
    plt.figlegend(
        *plt.gca().get_legend_handles_labels(),
        ncol=4,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.03),
        frameon=False,
    )
    # add metrics to subplot titles
    for ax in fig.axes:
        model_name = ax.get_title()
        assert model_name in models
        if not show_metrics:
            continue
        F1, FPR, FNR, DAF = (
            df_metrics[model_name][x] for x in "F1 FPR FNR DAF".split()
        )
        ax.set(title=f"{model_name} · {F1=:.2f} · {FPR=:.2f} · {FNR=:.2f} · {DAF=:.2f}")
else:
    for anno in fig.layout.annotations:
        model_name = anno.text = anno.text.split("=").pop()
        if model_name not in models or not show_metrics:
            continue
        F1, FPR, FNR, DAF = (
            df_metrics[model_name][x] for x in "F1 FPR FNR DAF".split()
        )
        anno.text = f"{model_name} · {F1=:.2f} · {FPR=:.2f} · {FNR=:.2f} · {DAF=:.2f}"

    # horizontal legend at the top
    legend = dict(yanchor="top", y=1, xanchor="right", x=1)
    fig.update_layout(legend=legend, margin=dict(t=50, b=30, l=40, r=0))
    fig.update_yaxes(range=[0, 3_000], title_text=None)

    # for trace in fig.data:
    #     # no need to store all 250k x values in plot, leads to 1.7 MB file,
    #     # subsample every 10th point is enough to see the distribution
    #     trace.x = trace.x[::10]

    # increase height of figure
    fig.show()


# %%
img_name = f"hist-{which_energy}-energy-vs-hull-dist-models"
# save_fig(fig, f"{FIGS}/{img_name}.svelte")
n_models = len(fig.layout.annotations)
save_fig(fig, f"{STATIC}/{img_name}.webp", scale=3, height=100 * n_models)
save_fig(fig, f"{ROOT}/tmp/figures/{img_name}.pdf", height=600, width=600)
