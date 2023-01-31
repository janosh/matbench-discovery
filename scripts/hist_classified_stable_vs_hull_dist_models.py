# %%
from pymatviz.utils import save_fig

from matbench_discovery import STATIC, today
from matbench_discovery.data import load_df_wbm_preds
from matbench_discovery.energy import stable_metrics
from matbench_discovery.plots import Backend, hist_classified_stable_vs_hull_dist, plt

__author__ = "Janosh Riebesell"
__date__ = "2022-12-01"

"""
Histogram of the energy difference (either according to DFT ground truth [default] or
model predicted energy) to the convex hull for materials in the WBM data set. The
histogram stacks true/false positives/negatives with different colors.
"""


# %%
models = sorted(
    "CGCNN, Voronoi Random Forest, Wrenformer, MEGNet, M3GNet, BOWSR MEGNet".split(", ")
)
df_wbm = load_df_wbm_preds(models).round(3)

e_form_col = "e_form_per_atom_mp2020_corrected"
each_true_col = "e_above_hull_mp2020_corrected_ppd_mp"


# %%
hover_cols = (df_wbm.index.name, e_form_col, each_true_col, "formula")
e_form_preds = "e_form_per_atom_pred"
each_pred_col = "e_above_hull_pred"
facet_col = "Model"

df_melt = df_wbm.melt(
    id_vars=hover_cols,
    value_vars=models,
    var_name=facet_col,
    value_name=e_form_preds,
)

df_melt[each_pred_col] = (
    df_melt[each_true_col] + df_melt[e_form_preds] - df_melt[e_form_col]
)


# %%
backend: Backend = "plotly"
rows, cols = len(models) // 2, 2
kwds = (
    dict(facet_col=facet_col, facet_col_wrap=cols, barmode="stack")
    if backend == "plotly"
    else dict(by=facet_col, figsize=(20, 20), layout=(rows, cols), bins=500)
)

fig = hist_classified_stable_vs_hull_dist(
    df=df_melt,
    each_true_col=each_true_col,
    each_pred_col=each_pred_col,
    which_energy=(which_energy := "true"),
    backend=backend,
    rolling_acc=None,
    **kwds,  # type: ignore[arg-type]
)


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
        df_model = df_melt[df_melt[facet_col] == model_name]
        metrics = stable_metrics(df_model[each_true_col], df_model[each_pred_col])

        DAF, acc, F1 = metrics["DAF"], metrics["Accuracy"], metrics["F1"]
        text = f" {model_name} · {DAF = :.2f} · {acc = :.2f} · {F1 = :.2f}"
        ax.set(title=text)
else:
    for anno in fig.layout.annotations:
        model_name = anno.text.split("=").pop()
        if model_name not in models:
            continue
        df_model = df_melt[df_melt[facet_col] == model_name]
        metrics = stable_metrics(df_model[each_true_col], df_model[each_pred_col])
        F1, FPR, FNR, DAF = (metrics[x] for x in "F1 FPR FNR DAF".split())
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
    fig.layout.height = 800
    fig.show()


# %%
img_path = f"{today}-hist-{which_energy}-energy-vs-hull-dist-models"
# save_fig(fig, f"{FIGS}/{img_path}.svelte")
save_fig(fig, f"{STATIC}/{img_path}.webp", scale=3, height=1000, width=1200)
# save_fig(fig, f"{STATIC}/{img_path}.webp", dpi=300)
