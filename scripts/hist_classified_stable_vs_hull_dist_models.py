# %%
from plotly.subplots import make_subplots
from pymatviz.utils import save_fig

from matbench_discovery import STATIC, today
from matbench_discovery.data import load_df_wbm_with_preds
from matbench_discovery.plots import (
    Backend,
    WhichEnergy,
    hist_classified_stable_vs_hull_dist,
    plt,
)

__author__ = "Janosh Riebesell"
__date__ = "2022-12-01"

"""
Histogram of the energy difference (either according to DFT ground truth [default] or
model predicted energy) to the convex hull for materials in the WBM data set. The
histogram stacks true/false positives/negatives with different colors.
"""


# %%
models = sorted(
    "CGCNN, Voronoi RF, Wrenformer, MEGNet, M3GNet, BOWSR MEGNet".split(", ")
)
df_wbm = load_df_wbm_with_preds(models=models).round(3)

e_form_col = "e_form_per_atom_mp2020_corrected"
e_above_hull_col = "e_above_hull_mp2020_corrected_ppd_mp"


# %%
which_energy: WhichEnergy = "true"
model_name = "Wrenformer"

backend: Backend = "matplotlib"
rows, cols = len(models) // 3, 3
if backend == "matplotlib":
    fig, axs = plt.subplots(nrows=rows, ncols=cols, figsize=(7 * cols, 5 * rows))
else:
    x_title = "distance to convex hull (eV/atom)"
    fig = make_subplots(
        rows=rows, cols=cols, y_title="Count", x_title=x_title, subplot_titles=models
    )


for idx, model_name in enumerate(models):
    ax, metrics = hist_classified_stable_vs_hull_dist(
        e_above_hull_true=df_wbm[e_above_hull_col],
        e_above_hull_pred=df_wbm[e_above_hull_col]
        + (df_wbm[model_name] - df_wbm[e_form_col]),
        which_energy=which_energy,
        ax=axs.flat[idx] if backend == "matplotlib" else None,
        backend=backend,
    )
    enrichment, acc, F1 = metrics["enrichment"], metrics["accuracy"], metrics["f1"]
    text = f"{enrichment = :.2f}\n{acc = :.2f}\n{F1 = :.2f}"

    # n_preds = len(df_wbm[model_name].dropna())
    # n_missing = len(df_wbm) - n_preds
    # title = f"{model_name} ({len(df_wbm[model_name].dropna()):,})"
    title = model_name
    if backend == "matplotlib":
        y_anno = 0.7 if model_name == "M3GNet" else 0.25
        ax.text(0.02, y_anno, text, fontsize=16, transform=ax.transAxes)
        ax.set(title=title)

    # no need to store all 250k x values in plot, leads to 1.7 MB file, subsample every 10th
    # point is enough to see the distribution

    else:
        for trace in ax.data:
            trace.x = trace.x[::10]
        fig.add_annotation(text=text, x=0.5, y=0.5, showarrow=False)
        row, col = idx % rows + 1, idx // rows + 1
        fig.add_traces(ax.data, rows=row, cols=col)
        # fig.update_xaxes(title_text=title, row=row, col=col)


if backend == "matplotlib":
    # fig.suptitle(f"{today} {which_energy=}", y=1.07, fontsize=16)
    plt.figlegend(
        *ax.get_legend_handles_labels(),
        ncol=10,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.05),
        frameon=False,
    )
else:
    fig.update_xaxes(range=[-0.4, 0.4])
    fig.update_layout(showlegend=False, barmode="stack")


fig.show()


# %%
img_path = f"{STATIC}/{today}-wbm-hull-dist-hist-models"
# save_fig(fig, f"{img_path}.html")
# save_fig(fig, f"{img_path}.png", scale=3, height=700, width=1000)
save_fig(fig, f"{img_path}.png", dpi=300)
# save_fig(fig, f"{img_path}.pdf")
