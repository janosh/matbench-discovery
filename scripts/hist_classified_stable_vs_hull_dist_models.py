# %%
from plotly.subplots import make_subplots

from matbench_discovery import FIGS, today
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
models = (
    "Wren, CGCNN, CGCNN IS2RE, CGCNN RS2RE, Voronoi RF, "
    "Wrenformer, MEGNet, M3GNet, BOWSR MEGNet"
).split(", ")
df_wbm = load_df_wbm_with_preds(models=models).round(3)

target_col = "e_form_per_atom_mp2020_corrected"
e_above_hull_col = "e_above_hull_mp2020_corrected_ppd_mp"


# %%
which_energy: WhichEnergy = "true"
model_name = "Wrenformer"

backend: Backend = "matplotlib"
if backend == "matplotlib":
    fig, axs = plt.subplots(nrows=3, ncols=3, figsize=(18, 12))
else:
    fig = make_subplots(rows=3, cols=3)


for idx, model_name in enumerate(models):
    ax, metrics = hist_classified_stable_vs_hull_dist(
        e_above_hull_true=df_wbm[e_above_hull_col],
        e_above_hull_pred=df_wbm[e_above_hull_col]
        + (df_wbm[model_name] - df_wbm[target_col]),
        which_energy=which_energy,
        ax=axs.flat[idx],
        backend=backend,
    )
    title = f"{model_name} ({len(df_wbm[model_name].dropna()):,})"
    text = f"Enrichment\nFactor = {metrics['enrichment']:.3}"

    if backend == "matplotlib":
        ax.text(0.02, 0.25, text, fontsize=16, transform=ax.transAxes)
        ax.set(title=title)

    else:
        ax.add_annotation(text=text, x=0.5, y=0.5, showarrow=False)
        ax.update_xaxes(title_text=title)

        for trace in ax.data:
            fig.append_trace(trace, row=idx % 3 + 1, col=idx // 3 + 1)

if backend == "matplotlib":
    fig.suptitle(f"{today} {which_energy=}", y=1.07, fontsize=16)
    plt.figlegend(
        *ax.get_legend_handles_labels(),
        ncol=10,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.05),
        frameon=False,
    )

fig.show()


# %%
img_path = f"{FIGS}/{today}-wbm-hull-dist-hist-models"
# if hasattr(fig, "write_image"):
#     fig.write_image(f"{img_path}.pdf")
#     fig.write_html(f"{img_path}.html", include_ploltyjs="cdn")
# else:
#     fig.savefig(f"{img_path}.pdf")
