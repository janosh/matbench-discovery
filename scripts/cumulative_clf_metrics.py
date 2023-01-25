# %%
import pandas as pd
from pymatviz.utils import save_fig

from matbench_discovery import FIGS, today
from matbench_discovery.data import load_df_wbm_preds
from matbench_discovery.plots import cumulative_precision_recall

__author__ = "Janosh Riebesell, Rhys Goodall"
__date__ = "2022-12-04"


# %%
models = (
    "CGCNN, Voronoi Random Forest, Wrenformer, MEGNet, M3GNet, BOWSR MEGNet"
).split(", ")

df_wbm = load_df_wbm_preds(models=models).round(3)

# df_wbm.columns = [f"{col}_e_form" if col in models else col for col in df_wbm]
e_form_col = "e_form_per_atom_mp2020_corrected"
e_above_hull_col = "e_above_hull_mp2020_corrected_ppd_mp"


# %%
df_e_above_hull_pred = pd.DataFrame()
for model in models:
    e_above_hul_pred = df_wbm[e_above_hull_col] + df_wbm[model] - df_wbm[e_form_col]
    df_e_above_hull_pred[model] = e_above_hul_pred

fig, df_metric = cumulative_precision_recall(
    e_above_hull_true=df_wbm[e_above_hull_col],
    df_preds=df_e_above_hull_pred,
    project_end_point="xy",
    backend=(backend := "plotly"),
)

title = f"{today} - Cumulative Precision, Recall, F1 scores for classifying stable materials"
# xlabel_cumulative = "Materials predicted stable sorted by hull distance"
if backend == "matplotlib":
    fig.suptitle(title)
    # fig.text(0.5, -0.08, xlabel_cumulative, ha="center", fontdict={"size": 16})
elif backend == "plotly":
    # place legend in lower right corner
    fig.update_layout(
        title=title,
        legend=dict(yanchor="bottom", y=0.02, xanchor="right", x=1),
    )
    fig.update_xaxes(matches=None, showticklabels=True)
    fig.update_yaxes(matches=None, showticklabels=True)

fig.show(responsive=True)


# %%
# file will be served by site
# so we round y floats to reduce file size since
for trace in fig.data:
    assert isinstance(trace.y[0], float)
    trace.y = [round(y, 3) for y in trace.y]

img_path = f"{FIGS}/{today}-cumulative-clf-metrics"
# save_fig(fig, f"{img_path}.pdf")
save_fig(fig, f"{img_path}.svelte")
# save_fig(fig, f"{img_path}.webp", scale=3)
