# %%
from sklearn.metrics import f1_score

from matbench_discovery import ROOT, today
from matbench_discovery.load_preds import load_df_wbm_with_preds
from matbench_discovery.plots import cumulative_precision_recall

__author__ = "Janosh Riebesell, Rhys Goodall"
__date__ = "2022-12-04"


# %%
models = (
    "Wren, CGCNN IS2RE, CGCNN RS2RE, Voronoi RF, "
    "Wrenformer, MEGNet, M3GNet, BOWSR MEGNet"
).split(", ")

df_wbm = load_df_wbm_with_preds(models=models).round(3)

# df_wbm.columns = [f"{col}_e_form" if col in models else col for col in df_wbm]
target_col = "e_form_per_atom_mp2020_corrected"
e_above_hull_col = "e_above_hull_mp2020_corrected_ppd_mp"


# %%
for model in models:
    pred_col = f"{model}_e_form"
    F1 = f1_score(df_wbm[e_above_hull_col] < 0, df_wbm[model] < 0)
    plot_label = f"{model} {F1=:.2}"
    df_wbm[plot_label] = df_wbm[e_above_hull_col] + df_wbm[model] - df_wbm[target_col]

fig, df_metric = cumulative_precision_recall(
    e_above_hull_true=df_wbm[e_above_hull_col],
    df_preds=df_wbm.filter(like="F1="),
    project_end_point="xy",
    backend=(backend := "plotly"),
    show_optimal=True,
)


title = f"{today} - Cumulative Precision and Recall for Stable Materials"
# xlabel_cumulative = "Materials predicted stable sorted by hull distance"
if backend == "matplotlib":
    fig.suptitle(title)
    # fig.text(0.5, -0.08, xlabel_cumulative, ha="center", fontdict={"size": 16})
elif backend == "plotly":
    fig.update_layout(title=title)


fig.show()


# %%
img_path = f"{ROOT}/figures/{today}-precision-recall-curves"
if hasattr(fig, "write_image"):
    # fig.write_image(f"{img_path}.pdf")
    fig.write_html(f"{img_path}.svelte", include_plotlyjs=False, full_html=False)
# else:
#     fig.savefig(f"{img_path}.pdf")
