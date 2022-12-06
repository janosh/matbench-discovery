# %%
from matbench_discovery import ROOT, today
from matbench_discovery.load_preds import load_df_wbm_with_preds
from matbench_discovery.plots import rolling_mae_vs_hull_dist

__author__ = "Rhys Goodall, Janosh Riebesell"
__date__ = "2022-06-18"


# %%
df_wbm = load_df_wbm_with_preds(models=["Wren", "Wrenformer"]).round(3)

e_above_hull_col = "e_above_hull_mp2020_corrected_ppd_mp"
target_col = "e_form_per_atom_mp2020_corrected"


# %%
model_name = "Wrenformer"
ax = rolling_mae_vs_hull_dist(
    e_above_hull_true=df_wbm[e_above_hull_col],
    e_above_hull_error=df_wbm[target_col] - df_wbm[model_name],
    label=model_name,
    backend=(backend := "plotly"),
)

title = f"{today} {model_name}"
if backend == "matplotlib":
    fig = ax.figure
    fig.set_size_inches(10, 9)
    ax.legend(loc="lower right", frameon=False)
    ax.set(title=title)
elif backend == "plotly":
    ax.update_layout(title=dict(text=title, x=0.5))
    ax.show()

img_path = f"{ROOT}/figures/{today}-rolling-mae-vs-hull-dist.pdf"
# fig.savefig(img_path)
