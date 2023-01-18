# %%
from matbench_discovery import FIGS, today
from matbench_discovery.data import load_df_wbm_with_preds
from matbench_discovery.plots import rolling_mae_vs_hull_dist

__author__ = "Rhys Goodall, Janosh Riebesell"
__date__ = "2022-06-18"


# %%
model = "Wrenformer"
df_wbm = load_df_wbm_with_preds([model]).round(3)

e_above_hull_col = "e_above_hull_mp2020_corrected_ppd_mp"
e_form_col = "e_form_per_atom_mp2020_corrected"


# %%
ax, df_err, df_std = rolling_mae_vs_hull_dist(
    e_above_hull_true=df_wbm[e_above_hull_col],
    e_above_hull_errors={model: df_wbm[e_form_col] - df_wbm[model]},
    # label=model,
    backend=(backend := "plotly"),
    # template="plotly_white",
)

title = f"{today} {model}"
if backend == "matplotlib":
    fig = ax.figure
    fig.set_size_inches(6, 5)
    ax.legend(loc="lower right", frameon=False)
    ax.set(title=title)
    for line in ax.lines:
        line._linewidth *= 2
elif backend == "plotly":
    ax.update_layout(title=dict(text=title, x=0.5))
    ax.show()

img_path = f"{FIGS}/{today}-rolling-mae-vs-hull-dist.pdf"
# fig.savefig(img_path)
