# %%
import matplotlib.pyplot as plt
import pandas as pd
from pymatgen.core import Structure
from pymatviz import plot_structure_2d, ptable_heatmap_plotly

from matbench_discovery import ROOT
from matbench_discovery.metrics import classify_stable
from matbench_discovery.preds import df_each_err, df_each_pred, df_wbm, each_true_col

__author__ = "Janosh Riebesell"
__date__ = "2023-02-15"

df_each_err[each_true_col] = df_wbm[each_true_col]
mean_ae_col = "All models mean absolute error (eV/atom)"
df_each_err[mean_ae_col] = df_wbm[mean_ae_col] = df_each_err.abs().mean(axis=1)


# %%
cse_path = f"{ROOT}/data/wbm/2022-10-19-wbm-computed-structure-entries.json.bz2"
df_cse = pd.read_json(cse_path).set_index("material_id")


# %%
n_rows, n_cols = 5, 4
for which in ("best", "worst"):
    fig, axs = plt.subplots(n_rows, n_cols, figsize=(3 * n_rows, 4 * n_cols))
    n_axs = len(axs.flat)

    errs = (
        df_each_err.mean_ae.nsmallest(n_axs)
        if which == "best"
        else df_each_err.mean_ae.nlargest(n_axs)
    )
    title = f"{which} {len(errs)} structures (across {len(list(df_each_pred))} models)"
    fig.suptitle(title, fontsize=16, fontweight="bold", y=0.95)

    for idx, (ax, (id, err)) in enumerate(zip(axs.flat, errs.items()), 1):
        struct = Structure.from_dict(
            df_cse.computed_structure_entry.loc[id]["structure"]
        )
        plot_structure_2d(struct, ax=ax)
        _, spg_num = struct.get_space_group_info()
        formula = struct.composition.reduced_formula
        ax.set_title(
            f"{idx}. {formula} (spg={spg_num})\n{id} {err=:.2f}", fontweight="bold"
        )

    fig.savefig(f"{ROOT}/tmp/figures/{which}-{len(errs)}-structures.webp", dpi=300)


# %% plotly scatter plot of largest model errors with points sized by mean error and
# colored by true stability
fig = df_wbm.nlargest(200, mean_ae_col).plot.scatter(
    x=each_true_col,
    y=mean_ae_col,
    color=each_true_col,
    size=mean_ae_col,
    backend="plotly",
)
fig.layout.coloraxis.colorbar.update(
    title="DFT distance to convex hull (eV/atom)",
    title_side="top",
    yanchor="bottom",
    y=1,
    xanchor="center",
    x=0.5,
    orientation="h",
    thickness=12,
)
fig.show()


# %% find materials that were misclassified by all models
for model in df_each_pred:
    true_pos, false_neg, false_pos, true_neg = classify_stable(
        df_each_pred[model], df_wbm[each_true_col]
    )
    df_wbm[f"{model}_true_pos"] = true_pos
    df_wbm[f"{model}_false_neg"] = false_neg
    df_wbm[f"{model}_false_pos"] = false_pos
    df_wbm[f"{model}_true_neg"] = true_neg


df_wbm["all_true_pos"] = df_wbm.filter(like="_true_pos").all(axis=1)
df_wbm["all_false_neg"] = df_wbm.filter(like="_false_neg").all(axis=1)
df_wbm["all_false_pos"] = df_wbm.filter(like="_false_pos").all(axis=1)
df_wbm["all_true_neg"] = df_wbm.filter(like="_true_neg").all(axis=1)

df_wbm.filter(like="all_").sum()


# %%
ptable_heatmap_plotly(df_wbm[df_wbm.all_false_pos].formula, colorscale="Viridis")
ptable_heatmap_plotly(df_wbm[df_wbm.all_false_neg].formula, colorscale="Viridis")


# %%
df_each_err.abs().mean().sort_values()
df_each_err.abs().mean(axis=1).nlargest(25)


# %% get mean distance to convex hull for each classification
df_wbm.query("all_true_pos").describe()
df_wbm.query("all_false_pos").describe()
