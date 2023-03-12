"""Analyze structures and composition with largest mean error across all models.
Maybe there's some chemistry/region of materials space that all models struggle with?
Might point to deficiencies in the data or models architecture.
"""


# %%
import itertools

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matminer.featurizers.site import CrystalNNFingerprint
from matminer.featurizers.structure import SiteStatsFingerprint
from pymatgen.core import Composition, Element, Structure
from pymatviz import count_elements, plot_structure_2d, ptable_heatmap_plotly
from tqdm import tqdm

from matbench_discovery import MODELS, ROOT
from matbench_discovery.data import DATA_FILES
from matbench_discovery.data import df_wbm as df_summary
from matbench_discovery.metrics import classify_stable
from matbench_discovery.preds import (
    df_each_err,
    df_each_pred,
    df_metrics,
    df_preds,
    each_true_col,
)

__author__ = "Janosh Riebesell"
__date__ = "2023-02-15"

df_each_err[each_true_col] = df_preds[each_true_col]
mean_ae_col = "All models MAE (eV/atom)"
df_each_err[mean_ae_col] = df_preds[mean_ae_col] = df_each_err.abs().mean(axis=1)


# %%
df_wbm = pd.read_json(DATA_FILES.wbm_cses_plus_init_structs).set_index("material_id")


# %%
n_rows, n_cols = 5, 4
for good_bad, init_final in itertools.product(("best", "worst"), ("initial", "final")):
    fig, axs = plt.subplots(n_rows, n_cols, figsize=(4 * n_cols, 3 * n_rows))
    n_structs = len(axs.flat)
    struct_col = {
        "initial": "initial_structure",
        "final": "computed_structure_entry",
    }[init_final]

    errs = {
        "best": df_each_err[mean_ae_col].nsmallest(n_structs),
        "worst": df_each_err[mean_ae_col].nlargest(n_structs),
    }[good_bad]
    title = (
        f"{good_bad.title()} {len(errs)} {init_final} structures (across "
        f"{len(list(df_each_pred))} models)\nErrors in (ev/atom)"
    )
    fig.suptitle(title, fontsize=20, fontweight="bold", y=1.05)

    for idx, (ax, (id, error)) in enumerate(zip(axs.flat, errs.items()), 1):
        struct = df_wbm[struct_col].loc[id]
        if init_final == "relaxed":
            struct = struct["structure"]
        struct = Structure.from_dict(struct)
        plot_structure_2d(struct, ax=ax)
        _, spg_num = struct.get_space_group_info()
        formula = struct.composition.reduced_formula
        ax.set_title(
            f"{idx}. {formula} (spg={spg_num})\n{id} {error=:.2f}", fontweight="bold"
        )
    out_path = f"{ROOT}/tmp/figures/{good_bad}-{len(errs)}-structures-{init_final}.webp"
    fig.savefig(out_path, dpi=300)


# %%
n_structs = 100
worst_ids = df_each_err[mean_ae_col].nlargest(n_structs).index.tolist()
best_ids = df_each_err[mean_ae_col].nsmallest(n_structs).index.tolist()

best_init_structs = df_wbm.initial_structure.loc[best_ids].map(Structure.from_dict)
worst_init_structs = df_wbm.initial_structure.loc[worst_ids].map(Structure.from_dict)
best_final_structs = df_wbm.computed_structure_entry.loc[best_ids].map(
    lambda cse: Structure.from_dict(cse["structure"])
)
worst_final_structs = df_wbm.computed_structure_entry.loc[worst_ids].map(
    lambda cse: Structure.from_dict(cse["structure"])
)


# %%
cnn_fp = CrystalNNFingerprint.from_preset("ops")
site_stats_fp = SiteStatsFingerprint(
    cnn_fp, stats=("mean", "std_dev", "minimum", "maximum")
)

worst_fp_diff_norms = (
    worst_final_structs.map(site_stats_fp.featurize).map(np.array)
    - worst_init_structs.map(site_stats_fp.featurize).map(np.array)
).map(np.linalg.norm)

best_fp_diff_norms = (
    best_final_structs.map(site_stats_fp.featurize).map(np.array)
    - best_init_structs.map(site_stats_fp.featurize).map(np.array)
).map(np.linalg.norm)

df_fp = pd.DataFrame(
    [worst_fp_diff_norms.values, best_fp_diff_norms.values],
    index=["highest-error structures", "lowest-error structures"],
).T


# %%
fig = df_fp.plot.hist(backend="plotly", nbins=50, barmode="overlay", opacity=0.8)
title = (
    f"SiteStatsFingerprint norm-diff between initial/final {n_structs}<br>"
    f"highest/lowest-error structures (mean over {len(list(df_each_pred))} models)"
)
fig.layout.title.update(text=title, font_size=20, xanchor="center", x=0.5)
fig.layout.legend.update(
    title="", yanchor="top", y=0.98, xanchor="right", x=0.98, font_size=16
)
fig.layout.xaxis.title = "|SSFP<sub>initial</sub> - SSFP<sub>final</sub>|"
fig.show()
fig.write_image(
    f"{ROOT}/tmp/figures/init-final-fp-diff-norms.webp", width=1000, scale=2
)


# %% plotly scatter plot of largest model errors with points sized by mean error and
# colored by true stability
fig = df_preds.nlargest(200, mean_ae_col).plot.scatter(
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
        df_each_pred[model], df_preds[each_true_col]
    )
    df_preds[f"{model}_true_pos"] = true_pos
    df_preds[f"{model}_false_neg"] = false_neg
    df_preds[f"{model}_false_pos"] = false_pos
    df_preds[f"{model}_true_neg"] = true_neg


df_preds["all_true_pos"] = df_preds.filter(like="_true_pos").all(axis=1)
df_preds["all_false_neg"] = df_preds.filter(like="_false_neg").all(axis=1)
df_preds["all_false_pos"] = df_preds.filter(like="_false_pos").all(axis=1)
df_preds["all_true_neg"] = df_preds.filter(like="_true_neg").all(axis=1)

df_preds.filter(like="all_").sum()


# %%
elem_counts: dict[str, pd.Series] = {}
for col in ("all_false_neg", "all_false_pos"):
    elem_counts[col] = elem_counts.get(col, count_elements(df_preds.query(col).formula))
    fig = ptable_heatmap_plotly(elem_counts[col], font_size=10)
    fig.layout.title = col
    fig.show()


# %% scatter plot error by element against prevalence in training set
df_mp = pd.read_csv(DATA_FILES.mp_energies, na_filter=False).set_index("material_id")
# compute number of samples per element in training set
# counting element occurrences not weighted by composition, assuming model don't learn
# much more about iron and oxygen from Fe2O3 than from FeO

count_col = "MP Occurrences"
df_elem_err = count_elements(df_mp.formula_pretty, count_mode="occurrence").to_frame(
    name=count_col
)

title = "Number of MP structures containing each element"
fig = df_elem_err[count_col].plot.bar(backend="plotly", title=title)
fig.update_layout(showlegend=False)
fig.show()

fig = ptable_heatmap_plotly(df_elem_err[count_col], font_size=10)
fig.layout.title.update(text=title, x=0.35, y=0.9, font_size=20)
fig.show()


# %% map average model error onto elements
df_summary["fractional_composition"] = [
    Composition(comp).fractional_composition for comp in tqdm(df_summary.formula)
]

df_frac_comp = pd.json_normalize(
    [comp.as_dict() for comp in df_summary["fractional_composition"]]
).set_index(df_summary.index)
assert all(
    df_frac_comp.sum(axis=1).round(6) == 1
), "composition fractions don't sum to 1"

(len(df_frac_comp) - df_frac_comp.isna().sum()).sort_values().plot.bar(backend="plotly")

# df_frac_comp = df_frac_comp.dropna(axis=1, thresh=100)  # remove Xe with only 1 entry


# %%
for model in (*df_metrics, mean_ae_col):
    df_elem_err[model] = (
        df_frac_comp * df_each_err[model].abs().values[:, None]
    ).mean()
    fig = ptable_heatmap_plotly(
        df_elem_err[model],
        precision=".2f",
        fill_value=None,
        cbar_max=0.2,
        colorscale="Turbo",
    )
    fig.layout.title.update(text=model, x=0.35, y=0.9, font_size=20)
    fig.show()


# %%
df_elem_err.to_json(f"{MODELS}/per-element/per-element-model-each-errors.json")


# %%
df_elem_err["elem_name"] = [Element(el).long_name for el in df_elem_err.index]
fig = df_elem_err.plot.scatter(
    x=count_col,
    y=mean_ae_col,
    backend="plotly",
    hover_name="elem_name",
    text=df_elem_err.index.where(
        (df_elem_err[mean_ae_col] > 0.04) | (df_elem_err[count_col] > 10_000)
    ),
    title="Correlation between element-error and element-occurrence in<br>training "
    f"set: {df_elem_err[mean_ae_col].corr(df_elem_err[count_col]):.2f}",
    hover_data={mean_ae_col: ":.2f", count_col: ":,.0f"},
)

fig.update_traces(textposition="top center")
fig.show()

# save_fig(fig, f"{ROOT}/tmp/figures/element-occu-vs-err.webp", scale=2)
# save_fig(fig, f"{ROOT}/tmp/figures/element-occu-vs-err.pdf")


# %%
df_each_err.abs().mean().sort_values()
df_each_err.abs().mean(axis=1).nlargest(25)


# %% get mean distance to convex hull for each classification
df_preds.query("all_true_pos").describe()
df_preds.query("all_false_pos").describe()
