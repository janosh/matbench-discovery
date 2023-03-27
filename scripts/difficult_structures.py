"""Analyze structures and composition with largest mean error across all models.
Maybe there's some chemistry/region of materials space that all models struggle with?
Might point to deficiencies in the data or models architecture.
"""


# %%
import itertools

import matplotlib.pyplot as plt
import pandas as pd
import plotly.express as px
import plotly.graph_objs as go
from pymatgen.core import Composition, Element, Structure
from pymatviz import count_elements, plot_structure_2d, ptable_heatmap_plotly
from pymatviz.utils import save_fig
from tqdm import tqdm

from matbench_discovery import FIGS, MODELS, ROOT
from matbench_discovery.data import DATA_FILES, df_wbm
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
m3ae_col = "Mean of models MAE"
df_each_err[m3ae_col] = df_preds[m3ae_col] = df_each_err.abs().mean(axis=1)
fp_diff_col = "site_stats_fingerprint_init_final_norm_diff"


# %%
df_cse = pd.read_json(DATA_FILES.wbm_cses_plus_init_structs).set_index("material_id")


# %% plot the highest and lowest error structures before and after relaxation
n_rows, n_cols = 5, 4
for good_or_bad, init_or_final in itertools.product(
    ("best", "worst"), ("initial", "final")
):
    fig, axs = plt.subplots(n_rows, n_cols, figsize=(4 * n_cols, 3 * n_rows))
    n_structs = len(axs.flat)
    struct_col = {
        "initial": "initial_structure",
        "final": "computed_structure_entry",
    }[init_or_final]

    errs = {
        "best": df_each_err[m3ae_col].nsmallest(n_structs),
        "worst": df_each_err[m3ae_col].nlargest(n_structs),
    }[good_or_bad]
    title = (
        f"{good_or_bad.title()} {len(errs)} {init_or_final} structures (across "
        f"{len(list(df_each_pred))} models)\nErrors in (ev/atom)"
    )
    fig.suptitle(title, fontsize=20, fontweight="bold", y=1.05)

    for idx, (ax, (mat_id, error)) in enumerate(zip(axs.flat, errs.items()), 1):
        struct = df_cse[struct_col].loc[mat_id]
        if "structure" in struct:
            struct = struct["structure"]
        struct = Structure.from_dict(struct)
        plot_structure_2d(struct, ax=ax)
        _, spg_num = struct.get_space_group_info()
        formula = struct.composition.reduced_formula
        ax.set_title(
            f"{idx}. {formula} (spg={spg_num})\n{mat_id} {error=:.2f}",
            fontweight="bold",
        )
    out_path = (
        f"{ROOT}/tmp/figures/{good_or_bad}-{len(errs)}-structures-{init_or_final}.webp"
    )
    # fig.savefig(out_path, dpi=300)


# %%
n_structs = 100
fig = go.Figure()
for idx, model in enumerate(df_metrics):
    errors = df_each_err[model].abs().nlargest(n_structs)
    model_mae = errors.mean().round(3)
    scatter = go.Scatter(
        x=df_wbm.loc[errors.index][fp_diff_col].values,
        y=errors.values,
        mode="markers",
        name=f"{model} · MAE={model_mae:.2f}",
        visible="legendonly" if idx else True,
        hovertemplate=(
            "material ID: %{customdata[0]}<br>"
            "formula: %{customdata[1]}<br>"
            "FP norm diff: %{x}<br>"
            "error: %{y} eV/atom"
        ),
        customdata=df_wbm.loc[errors.index][["material_id", "formula"]].values,
        legendrank=model_mae,
    )
    fig.add_trace(scatter)

title = (
    f"Norm-diff between initial/final SiteStatsFingerprint<br>"
    f"of the {n_structs} highest-error structures for each model"
)
fig.layout.title.update(text=title, xanchor="center", x=0.5)
fig.layout.legend.update(
    title="", yanchor="top", y=0.98, xanchor="right", x=0.98, font_size=16
)
fig.layout.xaxis.title = "|SSFP<sub>initial</sub> - SSFP<sub>final</sub>|"
fig.layout.yaxis.title = "Absolute error (eV/atom)"

fig.show()

save_fig(fig, f"{FIGS}/largest-each-errors-fp-diff-models.svelte")


# %% plotly scatter plot of largest model errors with points sized by mean error and
# colored by true stability
fig = df_preds.nlargest(200, m3ae_col).plot.scatter(
    x=each_true_col,
    y=m3ae_col,
    color=each_true_col,
    size=m3ae_col,
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
df_wbm["fractional_composition"] = [
    Composition(comp).fractional_composition for comp in tqdm(df_wbm.formula)
]

df_frac_comp = pd.json_normalize(
    [comp.as_dict() for comp in df_wbm["fractional_composition"]]
).set_index(df_wbm.index)
assert all(
    df_frac_comp.sum(axis=1).round(6) == 1
), "composition fractions don't sum to 1"

(len(df_frac_comp) - df_frac_comp.isna().sum()).sort_values().plot.bar(backend="plotly")

# df_frac_comp = df_frac_comp.dropna(axis=1, thresh=100)  # remove Xe with only 1 entry


# %%
for model in (*df_metrics, m3ae_col):
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
    y=m3ae_col,
    backend="plotly",
    hover_name="elem_name",
    text=df_elem_err.index.where(
        (df_elem_err[m3ae_col] > 0.04) | (df_elem_err[count_col] > 10_000)
    ),
    title="Correlation between element-error and element-occurrence in<br>training "
    f"set: {df_elem_err[m3ae_col].corr(df_elem_err[count_col]):.2f}",
    hover_data={m3ae_col: ":.2f", count_col: ":,.0f"},
)

fig.update_traces(textposition="top center")
fig.show()

# save_fig(fig, f"{ROOT}/tmp/figures/element-occu-vs-err.webp", scale=2)
# save_fig(fig, f"{ROOT}/tmp/figures/element-occu-vs-err.pdf")


# %% TODO investigate if structures with largest mean of models error can be attributed
# to SFT gone wrong. would be cool if models can be run across large databases as
# correctness checkers
df_each_err.abs().mean().sort_values()
df_each_err.abs().mean(axis=1).nlargest(25)


# %%
n_points = 1000
largest_fp_diff = df_wbm[fp_diff_col].nlargest(n_points)

fig = go.Figure()
colors = px.colors.qualitative.Plotly

for idx, model in enumerate(df_metrics):
    color = colors[idx]
    model_mae = df_each_err[model].loc[largest_fp_diff.index].abs().mean()

    visible = "legendonly" if idx else True
    scatter = go.Scatter(
        x=largest_fp_diff.values,
        y=df_each_err[model].loc[largest_fp_diff.index].abs(),
        mode="markers",
        name=f"{model} · MAE={model_mae:.2f}",
        visible=visible,
        hovertemplate=(
            "ID: %{customdata[0]}<br>"
            "formula: %{customdata[1]}<br>"
            "FP diff: %{x}<br>"
            "error: %{y}<extra></extra>"
        ),
        customdata=df_preds[["material_id", "formula"]]
        .loc[largest_fp_diff.index]
        .values,
        legendgroup=model,
        marker=dict(color=color),
        legendrank=model_mae,
    )
    fig.add_trace(scatter)
    # add dashed mean line for each model that toggles with the scatter plot
    # fig.add_hline(
    #     y=model_mae,
    #     line=dict(dash="dash"),
    #     annotation=dict(text=f"{model} mean", x=0, xanchor="left", font_size=10),
    # )
    # get color from scatter plot
    fig.add_scatter(
        x=[largest_fp_diff.min(), largest_fp_diff.max()],
        y=[model_mae, model_mae],
        line=dict(dash="dash", width=2, color=color),
        legendgroup=model,
        showlegend=False,
        visible=visible,
    )


title = (
    f"Absolute errors in model-predicted E<sub>above hull</sub> for {n_points} "
    "structures<br>with largest norm-diff of initial/final SiteStatsFingerprint in WBM "
    "test set"
)
fig.layout.title.update(text=title, xanchor="center", x=0.5)
fig.layout.legend.update(
    title="",
    yanchor="top",
    y=0.98,
    xanchor="right",
    x=0.98,
    font_size=12,
    tracegroupgap=0,
)
fig.layout.xaxis.title = "|SSFP<sub>initial</sub> - SSFP<sub>final</sub>|"
fig.layout.yaxis.title = "|E<sub>above hull</sub> error| (eV/atom)"

fig.show()

save_fig(fig, f"{FIGS}/largest-fp-diff-each-error-models.svelte")
save_fig(fig, f"{ROOT}/tmp/figures/large-fp-diff-vs-each-error.webp", scale=2)
