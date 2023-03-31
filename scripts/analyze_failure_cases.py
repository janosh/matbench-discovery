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
from sklearn.metrics import r2_score
from tqdm import tqdm

from matbench_discovery import FIGS, MODELS, ROOT
from matbench_discovery.data import DATA_FILES, df_wbm
from matbench_discovery.metrics import classify_stable
from matbench_discovery.preds import (
    df_each_err,
    df_each_pred,
    df_metrics,
    df_preds,
    each_pred_col,
    each_true_col,
)

__author__ = "Janosh Riebesell"
__date__ = "2023-02-15"

m3ae_col = "Mean of models MAE"
df_each_err[m3ae_col] = df_each_err.abs().mean(axis=1)
fp_diff_col = "site_stats_fingerprint_init_final_norm_diff"


# %%
df_mp = pd.read_csv(DATA_FILES.mp_energies, na_filter=False).set_index("material_id")
# compute number of samples per element in training set
# counting element occurrences not weighted by composition, assuming model don't learn
# much more about iron and oxygen from Fe2O3 than from FeO

train_count_col = "MP Occurrences"
df_elem_err = count_elements(df_mp.formula_pretty, count_mode="occurrence").to_frame(
    name=train_count_col
)


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
n_structs = 1000
fig = go.Figure()
for idx, model in enumerate((m3ae_col, *df_metrics)):
    large_errors = df_each_err[model].abs().nlargest(n_structs)
    small_errors = df_each_err[model].abs().nsmallest(n_structs)
    for label, errors in zip(("min", "max"), (large_errors, small_errors)):
        scatter = go.Histogram(
            x=df_wbm.loc[errors.index][fp_diff_col].values,
            name=f"{model} err<sub>{label}</sub>",
            visible="legendonly" if idx else True,
            legendgroup=model,
            hovertemplate=("SSFP diff: %{x:.2f}<br>Count: %{y}"),
        )
        fig.add_trace(scatter)

title = (
    f"Norm-diff between initial/final SiteStatsFingerprint<br>"
    f"of the {n_structs} highest and lowest error structures for each model"
)
fig.layout.title.update(text=title, xanchor="center", x=0.5)
fig.layout.legend.update(
    title="",
    yanchor="top",
    y=0.98,
    xanchor="right",
    x=0.98,
    font_size=12,
    orientation="h",
)
fig.layout.xaxis.title = "|SSFP<sub>initial</sub> - SSFP<sub>final</sub>|"
fig.layout.yaxis.title = "Count"

fig.show()

# save_fig(fig, f"{FIGS}/hist-largest-each-errors-fp-diff-models.svelte")


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
    title="", yanchor="top", y=0.98, xanchor="right", x=0.98, font_size=14
)
fig.layout.xaxis.title = "|SSFP<sub>initial</sub> - SSFP<sub>final</sub>|"
fig.layout.yaxis.title = "Absolute error (eV/atom)"

fig.show()

save_fig(fig, f"{FIGS}/scatter-largest-each-errors-fp-diff-models.svelte")


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
title = "Number of MP structures containing each element"
fig = df_elem_err[train_count_col].plot.bar(backend="plotly", title=title)
fig.update_layout(showlegend=False)
fig.show()

fig = ptable_heatmap_plotly(df_elem_err[train_count_col], font_size=10)
fig.layout.title.update(text=title, x=0.35, y=0.9, font_size=20)
fig.show()


# %% map average model error onto elements
frac_comp_col = "fractional composition"
df_wbm[frac_comp_col] = [
    Composition(comp).fractional_composition for comp in tqdm(df_wbm.formula)
]

df_frac_comp = pd.json_normalize(
    [comp.as_dict() for comp in df_wbm[frac_comp_col]]
).set_index(df_wbm.index)
assert all(
    df_frac_comp.sum(axis=1).round(6) == 1
), "composition fractions don't sum to 1"

# bar plot showing number of structures in MP containing each element
(len(df_frac_comp) - df_frac_comp.isna().sum()).sort_values().plot.bar(backend="plotly")

# df_frac_comp = df_frac_comp.dropna(axis=1, thresh=100)  # remove Xe with only 1 entry


# %%
comp_col = "composition"
df_wbm[comp_col] = [Composition(comp) for comp in tqdm(df_wbm.formula)]

test_set_std_col = "Test set standard deviation (eV/atom)"
for elem in tqdm(df_elem_err.index):
    mask = df_wbm[comp_col].map(lambda comp: elem in comp)  # noqa: B023
    elem_test_set_std = df_wbm[mask][each_true_col].std()
    df_elem_err.loc[elem, test_set_std_col] = elem_test_set_std


# %%
fig = ptable_heatmap_plotly(
    df_elem_err[test_set_std_col], precision=".2f", colorscale="Inferno"
)
fig.show()


# %%
normalized = True
# cs_range = (None, 0.5)
for model in (*df_metrics, m3ae_col):
    df_elem_err[model] = (
        df_frac_comp * df_each_err[model].abs().values[:, None]
    ).mean()
    per_elem_err = df_elem_err[model]
    per_elem_err.name = f"{model} (eV/atom)"
    if normalized:
        per_elem_err /= df_elem_err[test_set_std_col]
        per_elem_err.name = f"{model} (normalized by test set std)"
    fig = ptable_heatmap_plotly(per_elem_err, precision=".2f", colorscale="Inferno")
    fig.show()


# %%
df_elem_err.round(4).to_json(f"{MODELS}/tmi/per-element-model-each-errors.json")


# %% check correlation and R2 of elemental prevalence in MP training data vs.
# model error
df_elem_err["elem_name"] = [Element(el).long_name for el in df_elem_err.index]
R2 = r2_score(*df_elem_err[[train_count_col, m3ae_col]].dropna().values.T)
r_P = df_elem_err[m3ae_col].corr(df_elem_err[train_count_col])

fig = df_elem_err.plot.scatter(
    x=train_count_col,
    y=m3ae_col,
    backend="plotly",
    hover_name="elem_name",
    text=df_elem_err.index.where(
        (df_elem_err[m3ae_col] > 0.04) | (df_elem_err[train_count_col] > 6_000)
    ),
    title="Per-element error vs element-occurrence in MP training "
    f"set: r<sub>Pearson</sub>={r_P:.2f}, R<sup>2</sup>={R2:.2f}",
    hover_data={m3ae_col: ":.2f", train_count_col: ":,.0f"},
)
fig.update_traces(textposition="top center")  # place text above scatter points
fig.layout.title.update(xanchor="center", x=0.5)
fig.show()

save_fig(fig, f"{FIGS}/element-prevalence-vs-error.svelte")
# save_fig(fig, f"{ROOT}/tmp/figures/element-prevalence-vs-error.pdf")


# %% TODO investigate if structures with largest mean of models error can be attributed
# to DFT gone wrong. would be cool if models can be run across large databases as
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
    font_size=14,
    tracegroupgap=0,
)
fig.layout.xaxis.title = "|SSFP<sub>initial</sub> - SSFP<sub>final</sub>|"
fig.layout.yaxis.title = "|E<sub>above hull</sub> error| (eV/atom)"

fig.show()

save_fig(fig, f"{FIGS}/largest-fp-diff-each-error-models.svelte")
# save_fig(fig, f"{ROOT}/tmp/figures/large-fp-diff-vs-each-error.webp", scale=2)


# %% plot EACH errors against least prevalent element in structure (by occurrence in
# MP training set)
n_examp_for_rarest_elem_col = "Examples for rarest element in structure"
df_wbm["composition"] = df_wbm.get("composition", df_wbm.formula.map(Composition))
df_elem_err.loc[list(map(str, df_wbm.composition[0]))][train_count_col].min()
df_wbm[n_examp_for_rarest_elem_col] = [
    df_elem_err.loc[list(map(str, Composition(formula)))][train_count_col].min()
    for formula in tqdm(df_wbm.formula)
]


# %%
n_bins = 100
df_melt = (
    df_each_err.abs()
    .reset_index()
    .melt(var_name="Model", value_name=each_pred_col, id_vars="material_id")
    .set_index("material_id")
)
df_melt[n_examp_for_rarest_elem_col] = df_wbm[n_examp_for_rarest_elem_col]
df_melt["x_bin"] = pd.cut(df_melt[n_examp_for_rarest_elem_col], bins=n_bins)
df_melt["y_bin"] = pd.cut(df_melt[each_pred_col], bins=n_bins)

df_plot = df_melt.reset_index().groupby(["x_bin", "y_bin", "Model"]).first().dropna()
df_plot = df_plot.reset_index().set_index("material_id")
df_plot["formula"] = df_wbm.formula
print(f"{len(df_plot)=:,} / {len(df_melt)=:,} = {len(df_plot)/len(df_melt):.1%}")


# %%
fig = px.scatter(
    df_plot.reset_index(),
    x=n_examp_for_rarest_elem_col,
    y=each_pred_col,
    color="Model",
    facet_col="Model",
    facet_col_wrap=3,
    hover_data=dict(material_id=True, formula=True, Model=False),
    title="Absolute errors in model-predicted E<sub>above hull</sub> vs. occurrence "
    "count in MP training set<br>of least prevalent element in structure",
)
fig.layout.update(showlegend=False)
fig.layout.title.update(x=0.5, xanchor="center", y=0.95)
fig.layout.margin.update(t=100)
# remove axis labels
fig.update_xaxes(title="")
fig.update_yaxes(title="")
for anno in fig.layout.annotations:
    anno.text = anno.text.split("=")[1]

fig.add_annotation(
    text="MP occurrence count of least prevalent element in structure",
    x=0.5,
    y=-0.18,
    xref="paper",
    yref="paper",
    showarrow=False,
)
fig.add_annotation(
    text="Absolute error in E<sub>above hull</sub>",
    x=-0.07,
    y=0.5,
    xref="paper",
    yref="paper",
    showarrow=False,
    textangle=-90,
)

fig.show()
save_fig(fig, f"{FIGS}/each-error-vs-occu-count-rarest-element-in-struct.svelte")


# %%
fig = go.Figure()
for idx, model in enumerate(df_metrics):
    scatter = go.Scatter(
        x=df_wbm[n_examp_for_rarest_elem_col],
        y=df_each_err[model].abs(),
        mode="markers",
        name=model,
        hovertemplate=(
            "ID: %{customdata[0]}<br>"
            "formula: %{customdata[1]}<br>"
            "least prevalent element: %{customdata[2]}<br>"
            "error: %{y}<extra></extra>"
        ),
        customdata=df_preds[["material_id", "formula"]]
        .loc[df_each_err[model].index]
        .values,
        legendgroup=model,
        visible="legendonly" if idx else True,
        # marker=dict(color=color),
        legendrank=df_each_err[model].abs().mean(),
    )
    fig.add_trace(scatter)
    break
fig.show()

# save_fig(fig, f"{FIGS}/least-prevalent-element-vs-error.svelte")
# save_fig(fig, f"{ROOT}/tmp/figures/least-prevalent-element-vs-error.pdf")


# %%
tsne_cols = ["t-SNE 1", "t-SNE 2"]
df_wbm[tsne_cols] = pd.read_csv(
    f"{ROOT}/data/wbm/tsne/onehot-112-composition-2d.csv", index_col=0
)

df_wbm["wbm_step"] = df_wbm.index.str.split("-").str[1]


# %%
fig = px.scatter(
    df_wbm,
    x=tsne_cols[0],
    y=tsne_cols[1],
    color=(df_wbm.bandgap_pbe > 1).map({True: "band gap", False: "no gap"}),
    hover_name="material_id",
    hover_data=("formula",),
)
fig.show()

save_fig(
    fig, f"{ROOT}/tmp/figures/tsne-2d-composition-by-wbm-step-bandgap.png", scale=3
)
