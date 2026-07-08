"""Histogram of the energy difference (either according to DFT ground truth [default] or
model predicted energy) to the convex hull for materials in the WBM data set. The
histogram stacks true/false positives/negatives with different colors.
"""

# %%
from typing import Any, Final

from pymatviz.enums import Key

from matbench_discovery import figs, plots
from matbench_discovery.cli import cli_args, complete_models
from matbench_discovery.data import load_df_wbm_with_preds
from matbench_discovery.enums import MbdKey, TestSubset
from matbench_discovery.metrics.discovery import classify_stable, dfs_metrics
from matbench_discovery.plots import hist_classified_stable_vs_hull_dist

__author__ = "Janosh Riebesell"
__date__ = "2022-12-01"

show_non_compliant = globals().get("show_non_compliant", cli_args.show_non_compliant)
models_to_plot = complete_models(show_non_compliant=show_non_compliant)
test_subset: TestSubset = globals().get("test_subset", TestSubset.uniq_protos)
models_to_plot = sorted(  # sort models by F1
    models_to_plot,
    key=lambda model: -dfs_metrics[test_subset][model.label][Key.f1.symbol],
)
df_preds = load_df_wbm_with_preds(models=models_to_plot, subset=test_subset)


# %%
n_cols = 3
models_to_plot, n_rows = plots.calc_tile_grid(
    models_to_plot, n_cols, use_full_rows=cli_args.use_full_rows
)

hover_cols = (df_preds.index.name, MbdKey.e_form_dft, MbdKey.each_true, Key.formula)
facet_col = "Model"

# 'true' or 'pred': whether to put DFT or model-predicted hull distances on the x-axis
which_energy: Final = "pred"
hist_clf_kwargs = dict(
    facet_col=facet_col,
    facet_col_wrap=n_cols,
    category_orders={facet_col: [m.label for m in models_to_plot]},
    facet_col_spacing=0.04,
    facet_row_spacing=0.04,
)

df_melt = df_preds.melt(
    id_vars=hover_cols,
    value_vars=[model.label for model in models_to_plot],
    var_name=facet_col,
    value_name=Key.e_form_pred,
)

df_melt[Key.each_pred] = (
    df_melt[MbdKey.each_true] + df_melt[Key.e_form_pred] - df_melt[MbdKey.e_form_dft]
)

fig = hist_classified_stable_vs_hull_dist(
    df=df_melt,
    each_true_col=MbdKey.each_true,
    each_pred_col=Key.each_pred,
    which_energy=which_energy,
    rolling_acc=None,
    **hist_clf_kwargs,  # ty: ignore[invalid-argument-type]
)

metrics_in_plot_titles = True
for anno in fig.layout.annotations:
    model_name = anno.text = anno.text.split("=", 1).pop()
    if (
        model_name not in [m.label for m in models_to_plot]
        or not metrics_in_plot_titles
    ):
        continue
    F1 = dfs_metrics[test_subset][model_name]["F1"]
    anno.text = f"{model_name} · {F1=:.2f}"

# set the figure size based on the number of rows and columns
fig.layout.height = 230 * n_rows
fig.layout.width = 280 * n_cols
fig.layout.paper_bgcolor = "rgba(0,0,0,0)"

# set the shared y and x axis ranges
fig.update_yaxes(range=[0, 9_000], title_text=None, matches=None, tickformat="s")
fig.update_xaxes(range=[-0.4, 0.4], title_text=None, matches=None)

# shared x/y axis titles and standardized margins
plots.style_tiled_fig(
    fig, MbdKey.each_true.label, Key.each_pred.label, n_rows=n_rows, n_cols=n_cols
)

# place the legend above the subplots
fig.layout.legend.update(
    y=1.03, xanchor="center", x=0.5, bgcolor="rgba(0,0,0,0)", orientation="h"
)

fig.show()


# %%
# site payload: per-model stability-classification counts on shared hull-dist bins
# (binned over the displayed x-range only; the old fig binned (-0.7, 0.7) but clipped
# the view to +/-0.4 eV/atom)
hist_clf_kwargs: dict[str, Any] = dict(bins=64, value_range=(-0.45, 0.45))
clf_models: list[dict[str, object]] = []
for model in models_to_plot:
    df_model = df_melt.query(f"{facet_col} == {model.label!r}")
    true_pos, false_neg, false_pos, true_neg = classify_stable(
        df_model[MbdKey.each_true], df_model[Key.each_pred]
    )
    each_pred = df_model[Key.each_pred]
    f1_score = dfs_metrics[test_subset][model.label]["F1"]
    clf_models.append(
        {"key": model.key, "label": model.label, "f1": round(float(f1_score), 4)}
        | {
            key: figs.histogram(each_pred[mask], **hist_clf_kwargs)["y"]
            for key, mask in (
                ("tp", true_pos),
                ("fn", false_neg),
                ("fp", false_pos),
                ("tn", true_neg),
            )
        }
    )
# bin centers depend only on hist_clf_kwargs (bins/range), not the per-model data
bin_centers = figs.histogram([], **hist_clf_kwargs)["x"]
if show_non_compliant:  # site payload = full model set, sorted by F1 (site renders
    # panels in payload order)
    figs.write_site_payload(
        f"hist-clf-{which_energy}-hull-dist",
        {"bin_centers": bin_centers, "models": clf_models},
    )
