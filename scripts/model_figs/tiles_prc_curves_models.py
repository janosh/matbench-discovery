"""Plot ROC and PR (precision-recall) curves for each model."""

# %%
import pandas as pd
import plotly.express as px
import pymatviz as pmv
from pymatviz.enums import Key
from sklearn.metrics import precision_recall_curve
from tqdm import tqdm

from matbench_discovery import PDF_FIGS, STABILITY_THRESHOLD, plots
from matbench_discovery.cli import cli_args
from matbench_discovery.enums import MbdKey, TestSubset
from matbench_discovery.preds.discovery import df_each_pred, df_preds

__author__ = "Janosh Riebesell"
__date__ = "2023-01-30"


line = dict(dash="dash", width=0.5)
facet_col = "Model"
color_col = "Stability Threshold"
show_non_compliant = globals().get("show_non_compliant", cli_args.show_non_compliant)
models_to_plot = [
    model.label
    for model in cli_args.models
    if model.is_complete and (show_non_compliant or model.is_compliant)
]


test_subset = globals().get("test_subset", TestSubset.uniq_protos)

if test_subset == TestSubset.uniq_protos:
    df_preds = df_preds.query(MbdKey.uniq_proto)
    df_each_pred = df_each_pred.loc[df_preds.index]

df_each_pred = df_each_pred[models_to_plot]


# %%
df_prc = pd.DataFrame()

n_cols = 3
models_to_plot, n_rows = plots.calc_tile_grid(
    models_to_plot, n_cols, use_full_rows=globals().get("use_full_rows", True)
)

for model in (pbar := tqdm(models_to_plot, desc="Calculating PRC curves")):
    pbar.set_postfix_str(model)
    na_mask = df_preds[MbdKey.each_true].isna() | df_each_pred[model].isna()
    y_true = (df_preds[~na_mask][MbdKey.each_true] <= STABILITY_THRESHOLD).astype(int)
    # stable is the positive class (consistent with stable_metrics and the rest of
    # the benchmark). negate predicted hull distance so higher score = more stable.
    y_score = -df_each_pred[model][~na_mask]
    prec, recall, thresholds = precision_recall_curve(y_true, y_score, pos_label=1)
    dct = {
        Key.precision: prec[:-1],
        Key.recall: recall[:-1],
        color_col: -thresholds,  # undo negation to recover hull-distance thresholds
        facet_col: model,
    }
    df_prc = pd.concat([df_prc, pd.DataFrame(dct).round(3)])


# %%
fig = px.scatter(
    df_prc.iloc[:: len(df_prc) // 500 or 1],
    x=Key.recall,
    y=Key.precision,
    facet_col=facet_col,
    facet_col_wrap=n_cols,
    facet_col_spacing=0.04,
    facet_row_spacing=0.04,
    width=280 * n_cols,
    height=230 * n_rows,
    color=color_col,
    range_color=(-0.5, 1),
    hover_name=facet_col,
    hover_data={facet_col: False},
    category_orders={facet_col: models_to_plot},
)

for anno in fig.layout.annotations:
    anno.text = anno.text.split("=", 1)[1]  # remove Model= from subplot titles

# shared x/y axis titles and standardized margins
plots.style_tiled_fig(
    fig, Key.recall.label, Key.precision.label, n_rows=n_rows, n_cols=n_cols
)

# place the colorbar above the subplots
fig.layout.coloraxis.colorbar.update(
    x=0.5, y=1.03, thickness=11, len=0.8, orientation="h"
)

# set the shared x and y axis ranges to (0, 1). stable-class precision can reach
# down to the ~16-33% stable prevalence, so don't truncate the y-axis
fig.update_xaxes(title="", range=[0, 1], matches=None)
fig.update_yaxes(title="", range=[0, 1], matches=None)

fig.layout.template = "pymatviz_white"

fig.show()


# %%
img_suffix = "" if show_non_compliant else "-only-compliant"
img_name = f"prc-models-{n_rows}x{n_cols}{img_suffix}"
pmv.save_fig(fig, f"{PDF_FIGS}/{img_name}.pdf")
fig.show()
