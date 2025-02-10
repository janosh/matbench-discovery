"""Plot ROC and PR (precision-recall) curves for each model."""

# %%
import pymatviz as pmv

from matbench_discovery import PDF_FIGS, SITE_FIGS, STABILITY_THRESHOLD
from matbench_discovery import plots as plots
from matbench_discovery.cli import cli_args
from matbench_discovery.enums import MbdKey, TestSubset
from matbench_discovery.models import MODEL_METADATA, model_is_compliant
from matbench_discovery.preds.discovery import df_each_pred, df_preds

__author__ = "Janosh Riebesell"
__date__ = "2023-01-30"


test_subset = globals().get("test_subset", TestSubset.uniq_protos)

if test_subset == TestSubset.uniq_protos:
    df_preds = df_preds.query(MbdKey.uniq_proto)
    df_each_pred = df_each_pred.loc[df_preds.index]

show_non_compliant = globals().get("show_non_compliant", False)
models_to_plot = [
    model.label
    for model in cli_args.models
    if show_non_compliant or model_is_compliant(MODEL_METADATA[model.label])
]
df_each_pred = df_each_pred[models_to_plot]


# %% Convert E_(hull dist) continuous targets to binary classification labels
binary_targets = (df_preds[MbdKey.each_true] > STABILITY_THRESHOLD).astype(int)

fig = pmv.roc_curve_plotly(
    targets=binary_targets,
    probs_positive={col: (df_each_pred[col].to_numpy()) for col in df_each_pred},
)
fig.show()

img_suffix = "" if show_non_compliant else "-only-compliant"
img_name = f"roc-models{img_suffix}"


# %%
pmv.save_fig(fig, f"{SITE_FIGS}/{img_name}.svelte")
pmv.save_fig(fig, f"{PDF_FIGS}/{img_name}.pdf")
