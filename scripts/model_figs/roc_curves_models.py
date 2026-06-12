"""Plot ROC and PR (precision-recall) curves for each model."""

# %%
import pymatviz as pmv

from matbench_discovery import PDF_FIGS, STABILITY_THRESHOLD, figs
from matbench_discovery.cli import cli_args, complete_models
from matbench_discovery.enums import MbdKey, Model, TestSubset
from matbench_discovery.metrics.discovery import df_metrics
from matbench_discovery.preds.discovery import df_each_pred, df_preds

__author__ = "Janosh Riebesell"
__date__ = "2023-01-30"


test_subset = globals().get("test_subset", TestSubset.uniq_protos)

if test_subset == TestSubset.uniq_protos:
    df_preds = df_preds.query(MbdKey.uniq_proto)
    df_each_pred = df_each_pred.loc[df_preds.index]

show_non_compliant = globals().get("show_non_compliant", cli_args.show_non_compliant)
models_to_plot = [
    model.label for model in complete_models(show_non_compliant=show_non_compliant)
]


# %% Convert E_(hull dist) continuous targets to binary classification labels
binary_targets = (df_preds[MbdKey.each_true] > STABILITY_THRESHOLD).astype(int)

fig = pmv.roc_curve(
    targets=binary_targets,
    probs_positive={
        col: (df_each_pred[col].to_numpy()) for col in df_each_pred[models_to_plot]
    },
)
# Show only the top N models by default
show_n_best_models = 10
best_models = (
    df_metrics[models_to_plot].loc["F1"].sort_values().index[-show_n_best_models:]
)
for trace in fig.data:
    trace.visible = (
        True
        if any(
            (trace.name or "").startswith((model_name, "No skill"))
            for model_name in best_models
        )
        else "legendonly"
    )

fig.show()


# %%
img_suffix = "" if show_non_compliant else "-only-compliant"
img_name = f"roc-models{img_suffix}"
roc_models = []
for trace in fig.data:
    # skip the "No skill" random-classifier diagonal (drawn inline on the site)
    if " · AUC=" not in (trace.name or ""):
        continue
    label, auc_str = trace.name.split(" · AUC=")
    # ROC staircases at full resolution are ~4x over-resolved for a 480px panel
    fpr, tpr = figs.lttb(*figs.trace_xy(trace), 200)
    roc_models.append(
        {
            "key": Model.from_label(label).key,
            "label": label,
            "auc": float(auc_str),
            "fpr": figs.round_list(fpr),
            "tpr": figs.round_list(tpr),
        }
    )
if show_non_compliant:  # site payload = full model set;
    # the compliant-only variant exists as PDF only (paper SI)
    figs.write_site_payload("roc-models", {"models": roc_models})
pmv.save_fig(fig, f"{PDF_FIGS}/{img_name}.pdf")
