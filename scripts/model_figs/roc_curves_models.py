"""Generate receiver-operating-characteristic curve payloads for each model."""

# %%
import sklearn.metrics as sk_metrics

from matbench_discovery import STABILITY_THRESHOLD, figs
from matbench_discovery.cli import complete_models, shared_payload_test_subset
from matbench_discovery.enums import MbdKey, TestSubset
from matbench_discovery.preds.discovery import df_each_pred, df_preds

test_subset = shared_payload_test_subset()
if test_subset == TestSubset.uniq_protos:
    df_preds = df_preds.query(MbdKey.uniq_proto)
    df_each_pred = df_each_pred.loc[df_preds.index]


# %% Convert E_(hull dist) continuous targets to binary classification labels
binary_targets = (df_preds[MbdKey.each_true] > STABILITY_THRESHOLD).astype(int)


# %%
roc_models: list[dict[str, object]] = []
for model in complete_models():
    model_scores = df_each_pred[model.label].dropna()
    targets = binary_targets.loc[model_scores.index]
    fpr, tpr, _thresholds = sk_metrics.roc_curve(targets, model_scores)
    auc = sk_metrics.roc_auc_score(targets, model_scores)
    # ROC staircases at full resolution are ~4x over-resolved for a 480px panel
    fpr, tpr = figs.lttb(fpr, tpr, 200)
    roc_models.append(
        {
            "key": model.key,
            "label": model.label,
            "auc": round(float(auc), 2),
            "fpr": figs.round_list(fpr),
            "tpr": figs.round_list(tpr),
        }
    )
figs.write_site_payload("roc-models", {"models": roc_models})
