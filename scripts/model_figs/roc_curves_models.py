"""Generate receiver-operating-characteristic curve payloads for each model."""

# %%
import sklearn.metrics as sk_metrics

from matbench_discovery import STABILITY_THRESHOLD, figs, payload_numerics
from matbench_discovery.cli import cli_args, complete_models
from matbench_discovery.data import load_discovery_predictions
from matbench_discovery.enums import MbdKey, TestSubset

df_preds, df_each_pred, _df_each_err = load_discovery_predictions()
test_subset = cli_args.test_subset
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
    fpr, tpr = payload_numerics.lttb(fpr, tpr, 200)
    model_data: dict[str, object] = {
        **figs.discovery_model_identity(model),
        "auc": round(float(auc), 2),
        "fpr": payload_numerics.round_list(fpr),
        "tpr": payload_numerics.round_list(tpr),
    }
    roc_models.append(model_data)
provenance = figs.build_discovery_payload_provenance(
    generator=__file__,
    test_subset=test_subset.value,
    source_files={"payload_numerics": payload_numerics.__file__},
    parameters={
        "auc_decimals": 2,
        "coordinate_decimals": payload_numerics.COORD_DECIMALS,
        "lttb_points": 200,
        "stability_threshold": STABILITY_THRESHOLD,
    },
    packages=("scikit-learn",),
)
figs.write_site_payload("roc-models", {"provenance": provenance, "models": roc_models})
