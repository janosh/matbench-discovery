"""Generate binned stability-classification payload data for each model.

The histograms use model-predicted energy to the convex hull for materials in the WBM
data set and separate true/false positives/negatives.
"""

# %%
from matbench_discovery import ROOT, STABILITY_THRESHOLD, figs, payload_numerics
from matbench_discovery.cli import cli_args, complete_models
from matbench_discovery.data import load_df_wbm_with_preds
from matbench_discovery.enums import MbdKey, TestSubset
from matbench_discovery.metrics.discovery import classify_stable, stable_metrics

test_subset = cli_args.test_subset
models_to_plot = complete_models()
load_subset = test_subset if test_subset == TestSubset.uniq_protos else None
df_preds = load_df_wbm_with_preds(models=models_to_plot, subset=load_subset)


# %%
# site payload: per-model stability-classification counts on shared hull-dist bins
# (binned over the displayed x-range only)
n_bins = 64
hist_range = (-0.45, 0.45)
clf_models: list[dict[str, object]] = []
for model in models_to_plot:
    each_pred = (
        df_preds[MbdKey.each_true] + df_preds[model.label] - df_preds[MbdKey.e_form_dft]
    )
    true_pos, false_neg, false_pos, true_neg = classify_stable(
        df_preds[MbdKey.each_true], each_pred
    )
    f1_score = stable_metrics(
        df_preds[MbdKey.each_true],
        each_pred,
        stability_threshold=STABILITY_THRESHOLD,
    )["F1"]
    clf_models.append(
        {
            **figs.discovery_model_identity(model),
            "f1": round(f1_score, 4),
        }
        | {
            key: payload_numerics.histogram(
                each_pred[mask], bins=n_bins, value_range=hist_range
            )["y"]
            for key, mask in (
                ("tp", true_pos),
                ("fn", false_neg),
                ("fp", false_pos),
                ("tn", true_neg),
            )
        }
    )
# bin centers depend only on the bin count/range, not the per-model data
bin_centers = payload_numerics.histogram([], bins=n_bins, value_range=hist_range)["x"]
provenance = figs.build_discovery_payload_provenance(
    generator=__file__,
    test_subset=test_subset.value,
    source_files={
        "payload_numerics": payload_numerics.__file__,
        "stability_metrics": f"{ROOT}/matbench_discovery/metrics/discovery.py",
    },
    parameters={
        "coordinate_decimals": payload_numerics.COORD_DECIMALS,
        "f1_decimals": 4,
        "histogram_bins": n_bins,
        "histogram_range": list(hist_range),
        "stability_threshold": STABILITY_THRESHOLD,
    },
    packages=("scikit-learn",),
    prediction_round_decimals=None,
)
figs.write_site_payload(
    "hist-clf-pred-hull-dist",
    {
        **provenance,
        "bin_centers": bin_centers,
        "models": clf_models,
    },
)
