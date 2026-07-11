"""Generate binned stability-classification payload data for each model.

The histograms use model-predicted energy to the convex hull for materials in the WBM
data set and separate true/false positives/negatives.
"""

from matbench_discovery import figs
from matbench_discovery.cli import complete_models
from matbench_discovery.data import load_df_wbm_with_preds
from matbench_discovery.enums import MbdKey, TestSubset
from matbench_discovery.metrics.discovery import classify_stable, dfs_metrics

__author__ = "Janosh Riebesell"
__date__ = "2022-12-01"

models_to_plot = complete_models()
test_subset: TestSubset = globals().get("test_subset", TestSubset.uniq_protos)
df_preds = load_df_wbm_with_preds(models=models_to_plot, subset=test_subset)


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
    f1_score = dfs_metrics[test_subset][model.label]["F1"]
    clf_models.append(
        {"key": model.key, "label": model.label, "f1": round(float(f1_score), 4)}
        | {
            key: figs.histogram(each_pred[mask], bins=n_bins, value_range=hist_range)[
                "y"
            ]
            for key, mask in (
                ("tp", true_pos),
                ("fn", false_neg),
                ("fp", false_pos),
                ("tn", true_neg),
            )
        }
    )
# bin centers depend only on the bin count/range, not the per-model data
bin_centers = figs.histogram([], bins=n_bins, value_range=hist_range)["x"]
figs.write_site_payload(
    "hist-clf-pred-hull-dist",
    {"bin_centers": bin_centers, "models": clf_models},
)
