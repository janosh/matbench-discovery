"""Generate cumulative precision and recall payload data for all models.

Cumulative here means descending the list of test-set materials ranked by
model-predicted stability, starting from the most stable, and updating each metric
after every new material. This simulates a materials screening campaign and shows the
expected hit rate for a given DFT calculation budget.
"""

# %%
import numpy as np

from matbench_discovery import STABILITY_THRESHOLD, figs
from matbench_discovery.cli import complete_models, shared_payload_test_subset
from matbench_discovery.data import load_discovery_predictions
from matbench_discovery.enums import MbdKey, TestSubset
from matbench_discovery.metrics.discovery import classify_stable

df_preds, df_each_pred, _df_each_err = load_discovery_predictions()
test_subset = shared_payload_test_subset()
if test_subset == TestSubset.uniq_protos:
    df_preds = df_preds.query(MbdKey.uniq_proto)
    df_each_pred = df_each_pred.loc[df_preds.index]


# %%
cum_pr_models = []
for model in complete_models():
    each_pred = (
        df_each_pred[model.label].sort_index(kind="stable").sort_values(kind="stable")
    )
    each_true = df_preds[MbdKey.each_true].loc[each_pred.index]
    true_pos, false_neg, false_pos, _true_neg = classify_stable(
        each_true, each_pred, stability_threshold=STABILITY_THRESHOLD
    )
    n_true_pos_cum = true_pos.cumsum()  # all pd.Series, cumsum stays a Series
    precision_cum = n_true_pos_cum / (n_true_pos_cum + false_pos.cumsum())
    recall_cum = n_true_pos_cum / (n_true_pos_cum + false_neg.cumsum()).iloc[-1]
    # number of materials the model predicts stable = where its curve ends
    n_pred_stable = int((each_pred <= STABILITY_THRESHOLD).sum())
    if n_pred_stable < 2:  # can't happen for real models (thousands stable)
        raise ValueError(f"{model.label} predicts {n_pred_stable} stable materials")
    # log2-spaced sampling for higher density at the start of the discovery campaign
    # where metrics fluctuate most. Rounded to ints since x counts screened materials
    # (also compresses better).
    log_xs = np.logspace(0, np.log2(n_pred_stable - 1), 100, base=2)
    xs = np.unique([*log_xs.round().astype(int), n_pred_stable])
    # xs are 1-based material counts, so the curves can be sampled exactly by indexing
    # (no interpolation: a spline fit here only added LAPACK float noise)
    precision, recall = (
        curve.to_numpy()[xs - 1] for curve in (precision_cum, recall_cum)
    )
    cum_pr_models.append(
        {
            "key": model.key,
            "label": model.label,
            "x": figs.round_list(xs),
            "precision": figs.round_list(precision),
            "recall": figs.round_list(recall),
            # [n materials predicted stable, precision there, recall there]
            "end": [
                n_pred_stable,
                round(float(precision_cum.iloc[n_pred_stable - 1]), 5),
                round(float(recall_cum.iloc[n_pred_stable - 1]), 5),
            ],
        }
    )
n_stable = int((df_preds[MbdKey.each_true] <= STABILITY_THRESHOLD).sum())
figs.write_site_payload(
    "cumulative-precision-recall",
    {"n_stable": n_stable, "models": cum_pr_models},
)
