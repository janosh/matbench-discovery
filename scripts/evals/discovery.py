"""Calculate discovery metrics for selected models and write them to model YAMLs.

The v1 paper's metric-table PDFs are frozen publication artifacts and intentionally
do not track the evolving online leaderboard.
"""

import traceback

from pymatviz.enums import Key

from matbench_discovery import STABILITY_THRESHOLD
from matbench_discovery.cli import cli_args
from matbench_discovery.data import df_wbm
from matbench_discovery.enums import MbdKey, TestSubset
from matbench_discovery.metrics import discovery

if __name__ == "__main__":
    for model in cli_args.models:
        if not model.is_complete:
            print(f"\nSkipping {model.label}: incomplete discovery metrics")
    cli_args.models = [model for model in cli_args.models if model.is_complete]
    if not cli_args.models:
        raise SystemExit(0)

    from matbench_discovery.preds.discovery import df_each_pred, df_preds

    uniq_protos_idx = df_wbm.query(MbdKey.uniq_proto).index
    # dummy discovery rate of stable crystals when selecting randomly from the unique
    # prototype subset, used to compute the discovery acceleration factor (DAF)
    uniq_proto_prevalence = (
        df_wbm.loc[uniq_protos_idx, MbdKey.each_true] <= STABILITY_THRESHOLD
    ).mean()

    for model in cli_args.models:
        try:
            print(f"\nProcessing {model.label}...")
            model_preds = df_preds[model.label]
            each_true = df_preds[MbdKey.each_true]
            each_pred = df_each_pred[model.label]
            each_pred_uniq_protos = each_pred.loc[uniq_protos_idx]
            # 10k most stable = lowest *predicted hull distance* (each_pred) within
            # the unique prototype subset, not lowest raw formation energy
            most_stable_10k = each_pred_uniq_protos.nsmallest(10_000)

            full_metrics = discovery.stable_metrics(each_true, each_pred, fillna=True)
            uniq_proto_metrics = discovery.stable_metrics(
                each_true.loc[uniq_protos_idx], each_pred_uniq_protos, fillna=True
            )
            stable_10k_metrics = discovery.stable_metrics(
                each_true.loc[most_stable_10k.index], most_stable_10k, fillna=True
            )
            # DAF on these subsets is relative to the uniq-proto prevalence
            for metrics in (uniq_proto_metrics, stable_10k_metrics):
                metrics[str(Key.daf.symbol)] = (
                    metrics["Precision"] / uniq_proto_prevalence
                )

            for test_subset, (metrics, subset_idx) in {
                TestSubset.full_test_set: (full_metrics, slice(None)),
                TestSubset.uniq_protos: (uniq_proto_metrics, uniq_protos_idx),
                TestSubset.most_stable_10k: (
                    stable_10k_metrics,
                    most_stable_10k.index,
                ),
            }.items():
                # cast all values to float (incl. TP/FP/TN/FN counts) to match the
                # number format of previous YAML writes
                rounded_metrics: dict[str, str | float] = {
                    key: round(float(value), 3) for key, value in metrics.items()
                }
                discovery.write_metrics_to_yaml(
                    model, rounded_metrics, model_preds.loc[subset_idx], test_subset
                )
                print(f"\tUpdated discovery metrics for {test_subset}")
        except (ValueError, OSError, KeyError):
            print(f"\tError processing {model.label}: {traceback.format_exc()}")
