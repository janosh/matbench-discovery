"""Calculate discovery metrics for selected models and write them to model YAMLs.

The v1 paper's metric-table PDFs are frozen publication artifacts and intentionally
do not track the evolving online leaderboard.
"""

import traceback

from matbench_discovery.cli import cli_args
from matbench_discovery.metrics import discovery

if __name__ == "__main__":
    for model in cli_args.models:
        if not model.is_complete:
            print(f"\nSkipping {model.label}: incomplete discovery metrics")
    cli_args.models = [model for model in cli_args.models if model.is_complete]
    if not cli_args.models:
        raise SystemExit(0)

    from matbench_discovery.preds.discovery import df_preds

    uniq_proto_prevalence = discovery.wbm_uniq_proto_prevalence()
    for model in cli_args.models:
        try:
            print(f"\nProcessing {model.label}...")
            model_preds = df_preds[model.label]
            subset_indices = discovery.discovery_subset_indices(df_preds, model_preds)
            metrics_by_subset = discovery.calc_discovery_metrics(
                df_preds,
                model_preds,
                subset_indices=subset_indices,
                uniq_proto_prevalence=uniq_proto_prevalence,
            )
            discovery.write_all_metrics_to_yaml(
                model,
                metrics_by_subset,
                df_preds,
                model_preds,
                subset_indices=subset_indices,
            )
            for test_subset in metrics_by_subset:
                print(f"\tUpdated discovery metrics for {test_subset}")
        except (ValueError, OSError, KeyError):
            print(f"\tError processing {model.label}: {traceback.format_exc()}")
