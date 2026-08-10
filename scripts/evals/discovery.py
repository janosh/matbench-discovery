"""Calculate discovery metrics for selected models and write them to model YAMLs.

The v1 paper's metric-table PDFs are frozen publication artifacts and intentionally
do not track the evolving online leaderboard.
"""

from matbench_discovery.cli import cli_args, complete_models
from matbench_discovery.enums import Model
from matbench_discovery.metrics import discovery
from scripts.evals import evaluate_models


def main() -> int:
    """Evaluate discovery metrics, returning 0 if any model succeeds and 1 otherwise."""
    models = complete_models()
    for model in cli_args.models:
        if model not in models:
            print(f"\nSkipping {model.label}: incomplete discovery metrics")
    if not models:  # nothing selected has predictions, which is not an error
        return 0

    # deferred: loading every model's WBM predictions costs several GB
    from matbench_discovery.data import load_discovery_predictions

    df_preds, _df_each_pred, _df_each_err = load_discovery_predictions()
    uniq_proto_prevalence = discovery.wbm_uniq_proto_prevalence()

    def evaluate_one(model: Model) -> None:
        print(f"\nProcessing {model.label}...")
        metric_reference, model_preds = discovery.prepare_model_predictions(
            df_preds, df_preds[model.label]
        )
        metrics_by_subset = discovery.calc_discovery_metrics(
            metric_reference, model_preds, uniq_proto_prevalence=uniq_proto_prevalence
        )
        discovery.write_all_metrics_to_yaml(
            model, metrics_by_subset, metric_reference, model_preds
        )
        for test_subset in metrics_by_subset:
            print(f"\tUpdated discovery metrics for {test_subset}")

    return evaluate_models("discovery", models, evaluate_one)


if __name__ == "__main__":
    raise SystemExit(main())
