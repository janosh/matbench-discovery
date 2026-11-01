"""Calculate discovery metrics for selected models and write them to model YAMLs.

The v1 paper's metric-table PDFs are frozen publication artifacts and intentionally
do not track the evolving online leaderboard.
"""

from matbench_discovery.cli import cli_args, complete_models
from matbench_discovery.discovery import write_discovery_metrics
from matbench_discovery.enums import Model
from scripts.evals import evaluate_models


def main() -> int:
    """Evaluate complete discovery models and fail if any evaluation fails."""
    models = complete_models()
    for model in cli_args.models:
        if model not in models:
            print(f"\nSkipping {model.label}: incomplete discovery metrics")
    if not models:  # nothing selected has predictions, which is not an error
        return 0

    # deferred: loading every model's WBM predictions costs several GB
    from matbench_discovery.data import load_discovery_predictions

    df_preds, _df_each_pred, _df_each_err = load_discovery_predictions()

    def evaluate_one(model: Model) -> None:
        print(f"\nProcessing {model.label}...")
        written = write_discovery_metrics(model, df_preds, df_preds[model.key])
        for test_subset in written:
            print(f"\tUpdated discovery metrics for {test_subset}")

    return evaluate_models("discovery", models, evaluate_one)


if __name__ == "__main__":
    raise SystemExit(main())
