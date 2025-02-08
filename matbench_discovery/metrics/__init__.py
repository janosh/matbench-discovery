"""Routines for calculating metrics for different modeling tasks."""

from collections.abc import Sequence

import pandas as pd

from matbench_discovery.models import MODEL_METADATA


def metrics_df_from_yaml(nested_keys: Sequence[str]) -> pd.DataFrame:
    """Extract metrics from MODEL_METADATA into a DataFrame. MODEL_METADATA in turn
    reads directly from the model YAML files.

    Returns a DataFrame with models as rows and metrics as columns. To calculate and
    write discovery metrics into the model YAML files in the first place, run
    python matbench_discovery/preds/discovery.py model1 model2 ...
    where the model names are Model enum values.
    """
    out_dict = {}
    for model_name, metadata in MODEL_METADATA.items():
        metrics = None
        try:
            combined_metrics: dict[str, float] = {}
            for nested_key in nested_keys:
                metrics = metadata.get("metrics", {})
                for sub_key in nested_key.split("."):
                    if not isinstance(metrics, dict):
                        break  # go straight to next model if sub_key returned non-dict
                    metrics = metrics.get(sub_key, {})
                if isinstance(metrics, dict):
                    combined_metrics |= metrics
            if combined_metrics:
                out_dict[model_name] = combined_metrics

        except Exception as exc:
            exc.add_note(f"{model_name=} with {metrics=}")
            raise

    # Return DataFrame with models as rows instead of columns
    return pd.DataFrame.from_dict(out_dict, orient="index")
