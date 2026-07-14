"""Routines for calculating metrics for different modeling tasks."""

from collections.abc import Sequence

import pandas as pd

from matbench_discovery.enums import Model


def metrics_df_from_yaml(nested_keys: Sequence[str]) -> pd.DataFrame:
    """Extract metrics from Model enum metadata into a DataFrame.

    Returns a DataFrame with models as rows and metrics as columns. To calculate and
    write discovery metrics into the model YAML files in the first place, run
    python scripts/evals/discovery.py --models model1 model2 ...
    where the model names are Model enum values.
    """
    out_dict = {}
    for model in Model:
        try:
            if not model.is_active:
                continue
            combined_metrics: dict[object, object] = {}
            for nested_key in nested_keys:
                metrics: object = model.metadata.get("metrics", {})
                for sub_key in nested_key.split("."):
                    if not isinstance(metrics, dict):
                        break
                    metrics = metrics.get(sub_key) or {}
                if isinstance(metrics, dict):
                    combined_metrics = {**combined_metrics, **metrics}
            if combined_metrics:
                out_dict[model.label] = combined_metrics

        except Exception as exc:
            exc.add_note(f"{model.label=}")
            raise

    # Return DataFrame with models as rows instead of columns
    return pd.DataFrame.from_dict(out_dict, orient="index")
