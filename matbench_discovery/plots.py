"""Small plotting-related helpers retained for model evaluation scripts."""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import pandas as pd
    import wandb


def stable_screening_sort(srs_pred: pd.Series) -> pd.Series:
    """Sort predictions by value with ties broken deterministically by material ID."""
    return srs_pred.sort_index(kind="stable").sort_values(kind="stable")


def wandb_scatter(table: wandb.Table, fields: dict[str, str], **kwargs: str) -> None:
    """Log a parity scatter plot using the project's custom WandB Vega spec.

    Args:
        table: WandB data table.
        fields: Map from table columns to Vega fields. Must specify ``x`` and ``y``.
        **kwargs: String fields passed to ``wandb.plot_table``.

    Raises:
        ValueError: If ``fields`` does not contain both ``x`` and ``y``.
    """
    import wandb

    if not {"x", "y"} <= fields.keys():
        raise ValueError(f"{fields=} must specify x=str and y=str column names")

    if "form" in fields["x"] and "form" in fields["y"]:
        kwargs.setdefault("x_label", "DFT formation energy (eV/atom)")
        kwargs.setdefault("y_label", "Predicted formation energy (eV/atom)")

    scatter_plot = wandb.plot_table(
        vega_spec_name="janosh/scatter-parity",
        data_table=table,
        fields=fields,
        string_fields=kwargs,
    )
    wandb.log({"true_pred_scatter": scatter_plot})
