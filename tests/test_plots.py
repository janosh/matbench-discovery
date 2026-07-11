"""Tests for the retained plotting helper functions."""

import pandas as pd
import pytest
import wandb

from matbench_discovery.plots import stable_screening_sort, wandb_scatter


def test_stable_screening_sort_breaks_ties_by_material_id() -> None:
    """Equal predictions sort deterministically by material ID."""
    predictions = pd.Series([0.2, 0.1, 0.1], index=["mat-c", "mat-b", "mat-a"])
    sorted_predictions = stable_screening_sort(predictions)
    assert list(sorted_predictions.index) == ["mat-a", "mat-b", "mat-c"]


@pytest.mark.parametrize(
    ("fields", "expected_labels"),
    [
        ({"x": "col_x"}, None),
        ({"y": "col_y"}, None),
        ({"x": "col_x", "z": "col_z"}, None),
        ({"x": "col_x", "y": "col_y"}, {}),
        (
            {"x": "e_form_true", "y": "e_form_pred"},
            {
                "x_label": "DFT formation energy (eV/atom)",
                "y_label": "Predicted formation energy (eV/atom)",
            },
        ),
    ],
)
def test_wandb_scatter(
    fields: dict[str, str],
    expected_labels: dict[str, str] | None,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """WandB scatter validates fields and supplies formation-energy labels."""
    if expected_labels is None:
        with pytest.raises(ValueError, match="must specify x=str and y=str"):
            wandb_scatter(object(), fields)  # ty: ignore[invalid-argument-type]
        return

    plot_kwargs: dict[str, object] = {}
    logged_payload: dict[str, object] = {}

    def plot_table(**kwargs: object) -> object:
        """Capture WandB plot arguments."""
        plot_kwargs.update(kwargs)
        return object()

    def log(payload: dict[str, object]) -> None:
        """Capture the logged scatter plot."""
        logged_payload.update(payload)

    monkeypatch.setattr(wandb, "plot_table", plot_table)
    monkeypatch.setattr(wandb, "log", log)
    wandb_scatter(object(), fields)  # ty: ignore[invalid-argument-type]

    assert plot_kwargs["string_fields"] == expected_labels
    assert set(logged_payload) == {"true_pred_scatter"}
