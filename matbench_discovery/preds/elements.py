"""Per-element prediction-error reductions."""

import pandas as pd


def mean_abs_error_by_element(
    compositions: pd.DataFrame, errors: pd.Series
) -> pd.Series:
    """Average absolute errors over predictions available for each element."""
    presence = compositions.notna()
    return (
        presence.mul(errors.abs(), axis="index").sum()
        / presence.mul(errors.notna(), axis="index").sum()
    )
