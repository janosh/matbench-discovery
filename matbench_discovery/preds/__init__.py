"""Load discovery predictions and derive errors, composition, and element statistics."""

from collections.abc import Sequence

import pandas as pd
from pymatgen.core import Composition
from pymatviz.enums import Key

from matbench_discovery import ROOT
from matbench_discovery.data import load_df_wbm_with_preds
from matbench_discovery.enums import MbdKey, Model, TestSubset

TEST_SET_STD_COL = "Test set standard deviation"
TRAIN_COUNT_COL = "MP Occurrences"
MP_COUNTS_PATH = f"{ROOT}/site/src/routes/data/mp-element-counts-by-occurrence.json"


def load_prediction_errors(
    models: Sequence[Model], subset: TestSubset | None = None
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load rounded WBM predictions and per-model formation-energy errors.

    Formation-energy and energy-above-hull errors are identical because every model
    uses the same fixed DFT convex hull.
    """
    predictions = load_df_wbm_with_preds(models=models, subset=subset).round(3)
    model_labels = [model.label for model in models]
    errors = predictions[model_labels].sub(predictions[MbdKey.e_form_dft], axis="index")
    return predictions, errors


def derive_element_data(
    predictions: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return element amounts and benchmark occurrence/statistics columns."""
    compositions = pd.DataFrame(
        Composition(formula).as_dict() for formula in predictions[Key.formula]
    ).set_index(predictions.index)
    element_data = (
        pd.read_json(MP_COUNTS_PATH, typ="series").rename(TRAIN_COUNT_COL).to_frame()
    )
    element_data.index.name = "symbol"
    element_data[TEST_SET_STD_COL] = (
        compositions.where(pd.isna, 1)
        * predictions[MbdKey.each_true].to_numpy()[:, None]
    ).std()
    return compositions, element_data
