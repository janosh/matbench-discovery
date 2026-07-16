"""Routines for loading/processing predictions of different modeling tasks."""

from collections.abc import Sequence

import pandas as pd
from pymatgen.core import Composition
from pymatviz.enums import Key

from matbench_discovery import ROOT
from matbench_discovery.data import load_df_wbm_with_preds
from matbench_discovery.enums import MbdKey, Model, TestSubset

test_set_std_col = "Test set standard deviation"
train_count_col = "MP Occurrences"


def load_per_element_errors(
    models: Sequence[Model], subset: TestSubset | None = None
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Load WBM predictions and derive per-element error analysis dataframes.

    Args:
        models (Sequence[Model]): Models whose predictions to load.
        subset (TestSubset | None): Test subset passed to load_df_wbm_with_preds.

    Returns:
        tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
            - df_preds: WBM summary with model predictions (rounded to 3 decimals)
            - df_each_err: per-model error in predicted energy above convex hull
              (eV/atom)
            - df_comp: element amounts per structure (NaN where element absent)
            - df_elem_err: per-element MP occurrence counts (train_count_col) and
              test-set hull-dist standard deviation (test_set_std_col)
    """
    df_preds = load_df_wbm_with_preds(models=models, subset=subset).round(3)

    # error in predicted energy above convex hull (EACH) for each model, computed as
    # the formation energy error since the two are identical: the convex hull is built
    # from fixed DFT reference energies, so each_pred - each_true = e_form_pred -
    # e_form_dft
    df_each_err = pd.DataFrame(index=df_preds.index)
    for model in models:
        df_each_err[model.label] = df_preds[model.label] - df_preds[MbdKey.e_form_dft]

    # element amounts per structure, used to project errors onto periodic table
    df_comp = pd.DataFrame(
        Composition(comp).as_dict() for comp in df_preds[Key.formula]
    ).set_index(df_preds.index)

    # number of samples per element in training set. we count raw element occurrence,
    # i.e. not weighted by composition, based on the hypothesis that models don't
    # learn more about iron and oxygen from Fe2O3 than from FeO
    counts_path = f"{ROOT}/site/src/routes/data/mp-element-counts-by-occurrence.json"
    df_elem_err = pd.read_json(counts_path, typ="series")
    df_elem_err = df_elem_err.reset_index(name=train_count_col).set_index("index")
    df_elem_err.index.name = "symbol"

    # compute std dev of DFT hull dist for each element in test set
    df_elem_err[test_set_std_col] = (
        df_comp.where(pd.isna, 1) * df_preds[MbdKey.each_true].to_numpy()[:, None]
    ).std()

    return df_preds, df_each_err, df_comp, df_elem_err
