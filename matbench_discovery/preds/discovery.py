"""Centralize data-loading and computing metrics for plotting scripts."""

import pandas as pd

from matbench_discovery.data import load_df_wbm_with_preds
from matbench_discovery.enums import MbdKey, Model

__author__ = "Janosh Riebesell"
__date__ = "2023-02-04"


# load WBM summary dataframe with all models' formation energy predictions (eV/atom)
df_preds = load_df_wbm_with_preds().round(3)

# dataframe of all models' energy above convex hull (EACH) predictions (eV/atom)
df_each_pred = pd.DataFrame()
for model in Model:
    df_each_pred[model.label] = (
        df_preds[MbdKey.each_true] + df_preds[model.label] - df_preds[MbdKey.e_form_dft]
    )

"""
To avoid confusion for anyone reading this code, df_each_err calculates the formation
energy MAE but reports it as the MAE for the energy above the convex hull prediction.
The former is more easily calculated but the two quantities are the same. The formation
energy of a material is the difference in energy between a material and its
constituent elements in their standard states. The distance to the convex hull is
defined as the difference between a material's formation energy and the minimum
formation energy of all possible stable materials made from the same elements. Since
the formation energy of a material is used to calculate the distance to the convex
hull, the error of a formation energy prediction directly determines the error in the
distance to the convex hull prediction.

A further point of clarification: whenever we say convex hull distance we mean
the signed distance that is positive for thermodynamically unstable materials above
the hull and negative for stable materials below it.
"""


# dataframe of all model prediction errors for energy above convex hull (EACH) (eV/atom)
df_each_err = pd.DataFrame()
for model in Model:
    df_each_err[model.label] = df_preds[model.label] - df_preds[MbdKey.e_form_dft]

df_each_err[MbdKey.each_err_models] = df_preds[MbdKey.each_err_models] = (
    df_each_err.abs().mean(axis=1)
)
