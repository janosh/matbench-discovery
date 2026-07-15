"""Centralize data-loading of model predictions for plotting scripts.

Note on df_each_err: it calculates the formation energy MAE but reports it as the MAE
for the energy above the convex hull (EACH) prediction. The former is more easily
calculated but the two quantities are identical: in Matbench Discovery, the convex hull
is fixed (constructed from DFT reference energies, not model predictions), so EACH and
formation energy differ only by a composition-dependent constant. Linear
transformations leave the MAE invariant, hence
MAE(e_form_pred - e_form_dft) = MAE(each_pred - each_true). If the hull were instead
built from model predictions, systematic model errors could partially cancel out and
the two MAEs would differ. Whenever we say convex hull distance we mean the signed
distance that is positive for thermodynamically unstable materials above the hull and
negative for stable materials below it.
"""

from matbench_discovery.cli import complete_models
from matbench_discovery.data import load_df_wbm_with_preds
from matbench_discovery.enums import MbdKey

__author__ = "Janosh Riebesell"
__date__ = "2023-02-04"


# load WBM summary dataframe with all models' formation energy predictions (eV/atom)
models_to_load = tuple(complete_models())
df_preds = load_df_wbm_with_preds(models=models_to_load).round(3)

model_labels = [model.label for model in models_to_load]
each_shift = df_preds[MbdKey.each_true] - df_preds[MbdKey.e_form_dft]
df_each_pred = df_preds[model_labels].add(each_shift, axis=0)
df_each_err = df_preds[model_labels].sub(df_preds[MbdKey.e_form_dft], axis=0)
df_each_err[MbdKey.each_err_models] = df_preds[MbdKey.each_err_models] = (
    df_each_err.abs().mean(axis=1)
)
