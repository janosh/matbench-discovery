"""Centralize data-loading and computing metrics for plotting scripts."""

import pandas as pd
from pymatviz.enums import Key

from matbench_discovery import STABILITY_THRESHOLD
from matbench_discovery.data import Model, df_wbm, load_df_wbm_with_preds
from matbench_discovery.enums import MbdKey
from matbench_discovery.metrics.discovery import stable_metrics

__author__ = "Janosh Riebesell"
__date__ = "2023-02-04"


# load WBM summary dataframe with all models' formation energy predictions (eV/atom)
df_preds = load_df_wbm_with_preds().round(3)


df_metrics = pd.DataFrame()
df_metrics_10k = pd.DataFrame()  # look only at each model's 10k most stable predictions
df_metrics_uniq_protos = pd.DataFrame(index=df_metrics.index)

for df, title in (
    (df_metrics, "Metrics for Full Test Set"),
    (df_metrics_10k, "Metrics for 10k Most Stable Predictions"),
    (df_metrics_uniq_protos, "Metrics for unique non-MP prototypes"),
):
    df.attrs["title"] = title
    df.index.name = "model"

full_prevalence = (df_wbm[MbdKey.each_true] <= STABILITY_THRESHOLD).mean()
uniq_proto_prevalence = (
    df_wbm.query(MbdKey.uniq_proto)[MbdKey.each_true] <= STABILITY_THRESHOLD
).mean()

for model in Model:
    if not model.is_complete:
        continue

    each_pred = (
        df_preds[MbdKey.each_true] + df_preds[model.label] - df_preds[MbdKey.e_form_dft]
    )
    df_metrics[model.label] = stable_metrics(
        df_preds[MbdKey.each_true], each_pred, fillna=True
    )

    df_uniq_proto_preds = df_preds[df_wbm[MbdKey.uniq_proto]]

    each_pred_uniq_proto = (
        df_uniq_proto_preds[MbdKey.each_true]
        + df_uniq_proto_preds[model.label]
        - df_uniq_proto_preds[MbdKey.e_form_dft]
    )
    df_metrics_uniq_protos[model.label] = stable_metrics(
        df_uniq_proto_preds[MbdKey.each_true], each_pred_uniq_proto, fillna=True
    )
    df_metrics_uniq_protos.loc[Key.daf.symbol, model.label] = (
        df_metrics_uniq_protos[model.label]["Precision"] / uniq_proto_prevalence
    )

    # look only at each model's 10k most stable predictions in the unique prototype set
    most_stable_10k = each_pred_uniq_proto.nsmallest(10_000)
    df_metrics_10k[model.label] = stable_metrics(
        df_preds[MbdKey.each_true].loc[most_stable_10k.index],
        most_stable_10k,
        fillna=True,
    )
    df_metrics_10k.loc[Key.daf.symbol, model.label] = (
        df_metrics_10k[model.label]["Precision"] / uniq_proto_prevalence
    )


# pick F1 as primary metric to sort by
df_metrics = df_metrics.round(3).sort_values("F1", axis=1, ascending=False)
df_metrics_10k = df_metrics_10k.round(3).sort_values("F1", axis=1, ascending=False)
df_metrics_uniq_protos = df_metrics_uniq_protos.round(3).sort_values(
    "F1", axis=1, ascending=False
)

# dataframe of all models' energy above convex hull (EACH) predictions (eV/atom)
df_each_pred = pd.DataFrame()
for model in Model:
    if not model.is_complete:
        continue

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
    if not model.is_complete:
        continue

    df_each_err[model.label] = df_preds[model.label] - df_preds[MbdKey.e_form_dft]

df_each_err[MbdKey.each_err_models] = df_preds[MbdKey.each_err_models] = (
    df_each_err.abs().mean(axis=1)
)
