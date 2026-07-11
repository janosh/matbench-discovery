"""Generate rolling MAE and hull-distance density payloads for all models."""

# %%
import numpy as np

from matbench_discovery import figs
from matbench_discovery.cli import complete_models
from matbench_discovery.enums import MbdKey, TestSubset
from matbench_discovery.preds.discovery import df_each_pred, df_preds

__author__ = "Rhys Goodall, Janosh Riebesell"
__date__ = "2022-06-18"


# %%
test_subset = globals().get("test_subset", TestSubset.uniq_protos)

if test_subset == TestSubset.uniq_protos:
    df_preds = df_preds.query(MbdKey.uniq_proto)
    df_each_pred = df_each_pred.loc[df_preds.index]

window = 0.04
rolling_x = np.arange(-0.2, 0.2, 0.005)
rolling_models: list[dict[str, object]] = []
for model in complete_models():
    each_pred = df_each_pred[model.label].dropna()
    each_true = df_preds[MbdKey.each_true].loc[each_pred.index]
    abs_error = (each_pred - each_true).abs()
    rolling_mae = [
        abs_error[
            (each_true <= bin_center + window / 2)
            & (each_true > bin_center - window / 2)
        ].mean()
        for bin_center in rolling_x
    ]
    rolling_models.append(
        {
            "key": model.key,
            "label": model.label,
            "y": figs.round_list(rolling_mae),
        }
    )

# add negligible noise to prevent strange binning artifacts in the marginal plot
small_noise = np.random.default_rng(seed=0).random(len(df_preds)) * 1e-12
counts, bins = np.histogram(
    df_preds[MbdKey.each_true] + small_noise,
    bins=200,  # match the histogram clf plots.
    range=(-0.7, 0.7),
)
figs.write_site_payload(
    "rolling-mae-vs-hull-dist",
    {
        "x": figs.round_list(rolling_x),
        "models": rolling_models,
        # rolling count of test-set structures per hull-dist bin (drawn on y2)
        "density": {
            "x": figs.round_list((bins[:-1] + bins[1:]) / 2),
            "y": counts.tolist(),
        },
    },
)
