"""Generate rolling MAE and hull-distance density payloads for all models."""

# %%
import numpy as np

from matbench_discovery import figs, payload_numerics
from matbench_discovery.cli import cli_args, complete_models
from matbench_discovery.data import load_discovery_predictions
from matbench_discovery.enums import MbdKey, TestSubset

df_preds, df_each_pred, _df_each_err = load_discovery_predictions()
test_subset = cli_args.test_subset
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
            **figs.discovery_model_identity(model),
            "y": payload_numerics.round_list(rolling_mae),
        }
    )

# add negligible noise to prevent strange binning artifacts in the marginal plot
small_noise = np.random.default_rng(seed=0).random(len(df_preds)) * 1e-12
counts, bins = np.histogram(
    df_preds[MbdKey.each_true] + small_noise,
    bins=200,  # match the histogram clf plots.
    range=(-0.7, 0.7),
)
# rolling count of test-set structures per hull-dist bin (drawn on y2)
density = {
    "x": payload_numerics.round_list((bins[:-1] + bins[1:]) / 2),
    "y": counts.tolist(),
}
provenance = figs.build_discovery_payload_provenance(
    generator=__file__,
    test_subset=test_subset.value,
    source_files={"payload_numerics": payload_numerics.__file__},
    parameters={
        "coordinate_decimals": payload_numerics.COORD_DECIMALS,
        "density_bins": 200,
        "density_range": [-0.7, 0.7],
        "noise_magnitude": 1e-12,
        "noise_seed": 0,
        "rolling_start": -0.2,
        "rolling_step": 0.005,
        "rolling_stop": 0.2,
        "rolling_window": window,
    },
)
figs.write_site_payload(
    "rolling-mae-vs-hull-dist",
    {
        **provenance,
        "x": payload_numerics.round_list(rolling_x),
        "models": rolling_models,
        "density": density,
    },
)
