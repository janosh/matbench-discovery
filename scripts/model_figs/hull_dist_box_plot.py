"""Generate per-model hull-distance error box statistics for the site."""

# %%
from matbench_discovery import figs, payload_numerics
from matbench_discovery.cli import cli_args, complete_models
from matbench_discovery.data import load_discovery_predictions
from matbench_discovery.enums import MbdKey, TestSubset

df_preds, _df_each_pred, df_each_err = load_discovery_predictions()
test_subset = cli_args.test_subset
if test_subset == TestSubset.uniq_protos:
    df_preds = df_preds.query(MbdKey.uniq_proto)
    df_each_err = df_each_err.loc[df_preds.index]


# %%
box_models: list[dict[str, object]] = [
    {
        **figs.discovery_model_identity(model),
        "quantiles": payload_numerics.round_list(
            df_each_err[model.key].quantile((0.05, 0.25, 0.5, 0.75, 0.95))
        ),
    }
    for model in complete_models()
]


# %%
provenance = figs.build_discovery_payload_provenance(
    generator=__file__,
    test_subset=test_subset.value,
    source_files={"payload_numerics": payload_numerics.__file__},
    parameters={
        "coordinate_decimals": payload_numerics.COORD_DECIMALS,
        "quantiles": [0.05, 0.25, 0.5, 0.75, 0.95],
    },
)
figs.write_site_payload("box-hull-dist-errors", {**provenance, "models": box_models})
