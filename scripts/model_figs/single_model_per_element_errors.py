"""Generate per-element discovery-error payloads."""

# %%
import json

import pandas as pd

from matbench_discovery import SITE_DIR, figs, payload_numerics, preds
from matbench_discovery.cli import cli_args

models_to_plot = cli_args.models

df_predictions, df_each_err = preds.load_prediction_errors(
    models_to_plot, subset=cli_args.test_subset
)
df_comp, df_elem_err = preds.derive_element_data(df_predictions)
element_sources = {"prediction_error_loader": preds.__file__}
model_identities = {
    model.key: figs.discovery_model_identity(model) for model in models_to_plot
}


# %%
# Mean model error for structures containing each element against MP prevalence.
df_elem_present = df_comp.notna()
elem_present_counts = df_elem_present.sum()
element_prevalence_errors = {
    model.key: (
        df_elem_present.multiply(df_each_err[model.label].abs(), axis=0).sum()
        / elem_present_counts
    ).reindex(df_elem_err.index)
    for model in models_to_plot
}

figs.write_site_payload(
    "element-prevalence-vs-error",
    {
        "elements": [str(symbol) for symbol in df_elem_err.index],
        "occurrences": payload_numerics.round_list(df_elem_err[preds.TRAIN_COUNT_COL]),
        "models": [
            model_identities[model.key]
            | {"y": payload_numerics.round_list(element_prevalence_errors[model.key])}
            for model in models_to_plot
        ],
        **figs.build_discovery_payload_provenance(
            generator=__file__,
            test_subset=cli_args.test_subset.value,
            benchmark_inputs={"mp_element_occurrences": preds.MP_COUNTS_PATH},
            source_files=element_sources
            | {"payload_numerics": payload_numerics.__file__},
            parameters={
                "analysis": "element_prevalence",
                "coordinate_decimals": payload_numerics.COORD_DECIMALS,
            },
            packages=("pymatgen",),
        ),
    },
)


# %%
model_element_errors = pd.DataFrame(
    {
        model.key: (df_comp * df_each_err[model.label].abs().to_numpy()[:, None]).mean()
        for model in models_to_plot
    }
)
df_elem_err = pd.concat([df_elem_err, model_element_errors], axis="columns")


# %%
expected_cols = {preds.TRAIN_COUNT_COL, preds.TEST_SET_STD_COL, *model_identities}
if missing_cols := expected_cols - {*df_elem_err}:
    raise ValueError(f"{missing_cols=} not in {df_elem_err.columns=}")
if any(df_elem_err.isna().sum() > 35):
    raise ValueError("Too many NaNs in df_elem_err")

elem_err_models = [
    {
        **model_identities[model.key],
        "values": json.loads(df_elem_err[model.key].round(4).to_json()),
    }
    for model in models_to_plot
]
figs.write_site_payload(
    "per-element-each-errors",
    {
        **figs.build_discovery_payload_provenance(
            generator=__file__,
            test_subset=cli_args.test_subset.value,
            benchmark_inputs={"mp_element_occurrences": preds.MP_COUNTS_PATH},
            source_files=element_sources,
            parameters={"output_decimals": 4},
            packages=("pymatgen",),
        ),
        "mp_occurrences": json.loads(df_elem_err[preds.TRAIN_COUNT_COL].to_json()),
        "test_set_standard_deviation": json.loads(
            df_elem_err[preds.TEST_SET_STD_COL].to_json()
        ),
        "models": elem_err_models,
    },
    directory=f"{SITE_DIR}/routes/models",
)
