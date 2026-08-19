"""Generate data payloads for the TMI (Too Much Information) pages.

Analyzes structures and compositions with largest mean error across all models.
Maybe there's some chemistry/region of materials space that all models struggle with?
Might point to deficiencies in the data or models architecture.
"""

# %%
from matbench_discovery import figs, payload_numerics, preds
from matbench_discovery.cli import cli_args
from matbench_discovery.data import df_wbm

fp_diff_col = "site_stats_fingerprint_init_final_norm_diff"


def payload_provenance(parameters: dict[str, object]) -> dict[str, object]:
    """Build shared provenance for one TMI payload."""
    return figs.build_discovery_payload_provenance(
        generator=__file__,
        test_subset=cli_args.test_subset.value,
        source_files={
            "payload_numerics": payload_numerics.__file__,
            "prediction_error_loader": preds.__file__,
        },
        parameters={
            "coordinate_decimals": payload_numerics.COORD_DECIMALS,
        }
        | parameters,
    )


# %%
models_to_plot = cli_args.models
model_identities = {
    model.key: figs.discovery_model_identity(model) for model in models_to_plot
}

_predictions, df_each_err = preds.load_prediction_errors(
    models_to_plot, subset=cli_args.test_subset
)


# %% --- Fingerprint-based analysis ---
# Analyze correlation between relaxation change (measured by SiteStatsFingerprint diff)
# and model errors


# %% histogram of FP diff for structures with largest/smallest errors
n_structs = 1000
hist_largest_models: list[dict[str, object]] = []
for model in models_to_plot:
    large_errors = df_each_err[model.key].abs().nlargest(n_structs)
    small_errors = df_each_err[model.key].abs().nsmallest(n_structs)
    hist_entry = dict(model_identities[model.key])
    for label, errors in (("min", small_errors), ("max", large_errors)):
        fp_diff_values = df_wbm.loc[errors.index][fp_diff_col].to_numpy()
        hist_entry[f"err_{label}"] = payload_numerics.histogram(
            fp_diff_values, bins=100
        )
    hist_largest_models.append(hist_entry)

figs.write_site_payload(
    "hist-largest-each-errors-fp-diff",
    {
        **payload_provenance(
            {
                "histogram_bins": 100,
                "histogram_range": "per_model_data_min_max",
                "selected_structures": n_structs,
            }
        ),
        "models": hist_largest_models,
    },
)


# %% FP diff vs error for highest-error structures
n_structs = 100
each_errors_models: list[dict[str, object]] = []
for model in models_to_plot:
    errors = df_each_err[model.key].abs().nlargest(n_structs)
    model_mae = errors.mean()
    each_errors_models.append(
        model_identities[model.key]
        | {
            "mae": round(float(model_mae), 4),
            "x": payload_numerics.round_list(
                df_wbm.loc[errors.index][fp_diff_col].values
            ),
            "y": payload_numerics.round_list(errors.values),
        }
    )

figs.write_site_payload(
    "scatter-largest-each-errors-fp-diff",
    {
        **payload_provenance({"mae_decimals": 4, "selected_structures": n_structs}),
        "models": each_errors_models,
    },
)


# %% Errors for structures with largest FP diff (most relaxation change)
n_points = 1000
# filter to only materials in the predictions subset
df_largest_fp_diff = df_wbm.loc[
    df_wbm.index.intersection(df_each_err.index), fp_diff_col
].nlargest(n_points)

fp_diff_models: list[dict[str, object]] = []
for model in models_to_plot:
    abs_errors = df_each_err[model.key].loc[df_largest_fp_diff.index].abs()
    model_mae = abs_errors.mean()
    fp_diff_models.append(
        model_identities[model.key]
        | {
            "mae": round(float(model_mae), 4),
            "y": payload_numerics.round_list(abs_errors.values),
        }
    )

figs.write_site_payload(
    "scatter-largest-fp-diff-each-error",
    {
        **payload_provenance({"mae_decimals": 4, "selected_structures": n_points}),
        "fp_diff": payload_numerics.round_list(df_largest_fp_diff.values),
        "models": fp_diff_models,
    },
)
