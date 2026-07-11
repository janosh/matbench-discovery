"""Generate data payloads for the TMI (Too Much Information) pages.

Analyzes structures and compositions with largest mean error across all models.
Maybe there's some chemistry/region of materials space that all models struggle with?
Might point to deficiencies in the data or models architecture.
"""

# %%
import pandas as pd

from matbench_discovery import figs
from matbench_discovery.cli import cli_args
from matbench_discovery.data import df_wbm
from matbench_discovery.preds import load_per_element_errors, train_count_col

fp_diff_col = "site_stats_fingerprint_init_final_norm_diff"


# %%
models_to_plot = cli_args.models

_df_preds, df_each_err, df_comp, df_elem_err = load_per_element_errors(
    models_to_plot, subset=cli_args.test_subset
)


# %%
# Error by element against prevalence in the training set
# for checking correlation and R2 of elemental prevalence in MP training data vs.
# model error

# compute mean absolute error per element for each model
# weight by element presence (1 if element in structure, 0 otherwise)
df_elem_present = df_comp.where(pd.isna, 1).fillna(0)
elem_present_counts = df_elem_present.sum()
for model in models_to_plot:
    # for each element, compute mean error across structures containing it
    df_elem_err[model.label] = (
        df_elem_present.multiply(df_each_err[model.label].abs(), axis=0).sum()
        / elem_present_counts
    )

figs.write_site_payload(
    "element-prevalence-vs-error",
    {
        # element symbols + x (occurrence count per element) shared across models
        "elements": [str(symbol) for symbol in df_elem_err.index],
        "occurrences": figs.round_list(df_elem_err[train_count_col]),
        "models": [
            {"label": model.label, "y": figs.round_list(df_elem_err[model.label])}
            for model in models_to_plot
        ],
    },
    id_field="label",
)


# %% --- Fingerprint-based analysis ---
# Analyze correlation between relaxation change (measured by SiteStatsFingerprint diff)
# and model errors

model_labels = [model.label for model in models_to_plot]


# %% histogram of FP diff for structures with largest/smallest errors
n_structs = 1000
hist_largest_models: list[dict[str, object]] = []
for model in model_labels:
    large_errors = df_each_err[model].abs().nlargest(n_structs)
    small_errors = df_each_err[model].abs().nsmallest(n_structs)
    hist_entry: dict[str, object] = {"label": model}
    for label, errors in (("min", small_errors), ("max", large_errors)):
        fp_diff_values = df_wbm.loc[errors.index][fp_diff_col].to_numpy()
        hist_entry[f"err_{label}"] = figs.histogram(fp_diff_values, bins=100)
    hist_largest_models.append(hist_entry)

figs.write_site_payload(
    "hist-largest-each-errors-fp-diff",
    {"models": hist_largest_models},
    id_field="label",
)


# %% FP diff vs error for highest-error structures
n_structs = 100
each_errors_models: list[dict[str, object]] = []
for model in model_labels:
    errors = df_each_err[model].abs().nlargest(n_structs)
    model_mae = errors.mean()
    each_errors_models.append(
        {
            "label": model,
            "mae": round(float(model_mae), 4),
            "x": figs.round_list(df_wbm.loc[errors.index][fp_diff_col].values),
            "y": figs.round_list(errors.values),
        }
    )

figs.write_site_payload(
    "scatter-largest-each-errors-fp-diff",
    {"models": each_errors_models},
    id_field="label",
)


# %% Errors for structures with largest FP diff (most relaxation change)
n_points = 1000
# filter to only materials in the predictions subset
df_largest_fp_diff = df_wbm.loc[
    df_wbm.index.intersection(df_each_err.index), fp_diff_col
].nlargest(n_points)

fp_diff_models: list[dict[str, object]] = []
for model in model_labels:
    abs_errors = df_each_err[model].loc[df_largest_fp_diff.index].abs()
    model_mae = abs_errors.mean()
    fp_diff_models.append(
        {
            "label": model,
            "mae": round(float(model_mae), 4),
            "y": figs.round_list(abs_errors.values),
        }
    )

figs.write_site_payload(
    "scatter-largest-fp-diff-each-error",
    {
        "fp_diff": figs.round_list(df_largest_fp_diff.values),
        "models": fp_diff_models,
    },
    id_field="label",
)
