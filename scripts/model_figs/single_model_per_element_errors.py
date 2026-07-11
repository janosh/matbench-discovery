"""Analyze structures and composition with largest mean error across all models.
Maybe there's some chemistry/region of materials space that all models struggle with?
Might point to deficiencies in the data or models architecture.
"""

# %%
import json
from typing import cast

from tqdm.auto import tqdm

from matbench_discovery import SITE_DIR, figs
from matbench_discovery.cli import cli_args, is_full_model_run
from matbench_discovery.preds import (
    load_per_element_errors,
    test_set_std_col,
    train_count_col,
)

models_to_plot = cli_args.models

_df_preds, df_each_err, df_comp, df_elem_err = load_per_element_errors(
    models_to_plot, subset=cli_args.test_subset
)


# %%
for model in tqdm(models_to_plot, desc="Processing models"):
    df_elem_err[model.key] = (
        df_comp * df_each_err[model.label].abs().to_numpy()[:, None]
    ).mean()


# %%
expected_cols = {
    *[model.key for model in models_to_plot],
    train_count_col,
    test_set_std_col,
}
if missing_cols := expected_cols - {*df_elem_err}:
    raise ValueError(f"{missing_cols=} not in {df_elem_err.columns=}")
if any(df_elem_err.isna().sum() > 35):
    raise ValueError("Too many NaNs in df_elem_err")

# Each column (model key or metadata col) is one JSONL line, so concurrent
# submissions add different lines that git merges cleanly instead of rewriting one
# un-mergeable blob. Subset runs splice only their columns; full runs rewrite the
# roster. round(4) trims size; to_json maps NaN -> null; json.loads -> plain dict.
payload_path = f"{SITE_DIR}/routes/models/per-element-each-errors.jsonl"
elem_err_models = [
    {"key": str(col), "values": json.loads(cast("str", series.round(4).to_json()))}
    for col, series in df_elem_err.items()
]
figs.write_jsonl_payload(
    payload_path, {"models": elem_err_models}, full_run=is_full_model_run()
)
