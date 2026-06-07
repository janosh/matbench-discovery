"""Analyze structures and composition with largest mean error across all models.
Maybe there's some chemistry/region of materials space that all models struggle with?
Might point to deficiencies in the data or models architecture.
"""

# %%
import os

import pandas as pd
import pymatviz as pmv
from pymatviz.utils import si_fmt
from tqdm.auto import tqdm

from matbench_discovery import SITE_DIR
from matbench_discovery.cli import cli_args
from matbench_discovery.enums import TestSubset
from matbench_discovery.preds import (
    load_per_element_errors,
    test_set_std_col,
    train_count_col,
)

elem_col, size_col = "Element", "group"


# %%
test_subset = globals().get("test_subset", TestSubset.uniq_protos)
models_to_plot = cli_args.models

df_preds, df_each_err, df_comp, df_elem_err = load_per_element_errors(
    models_to_plot, subset=cli_args.test_subset
)


# %%
normalized = True
cs_range = (0, 0.5)  # same range for all plots
# cs_range = (None, None)  # different range for each plot
for model in tqdm(models_to_plot, desc="Processing models"):
    df_elem_err[model.key] = (
        df_comp * df_each_err[model.label].abs().to_numpy()[:, None]
    ).mean()
    # don't change series values in place, would change the df
    per_elem_err = df_elem_err[model.key].copy(deep=True)
    per_elem_err.name = f"{model.label} (eV/atom)"
    if normalized:
        per_elem_err /= df_elem_err[test_set_std_col]
        per_elem_err.name = f"{model.label} (normalized by test set std)"
    fig = pmv.ptable_heatmap(
        per_elem_err, fmt=".2f", colorscale="Inferno", cscale_range=cs_range
    )
    fig.show()

    # plot histogram of model errors for each element
    model_name = model.label
    df_each_err_elems = df_comp * (df_each_err[model_name].to_numpy()[:, None])
    df_each_err_elems = df_each_err_elems.drop(columns="Xe")
    n_samples_per_elem = (len(df_comp) - df_each_err_elems.isna().sum()).map(
        lambda x: si_fmt(x, fmt=".0f")
    )

    fig_ptable_each_errors = pmv.ptable_hists(
        df_each_err_elems,
        annotations=n_samples_per_elem,
        bins=100,
        subplot_kwargs=dict(shared_yaxes=False),
        log=False,
        colorbar=dict(title=f"{model_name} convex hull distance errors (eV/atom)"),
        x_range=(-0.5, 0.5),
        colorscale="Viridis",
    )

    fig_ptable_each_errors.show()


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

# Merge with existing data to preserve other models when running with --models subset.
# .json.gz: pandas infers gzip from the extension (read + write, compression=infer)
json_path = f"{SITE_DIR}/routes/models/per-element-each-errors.json.gz"
if os.path.isfile(json_path):
    df_existing = pd.read_json(json_path)
    # Update existing with new model columns (and refresh metadata cols)
    for col in df_elem_err.columns:
        df_existing[col] = df_elem_err[col]
    df_elem_err = df_existing

df_elem_err.round(4).to_json(json_path)
