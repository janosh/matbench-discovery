"""Analyze structures and composition with largest mean error across all models.
Maybe there's some chemistry/region of materials space that all models struggle with?
Might point to deficiencies in the data or models architecture.
"""

# %%
import os

import pandas as pd
import pymatviz as pmv
from pymatgen.core import Composition
from pymatviz.enums import Key
from pymatviz.utils import si_fmt
from tqdm.auto import tqdm

from matbench_discovery import ROOT, SITE_DIR
from matbench_discovery.cli import cli_args
from matbench_discovery.data import load_df_wbm_with_preds
from matbench_discovery.enums import MbdKey, TestSubset

test_set_std_col = "Test set standard deviation"
elem_col, size_col = "Element", "group"


# %%
test_subset = globals().get("test_subset", TestSubset.uniq_protos)
models_to_plot = cli_args.models

# Load predictions for specified models
df_preds = load_df_wbm_with_preds(
    models=models_to_plot, subset=cli_args.test_subset
).round(3)

# Calculate errors for each model
df_each_err = pd.DataFrame(index=df_preds.index)
for model in models_to_plot:
    df_each_err[model.label] = df_preds[model.label] - df_preds[MbdKey.e_form_dft]

# project average model error onto periodic table
df_comp = pd.DataFrame(
    Composition(comp).as_dict() for comp in df_preds[Key.formula]
).set_index(df_preds.index)


# %% compute number of samples per element in training set
# we count raw element occurrence, i.e. not weighted by composition, based on the
# hypothesis that models don't learn more about iron and oxygen from Fe2O3 than from FeO
counts_path = f"{ROOT}/site/src/routes/data/mp-element-counts-by-occurrence.json"
df_elem_err = pd.read_json(counts_path, typ="series")
train_count_col = "MP Occurrences"
df_elem_err = df_elem_err.reset_index(name=train_count_col).set_index("index")
df_elem_err.index.name = "symbol"

# compute std dev of DFT hull dist for each element in test set
df_elem_err[test_set_std_col] = (
    df_comp.where(pd.isna, 1) * df_preds[MbdKey.each_true].to_numpy()[:, None]
).std()


# %%
normalized = True
cs_range = (0, 0.5)  # same range for all plots
# cs_range = (None, None)  # different range for each plot
for model in tqdm(models_to_plot, desc="Processing models"):
    df_elem_err[model.label] = (
        df_comp * df_each_err[model.label].abs().to_numpy()[:, None]
    ).mean()
    # don't change series values in place, would change the df
    per_elem_err = df_elem_err[model.label].copy(deep=True)
    per_elem_err.name = f"{model.label} (eV/atom)"
    if normalized:
        per_elem_err /= df_elem_err[test_set_std_col]
        per_elem_err.name = f"{model.label} (normalized by test set std)"
    fig = pmv.ptable_heatmap_plotly(
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

    fig_ptable_each_errors = pmv.ptable_hists_plotly(
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
    *[model.label for model in models_to_plot],
    train_count_col,
    test_set_std_col,
}
if missing_cols := expected_cols - {*df_elem_err}:
    raise ValueError(f"{missing_cols=} not in {df_elem_err.columns=}")
if any(df_elem_err.isna().sum() > 35):
    raise ValueError("Too many NaNs in df_elem_err")

# Merge with existing data to preserve other models when running with --models subset
json_path = f"{SITE_DIR}/routes/models/per-element-each-errors.json"
if os.path.isfile(json_path):
    df_existing = pd.read_json(json_path)
    # Update existing with new model columns (and refresh metadata cols)
    for col in df_elem_err.columns:
        df_existing[col] = df_elem_err[col]
    df_elem_err = df_existing

df_elem_err.round(4).to_json(json_path)
