"""Analyze structures and composition with largest mean error across all models.
Maybe there's some chemistry/region of materials space that all models struggle with?
Might point to deficiencies in the data or models architecture.
"""


# %%
import pandas as pd
from pymatgen.core import Composition
from pymatviz import ptable_heatmap_plotly
from tqdm import tqdm

from matbench_discovery import MODELS, ROOT
from matbench_discovery.data import df_wbm
from matbench_discovery.preds import (
    df_each_err,
    df_metrics,
    df_preds,
    each_true_col,
    model_mean_err_col,
)

__author__ = "Janosh Riebesell"
__date__ = "2023-02-15"

df_each_err[model_mean_err_col] = df_preds[model_mean_err_col] = df_each_err.abs().mean(
    axis=1
)


# %% map average model error onto elements
frac_comp_col = "fractional composition"
df_wbm[frac_comp_col] = [
    Composition(comp).fractional_composition for comp in tqdm(df_wbm.formula)
]

df_frac_comp = pd.DataFrame(comp.as_dict() for comp in df_wbm[frac_comp_col]).set_index(
    df_wbm.index
)
assert all(
    df_frac_comp.sum(axis=1).round(6) == 1
), "composition fractions don't sum to 1"

# df_frac_comp = df_frac_comp.dropna(axis=1, thresh=100)  # remove Xe with only 1 entry


# %% compute number of samples per element in training set
# counting element occurrences not weighted by composition, assuming model don't learn
# much more about iron and oxygen from Fe2O3 than from FeO
df_elem_err = pd.read_json(
    f"{ROOT}/site/src/routes/about-the-data/mp-element-counts-occurrence.json",
    typ="series",
)
train_count_col = "MP Occurrences"
df_elem_err = df_elem_err.reset_index(name=train_count_col).set_index("index")


# %%
test_set_std_col = "Test set standard deviation"
df_elem_err[test_set_std_col] = (
    df_frac_comp.where(pd.isna, 1) * df_wbm[each_true_col].to_numpy()[:, None]
).std()


# %%
normalized = True
cs_range = (0, 0.5)  # same range for all plots
# cs_range = (None, None)  # different range for each plot
for model in (*df_metrics, model_mean_err_col):
    df_elem_err[model] = (
        df_frac_comp * df_each_err[model].abs().to_numpy()[:, None]
    ).mean()
    # don't change series values in place, would change the df
    per_elem_err = df_elem_err[model].copy(deep=True)
    per_elem_err.name = f"{model} (eV/atom)"
    if normalized:
        per_elem_err /= df_elem_err[test_set_std_col]
        per_elem_err.name = f"{model} (normalized by test set std)"
    fig = ptable_heatmap_plotly(
        per_elem_err, precision=".2f", colorscale="Inferno", cscale_range=cs_range
    )
    fig.show()


# %%
expected_cols = {
    *"ALIGNN, BOWSR, CGCNN, CGCNN+P, CHGNet, MACE, M3GNet, MEGNet, "
    "MP Occurrences, Mean error all models, Test set standard deviation, Voronoi RF, "
    "Wrenformer".split(", ")
}
assert {*df_elem_err} >= expected_cols
assert (df_elem_err.isna().sum() < 35).all()
df_elem_err.round(4).to_json(f"{MODELS}/per-element-each-errors.json")
