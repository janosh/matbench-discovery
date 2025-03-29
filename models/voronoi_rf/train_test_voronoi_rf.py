"""Train and test a Voronoi RandomForestRegressor model."""

# %%
import os
import sys
from importlib.metadata import version

import joblib
import pandas as pd
import wandb
from pymatviz.enums import Key
from sklearn.ensemble import RandomForestRegressor
from sklearn.impute import SimpleImputer
from sklearn.metrics import r2_score
from sklearn.pipeline import Pipeline

from matbench_discovery import ROOT, today
from matbench_discovery.data import df_wbm, glob_to_df
from matbench_discovery.enums import DataFiles, MbdKey, Task
from matbench_discovery.hpc import slurm_submit
from matbench_discovery.plots import wandb_scatter

sys.path.append(f"{ROOT}/models")

from voronoi_rf import featurizer

__author__ = "Janosh Riebesell"
__date__ = "2022-11-26"


# %%
module_dir = os.path.dirname(__file__)
task_type = Task.IS2RE
print(f"{task_type=}")

out_dir = f"{module_dir}/{today}-train-test"
os.makedirs(out_dir, exist_ok=True)
out_path = f"{out_dir}/e-form-preds-{task_type}.csv.gz"
model_path = f"{out_dir}/voronoi-rf-model.joblib"  # Path to save the model
if os.path.isfile(out_path):
    raise SystemExit(f"{out_path=} already exists, exciting early")

job_name = "train-test-voronoi-rf"

slurm_vars = slurm_submit(
    job_name=job_name,
    out_dir=out_dir,
    account="matgen",
    time="6:0:0",
)


# %%
train_path = f"{module_dir}/2022-11-25-features-mp/voronoi-features-mp-*.csv.bz2"
df_train = glob_to_df(train_path).set_index(Key.mat_id)
print(f"{df_train.shape=}")

df_mp = pd.read_csv(DataFiles.mp_energies.path, na_filter=False).set_index(Key.mat_id)

test_path = f"{module_dir}/2022-11-18-features-wbm-{task_type}.csv.bz2"
df_test = pd.read_csv(test_path).set_index(Key.mat_id)
print(f"{df_test.shape=}")


for df, df_targets, e_form_col in (
    (df_train, df_mp, Key.formation_energy_per_atom),
    (df_test, df_wbm, MbdKey.e_form_dft),
):
    df[Key.form_energy] = df_targets[e_form_col]
    e_form_col = str(e_form_col)
    n_nans = df_targets[e_form_col].isna().sum()
    if n_nans != 0:
        raise ValueError(f"{e_form_col=} has {n_nans} NaNs")

model_name = "Voronoi RandomForestRegressor"
n_estimators = 150
random_forest = RandomForestRegressor(n_estimators=n_estimators, n_jobs=-1, verbose=1)


run_params = dict(
    train_path=train_path,
    test_path=test_path,
    mp_energies_path=DataFiles.mp_energies.path,
    versions={dep: version(dep) for dep in ("scikit-learn", "matminer", "numpy")},
    model_name=model_name,
    train_target_col=Key.form_energy,
    test_target_col=MbdKey.e_form_dft,
    df_train=dict(shape=str(df_train.shape)),
    df_test=dict(shape=str(df_test.shape)),
    slurm_vars=slurm_vars,
    n_estimators=n_estimators,
)

wandb.init(project="matbench-discovery", name=job_name, config=run_params)


# %%
feature_names = featurizer.feature_labels()
n_nans = df_train[feature_names].isna().any(axis=1).sum()

print(f"train set NaNs: {n_nans:,} / {len(df_train):,} = {n_nans / len(df_train):.2%}")

df_train = df_train.dropna(subset=[*feature_names, Key.form_energy])


# %% SimpleImputer for the failed structures
model = Pipeline([("imputer", SimpleImputer()), ("model", random_forest)])


# %% model needs to be fit to create random_forest.estimators_
model.fit(df_train[feature_names], df_train[Key.form_energy])

# counting only the threshold values for node splitting, see https://stackoverflow.com/a/75995417
n_trainable_params = sum(tree.tree_.node_count for tree in random_forest.estimators_)
print(f"{n_trainable_params=:,}")
run_params["n_trainable_params"] = n_trainable_params

# Save the trained model to disk
print(f"Saving model to {model_path}")
joblib.dump(model, model_path)
model_desc = f"Voronoi RandomForestRegressor model trained for {task_type} task"
artifact = wandb.Artifact(
    name=f"voronoi_rf_model_{task_type}", type="model", description=model_desc
)
artifact.add_file(model_path)
wandb.run.log_artifact(artifact)


# %%
n_nans = df_test[feature_names].isna().any(axis=1).sum()
print(f"test set NaNs: {n_nans:,} / {len(df_test):,} = {n_nans / len(df_test):.2%}")

df_test = df_test.dropna(subset=feature_names)

pred_col = "e_form_per_atom_voronoi_rf"
df_test[pred_col] = model.predict(df_test[feature_names])
# saving preds first to df_test, then df_wbm avoids length mismatch errors between
# output array and df_wbm
df_wbm[pred_col] = df_test[pred_col]

df_wbm[pred_col].round(4).to_csv(out_path)

table = wandb.Table(
    dataframe=df_wbm[[Key.formula, MbdKey.e_form_dft, pred_col]].reset_index()
)

df_wbm[pred_col].isna().sum()
MAE = (df_wbm[MbdKey.e_form_dft] - df_wbm[pred_col]).abs().mean()
R2 = r2_score(*df_wbm[[MbdKey.e_form_dft, pred_col]].dropna().to_numpy().T)
title = f"{model_name} {task_type} {MAE=:.3} {R2=:.3}"
print(title)

wandb_scatter(table, fields=dict(x=MbdKey.e_form_dft, y=pred_col), title=title)
