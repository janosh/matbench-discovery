# %%
import os
from datetime import datetime

import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.impute import SimpleImputer
from sklearn.pipeline import Pipeline

from matbench_discovery.plot_scripts import df_wbm
from models.voronoi import featurizer

today = f"{datetime.now():%Y-%m-%d}"
module_dir = os.path.dirname(__file__)


# %%
train_path = f"{module_dir}/2022-11-18-voronoi-features-wbm.csv.bz2"
print(f"{train_path=}")
df_train = pd.read_csv(train_path)
print(f"{df_train.shape=}")

target_col = "e_form_per_atom_mp2020_corrected"
df_train[target_col] = df_wbm[target_col]

feature_names = featurizer.feature_labels()


# %%
failed_train = df_train[feature_names].isna().any(axis=1)
print(f"{failed_train:,} / {len(df_train):,} = {failed_train/len(df_train):.1%}")


# %%
model = Pipeline(
    [
        ("imputer", SimpleImputer()),  # For the failed structures
        ("model", RandomForestRegressor(n_estimators=150, n_jobs=-1, verbose=1)),
    ]
)


# %%
model.fit(df_train[feature_names], df_train[target_col])


# %%
task_type = "IS2RE"
print(f"{task_type=}")
test_path = f"{module_dir}/2022-11-18-voronoi-features-wbm.csv.bz2"
print(f"{test_path=}")
df_test = pd.read_csv(test_path)
print(f"{df_test.shape=}")


print("Train data shape:", df_test[feature_names].shape)

failed_test = df_test[feature_names].isna().any(axis=1)
print(f"{failed_test:,} / {len(df_train):,} = {failed_test/len(df_train):.1%}")


pred_col = "e_form_per_atom_voronoi_rf"
df_test[pred_col] = model.predict(df_test)


print(f"MAE = {(df_test[target_col] - df_test[pred_col]).abs().mean():=.3}")

cols = ["material_id", target_col, pred_col]

df_test[cols].to_csv(f"{module_dir}/voronoi-rf-e-form-preds-{task_type}.csv")


# %%
task_type = "RS2RE"
print(f"{task_type=}")
data_path = f"{module_dir}/wbm-{task_type}-vt.json.gz"
print(f"{data_path=}")
df_test = pd.read_csv(data_path)
print(f"{df_test.shape=}")

failed_test = df_test[feature_names].isna().any(axis=1)
print(f"{failed_test:,} / {len(df_train):,} = {failed_test/len(df_train):.1%}")

df_test[pred_col] = model.predict(df_test)

print(f"MAE = {(df_test[target_col] - df_test[pred_col]).abs().mean():=.3}")

df_test[cols].to_csv(f"{module_dir}/voronoi-rf-e-form-preds-{task_type}.csv")
