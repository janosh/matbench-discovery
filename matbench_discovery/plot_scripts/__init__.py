from collections.abc import Sequence
from typing import Any

import pandas as pd
from tqdm import tqdm

from matbench_discovery import ROOT

df_wbm = pd.read_csv(f"{ROOT}/data/wbm/2022-10-19-wbm-summary.csv")
df_wbm.index = df_wbm.material_id


data_paths = {
    "Wren": "data/2022-06-11-from-rhys/wren-mp-initial-structures.csv",
    "CGCNN IS2RE": "models/cgcnn/2022-11-23-test-cgcnn-wbm-IS2RE/"
    "cgcnn-ensemble-preds.csv",
    # "CGCNN IS2RE": "data/2022-06-11-from-rhys/cgcnn-mp-initial-structures.csv",
    "CGCNN RS2RE": "data/2022-06-11-from-rhys/cgcnn-mp-cse.csv",
    "Voronoi IS2RE": "data/2022-06-11-from-rhys/voronoi-mp-initial-structures.csv",
    "Voronoi RS2RE": "data/2022-06-11-from-rhys/voronoi-mp-cse.csv",
    "Wrenformer": "models/wrenformer/2022-11-15-wrenformer-IS2RE-preds.csv",
    "MEGNet": "models/megnet/2022-11-18-megnet-wbm-IS2RE/megnet-e-form-preds.csv",
    "M3GNet": "models/m3gnet/2022-10-31-m3gnet-wbm-IS2RE.json.gz",
    "Bowsr MEGNet": "models/bowsr/2022-11-22-bowsr-megnet-wbm-IS2RE.json.gz",
}


def load_model_preds(
    models: Sequence[str], pbar: bool = True, id_col: str = "material_id"
) -> dict[str, pd.DataFrame]:

    if mismatch := set(models) - set(data_paths):
        raise ValueError(f"Unknown models: {mismatch}")

    dfs: dict[str, pd.DataFrame] = {}

    for model_name in (bar := tqdm(models, disable=not pbar)):
        bar.set_description(model_name)
        data_path = data_paths[model_name]
        reader = pd.read_csv if ".csv" in data_path else pd.read_json
        df = reader(f"{ROOT}/{data_path}").set_index(id_col)
        dfs[model_name] = df

    return dfs


def load_df_wbm_with_preds(**kwargs: Any) -> pd.DataFrame:
    dfs = load_model_preds(**kwargs)
    df_out = df_wbm.copy()
    for model_name, df in dfs.items():
        if f"e_form_per_atom_{model_name.lower()}" in df:
            df_out[model_name] = df[f"e_form_per_atom_{model_name.lower()}"]
        elif len(pred_cols := df.filter(like="_pred_ens").columns) > 0:
            assert len(pred_cols) == 1
            df_out[model_name] = df[pred_cols[0]]
        elif len(pred_cols := df.filter(like=r"_pred_").columns) > 1:
            # make sure we average the expected number of ensemble member predictions
            assert len(pred_cols) == 10, f"{len(pred_cols) = }, expected 10"
            df_out[model_name] = df[pred_cols].mean(axis=1)
        elif "e_form_pred" in df:  # voronoi
            df_out[model_name] = df.e_form_pred
        else:
            raise ValueError(
                f"No pred col for {model_name=}, available cols={list(df)}"
            )
    return df_out
