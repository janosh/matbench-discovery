from collections.abc import Sequence
from glob import glob
from typing import Any, Callable

import pandas as pd
from tqdm import tqdm

from matbench_discovery import ROOT

df_wbm = pd.read_csv(f"{ROOT}/data/wbm/2022-10-19-wbm-summary.csv")
df_wbm.index = df_wbm.material_id


data_paths = {
    "Wren": "data/2022-06-11-from-rhys/wren-mp-initial-structures.csv",
    "CGCNN": "models/cgcnn/2022-11-23-test-cgcnn-wbm-IS2RE/cgcnn-ensemble-preds.csv",
    "CGCNN IS2RE": "data/2022-06-11-from-rhys/cgcnn-mp-initial-structures.csv",
    "CGCNN RS2RE": "data/2022-06-11-from-rhys/cgcnn-mp-cse.csv",
    # "Voronoi IS2RE": "data/2022-06-11-from-rhys/voronoi-mp-initial-structures.csv",
    "Voronoi RF": "models/voronoi/2022-11-27-train-test/e-form-preds-IS2RE.csv",
    # "Voronoi RS2RE": "data/2022-06-11-from-rhys/voronoi-mp-cse.csv",
    "Wrenformer": "models/wrenformer/2022-11-15-wrenformer-IS2RE-preds.csv",
    "MEGNet": "models/megnet/2022-11-18-megnet-wbm-IS2RE/megnet-e-form-preds.csv",
    "M3GNet": "models/m3gnet/2022-10-31-m3gnet-wbm-IS2RE.csv",
    "BOWSR MEGNet": "models/bowsr/2022-11-22-bowsr-megnet-wbm-IS2RE.csv",
}


def glob_to_df(
    pattern: str, reader: Callable[[Any], pd.DataFrame] = None, pbar: bool = True
) -> pd.DataFrame:
    """Combine data files matching a glob pattern into a single dataframe.

    Args:
        pattern (str): Glob file pattern.
        reader (Callable[[Any], pd.DataFrame], optional): Function that loads data from
            disk. Defaults to pd.read_csv if ".csv" in pattern else pd.read_json.
        pbar (bool, optional): Whether to show progress bar. Defaults to True.

    Returns:
        pd.DataFrame: Combined dataframe.
    """
    reader = reader or pd.read_csv if ".csv" in pattern else pd.read_json

    # prefix pattern with ROOT if not absolute path
    files = glob(pattern if pattern.startswith("/") else f"{ROOT}/{pattern}")
    if len(files) == 0:
        raise FileNotFoundError(f"No files matching glob {pattern=}")

    sub_dfs = {}  # used to join slurm job array results into single df
    for file in tqdm(files, disable=not pbar):
        df = reader(file)
        sub_dfs[file] = df

    return pd.concat(sub_dfs.values())


def load_model_preds(
    models: Sequence[str], pbar: bool = True, id_col: str = "material_id"
) -> dict[str, pd.DataFrame]:
    """Load model predictions from disk into dictionary of dataframes.

    Args:
        models (Sequence[str]): Model names must be keys of data_paths dict.
        pbar (bool, optional): Whether to show progress bar. Defaults to True.
        id_col (str, optional): Column to set as df.index. Defaults to "material_id".

    Raises:
        ValueError: On unknown model names.

    Returns:
        dict[str, pd.DataFrame]: Dictionary of dataframes, one for each model.
    """
    if mismatch := ", ".join(set(models) - set(data_paths)):
        raise ValueError(f"Unknown models: {mismatch}")

    dfs: dict[str, pd.DataFrame] = {}

    for model_name in (bar := tqdm(models, disable=not pbar)):
        bar.set_description(model_name)
        pattern = data_paths[model_name]
        df = glob_to_df(pattern, pbar=False).set_index(id_col)
        dfs[model_name] = df

    return dfs


def load_df_wbm_with_preds(**kwargs: Any) -> pd.DataFrame:
    dfs = load_model_preds(**kwargs)
    df_out = df_wbm.copy()
    for model_name, df in dfs.items():
        model_key = model_name.lower().replace(" ", "_")
        if f"e_form_per_atom_{model_key}" in df:
            df_out[model_name] = df[f"e_form_per_atom_{model_key}"]
        elif len(pred_cols := df.filter(like="_pred_ens").columns) > 0:
            assert len(pred_cols) == 1
            df_out[model_name] = df[pred_cols[0]]
        elif len(pred_cols := df.filter(like=r"_pred_").columns) > 1:
            # make sure we average the expected number of ensemble member predictions
            assert len(pred_cols) == 10, f"{len(pred_cols) = }, expected 10"
            df_out[model_name] = df[pred_cols].mean(axis=1)
        elif "e_form_per_atom_voronoi_rf" in df:  # new voronoi
            df_out[model_name] = df.e_form_per_atom_voronoi_rf
        elif "e_form_pred" in df:  # old voronoi
            df_out[model_name] = df.e_form_pred
        else:
            raise ValueError(
                f"No pred col for {model_name=}, available cols={list(df)}"
            )
    return df_out
