from __future__ import annotations

import os
import urllib.error
from collections.abc import Sequence
from glob import glob
from pathlib import Path
from typing import Any, Callable

import pandas as pd
from pymatgen.core import Structure
from pymatgen.entries.computed_entries import ComputedStructureEntry
from tqdm import tqdm

from matbench_discovery import ROOT

df_wbm = pd.read_csv(f"{ROOT}/data/wbm/2022-10-19-wbm-summary.csv")
df_wbm.index = df_wbm.material_id

# repo URL to raw files on GitHub
RAW_REPO_URL = "https://raw.githubusercontent.com/janosh/matbench-discovery"
# directory to cache downloaded data files
default_cache_dir = os.path.expanduser("~/.cache/matbench-discovery")

DATA_FILENAMES = {
    "mp-computed-structure-entries": "mp/2022-09-16-mp-computed-structure-entries.json.gz",
    "mp-elemental-ref-energies": "mp/2022-09-19-mp-elemental-ref-energies.json",
    "mp-energies": "mp/2022-08-13-mp-energies.json.gz",
    "mp-patched-phase-diagram": "mp/2022-09-18-ppd-mp.pkl.gz",
    "wbm-computed-structure-entries": "wbm/2022-10-19-wbm-computed-structure-entries.json.bz2",
    "wbm-initial-structures": "wbm/2022-10-19-wbm-init-structs.json.bz2",
    "wbm-summary": "wbm/2022-10-19-wbm-summary.csv",
}


def as_dict_handler(obj: Any) -> dict[str, Any] | None:
    """Pass this to json.dump(default=) or as pandas.to_json(default_handler=) to
    convert Python classes with a as_dict() method to dictionaries on serialization.
    Objects without a as_dict() method are replaced with None in the serialized data.
    """
    try:
        return obj.as_dict()  # all MSONable objects implement as_dict()
    except AttributeError:
        return None  # replace unhandled objects with None in serialized data
        # removes e.g. non-serializable AseAtoms from M3GNet relaxation trajectories


def load_train_test(
    data_names: str | Sequence[str] = ("summary",),
    version: str = "1.0.0",
    cache_dir: str | Path = default_cache_dir,
    hydrate: bool = False,
    **kwargs: Any,
) -> pd.DataFrame:
    """Download parts of or the full MP training data and WBM test data as pandas
    DataFrames. The full training and test sets are each about ~500 MB as compressed
    JSON which will be cached locally to cache_dir for faster re-loading unless
    cache_dir is set to None.

    Recognized data keys are mp-computed-structure-entries, mp-elemental-ref-energies,
    mp-energies, mp-patched-phase-diagram, wbm-computed-structure-entries,
    wbm-initial-structures, wbm-summary. See
    https://janosh.github.io/matbench-discovery/how-to-contribute for brief data descriptions.

    Args:
        data_names (str | list[str], optional): Which parts of the MP/WBM dataset to load.
            Can be any subset of the above data names or 'all'. Defaults to ["summary"].
        version (str, optional): Which version of the dataset to load. Defaults to
            '1.0.0'. Can be any git tag, branch or commit hash.
        cache_dir (str, optional): Where to cache data files on local drive. Defaults to
            '~/.cache/matbench-discovery'. Set to None to disable caching.
        hydrate (bool, optional): Whether to hydrate pymatgen objects. If False,
            Structures and ComputedStructureEntries are returned as dictionaries which
            can be hydrated on-demand with df.col.map(Structure.from_dict). Defaults to
            False as it noticeably increases load time.
        **kwargs: Additional keyword arguments passed to pandas.read_json or read_csv,
            depending on which file is loaded.

    Raises:
        ValueError: On bad version number or bad data names.

    Returns:
        pd.DataFrame: Single dataframe or dictionary of dfs if
        multiple data were requested.
    """
    if data_names == "all":
        data_names = list(DATA_FILENAMES)
    elif isinstance(data_names, str):
        data_names = [data_names]

    if missing := set(data_names) - set(DATA_FILENAMES):
        raise ValueError(f"{missing} must be subset of {set(DATA_FILENAMES)}")

    dfs = {}
    for key in data_names:
        file = DATA_FILENAMES[key]
        reader = pd.read_csv if file.endswith(".csv") else pd.read_json

        cache_path = f"{cache_dir}/{version}/{file}"
        if os.path.isfile(cache_path):
            print(f"Loading {key!r} from cached file at {cache_path!r}")
            df = reader(cache_path, **kwargs)
        else:
            url = f"{RAW_REPO_URL}/{version}/data/{file}"
            print(f"Downloading {key!r} from {url}")
            try:
                df = reader(url)
            except urllib.error.HTTPError as exc:
                raise ValueError(f"Bad {url=}") from exc
            if cache_dir and not os.path.isfile(cache_path):
                os.makedirs(os.path.dirname(cache_path), exist_ok=True)
                if ".csv" in file:
                    df.to_csv(cache_path)
                elif ".json" in file:
                    df.reset_index().to_json(
                        cache_path, default_handler=as_dict_handler
                    )
                else:
                    raise ValueError(f"Unexpected file type {file}")

        df = df.set_index("material_id")
        if hydrate:
            for col in df:
                if not isinstance(df[col].iloc[0], dict):
                    continue
                try:
                    df[col] = [
                        ComputedStructureEntry.from_dict(d)
                        for d in tqdm(df[col], desc=col)
                    ]
                except Exception:
                    df[col] = [Structure.from_dict(d) for d in tqdm(df[col], desc=col)]

        dfs[key] = df

    if len(data_names) == 1:
        return dfs[data_names[0]]
    return dfs


PRED_FILENAMES = {
    "CGCNN": "cgcnn/2023-01-26-test-cgcnn-wbm-IS2RE/cgcnn-ensemble-preds.csv",
    "Voronoi Random Forest": "voronoi/2022-11-27-train-test/e-form-preds-IS2RE.csv",
    "Wrenformer": "wrenformer/2022-11-15-wrenformer-IS2RE-preds.csv",
    "MEGNet": "megnet/2022-11-18-megnet-wbm-IS2RE/megnet-e-form-preds.csv",
    "M3GNet": "m3gnet/2022-10-31-m3gnet-wbm-IS2RE.csv",
    "BOWSR MEGNet": "bowsr/2023-01-23-bowsr-megnet-wbm-IS2RE.csv",
}


def glob_to_df(
    pattern: str,
    reader: Callable[[Any], pd.DataFrame] = None,
    pbar: bool = True,
    **kwargs: Any,
) -> pd.DataFrame:
    """Combine data files matching a glob pattern into a single dataframe.

    Args:
        pattern (str): Glob file pattern.
        reader (Callable[[Any], pd.DataFrame], optional): Function that loads data from
            disk. Defaults to pd.read_csv if ".csv" in pattern else pd.read_json.
        pbar (bool, optional): Whether to show progress bar. Defaults to True.
        **kwargs: Keyword arguments passed to reader (i.e. pd.read_csv or pd.read_json).

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
        df = reader(file, **kwargs)
        sub_dfs[file] = df

    return pd.concat(sub_dfs.values())


def load_df_wbm_preds(
    models: Sequence[str],
    pbar: bool = True,
    id_col: str = "material_id",
    return_model_dfs: bool = False,
    **kwargs: Any,
) -> pd.DataFrame:
    """Load WBM summary dataframe with model predictions from disk.

    Args:
        models (Sequence[str]): Model names must be keys of the dict
            matbench_discovery.data.PRED_FILENAMES.
        pbar (bool, optional): Whether to show progress bar. Defaults to True.
        id_col (str, optional): Column to set as df.index. Defaults to "material_id".
        return_model_dfs (bool, optional): Whether to return dict of dataframes for each
            model dfs. Defaults to False.
        **kwargs: Keyword arguments passed to glob_to_df().

    Raises:
        ValueError: On unknown model names.

    Returns:
        pd.DataFrame: WBM summary dataframe with model predictions.
    """
    if mismatch := ", ".join(set(models) - set(PRED_FILENAMES)):
        raise ValueError(f"Unknown models: {mismatch}")

    dfs: dict[str, pd.DataFrame] = {}

    for model_name in (bar := tqdm(models, disable=not pbar)):
        bar.set_description(model_name)
        pattern = f"models/{PRED_FILENAMES[model_name]}"
        df = glob_to_df(pattern, pbar=False, **kwargs).set_index(id_col)
        dfs[model_name] = df

    if return_model_dfs:
        return dfs

    df_out = df_wbm.copy()
    for model_name, df in dfs.items():
        model_key = model_name.lower().replace(" ", "_")
        if f"e_form_per_atom_{model_key}" in df:
            df_out[model_name] = df[f"e_form_per_atom_{model_key}"]

        elif len(pred_cols := df.filter(like="_pred_ens").columns) > 0:
            assert len(pred_cols) == 1
            df_out[model_name] = df[pred_cols[0]]
            if len(std_cols := df.filter(like="_std_ens").columns) > 0:
                df_out[f"{model_name}_std"] = df[std_cols[0]]

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
