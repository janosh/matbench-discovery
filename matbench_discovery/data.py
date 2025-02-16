"""Download, cache and hydrate data files from the Matbench Discovery Figshare article.

https://figshare.com/articles/dataset/22715158

Environment Variables:
    MBD_AUTO_DOWNLOAD_FILES: Controls whether to auto-download missing data files.
        Defaults to "true". Set to "false" to be prompted before downloading.
        This affects both model prediction files and dataset files.
    MBD_CACHE_DIR: Directory to cache downloaded data files.
        Defaults to DATA_DIR if the full repo was cloned, otherwise
        ~/.cache/matbench-discovery.
"""

import io
import os
import re
import sys
import zipfile
from collections import defaultdict
from collections.abc import Callable, Sequence
from glob import glob
from pathlib import Path
from typing import Any

import ase.io
import pandas as pd
import yaml
from ase import Atoms
from pymatviz.enums import Key
from ruamel.yaml import YAML
from tqdm import tqdm

from matbench_discovery import TEST_FILES
from matbench_discovery.enums import DataFiles, MbdKey, Model, TestSubset

round_trip_yaml = YAML()  # round-trippable YAML for updating model metadata files
round_trip_yaml.preserve_quotes = True
round_trip_yaml.width = 1000  # avoid changing line wrapping
round_trip_yaml.indent(mapping=2, sequence=4, offset=2)


def as_dict_handler(obj: Any) -> dict[str, Any] | None:
    """Pass this to json.dump(default=) or as pandas.to_json(default_handler=) to
    serialize Python classes with as_dict(). Warning: Objects without a as_dict() method
    are replaced with None in the serialized data.
    """
    try:
        return obj.as_dict()  # all MSONable objects implement as_dict()
    except AttributeError:
        return None  # replace unhandled objects with None in serialized data
        # removes e.g. non-serializable AseAtoms from M3GNet relaxation trajectories


def glob_to_df(
    pattern: str,
    *,
    reader: Callable[[Any], pd.DataFrame] | None = None,
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

    Raises:
        FileNotFoundError: If no files match the glob pattern.
        ValueError: If reader is None and the file extension is unrecognized.
    """
    if reader is None:
        if ".csv" in pattern.lower():
            reader = pd.read_csv
        elif ".json" in pattern.lower():
            reader = pd.read_json
        else:
            raise ValueError(f"Unsupported file extension in {pattern=}")

    files = glob(pattern)

    if len(files) == 0:
        # load mocked model predictions when running pytest (just first 500 lines
        # from MACE-MPA-0 WBM energy preds)
        if "pytest" in sys.modules or "CI" in os.environ:
            df_mock = pd.read_csv(f"{TEST_FILES}/mock-wbm-energy-preds.csv.gz")
            # .set_index( "material_id" )
            # make sure pred_cols for all models are present in df_mock
            for model in Model:
                with open(model.yaml_path) as file:
                    model_data = yaml.safe_load(file)

                pred_col = (
                    model_data.get("metrics", {}).get("discovery", {}).get("pred_col")
                )
                df_mock[pred_col] = df_mock["e_form_per_atom"]
            return df_mock
        raise FileNotFoundError(f"No files matching glob {pattern=}")

    sub_dfs = {}  # used to join slurm job array results into single df
    for file in tqdm(files, disable=not pbar):
        df_i = reader(file, **kwargs)
        sub_dfs[file] = df_i

    return pd.concat(sub_dfs.values())


def ase_atoms_from_zip(
    zip_filename: str | Path,
    *,
    file_filter: Callable[[str], bool] = lambda fname: fname.endswith(".extxyz"),
    filename_to_info: bool = False,
    limit: int | None = None,
) -> list[Atoms]:
    """Read ASE Atoms objects from a ZIP file containing extXYZ files.

    Args:
        zip_filename (str): Path to the ZIP file.
        file_filter (Callable[[str], bool], optional): Function to check if a file
            should be read. Defaults to lambda fname: fname.endswith(".extxyz").
        filename_to_info (bool, optional): If True, assign filename to Atoms.info.
            Defaults to False.
        limit (int, optional): Maximum number of files to read. Defaults to None.
            Use a small number to speed up debugging runs.

    Returns:
        list[Atoms]: ASE Atoms objects.
    """
    atoms_list = []
    with zipfile.ZipFile(zip_filename) as zip_file:
        desc = f"Reading ASE Atoms from {zip_filename=}"
        for filename in tqdm(zip_file.namelist()[:limit], desc=desc, mininterval=5):
            if not file_filter(filename):
                continue
            with zip_file.open(filename) as file:
                content = io.TextIOWrapper(file, encoding="utf-8").read()
                atoms = ase.io.read(
                    io.StringIO(content), format="extxyz", index=slice(None)
                )  # reads multiple Atoms objects as frames if file contains trajectory
                if isinstance(atoms, Atoms):
                    atoms = [atoms]  # Wrap single Atoms object in a list
                if filename_to_info:
                    for atom in atoms:
                        atom.info["filename"] = filename
                atoms_list.extend(atoms)
    return atoms_list


def ase_atoms_to_zip(
    atoms_set: list[Atoms] | dict[str, Atoms], zip_filename: str | Path
) -> None:
    """Write ASE Atoms objects to a ZIP archive with each Atoms object as a separate
    extXYZ file, grouped by mat_id.

    Args:
        atoms_set (list[Atoms] | dict[str, Atoms]): Either a list of ASE Atoms objects
            (which should have a 'material_id' in their Atoms.info dictionary) or a
            dictionary mapping material IDs to Atoms objects.
        zip_filename (str | Path): Path to the ZIP file to write.
    """
    # Group atoms by mat_id to avoid overwriting files with the same name

    if isinstance(atoms_set, dict):
        atoms_dict = atoms_set
    else:
        atoms_dict = defaultdict(list)

        # If input is a list, get material ID from atoms.info falling back to formula
        # if missing
        for atoms in atoms_set:
            mat_id = atoms.info.get(Key.mat_id, f"no-id-{atoms.get_chemical_formula()}")
            atoms_dict[mat_id] += [atoms]

    # Write grouped atoms to the ZIP archive
    with zipfile.ZipFile(
        zip_filename, mode="w", compression=zipfile.ZIP_DEFLATED
    ) as zip_file:
        for mat_id, atoms_or_list in tqdm(
            atoms_dict.items(), desc=f"Writing ASE Atoms to {zip_filename=}"
        ):
            buffer = io.StringIO()  # string buffer to write the extxyz content
            for atoms in (
                atoms_or_list if isinstance(atoms_or_list, list) else [atoms_or_list]
            ):
                ase.io.write(
                    buffer, atoms, format="extxyz", append=True, write_info=True
                )

            # Write the combined buffer content to the ZIP file
            zip_file.writestr(f"{mat_id}.extxyz", buffer.getvalue())


df_wbm = pd.read_csv(DataFiles.wbm_summary.path)
# str() around Key.mat_id added for https://github.com/janosh/matbench-discovery/issues/81
df_wbm.index = df_wbm[str(Key.mat_id)]


def load_df_wbm_with_preds(
    *,
    models: Sequence[str | Model] = (),
    pbar: bool = True,
    id_col: str = Key.mat_id,
    subset: pd.Index | Sequence[str] | TestSubset | None = None,
    max_error_threshold: float | None = 5.0,
    **kwargs: Any,
) -> pd.DataFrame:
    """Load WBM summary dataframe with model predictions from disk.

    Args:
        models (Sequence[str], optional): Model names must be keys of
            matbench_discovery.data.Model. Defaults to all models.
        pbar (bool, optional): Whether to show progress bar. Defaults to True.
        id_col (str, optional): Column to set as df.index. Defaults to "material_id".
        subset (pd.Index | Sequence[str] | 'uniq_protos' | None, optional):
            Subset of material IDs to keep. Defaults to None, which loads all materials.
            'uniq_protos' drops WBM structures with matching prototype in MP
            training set and duplicate prototypes in WBM test set (keeping only the most
            stable structure per prototype). This increases the 'OOD-ness' of WBM.
        max_error_threshold (float, optional): Maximum absolute error between predicted
            and DFT formation energies before a prediction is filtered out as
            unrealistic. Doing this filtering is acceptable as it could also be done by
            a practitioner doing a prospective discovery effort. Predictions exceeding
            this threshold will be ignored in all downstream calculations of metrics.
            Defaults to 5 eV/atom.
        **kwargs: Keyword arguments passed to glob_to_df().

    Raises:
        ValueError: On unknown model names.

    Returns:
        pd.DataFrame: WBM summary dataframe with model predictions.
    """
    valid_models = {model.name for model in Model}
    if models == ():
        models = tuple(valid_models)
    inv_label_map = {key.label: key.name for key in Model}
    # map pretty model names back to Model enum keys
    models = [inv_label_map.get(model, model) for model in models]
    if unknown_models := ", ".join(set(models) - valid_models):
        raise ValueError(f"{unknown_models=}, expected subset of {valid_models}")

    model_name: str = ""
    from matbench_discovery.data import df_wbm

    df_out = df_wbm.copy()

    try:
        prog_bar = tqdm(models, disable=not pbar, desc="Loading preds")
        for model_name in prog_bar:
            prog_bar.set_postfix_str(model_name)

            # use getattr(name) in case model_name is already a Model enum
            model = Model[getattr(model_name, "name", model_name)]

            df_preds = glob_to_df(model.discovery_path, pbar=False, **kwargs)

            with open(model.yaml_path) as file:
                model_data = yaml.safe_load(file)

            pred_col = (
                model_data.get("metrics", {}).get("discovery", {}).get("pred_col")
            )
            if not pred_col:
                raise ValueError(
                    f"pred_col not specified for {model_name} in {model.yaml_path!r}"
                )

            if pred_col not in df_preds:
                raise ValueError(
                    f"{pred_col=} set in {model.yaml_path!r}:metrics.discovery."
                    f"pred_col not found in {model.discovery_path}"
                )

            df_out[model.label] = df_preds.set_index(id_col)[pred_col]
            if max_error_threshold is not None:
                if max_error_threshold < 0:
                    raise ValueError("max_error_threshold must be a positive number")
                # Apply centralized model prediction cleaning criterion (see doc string)
                bad_mask = (
                    abs(df_out[model.label] - df_out[MbdKey.e_form_dft])
                ) > max_error_threshold
                df_out.loc[bad_mask, model.label] = pd.NA
                n_preds, n_bad = len(df_out[model.label].dropna()), sum(bad_mask)
                if n_bad > 0:
                    print(
                        f"{n_bad:,} of {n_preds:,} unrealistic preds for {model_name}"
                    )
    except Exception as exc:
        exc.add_note(f"Failed to load {model_name=}")
        raise

    if subset == TestSubset.uniq_protos:
        df_out = df_out.query(MbdKey.uniq_proto)
    elif subset is not None:
        df_out = df_out.loc[subset]

    return df_out


def update_yaml_at_path(
    file_path: str | Path,
    dotted_path: str,
    data: dict[str, Any],
) -> dict[str, Any]:
    """Update a YAML file at a specific dotted path with new data.

    Args:
        file_path (str | Path): Path to YAML file to update
        dotted_path (str): Dotted path to update (e.g. 'metrics.discovery')
        data (dict[str, Any]): Data to write at the specified path

    Returns:
        dict[str, Any]: The complete updated YAML data written to file.

    Example:
        update_yaml_at_path(
            "models/mace/mace-mp-0.yml",
            "metrics.discovery",
            dict(mae=0.1, rmse=0.2),
        )
    """
    # raise on repeated or trailing dots in dotted path
    if not re.match(r"^[a-zA-Z0-9-+=_]+(\.[a-zA-Z0-9-+=_]+)*$", dotted_path):
        raise ValueError(f"Invalid dotted path: {dotted_path}")

    with open(file_path) as file:
        yaml_data = round_trip_yaml.load(file)

    # Navigate to the correct nested level
    current = yaml_data
    *parts, last = dotted_path.split(".")

    for part in parts:
        if part not in current:
            current[part] = {}
        current = current[part]

    # Update the data at the final level
    if last not in current:
        current[last] = {}
    # Replace the entire section to preserve comments
    current[last] = data

    # Write back to file
    with open(file_path, mode="w") as file:
        round_trip_yaml.dump(yaml_data, file)

    return yaml_data
