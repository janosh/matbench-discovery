"""Download, cache and hydrate data files from the Matbench Discovery Figshare article.

https://figshare.com/articles/dataset/22715158
"""

import io
import os
import sys
import zipfile
from collections import defaultdict
from collections.abc import Callable, Sequence
from enum import EnumMeta, StrEnum, _EnumDict
from glob import glob
from pathlib import Path
from typing import Any, Literal, Self, TypeVar

import ase.io
import pandas as pd
import plotly.express as px
import requests
import yaml
from ase import Atoms
from pymatviz.enums import Key
from ruamel.yaml import YAML
from tqdm import tqdm

from matbench_discovery import DATA_DIR, ROOT, pkg_is_editable
from matbench_discovery.enums import MbdKey, TestSubset

# ruff: noqa: T201
T = TypeVar("T", bound="Files")

# repo URL to raw files on GitHub
RAW_REPO_URL = "https://github.com/janosh/matbench-discovery/raw"
# directory to cache downloaded data files
DEFAULT_CACHE_DIR = os.getenv(
    "MATBENCH_DISCOVERY_CACHE_DIR",
    DATA_DIR if pkg_is_editable else os.path.expanduser("~/.cache/matbench-discovery"),
)


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
    """
    reader = reader or pd.read_csv if ".csv" in pattern.lower() else pd.read_json

    # prefix pattern with ROOT if not absolute path
    files = glob(pattern)
    if len(files) == 0:
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
        for filename in tqdm(zip_file.namelist()[:limit], desc=desc):
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

        # If input is a list, get material ID from atoms.info falling back to formula if missing
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


def download_file(file_path: str, url: str) -> None:
    """Download the file from the given URL to the given file path."""
    file_dir = os.path.dirname(file_path)
    os.makedirs(file_dir, exist_ok=True)
    try:
        response = requests.get(url, timeout=5)

        response.raise_for_status()

        with open(file_path, "wb") as file:
            file.write(response.content)
    except requests.exceptions.RequestException as exc:
        print(f"Error downloading {url=}\nto {file_path=}.\n{exc!s}")


class MetaFiles(EnumMeta):
    """Metaclass of Files enum that adds base_dir and (member|label)_map class
    properties.
    """

    _base_dir: str

    def __new__(
        cls,
        name: str,
        bases: tuple[type, ...],
        namespace: _EnumDict,
        base_dir: str = DEFAULT_CACHE_DIR,
        **kwargs: Any,
    ) -> "MetaFiles":
        """Create new Files enum with given base directory."""
        obj = super().__new__(cls, name, bases, namespace, **kwargs)
        obj._base_dir = base_dir  # noqa: SLF001
        return obj

    # Improvement 2: Add type annotation for cls
    @property
    def member_map(cls: type[T]) -> dict[str, "Files"]:  # type: ignore[misc]
        """Map of member names to member objects."""
        return cls._member_map_  # type: ignore[return-value]

    @property
    def label_map(cls: type[T]) -> dict[str, str]:  # type: ignore[misc]
        """Map of member names to member labels."""
        return {k: v.label for k, v in cls.member_map.items()}

    @property
    def base_dir(cls) -> str:
        """Return the base directory of the file URL."""
        return cls._base_dir


class Files(StrEnum, metaclass=MetaFiles):
    """Enum of data files with associated file directories and URLs."""

    def __new__(
        cls, file_path: str, url: str | None = None, label: str | None = None
    ) -> Self:
        """Create a new member of the FileUrls enum with a given URL where to load the
        file from and directory where to save it to.
        """
        obj = str.__new__(cls)
        obj._value_ = file_path.split("/")[-1]  # use file name as enum value

        obj._rel_path = file_path  # type: ignore[attr-defined] # noqa: SLF001
        obj._url = url  # type: ignore[attr-defined] # noqa: SLF001
        obj._label = label  # type: ignore[attr-defined] # noqa: SLF001

        return obj

    def __str__(self) -> str:
        """File path associated with the file URL. Use str(DataFiles.some_key) if you
        want the absolute file path without auto-downloading the file if it doesn't
        exist yet, e.g. for use in script that generates the file in the first place.
        """
        return f"{type(self).base_dir}/{self._rel_path}"  # type: ignore[attr-defined]

    def __repr__(self) -> str:
        """Return enum attribute's string representation."""
        return f"{type(self).__name__}.{self.name}"

    @property
    def path(self) -> str:
        """Return the file path associated with the file URL if it exists, otherwise
        download the file first, then return the path.
        """
        key, url, rel_path = self.name, self._url, self._rel_path  # type: ignore[attr-defined]
        abs_path = f"{type(self).base_dir}/{rel_path}"
        if not os.path.isfile(abs_path):
            is_ipython = hasattr(__builtins__, "__IPYTHON__")
            # default to 'y' if not in interactive session, and user can't answer
            answer = (
                input(
                    f"{abs_path!r} associated with {key=} does not exist. Download it "
                    "now? This will cache the file for future use. [y/n] "
                )
                if is_ipython or sys.stdin.isatty()
                else "y"
            )
            if answer.lower().strip() == "y":
                if not is_ipython:
                    print(f"Downloading {key!r} from {url} to {abs_path} for caching")
                download_file(abs_path, url)
        return abs_path

    @property
    def url(self) -> str:
        """Return the URL associated with the file URL."""
        return self._url  # type: ignore[attr-defined]

    @property
    def rel_path(self) -> str:
        """Return the relative path of the file associated with the file URL."""
        return self._rel_path  # type: ignore[attr-defined]

    @property
    def label(self) -> str:
        """Return the label associated with the file URL."""
        return self._label  # type: ignore[attr-defined]


class DataFiles(Files):
    """Enum of data files with associated file directories and URLs."""

    mp_computed_structure_entries = (
        "mp/2023-02-07-mp-computed-structure-entries.json.gz",
        "https://figshare.com/ndownloader/files/40344436",
    )
    mp_elemental_ref_entries = (
        "mp/2023-02-07-mp-elemental-reference-entries.json.gz",
        "https://figshare.com/ndownloader/files/40387775",
    )
    mp_energies = (
        "mp/2023-01-10-mp-energies.csv.gz",
        "https://figshare.com/ndownloader/files/49083124",
    )
    mp_patched_phase_diagram = (
        "mp/2023-02-07-ppd-mp.pkl.gz",
        "https://figshare.com/ndownloader/files/48241624",
    )
    mp_trj_extxyz = (
        "mp/2024-09-03-mp-trj.extxyz.zip",
        "https://figshare.com/ndownloader/files/49034296",
    )
    # snapshot of every task (calculation) in MP as of 2023-03-16 (14 GB)
    all_mp_tasks = (
        "mp/2023-03-16-all-mp-tasks.zip",
        "https://figshare.com/ndownloader/files/43350447",
    )

    wbm_computed_structure_entries = (
        "wbm/2022-10-19-wbm-computed-structure-entries.json.bz2",
        "https://figshare.com/ndownloader/files/40344463",
    )
    wbm_relaxed_atoms = (
        "wbm/2024-08-04-wbm-relaxed-atoms.extxyz.zip",
        "https://figshare.com/ndownloader/files/48169600",
    )
    wbm_initial_structures = (
        "wbm/2022-10-19-wbm-init-structs.json.bz2",
        "https://figshare.com/ndownloader/files/40344466",
    )
    wbm_initial_atoms = (
        "wbm/2024-08-04-wbm-initial-atoms.extxyz.zip",
        "https://figshare.com/ndownloader/files/48169597",
    )
    wbm_cses_plus_init_structs = (
        "wbm/2022-10-19-wbm-computed-structure-entries+init-structs.json.bz2",
        "https://figshare.com/ndownloader/files/40344469",
    )
    wbm_summary = (
        "wbm/2023-12-13-wbm-summary.csv.gz",
        "https://figshare.com/ndownloader/files/44225498",
    )
    alignn_checkpoint = (
        "2023-06-02-pbenner-best-alignn-model.pth.zip",
        "https://figshare.com/ndownloader/files/41233560",
    )
    mp_trj = (
        "mp/2022-09-16-mp-trj.json",
        "https://figshare.com/ndownloader/files/41619375",
    )


df_wbm = pd.read_csv(DataFiles.wbm_summary.path)
# str() around Key.mat_id added for https://github.com/janosh/matbench-discovery/issues/81
df_wbm.index = df_wbm[str(Key.mat_id)]


# ruff: noqa: E501, ERA001 (ignore long lines in class Model)
class Model(Files, base_dir=f"{ROOT}/models"):
    """Data files provided by Matbench Discovery.
    See https://janosh.github.io/matbench-discovery/contribute for data descriptions.
    """

    alignn = "alignn/2023-06-02-alignn-wbm-IS2RE.csv.gz", "alignn/alignn.yml", "ALIGNN"
    # alignn_pretrained = "alignn/2023-06-03-mp-e-form-alignn-wbm-IS2RE.csv.gz", "alignn/alignn.yml", "ALIGNN Pretrained"
    # alignn_ff = "alignn_ff/2023-07-11-alignn-ff-wbm-IS2RE.csv.gz", "alignn/alignn-ff.yml", "ALIGNN FF"

    # BOWSR optimizer coupled with original megnet
    bowsr_megnet = "bowsr/2023-01-23-bowsr-megnet-wbm-IS2RE.csv.gz", "bowsr/bowsr.yml", "BOWSR"  # fmt: skip

    # default CHGNet model from publication with 400,438 params
    chgnet = "chgnet/2023-12-21-chgnet-0.3.0-wbm-IS2RE.csv.gz", "chgnet/chgnet.yml", "CHGNet"  # fmt: skip
    # chgnet_no_relax = "chgnet/2023-12-05-chgnet-0.3.0-wbm-IS2RE-no-relax.csv.gz", None, "CHGNet No Relax"

    # CGCNN 10-member ensemble
    cgcnn = "cgcnn/2023-01-26-cgcnn-ens=10-wbm-IS2RE.csv.gz", "cgcnn/cgcnn.yml", "CGCNN"

    # CGCNN 10-member ensemble with 5-fold training set perturbations
    cgcnn_p = "cgcnn/2023-02-05-cgcnn-perturb=5-wbm-IS2RE.csv.gz", "cgcnn/cgcnn+p.yml", "CGCNN+P"  # fmt: skip

    # original M3GNet straight from publication, not re-trained
    m3gnet = "m3gnet/2023-12-28-m3gnet-wbm-IS2RE.csv.gz", "m3gnet/m3gnet.yml", "M3GNet"
    # m3gnet_direct = "m3gnet/2023-05-30-m3gnet-direct-wbm-IS2RE.csv.gz", None, "M3GNet DIRECT"
    # m3gnet_ms = "m3gnet/2023-06-01-m3gnet-manual-sampling-wbm-IS2RE.csv.gz", None, "M3GNet MS"

    # MACE-MP as published in https://arxiv.org/abs/2401.00096 trained on MPtrj
    mace = "mace/2023-12-11-mace-wbm-IS2RE-FIRE.csv.gz", "mace/mace.yml", "MACE"
    # mace_alex = "mace/2024-08-09-mace-wbm-IS2RE-FIRE.csv.gz", None, "MACE Alex"
    # https://github.com/ACEsuit/mace-mp/releases/tag/mace_mp_0b
    # mace_0b = "mace/2024-07-20-mace-wbm-IS2RE-FIRE.csv.gz", None, "MACE 0b"

    # original MEGNet straight from publication, not re-trained
    megnet = "megnet/2022-11-18-megnet-wbm-IS2RE.csv.gz", "megnet/megnet.yml", "MEGNet"

    # SevenNet trained on MPtrj
    sevennet = "sevennet/2024-07-11-sevennet-preds.csv.gz", "sevennet/sevennet.yml", "SevenNet"  # fmt: skip

    # Magpie composition+Voronoi tessellation structure features + sklearn random forest
    voronoi_rf = "voronoi_rf/2022-11-27-train-test/e-form-preds-IS2RE.csv.gz", "voronoi_rf/voronoi-rf.yml", "Voronoi RF"  # fmt: skip

    # wrenformer 10-member ensemble
    wrenformer = "wrenformer/2022-11-15-wrenformer-ens=10-IS2RE-preds.csv.gz", "wrenformer/wrenformer.yml", "Wrenformer"  # fmt: skip

    # --- Proprietary Models
    # GNoME
    gnome = "gnome/2023-11-01-gnome-preds-50076332.csv.gz", "gnome/gnome.yml", "GNoME"

    # MatterSim
    mattersim = "mattersim/2024-06-16-mattersim-wbm-IS2RE.csv.gz", "mattersim/mattersim.yml", "MatterSim"  # fmt: skip

    # ORB
    orb = "orb/orbff-v2-20241011.csv.gz", "orb/orb.yml", "ORB"
    orb_mptrj = "orb/orbff-mptrj-only-v2-20241014.csv.gz", "orb/orb-mptrj.yml", "ORB MPtrj"  # fmt: skip

    # fairchem
    eqv2_s_dens = "eqV2/eqV2-s-dens-mp.csv.gz", "eqV2/eqV2-s-dens-mp.yml", "eqV2 S DeNS"  # fmt: skip
    eqv2_m = "eqV2/eqV2-m-omat-mp-salex.csv.gz", "eqV2/eqV2-m-omat-mp-salex.yml", "eqV2 M"  # fmt: skip

    # --- Model Combos
    # # CHGNet-relaxed structures fed into MEGNet for formation energy prediction
    # chgnet_megnet = "chgnet/2023-03-06-chgnet-0.2.0-wbm-IS2RE.csv.gz", None, "CHGNet→MEGNet"
    # # M3GNet-relaxed structures fed into MEGNet for formation energy prediction
    # m3gnet_megnet = "m3gnet/2022-10-31-m3gnet-wbm-IS2RE.csv.gz", None, "M3GNet→MEGNet"
    # megnet_rs2re = "megnet/2023-08-23-megnet-wbm-RS2RE.csv.gz", None, "MEGNet RS2RE"


px.defaults.labels |= Model.label_map


def load_df_wbm_with_preds(
    *,
    models: Sequence[str] = (),
    pbar: bool = True,
    id_col: str = Key.mat_id,
    subset: pd.Index | Sequence[str] | Literal[TestSubset.uniq_protos] | None = None,
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
    inv_label_map = {v: k for k, v in Model.label_map.items()}
    # map pretty model names back to Model enum keys
    models = {inv_label_map.get(model, model) for model in models}
    if unknown_models := ", ".join(models - valid_models):
        raise ValueError(f"{unknown_models=}, expected subset of {valid_models}")

    model_name: str = ""
    from matbench_discovery.data import df_wbm

    df_out = df_wbm.copy()

    try:
        prog_bar = tqdm(models, disable=not pbar, desc="Loading preds")
        for model_name in prog_bar:
            prog_bar.set_postfix_str(model_name)
            pred_file = Model[model_name]
            df_preds = glob_to_df(pred_file.path, pbar=False, **kwargs)

            # Get prediction column name from metadata
            model_key = getattr(Model, model_name)
            model_label = model_key.label
            model_yaml_path = f"{Model.base_dir}/{model_key.url}"
            with open(model_yaml_path) as file:
                model_data = yaml.safe_load(file)

            pred_col = (
                model_data.get("metrics", {}).get("discovery", {}).get("pred_col")
            )
            if not pred_col:
                raise ValueError(
                    f"pred_col not specified for {model_name} in {model_yaml_path!r}"
                )

            if pred_col not in df_preds:
                raise ValueError(f"{pred_col=} not found in {pred_file.path}")

            df_out[model_label] = df_preds.set_index(id_col)[pred_col]
            if max_error_threshold is not None:
                if max_error_threshold < 0:
                    raise ValueError("max_error_threshold must be a positive number")
                # Apply centralized model prediction cleaning criterion (see doc string)
                bad_mask = (
                    abs(df_out[model_label] - df_out[MbdKey.e_form_dft])
                ) > max_error_threshold
                df_out.loc[bad_mask, model_label] = pd.NA
                n_preds, n_bad = len(df_out[model_label].dropna()), sum(bad_mask)
                if n_bad > 0:
                    print(
                        f"{n_bad:,} of {n_preds:,} unrealistic preds for {model_name}"
                    )
    except Exception as exc:
        exc.add_note(f"Failed to load {model_name=}")
        raise

    if subset == TestSubset.uniq_protos:
        df_out = df_out.query(Key.uniq_proto)
    elif subset is not None:
        df_out = df_out.loc[subset]

    return df_out
