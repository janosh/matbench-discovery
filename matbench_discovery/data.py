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
from datetime import date
from decimal import Decimal, InvalidOperation
from glob import glob
from pathlib import Path
from typing import Any, Final, TypedDict

import ase.io
import pandas as pd
import yaml
from ase import Atoms
from filelock import FileLock
from pymatviz.enums import Key
from ruamel.yaml import YAML
from tqdm import tqdm

from matbench_discovery import DATA_DIR, TEST_FILES
from matbench_discovery.enums import DataFiles, MbdKey, Model, TestSubset

round_trip_yaml = YAML()  # round-trippable YAML for updating model metadata files
round_trip_yaml.preserve_quotes = True
round_trip_yaml.width = 1000  # avoid changing line wrapping
round_trip_yaml.indent(mapping=2, sequence=4, offset=2)

ISO_DATE_PATTERN: Final = re.compile(r"^\d{4}-\d{2}-\d{2}$")
MOYO_VERSION_PATTERN: Final = re.compile(
    r"^[0-9]+(?:\.[0-9]+)*(?:[-+][0-9A-Za-z.-]+)?$"
)
_GEO_OPT_ANALYSIS_SUFFIX: Final = re.compile(
    r"^geo-opt-symprec=([^=]+)-moyo=([^=]+)\.csv\.gz$"
)

# role -> suffix after ``YYYY-MM-DD-``
ARTIFACT_SUFFIXES: Final[dict[str, str]] = {
    "discovery": "discovery.csv.gz",
    "geo_opt": "geo-opt.jsonl.gz",
    "phonons_kappa_103": "phonons-kappa-103.json.gz",
    "phonons_kappa_103_forces": "phonons-kappa-103-forces.json.gz",
    "phonons_kappa_103_run_info": "phonons-kappa-103-run-info.json",
    "md_metrics": "md-metrics.csv.gz",
    "diatomics": "diatomics.json.gz",
}

FILE_REF_KEYS: Final = frozenset(
    {"pred_file", "analysis_file", "force_file", "run_info_file"}
)
# Tasks that require forces (not applicable when targets == "E")
FORCE_TASKS: Final = frozenset({"geo_opt", "phonons", "md", "diatomics"})
BENCHMARK_TASKS: Final = ("discovery", "geo_opt", "phonons", "md", "diatomics")
_COVERAGE_META_KEYS: Final = frozenset({"status", "reason"})


class FileRef(TypedDict, total=False):
    """Local path plus optional download URL and content checksums."""

    name: str
    url: str
    size: int
    md5: str


def task_coverage(metadata: dict[str, Any], task: str) -> tuple[str, str | None]:
    """Derive task coverage from metrics, targets, and lifecycle.

    Precedence:
    1. Explicit ``metrics.<task>.status`` (optional exception note)
    2. Metrics results present → complete
    3. ``lifecycle: aborted`` → not_available
    4. ``targets: E`` for force tasks → not_applicable
    5. Else → pending
    """
    metrics_root = metadata.get("metrics")
    task_metrics = metrics_root.get(task) if isinstance(metrics_root, dict) else None
    if isinstance(task_metrics, dict):
        status = task_metrics.get("status")
        if isinstance(status, str):
            reason = task_metrics.get("reason")
            return status, reason if isinstance(reason, str) else None
        if any(key not in _COVERAGE_META_KEYS for key in task_metrics):
            return "complete", None

    if metadata.get("lifecycle") == "aborted":
        return "not_available", None
    if metadata.get("targets") == "E" and task in FORCE_TASKS:
        return "not_applicable", "energy-only (no forces)"
    return "pending", None


def make_file_ref(
    name: str,
    *,
    url: str | None = None,
    size: int | None = None,
    md5: str | None = None,
) -> FileRef:
    """Build a file reference; ``size``/``md5`` must both be set or both omitted."""
    if (size is None) ^ (md5 is None):
        raise ValueError("size and md5 must both be set or both omitted")
    if name.startswith("models/"):
        parse_artifact_filename(name)
    ref: FileRef = {"name": name}
    if url is not None:
        ref["url"] = url
    if size is not None and md5 is not None:
        ref["size"], ref["md5"] = size, md5
    return ref


def file_ref_name(ref: object) -> str | None:
    """Return the local path from a nested file ref or legacy string path."""
    if isinstance(ref, str):
        return ref
    if isinstance(ref, dict):
        name = ref.get("name")
        return name if isinstance(name, str) else None
    return None


def file_ref_url(ref: object) -> str | None:
    """Return the download URL from a nested file ref, if present."""
    if isinstance(ref, dict):
        url = ref.get("url")
        return url if isinstance(url, str) else None
    return None


def canonical_scientific_notation(value: float | str | Decimal) -> str:
    """Format a positive finite number without exponent padding (e.g. ``1e-5``)."""
    try:
        decimal_value = Decimal(str(value))
    except InvalidOperation as exc:
        raise ValueError(f"Invalid numeric value {value!r}") from exc
    if not decimal_value.is_finite() or decimal_value <= 0:
        raise ValueError(f"Expected a positive finite number, got {value!r}")

    normalized = decimal_value.normalize()
    sign, digits, exponent = normalized.as_tuple()
    if not isinstance(exponent, int):
        raise TypeError(f"Expected a finite decimal exponent, got {value!r}")
    coefficient = "".join(map(str, digits))
    scientific_exponent = exponent + len(coefficient) - 1
    mantissa = coefficient[0]
    if len(coefficient) > 1:
        mantissa += f".{coefficient[1:].rstrip('0')}"
        mantissa = mantissa.rstrip(".")
    if sign:
        mantissa = f"-{mantissa}"
    return f"{mantissa}e{scientific_exponent}"


def _iso_date(value: date | str) -> str:
    """Return a validated ISO calendar date."""
    iso_date = value.isoformat() if isinstance(value, date) else value
    if not ISO_DATE_PATTERN.fullmatch(iso_date):
        raise ValueError(f"Expected an ISO date, got {value!r}")
    try:
        date.fromisoformat(iso_date)
    except ValueError as exc:
        raise ValueError(f"Invalid ISO date {value!r}") from exc
    return iso_date


def artifact_date_from_prefix(prefix: str, *, fallback: str) -> str:
    """Return the leading ISO date from an artifact prefix, else ``fallback``."""
    try:
        return _iso_date(os.path.basename(prefix)[:10])
    except ValueError:
        return fallback


def artifact_filename(
    artifact_date: date | str,
    role: str,
    *,
    symprec: float | str | Decimal | None = None,
    moyo_version: str | None = None,
) -> str:
    """Return a canonical dated artifact basename for ``role``."""
    iso_date = _iso_date(artifact_date)
    if role == "geo_opt_analysis":
        if symprec is None or moyo_version is None:
            raise ValueError(
                "symprec and moyo_version are required for geo_opt_analysis"
            )
        if not MOYO_VERSION_PATTERN.fullmatch(moyo_version):
            raise ValueError(f"Invalid moyo version {moyo_version!r}")
        suffix = (
            f"geo-opt-symprec={canonical_scientific_notation(symprec)}"
            f"-moyo={moyo_version}.csv.gz"
        )
    else:
        if symprec is not None or moyo_version is not None:
            raise ValueError("symprec and moyo_version are only for geo_opt_analysis")
        try:
            suffix = ARTIFACT_SUFFIXES[role]
        except KeyError as exc:
            raise ValueError(f"Unknown artifact role {role!r}") from exc
    return f"{iso_date}-{suffix}"


def parse_artifact_filename(filename: str) -> str:
    """Validate a canonical artifact filename and return its role key."""
    basename = filename.rsplit("/", maxsplit=1)[-1]
    if not ISO_DATE_PATTERN.match(basename[:10]) or basename[10:11] != "-":
        raise ValueError(f"Not a canonical model artifact filename: {filename!r}")
    artifact_date, suffix = basename[:10], basename[11:]
    _iso_date(artifact_date)
    for role, expected_suffix in ARTIFACT_SUFFIXES.items():
        if suffix == expected_suffix:
            return role
    if match := _GEO_OPT_ANALYSIS_SUFFIX.fullmatch(suffix):
        symprec, moyo_version = match.groups()
        if (
            artifact_filename(
                artifact_date,
                "geo_opt_analysis",
                symprec=symprec,
                moyo_version=moyo_version,
            )
            == basename
        ):
            return "geo_opt_analysis"
    raise ValueError(f"Not a canonical model artifact filename: {filename!r}")


with open(f"{DATA_DIR}/datasets.yml", encoding="utf-8") as file:
    DATASETS = yaml.safe_load(stream=file)


def as_dict_handler(obj: object) -> dict[str, Any] | None:
    """Pass this to json.dump(default=) or as pandas.to_json(default_handler=) to
    serialize Python classes with as_dict(). Warning: Objects without a as_dict() method
    are replaced with None in the serialized data.
    """
    as_dict = getattr(obj, "as_dict", None)  # all MSONable objects implement as_dict()
    # objects without as_dict() serialize to None, which drops e.g. non-serializable
    # AseAtoms from M3GNet relaxation trajectories
    return as_dict() if callable(as_dict) else None


def glob_to_df(
    pattern: str,
    *,
    reader: Callable[..., pd.DataFrame] | None = None,
    pbar: bool = True,
    **kwargs: object,
) -> pd.DataFrame:
    """Combine data files matching a glob pattern into a single dataframe.

    Args:
        pattern (str): Glob file pattern.
        reader (Callable[..., pd.DataFrame], optional): Function that loads data from
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
        # dict value type erases pandas' concrete read_csv/read_json overloads, keeping
        # reader as the broad Callable[..., DataFrame] so it can forward **kwargs
        ext_to_reader: dict[str, Callable[..., pd.DataFrame]] = {
            ".csv": pd.read_csv,
            ".json": pd.read_json,
        }
        reader = next(
            (fn for ext, fn in ext_to_reader.items() if ext in pattern.lower()), None
        )
        if reader is None:
            raise ValueError(f"Unsupported file extension in {pattern=}")

    files = glob(pattern)

    if len(files) == 0:
        # load mocked model predictions when running pytest (just first 500 lines
        # from MACE-MPA-0 WBM energy preds)
        if "pytest" in sys.modules or "CI" in os.environ:
            df_mock = pd.read_csv(f"{TEST_FILES}/mock-wbm-energy-preds.csv.gz")
            if (
                "e_form_per_atom" not in df_mock
                and str(Key.formation_energy_per_atom) in df_mock
            ):
                df_mock["e_form_per_atom"] = df_mock[str(Key.formation_energy_per_atom)]
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
    file_filter: Callable[[str, int], bool] = lambda filename, _idx: filename.endswith(
        ".extxyz"
    ),
    filename_to_info: bool = False,
    limit: int | slice | None = None,
) -> list[Atoms]:
    """Read ASE Atoms objects from a ZIP file containing extXYZ files.

    Args:
        zip_filename (str): Path to the ZIP file.
        file_filter (Callable[[str, int], bool], optional): Function to check if a file
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
        filenames = zip_file.namelist()
        if limit is not None:
            slice_lim = slice(limit) if isinstance(limit, int) else limit
            filenames = filenames[slice_lim]

        for idx, filename in tqdm(enumerate(filenames), desc=desc, mininterval=5):
            if not file_filter(filename, idx):
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

# formation-energy predictions further than this from DFT (eV/atom) are treated as
# unrealistic outliers and masked out of all downstream metrics
MAX_E_FORM_ERROR_THRESHOLD = 5.0


def resolve_discovery_pred_col(
    df_preds: pd.DataFrame,
    *,
    path: str,
    preferred: str | None = None,
) -> str:
    """Return the formation-energy prediction column from a discovery CSV.

    Resolution order:
    1. ``preferred`` from YAML ``metrics.discovery.pred_col`` when present in the
       frame (needed for multi-column ensemble artifacts on Figshare)
    2. Canonical ``e_form_per_atom`` (current runner output / rewritten locals)
    3. A unique legacy ``e_form_per_atom_*`` column
    """
    canonical = "e_form_per_atom"
    if preferred and preferred in df_preds.columns:
        return preferred
    if canonical in df_preds.columns:
        return canonical
    legacy = [col for col in df_preds.columns if str(col).startswith(f"{canonical}_")]
    if len(legacy) == 1:
        return legacy[0]
    raise ValueError(f"{canonical} column not found in {path}")


def load_df_wbm_with_preds(
    *,
    models: Sequence[str | Model] | None = None,
    pbar: bool = True,
    id_col: str = Key.mat_id,
    subset: pd.Index | Sequence[str] | TestSubset | None = None,
    max_error_threshold: float | None = MAX_E_FORM_ERROR_THRESHOLD,
    nrows: int | None = None,
) -> pd.DataFrame:
    """Load WBM summary dataframe with model predictions from disk.

    Args:
        models (Sequence[str | Model] | None, optional): Models to load, given as
            Model members, enum names, or display labels. Defaults to active models
            when None. Pass an empty sequence to load no model predictions.
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
        nrows (int | None, optional): Only read the first nrows of each model's
            prediction file (passed to the pandas reader). Useful to speed up tests.
            Defaults to None, i.e. read all rows.

    Raises:
        ValueError: On unknown model names.

    Returns:
        pd.DataFrame: WBM summary dataframe with model predictions.
    """
    if models is None:
        models_to_load = Model.active()
    else:
        models_to_load = tuple(map(Model.from_ref, models))

    if max_error_threshold is not None and max_error_threshold < 0:
        raise ValueError(f"{max_error_threshold=} must be a positive number")

    model_name = ""
    df_out = df_wbm.copy()

    try:
        prog_bar = tqdm(models_to_load, disable=not pbar, desc="Loading preds")
        for model in prog_bar:
            model_name = model.name
            prog_bar.set_postfix_str(model_name)

            df_preds = glob_to_df(model.discovery_path, pbar=False, nrows=nrows)
            discovery_metrics = model.metrics.get("discovery", {})
            preferred = (
                discovery_metrics.get("pred_col")
                if isinstance(discovery_metrics, dict)
                else None
            )
            pred_col = resolve_discovery_pred_col(
                df_preds,
                path=model.discovery_path,
                preferred=preferred if isinstance(preferred, str) else None,
            )

            index_column = "material_id" if "material_id" in df_preds else id_col
            df_out[model.label] = df_preds.set_index(index_column)[pred_col]
            if max_error_threshold is not None:
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


def update_yaml_file(
    file_path: str | Path,
    dotted_path: str,
    data: dict[str, Any] | Callable[[dict[str, Any]], dict[str, Any]],
    *,
    preserve_existing: bool = True,
) -> dict[str, Any]:
    """Update a YAML file at a specific dotted path with new data.

    Uses file locking to prevent race conditions when multiple processes
    try to update the same file simultaneously.

    Args:
        file_path (str | Path): Path to YAML file to update
        dotted_path (str): Dotted path to update (e.g. 'metrics.discovery')
        data (dict | Callable): Data to write, or a locked section transformer.
        preserve_existing (bool): If True (default), keep existing keys in the target
            dict section that aren't in `data`. If False, replace the section entirely.

    Returns:
        dict[str, Any]: The complete updated YAML data written to file.

    Example:
        update_yaml_file(
            "models/mace/mace-mp-0.yml",
            "metrics.discovery",
            dict(mae=0.1, rmse=0.2),
        )
    """
    # raise on repeated or trailing dots in dotted path
    if not re.match(r"^[a-zA-Z0-9-+=_]+(\.[a-zA-Z0-9-+=_]+)*$", dotted_path):
        raise ValueError(f"Invalid {dotted_path=}")

    # Use a lock file to prevent race conditions
    lock_path = f"{file_path}.lock"
    with FileLock(lock_path):
        with open(file_path, encoding="utf-8") as file:
            yaml_data = round_trip_yaml.load(file)

        # Navigate to the correct nested level
        current = yaml_data
        *parts, last = dotted_path.split(".")

        for part in parts:
            # replace missing or scalar intermediate nodes (e.g. 'not available'
            # placeholders) with dicts so traversal can continue
            if not isinstance(current.get(part), dict):
                current[part] = {}
            current = current[part]

        # Update the data at the final level. By default, preserve existing keys when
        # replacing a dict section (and replace placeholder strings like 'not
        # available'). Pass preserve_existing=False to fully replace the section, so a
        # recompute drops keys that are no longer emitted (deprecated metrics).
        previous = current.get(last)
        # Callables own the merge (they receive a copy of the prior section). Plain
        # dict updates optionally keep unspecified prior keys via preserve_existing.
        if isinstance(data, dict):
            updated_data = data.copy()
            if preserve_existing and isinstance(previous, dict):
                for key, val in previous.items():
                    updated_data.setdefault(key, val)
        else:
            updated_data = data(dict(previous) if isinstance(previous, dict) else {})
        current[last] = updated_data

        # Write back to file
        with open(file_path, mode="w", encoding="utf-8") as file:
            round_trip_yaml.dump(yaml_data, file)

        return yaml_data
