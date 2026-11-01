"""Tests for data loading and IO helpers."""

import io
import os
import tempfile
import zipfile
from collections.abc import Iterator
from datetime import date
from pathlib import Path
from typing import Any, Self
from unittest.mock import Mock, patch

import ase.io
import numpy as np
import pandas as pd
import pytest
from ase import Atoms
from pymatviz.enums import Key
from ruamel.yaml.comments import CommentedMap

from matbench_discovery import TEST_FILES, preds
from matbench_discovery import data as data_module
from matbench_discovery.data import (
    artifact_filename,
    ase_atoms_from_zip,
    canonical_scientific_notation,
    df_wbm,
    file_ref_name,
    file_ref_url,
    glob_to_df,
    iter_file_refs,
    load_df_wbm_with_preds,
    make_file_ref,
    merge_file_ref,
    parse_artifact_filename,
    round_trip_yaml,
    task_coverage,
    update_yaml_file,
)
from matbench_discovery.enums import MbdKey, Model, TestSubset
from matbench_discovery.preds import elements

atoms1 = Atoms("H2O", positions=[[0, 0, 0], [0, 0, 1], [0, 1, 0]], cell=[5, 5, 5])
atoms1.info[Key.mat_id] = "structure1"
atoms1.info["formation_energy"] = 1.23
atoms1.info["forces"] = np.array([[0, 0, 1], [0, 1, 0], [1, 0, 0]])
atoms1.info["dict_data"] = {"a": 1, "b": 2}
atoms1.info["list_data"] = [1, 2, 3]

atoms2 = Atoms("CO2", positions=[[0, 0, 0], [0, 0, 1], [0, 0, -1]], cell=[6, 6, 6])
atoms2.info[Key.mat_id] = "structure2"
atoms2.info["bandgap"] = 2.34
atoms2.info["magmoms"] = np.array([0.1, -0.1, 0.1])

dummy_atoms = [atoms1, atoms2]


def write_atoms_zip(atoms_list: list[Atoms], zip_path: Path) -> None:
    """Write each Atoms object to its own extXYZ member named after its material ID."""
    with zipfile.ZipFile(zip_path, mode="w") as zip_file:
        for atoms in atoms_list:
            buffer = io.StringIO()
            ase.io.write(buffer, atoms, format="extxyz", write_info=True)
            zip_file.writestr(f"{atoms.info[Key.mat_id]}.extxyz", buffer.getvalue())


@pytest.fixture
def dummy_atoms_zip(tmp_path: Path) -> Path:
    """ZIP containing the module-level dummy_atoms structures."""
    zip_path = tmp_path / "test_structures.zip"
    write_atoms_zip(dummy_atoms, zip_path)
    return zip_path


def test_df_wbm() -> None:
    """The canonical WBM dataframe has its expected shape and core columns."""
    assert df_wbm.shape == (256_963, 18)
    assert df_wbm.index.name == Key.mat_id
    assert set(df_wbm) > {Key.formula, Key.mat_id, Key.bandgap_pbe}

    for col in (
        MbdKey.e_form_dft,
        MbdKey.each_true,
        MbdKey.init_protostructure_spglib,
        MbdKey.protostructure_spglib,
    ):
        assert col in df_wbm, f"{col=} not in {list(df_wbm)=}"


@pytest.mark.parametrize("pattern", ["*df.csv", "*df.json"])
def test_glob_to_df(
    pattern: str,
    tmp_path: Path,
    df_mixed: pd.DataFrame,
) -> None:
    """CSV and JSON glob patterns load dataframes with matching columns."""
    df_mixed.to_csv(f"{tmp_path}/dummy_df.csv", index=False)
    df_mixed.to_json(f"{tmp_path}/dummy_df.json")

    df_loaded = glob_to_df(f"{tmp_path}/{pattern}")
    assert df_loaded.shape == df_mixed.shape
    assert list(df_loaded) == list(df_mixed)


def test_glob_to_df_errors(tmp_path: Path) -> None:
    """Unsupported patterns and absent files raise their specific errors."""
    with pytest.raises(ValueError, match="Unsupported file extension in pattern='foo'"):
        glob_to_df("foo")
    with pytest.raises(FileNotFoundError, match="No files matching glob pattern="):
        glob_to_df(f"{tmp_path}/foo.csv")


@pytest.fixture
def mock_wbm_preds() -> Iterator[Mock]:
    """Serve 499 MACE-MPA-0 WBM preds for every model without resolving real files.

    Yields the resolve_verified_file mock so tests can count path resolutions.
    """
    df_mock = pd.read_csv(f"{TEST_FILES}/mock-wbm-energy-preds.csv.gz").rename(
        columns={str(Key.formation_energy_per_atom): "e_form_per_atom"}
    )
    with (
        patch(
            "matbench_discovery.enums.resolve_verified_file",
            side_effect=lambda abs_path, **_kwargs: abs_path,
        ) as resolve_file,
        patch.object(
            data_module.pd, "read_csv", side_effect=lambda *_a, **_kw: df_mock.copy()
        ),
    ):
        yield resolve_file


def test_atoms_zip_round_trip(dummy_atoms_zip: Path) -> None:
    """Zipped extXYZ structures load back without Atoms.info drift."""
    read_atoms = ase_atoms_from_zip(str(dummy_atoms_zip))
    assert len(read_atoms) == len(dummy_atoms)

    for original, read in zip(dummy_atoms, read_atoms, strict=True):
        assert original.get_chemical_formula() == read.get_chemical_formula()
        assert np.allclose(original.get_positions(), read.get_positions())
        assert np.allclose(original.get_cell(), read.get_cell())
        assert np.all(original.pbc == read.pbc)
        assert set(original.info) == set(read.info)
        for key, value in original.info.items():
            read_val = read.info[key]
            if isinstance(value, np.ndarray) or isinstance(read_val, np.ndarray):
                assert np.allclose(value, read_val), f"Mismatch in {key}"
            else:
                assert value == read_val, f"Mismatch in {key}"


def test_ase_atoms_from_zip_skips_non_extxyz_members(dummy_atoms_zip: Path) -> None:
    """Archive members without an .extxyz suffix are ignored."""
    with zipfile.ZipFile(dummy_atoms_zip, mode="a") as zip_file:
        zip_file.writestr("readme.txt", "not a structure")
    read_atoms = ase_atoms_from_zip(str(dummy_atoms_zip))
    assert [atoms.get_chemical_formula() for atoms in read_atoms] == ["H2O", "CO2"]


def test_ase_atoms_from_zip_bad_inputs(tmp_path: Path) -> None:
    """Empty zips return []; corrupt zips raise BadZipFile."""
    empty_zip = tmp_path / "empty.zip"
    with zipfile.ZipFile(empty_zip, mode="w"):
        pass
    assert ase_atoms_from_zip(str(empty_zip)) == []

    invalid_zip = tmp_path / "invalid.zip"
    invalid_zip.write_text("This is not a zip file")
    with pytest.raises(zipfile.BadZipFile):
        ase_atoms_from_zip(str(invalid_zip))


@pytest.mark.parametrize(
    ("limit", "expected_formulas"),
    [
        (1, ["H2"]),
        (2, ["H2", "N2"]),
        (None, ["H2", "N2", "O2", "F2", "Cl2"]),
        (10, ["H2", "N2", "O2", "F2", "Cl2"]),
        (slice(1, 3), ["N2", "O2"]),
        (slice(3, None), ["F2", "Cl2"]),
        (slice(0, None, 2), ["H2", "O2", "Cl2"]),
        (slice(None, None, -1), ["Cl2", "F2", "O2", "N2", "H2"]),
        (slice(5, 5), []),
        (slice(10, 20), []),
    ],
)
def test_ase_atoms_from_zip_with_limit(
    tmp_path: Path, limit: int | slice | None, expected_formulas: list[str]
) -> None:
    """Integer and slice limits select zip members by archive order."""
    bond_lengths = {"H2": 1.0, "N2": 1.1, "O2": 1.2, "F2": 1.4, "Cl2": 2.0}
    more_atoms = []
    for idx, (formula, bond_length) in enumerate(bond_lengths.items()):
        atoms = Atoms(formula, positions=[[0, 0, 0], [0, 0, bond_length]], cell=[5] * 3)
        atoms.info[Key.mat_id] = f"limit_test_{idx}"
        more_atoms.append(atoms)

    zip_path = tmp_path / "limit_test_structures.zip"
    write_atoms_zip(more_atoms, zip_path)
    read_atoms = ase_atoms_from_zip(str(zip_path), limit=limit)
    assert [atoms.get_chemical_formula() for atoms in read_atoms] == expected_formulas


@pytest.mark.parametrize(
    ("models", "max_error_threshold"),
    [
        ([], None),
        (["wrenformer"], None),
        (["wrenformer"], 5.0),
        (["wrenformer"], 0.1),
    ],
)
def test_load_df_wbm_with_preds(
    models: list[str], max_error_threshold: float | None, mock_wbm_preds: Mock
) -> None:
    """Model selection and outlier thresholds preserve the canonical WBM rows."""
    df_wbm_with_preds = load_df_wbm_with_preds(
        models=models, max_error_threshold=max_error_threshold, pbar=False
    )
    assert mock_wbm_preds.call_count == len(models)
    assert len(df_wbm_with_preds) == len(df_wbm)

    assert list(df_wbm_with_preds) == list(df_wbm) + [
        Model[model].key for model in models
    ]
    assert df_wbm_with_preds.index.name == Key.mat_id

    for model_name in models:
        model = Model[model_name]
        if max_error_threshold is not None:
            error = abs(
                df_wbm_with_preds[model.key] - df_wbm_with_preds[MbdKey.e_form_dft]
            )
            n_kept = df_wbm_with_preds[model.key].notna().sum()
            assert np.all(error[~error.isna()] <= max_error_threshold)
            assert 0 < n_kept < 499 if max_error_threshold < 1 else n_kept == 499
        else:
            assert df_wbm_with_preds[model.key].notna().sum() == 499


def test_load_df_wbm_with_preds_mock_data_models(mock_wbm_preds: Mock) -> None:
    """Default loading selects active models; inactive ones load by any ref."""
    inactive_model = Model["alphanet_v1_mptrj"]
    inactive_model_refs = (inactive_model, inactive_model.name, inactive_model.label)
    df_default = load_df_wbm_with_preds(pbar=False)
    inactive_cols = [
        list(load_df_wbm_with_preds(models=[model_ref], pbar=False))
        for model_ref in inactive_model_refs
    ]

    default_cols = list(df_default)
    assert default_cols == [*df_wbm, *(model.key for model in Model.active())]
    assert set(default_cols).isdisjoint(
        model.key for model in Model if not model.is_active
    )
    assert inactive_cols == [[*df_wbm, inactive_model.key]] * len(inactive_model_refs)
    n_loads = len(Model.active()) + len(inactive_model_refs)
    assert mock_wbm_preds.call_count == n_loads
    assert df_default[Model["mace_mpa_0"].key].notna().any()


@pytest.mark.network
def test_load_df_wbm_max_error_threshold() -> None:
    """Tighter outlier thresholds monotonically increase missing predictions."""
    model = Model["mace_mp_0"]
    df_no_threshold = load_df_wbm_with_preds(models=[model], max_error_threshold=None)
    df_high_threshold = load_df_wbm_with_preds(models=[model], max_error_threshold=10)
    df_low_threshold = load_df_wbm_with_preds(models=[model], max_error_threshold=0.1)

    missing_counts = [
        frame[model.key].isna().sum()
        for frame in (df_no_threshold, df_high_threshold, df_low_threshold)
    ]
    assert missing_counts[0] == 38
    assert missing_counts == sorted(missing_counts)


def test_load_df_wbm_with_preds_errors(df_float: pd.DataFrame) -> None:
    """Test error handling in load_df_wbm_with_preds function."""
    with pytest.raises(ValueError, match="not found in Model"):
        load_df_wbm_with_preds(models=["InvalidModel"])

    with pytest.raises(
        ValueError, match="max_error_threshold=-1 must be a positive number"
    ):
        load_df_wbm_with_preds(max_error_threshold=-1)

    with (
        patch(
            "matbench_discovery.enums.resolve_verified_file",
            side_effect=lambda abs_path, **_kwargs: abs_path,
        ),
        patch.object(data_module.pd, "read_csv", return_value=df_float),
        pytest.raises(ValueError, match=r"e_form_per_atom column not found in"),
    ):
        load_df_wbm_with_preds(models=["alignn"])


@pytest.mark.parametrize(
    "subset",
    [
        "unique_prototypes",
        TestSubset.uniq_protos,
        TestSubset.full_test_set,
        "full_test_set",
        ["wbm-1-1", "wbm-1-2"],
        pd.Index(["wbm-1-1", "wbm-1-2"]),
        None,
    ],
)
def test_load_df_wbm_with_preds_subset(
    subset: str | TestSubset | list[str] | pd.Index | None, mock_wbm_preds: Mock
) -> None:
    """Subset selectors return exactly the requested WBM material IDs."""
    # Default model selection is covered above; one prediction column exercises
    # subsetting without loading every model for each selector.
    model = Model["mace_mpa_0"]
    df_subset = load_df_wbm_with_preds(models=[model], subset=subset, pbar=False)
    assert mock_wbm_preds.call_count == 1
    if subset is None or (isinstance(subset, str) and subset == "full_test_set"):
        expected_index = df_wbm.index
    elif isinstance(subset, list | pd.Index):
        expected_index = df_wbm.loc[subset].index
    else:
        expected_index = df_wbm.query(MbdKey.uniq_proto).index
    assert df_subset.index.equals(expected_index)
    assert list(df_subset) == [*df_wbm, model.key]
    assert df_subset[model.key].notna().any()


def test_prediction_errors_and_element_enrichment(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Prediction errors stay independent from optional per-element enrichment."""
    models = [Model["mace_mp_0"], Model["mace_mpa_0"]]
    raw = pd.DataFrame(
        {
            Key.formula: ["Fe2O3", "FeO"],
            MbdKey.each_true: [0.1234, 0.4567],
            MbdKey.e_form_dft: [-1.2344, -2.3456],
            models[0].key: [-1.1111, -2.1111],
            models[1].key: [-1.0111, -2.0111],
        },
        index=pd.Index(["wbm-1", "wbm-2"], name=Key.mat_id),
    )
    monkeypatch.setattr(preds, "load_df_wbm_with_preds", lambda **_kwargs: raw)

    predictions, errors = preds.load_prediction_errors(models)
    expected_predictions = raw.round(3)
    expected_errors = expected_predictions[[model.key for model in models]].sub(
        expected_predictions[MbdKey.e_form_dft], axis="index"
    )
    pd.testing.assert_frame_equal(predictions, expected_predictions)
    pd.testing.assert_frame_equal(errors, expected_errors)

    monkeypatch.setattr(
        preds.pd,
        "read_json",
        lambda *_args, **_kwargs: pd.Series({"Fe": 10, "O": 20}),
    )
    compositions, element_data = preds.derive_element_data(predictions)
    assert compositions.to_dict("index") == {
        "wbm-1": {"Fe": 2.0, "O": 3.0},
        "wbm-2": {"Fe": 1.0, "O": 1.0},
    }
    assert element_data[preds.TRAIN_COUNT_COL].to_dict() == {"Fe": 10, "O": 20}
    assert element_data[preds.TEST_SET_STD_COL].to_dict() == pytest.approx(
        {"Fe": 0.2361736649163069, "O": 0.2361736649163069}
    )
    assert elements.mean_abs_error_by_element(
        compositions, pd.Series([2.0, float("nan")], index=raw.index)
    ).to_dict() == {"Fe": 2.0, "O": 2.0}


def test_update_yaml_file(tmp_path: Path) -> None:
    """Update YAML at dotted paths; preserve comments; callables own the merge."""
    test_file = tmp_path / "test.yml"
    initial_yaml = "metrics:\n  discovery: {mae: 0.1, pred_file: old.csv}\n"
    test_file.write_text(initial_yaml, encoding="utf-8")

    update_data = {"mae": 0.2, "rmse": 0.3}
    updated_yaml = update_yaml_file(str(test_file), "metrics.discovery", update_data)
    assert updated_yaml["metrics"]["discovery"] == {
        **update_data,
        "pred_file": "old.csv",
    }
    assert update_data == {"mae": 0.2, "rmse": 0.3}

    updated_yaml = update_yaml_file(
        str(test_file), "metrics.new.nested.path", {"value": 42}
    )
    assert updated_yaml["metrics"]["new"]["nested"]["path"] == {"value": 42}

    # null intermediate sections (bare `md:`) are replaced, not indexed into
    test_file.write_text("metrics:\n  md:\n", encoding="utf-8")
    updated_yaml = update_yaml_file(str(test_file), "metrics.md.nested", {"value": 1})
    assert updated_yaml["metrics"] == {"md": {"nested": {"value": 1}}}
    test_file.write_text(initial_yaml, encoding="utf-8")

    # the write lock must not be left next to the target, which for real model YAMLs
    # would leave an untracked .yml.lock sitting in the repo
    assert sorted(os.listdir(tmp_path)) == ["test.yml"]

    yaml_with_comments = """
metrics:
  discovery:  # Discovery metrics
    mae: 0.1  # Mean absolute error
    rmse: 0.2  # Root mean squared error
"""
    test_file.write_text(yaml_with_comments, encoding="utf-8")

    updated_yaml = update_yaml_file(
        str(test_file), "metrics.discovery", {"mae": 0.3, "rmse": 0.4}
    )
    assert updated_yaml["metrics"]["discovery"] == {"mae": 0.3, "rmse": 0.4}

    content = test_file.read_text(encoding="utf-8")
    assert "discovery:  # Discovery metrics\n    mae: 0.3" in content, f"{content=}"
    assert "mae: 0.3  # Mean absolute error" in content
    assert "rmse: 0.4  # Root mean squared error" in content

    commented_data = CommentedMap({"value": 1})
    commented_data.yaml_add_eol_comment("κ (W/mK)", "value")
    updated_yaml = update_yaml_file(str(test_file), "metrics.discovery", commented_data)
    assert updated_yaml["metrics"]["discovery"]["value"] == 1
    content = test_file.read_text(encoding="utf-8")
    assert "discovery:  # Discovery metrics\n    mae: 0.3" in content, f"{content=}"
    assert "value: 1  # κ (W/mK)" in content, f"{content=}"
    assert "mae: 0.3  # Mean absolute error" in content

    # Callables ignore preserve_existing and may drop unspecified prior keys.
    test_file.write_text(initial_yaml, encoding="utf-8")

    def replace_mae_only(section: dict[str, object]) -> dict[str, object]:
        """Return only the updated field (drops unspecified prior keys)."""
        assert section["pred_file"] == "old.csv"
        return {"mae": 0.9}

    updated = update_yaml_file(
        str(test_file), "metrics.discovery", replace_mae_only, preserve_existing=True
    )
    assert updated["metrics"]["discovery"] == {"mae": 0.9}

    # Dict updates honor preserve_existing.
    test_file.write_text(initial_yaml, encoding="utf-8")
    updated = update_yaml_file(
        str(test_file), "metrics.discovery", {"mae": 0.5}, preserve_existing=True
    )
    assert updated["metrics"]["discovery"] == {"mae": 0.5, "pred_file": "old.csv"}
    updated = update_yaml_file(
        str(test_file), "metrics.discovery", {"mae": 0.6}, preserve_existing=False
    )
    assert updated["metrics"]["discovery"] == {"mae": 0.6}

    with pytest.raises(FileNotFoundError):
        update_yaml_file(f"{tmp_path}/non-existent.yml", "path", {"data": 1})

    for path in ("metrics..discovery", "metrics..", "metrics.discovery..", "."):
        with pytest.raises(ValueError, match="Invalid dotted_path="):
            update_yaml_file(str(test_file), path, {"data": 1})


def test_update_yaml_file_lock_survives_differing_tmpdir(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Concurrent writers must share one lock even under different TMPDIR.

    Slurm sets TMPDIR per job and macOS per user, so a tempdir-derived lock path
    hands them separate locks and lets their read-modify-write cycles interleave
    into a corrupt YAML file.
    """
    test_file = f"{tmp_path}/test.yml"
    with open(test_file, mode="w", encoding="utf-8") as file:
        round_trip_yaml.dump({"metrics": {}}, file)

    lock_paths: list[str] = []

    class RecordingLock:
        """Stand-in for FileLock that records the lock path it would take."""

        def __init__(self, lock_file: str, *_args: object, **_kwargs: object) -> None:
            lock_paths.append(str(lock_file))

        def __enter__(self) -> Self:
            return self

        def __exit__(self, *_exc_info: object) -> None:
            return

    monkeypatch.setattr(data_module, "FileLock", RecordingLock)
    for job_dir in ("slurm-job-a", "slurm-job-b"):
        # gettempdir() silently falls back when TMPDIR is missing, so it must exist
        tmpdir = tmp_path / job_dir
        tmpdir.mkdir()
        monkeypatch.setenv("TMPDIR", str(tmpdir))
        monkeypatch.setattr(tempfile, "tempdir", None)  # drop gettempdir()'s cache
        update_yaml_file(test_file, "metrics.discovery", {"mae": 0.1})

    assert lock_paths[0] == lock_paths[1]
    assert tmp_path.as_posix() not in lock_paths[0]  # never beside the target


@pytest.mark.parametrize("use_symlink", [False, True])
@pytest.mark.parametrize("failure_stage", ["dump", "replace"])
def test_update_yaml_file_preserves_original_on_write_failure(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    use_symlink: bool,
    failure_stage: str,
) -> None:
    """Failed YAML writes preserve original bytes, symlinks, and file permissions."""
    target = tmp_path / "model.yml"
    original = "metrics:\n  discovery: {mae: 0.1}\n"
    target.write_text(original, encoding="utf-8")
    target.chmod(0o640)
    original_mode = target.stat().st_mode & 0o777
    path = target
    if use_symlink:
        path = tmp_path / "link.yml"
        path.symlink_to(target)
    expected_files = set(tmp_path.iterdir())

    def fail_dump(_data: object, stream: io.TextIOBase) -> None:
        """Simulate a partially written serialization failure."""
        stream.write("partial YAML")
        raise OSError("dump failed")

    def fail_replace(_source: str, _destination: str) -> None:
        """Simulate an unsuccessful atomic replacement."""
        raise OSError("replace failed")

    with monkeypatch.context() as patch_context:
        if failure_stage == "dump":
            patch_context.setattr(round_trip_yaml, "dump", fail_dump)
        else:
            patch_context.setattr(data_module.os, "replace", fail_replace)
        with pytest.raises(OSError, match=f"{failure_stage} failed"):
            update_yaml_file(str(path), "metrics.discovery", {"mae": 0.2})
    assert target.read_text(encoding="utf-8") == original
    assert path.is_symlink() is use_symlink
    assert set(tmp_path.iterdir()) == expected_files
    update_yaml_file(str(path), "metrics.discovery", {"mae": 0.3})
    with target.open(encoding="utf-8") as file:
        assert round_trip_yaml.load(file)["metrics"]["discovery"]["mae"] == 0.3
    assert target.stat().st_mode & 0o777 == original_mode
    assert path.is_symlink() is use_symlink
    assert set(tmp_path.iterdir()) == expected_files


# --- model artifact filenames / FileRef helpers ---

FILE_REF_PATH = "models/mace/mace-mp-0/2026-07-01-discovery.csv.gz"
UPDATED_FILE_REF_PATH = "models/mace/mace-mp-0/2026-07-02-discovery.csv.gz"
FILE_REF_URL = "https://figshare.com/files/1"


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        (1e-5, "1e-5"),
        ("0.00001", "1e-5"),
        (1e-2, "1e-2"),
        (2.5e-3, "2.5e-3"),
    ],
)
def test_canonical_scientific_notation(value: float | str, expected: str) -> None:
    """Positive finite values render without exponent padding."""
    assert canonical_scientific_notation(value) == expected


@pytest.mark.parametrize("value", [0, -1.0, "not-a-number"])
def test_canonical_scientific_notation_rejects_invalid(value: object) -> None:
    """Non-positive or non-numeric symprec values are rejected."""
    with pytest.raises(
        ValueError,
        match=r"Expected a positive finite number|Invalid numeric",
    ):
        canonical_scientific_notation(value)  # ty: ignore[invalid-argument-type]


@pytest.mark.parametrize(
    ("role", "expected_suffix"),
    [
        ("discovery", "discovery.csv.gz"),
        ("geo_opt", "geo-opt.jsonl.gz"),
        ("md_metrics", "md-metrics.csv.gz"),
        ("diatomics", "diatomics.json.gz"),
    ],
)
def test_artifact_filename_static_roles(role: str, expected_suffix: str) -> None:
    """Static artifact roles render dated canonical filenames."""
    assert artifact_filename("2026-07-01", role) == f"2026-07-01-{expected_suffix}"
    assert parse_artifact_filename(f"2026-07-01-{expected_suffix}") == role


def test_artifact_filename_geo_opt_analysis_round_trip() -> None:
    """Geo-opt analysis filenames round-trip; invalid calendar dates are rejected."""
    filename = artifact_filename(
        date(2026, 7, 1),
        "geo_opt_analysis",
        symprec=1e-5,
        moyo_version="0.4.2",
    )
    assert filename == "2026-07-01-geo-opt-symprec=1e-5-moyo=0.4.2.csv.gz"
    assert parse_artifact_filename(filename) == "geo_opt_analysis"
    with pytest.raises(ValueError, match="Invalid ISO date"):
        parse_artifact_filename("2026-02-30-discovery.csv.gz")


def test_make_file_ref() -> None:
    """Omit unset optionals; reject unpaired size/md5."""
    assert make_file_ref(FILE_REF_PATH) == {"name": FILE_REF_PATH}
    assert make_file_ref(FILE_REF_PATH, url=FILE_REF_URL, size=10, md5="a" * 32) == {
        "name": FILE_REF_PATH,
        "url": FILE_REF_URL,
        "size": 10,
        "md5": "a" * 32,
    }
    for size, md5 in ((10, None), (None, "a" * 32)):
        with pytest.raises(ValueError, match="size and md5"):
            make_file_ref(FILE_REF_PATH, size=size, md5=md5)


@pytest.mark.parametrize(
    ("same_name", "url", "replace", "expected_url", "keep_checksum"),
    [
        (True, None, False, FILE_REF_URL, True),
        (True, FILE_REF_URL, False, FILE_REF_URL, True),
        (
            True,
            "https://figshare.com/files/2",
            False,
            "https://figshare.com/files/2",
            False,
        ),
        (False, None, False, None, False),
        (True, None, True, None, False),
    ],
    ids=("unchanged", "same-url", "new-url", "new-name", "replace"),
)
def test_merge_file_ref(
    same_name: bool,
    url: str | None,
    replace: bool,
    expected_url: str | None,
    keep_checksum: bool,
) -> None:
    """Checksums survive only for an unchanged artifact and upload URL."""
    name = FILE_REF_PATH if same_name else UPDATED_FILE_REF_PATH
    existing = make_file_ref(
        FILE_REF_PATH,
        url=FILE_REF_URL,
        size=10,
        md5="a" * 32,
    )
    expected = make_file_ref(
        name,
        url=expected_url,
        size=10 if keep_checksum else None,
        md5="a" * 32 if keep_checksum else None,
    )
    assert merge_file_ref(existing, name, url=url, replace=replace) == expected


def test_file_ref_accessors() -> None:
    """Nested file refs expose names/URLs and can be traversed."""
    nested = {"name": FILE_REF_PATH, "url": FILE_REF_URL}
    assert file_ref_name(nested) == nested["name"]
    assert file_ref_url(nested) == nested["url"]
    assert file_ref_name("legacy/path.csv.gz") is None
    assert list(iter_file_refs({"metrics": {"pred_file": nested}})) == [
        (("metrics", "pred_file"), nested["name"])
    ]
    # reproducibility run files sit in a list: key paths carry the list index
    runs: list[dict[str, Any]] = [{"label": "a", "file": nested}, {"file": "b.gz"}]
    md_block = {"md": {"pred_file": None, "reproducibility": {"runs": runs}}}
    with pytest.raises(ValueError, match=r"at md\.reproducibility\.runs\.1\.file"):
        list(iter_file_refs(md_block))  # bare string instead of a FileRef
    runs[1]["file"] = {"name": "b.gz"}
    assert list(iter_file_refs(md_block)) == [
        (("md", "reproducibility", "runs", 0, "file"), nested["name"]),
        (("md", "reproducibility", "runs", 1, "file"), "b.gz"),
    ]


@pytest.mark.parametrize(
    ("metadata", "task", "expected_status"),
    [
        (
            {"metrics": {"discovery": {"pred_file": {"name": "x.csv.gz"}}}},
            "discovery",
            "complete",
        ),
        ({"targets": "E"}, "md", "not_applicable"),
        ({"lifecycle": "aborted"}, "discovery", "not_available"),
        ({}, "md", "pending"),
        # bare `metrics:` / `md:` keys load as None
        ({"metrics": None}, "md", "pending"),
        ({"metrics": {"md": None}}, "md", "pending"),
        (
            {"metrics": {"md": {"status": "pending", "reason": "not run yet"}}},
            "md",
            "pending",
        ),
    ],
)
def test_task_coverage_derivation(
    metadata: dict[str, object], task: str, expected_status: str
) -> None:
    """Coverage is derived from metrics, targets, lifecycle, or explicit status."""
    status, _reason = task_coverage(metadata, task)  # type: ignore[arg-type]
    assert status == expected_status
