"""Tests for data loading and IO helpers."""

import os
import sys
import zipfile
from datetime import date
from pathlib import Path
from typing import Any, cast
from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest
from ase import Atoms
from pymatviz.enums import Key
from ruamel.yaml.comments import CommentedMap

from matbench_discovery.data import (
    artifact_filename,
    as_dict_handler,
    ase_atoms_from_zip,
    ase_atoms_to_zip,
    canonical_scientific_notation,
    df_wbm,
    file_ref_name,
    file_ref_url,
    glob_to_df,
    iter_file_refs,
    load_df_wbm_with_preds,
    make_file_ref,
    parse_artifact_filename,
    round_trip_yaml,
    task_coverage,
    update_yaml_file,
)
from matbench_discovery.enums import MbdKey, Model, TestSubset

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


@pytest.fixture
def dummy_atoms_zip(tmp_path: Path) -> Path:
    """ZIP containing the module-level dummy_atoms structures."""
    zip_path = tmp_path / "test_structures.zip"
    ase_atoms_to_zip(dummy_atoms, zip_path)
    return zip_path


def test_as_dict_handler() -> None:
    class C:
        @staticmethod
        def as_dict() -> dict[str, Any]:
            return {"foo": "bar"}

    assert as_dict_handler(C()) == {"foo": "bar"}
    assert as_dict_handler(1) is None
    assert as_dict_handler("foo") is None
    assert as_dict_handler([1, 2, 3]) is None
    assert as_dict_handler({"foo": "bar"}) is None


def test_df_wbm() -> None:
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
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    os.makedirs(tmp_path, exist_ok=True)
    df_mixed.to_csv(f"{tmp_path}/dummy_df.csv", index=False)
    df_mixed.to_json(f"{tmp_path}/dummy_df.json")

    df_out = glob_to_df(f"{tmp_path}/{pattern}")
    assert df_out.shape == df_mixed.shape
    assert list(df_out) == list(df_mixed)

    with pytest.raises(ValueError, match="Unsupported file extension in pattern='foo'"):
        glob_to_df("foo")

    # Mock sys.modules without pytest to test file not found error
    mock_modules = dict(sys.modules)
    mock_modules.pop("pytest", None)  # remove pytest since glob_to_df returns mock data
    # if if finds pytest imported
    # also remove CI from os.environ
    monkeypatch.delenv("CI", raising=False)
    with (
        patch("sys.modules", mock_modules),
        pytest.raises(FileNotFoundError, match="No files matching glob pattern="),
    ):
        glob_to_df("foo.csv")


@pytest.mark.parametrize(
    "atoms_payload",
    [dummy_atoms, {"atoms1": atoms1, "atoms2": atoms2}],
)
def test_atoms_zip_round_trip(
    tmp_path: Path, atoms_payload: dict[str, Atoms] | list[Atoms]
) -> None:
    """List and dict ASE payloads round-trip through zip without info drift."""
    zip_path = tmp_path / "test_structures.zip"
    ase_atoms_to_zip(atoms_payload, zip_path)
    read_atoms = ase_atoms_from_zip(zip_path)
    assert len(read_atoms) == len(atoms_payload)

    orig_atoms = (
        atoms_payload.values() if isinstance(atoms_payload, dict) else atoms_payload
    )
    for original, read in zip(orig_atoms, read_atoms, strict=True):
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


def test_ase_atoms_from_zip_read_options(dummy_atoms_zip: Path) -> None:
    """file_filter and filename_to_info options select / annotate zip members."""
    filtered = ase_atoms_from_zip(
        dummy_atoms_zip, file_filter=lambda name, _idx: "1" in name
    )
    assert len(filtered) == 1
    assert filtered[0].get_chemical_formula() == "H2O"

    with_names = ase_atoms_from_zip(dummy_atoms_zip, filename_to_info=True)
    assert [atoms.info["filename"] for atoms in with_names] == [
        "structure1.extxyz",
        "structure2.extxyz",
    ]


def test_ase_atoms_from_zip_bad_inputs(tmp_path: Path) -> None:
    """Empty zips return []; corrupt zips raise BadZipFile."""
    empty_zip = tmp_path / "empty.zip"
    with zipfile.ZipFile(empty_zip, mode="w"):
        pass
    assert ase_atoms_from_zip(empty_zip) == []

    invalid_zip = tmp_path / "invalid.zip"
    invalid_zip.write_text("This is not a zip file")
    with pytest.raises(zipfile.BadZipFile):
        ase_atoms_from_zip(invalid_zip)


@pytest.mark.parametrize(
    ("limit", "expected_formulas"),
    [
        (1, ["H2O"]),
        (2, ["H2O", "CO2"]),
        (None, ["H2O", "CO2"]),
        (10, ["H2O", "CO2"]),
    ],
)
def test_ase_atoms_from_zip_with_limit(
    dummy_atoms_zip: Path, limit: int | None, expected_formulas: list[str]
) -> None:
    """Integer/None limits truncate or return the full zip contents."""
    read_atoms = ase_atoms_from_zip(dummy_atoms_zip, limit=limit)
    assert [atoms.get_chemical_formula() for atoms in read_atoms] == expected_formulas


@pytest.mark.parametrize(
    ("slice_limit", "expected_formulas"),
    [
        (slice(1, 3), ["N2", "O2"]),
        (slice(None, 2), ["H2", "N2"]),
        (slice(3, None), ["F2", "Cl2"]),
        (slice(0, None, 2), ["H2", "O2", "Cl2"]),
        (slice(None, None, -1), ["Cl2", "F2", "O2", "N2", "H2"]),
        (slice(5, 5), []),
        (slice(10, 20), []),
    ],
)
def test_ase_atoms_from_zip_with_slice_limit(
    tmp_path: Path, slice_limit: slice, expected_formulas: list[str]
) -> None:
    """Slice limits select zip members by archive order."""
    more_atoms = [
        Atoms("H2", positions=[[0, 0, 0], [0, 0, 1]], cell=[4, 4, 4]),
        Atoms("N2", positions=[[0, 0, 0], [0, 0, 1.1]], cell=[4, 4, 4]),
        Atoms("O2", positions=[[0, 0, 0], [0, 0, 1.2]], cell=[4, 4, 4]),
        Atoms("F2", positions=[[0, 0, 0], [0, 0, 1.4]], cell=[4, 4, 4]),
        Atoms("Cl2", positions=[[0, 0, 0], [0, 0, 2.0]], cell=[5, 5, 5]),
    ]
    for idx, atoms in enumerate(more_atoms):
        atoms.info[Key.mat_id] = f"slice_test_{idx}"

    zip_path = tmp_path / "slice_test_structures.zip"
    ase_atoms_to_zip(more_atoms, zip_path)
    read_atoms = ase_atoms_from_zip(zip_path, limit=slice_limit)
    assert [atoms.get_chemical_formula() for atoms in read_atoms] == expected_formulas


@pytest.mark.skipif(
    "CI" in os.environ,
    reason="CI uses mock data so don't check length against on-the-fly "
    "downloaded actual df_wbm",
)
@pytest.mark.parametrize("models", [[], ["wrenformer"]])
@pytest.mark.parametrize("max_error_threshold", [None, 5.0, 1.0])
def test_load_df_wbm_with_preds(
    models: list[str], max_error_threshold: float | None
) -> None:
    df_wbm_with_preds = load_df_wbm_with_preds(
        models=models, max_error_threshold=max_error_threshold
    )
    assert len(df_wbm_with_preds) == len(df_wbm)

    assert list(df_wbm_with_preds) == list(df_wbm) + [
        Model[model].label for model in models
    ]
    assert df_wbm_with_preds.index.name == Key.mat_id

    for model_name in models:
        model = Model[model_name]
        assert model.label in df_wbm_with_preds
        if max_error_threshold is not None:
            # Check if predictions exceeding the threshold are filtered out
            error = abs(
                df_wbm_with_preds[model.label] - df_wbm_with_preds[MbdKey.e_form_dft]
            )
            assert np.all(error[~error.isna()] <= max_error_threshold)
        else:
            # If no threshold is set, all predictions should be present
            assert df_wbm_with_preds[model.label].isna().sum() == 0


def test_load_df_wbm_with_preds_mock_data_models() -> None:
    """Test default and explicit model loading with pytest mock data."""
    inactive_model = Model.alphanet_v1_mptrj
    inactive_model_refs = (inactive_model, inactive_model.name, inactive_model.label)
    with patch("matbench_discovery.data.glob", return_value=[]):
        df_default = load_df_wbm_with_preds(pbar=False)
        inactive_cols = [
            list(load_df_wbm_with_preds(models=[model_ref], pbar=False))
            for model_ref in inactive_model_refs
        ]

    default_cols = list(df_default)
    assert default_cols == [*df_wbm, *(model.label for model in Model.active())]
    assert set(default_cols).isdisjoint(
        model.label for model in Model if not model.is_active
    )
    assert inactive_cols == [[*df_wbm, inactive_model.label]] * len(inactive_model_refs)


@pytest.mark.skipif(
    "CI" in os.environ, reason="CI uses mock data where other error thresholds apply"
)
def test_load_df_wbm_max_error_threshold() -> None:
    models: dict[str, int] = {  # map model to number of max allowed missing preds
        Model.mace_mp_0.label: 38  # before error is raised
    }
    df_no_thresh = load_df_wbm_with_preds(models=list(models))
    df_high_thresh = load_df_wbm_with_preds(models=list(models), max_error_threshold=10)
    df_low_thresh = load_df_wbm_with_preds(models=list(models), max_error_threshold=0.1)

    for model, n_missing in models.items():
        assert df_no_thresh[model].isna().sum() == n_missing
        assert df_high_thresh[model].isna().sum() <= df_no_thresh[model].isna().sum()
        assert df_high_thresh[model].isna().sum() <= df_low_thresh[model].isna().sum()


def test_load_df_wbm_with_preds_errors(df_float: pd.DataFrame) -> None:
    """Test error handling in load_df_wbm_with_preds function."""
    # Test invalid model name
    with pytest.raises(ValueError, match="not found in Model"):
        load_df_wbm_with_preds(models=["InvalidModel"])

    # Test negative error threshold
    with pytest.raises(
        ValueError, match="max_error_threshold=-1 must be a positive number"
    ):
        load_df_wbm_with_preds(max_error_threshold=-1)

    # Test missing canonical prediction column
    with (
        # Make glob return a non-empty list to skip the mock data loading path
        patch("matbench_discovery.data.glob", return_value=["dummy_file.csv"]),
        # Patch only the specific read_csv call in glob_to_df that loads predictions,
        # not the one that loads mock data
        patch("pandas.read_csv", return_value=df_float),
        pytest.raises(ValueError, match=r"e_form_per_atom column not found in"),
    ):
        load_df_wbm_with_preds(models=["alignn"])


@pytest.mark.parametrize(
    "subset",
    ["unique_prototypes", TestSubset.uniq_protos, ["wbm-1-1", "wbm-1-2"], None],
)
def test_load_df_wbm_with_preds_subset(
    subset: str | TestSubset | list[str] | None,
) -> None:
    """Test subset handling in load_df_wbm_with_preds."""
    df_wbm = load_df_wbm_with_preds(subset=subset)
    assert isinstance(df_wbm, pd.DataFrame)


def test_update_yaml_file(tmp_path: Path) -> None:
    """Update YAML at dotted paths; preserve comments; callables own the merge."""
    test_file = f"{tmp_path}/test.yml"

    initial_data = {"metrics": {"discovery": {"mae": 0.1, "pred_file": "old.csv"}}}
    with open(test_file, mode="w") as file:
        round_trip_yaml.dump(initial_data, file)

    update_data = {"mae": 0.2, "rmse": 0.3}
    updated_yaml = update_yaml_file(test_file, "metrics.discovery", update_data)
    result = updated_yaml["metrics"]["discovery"]
    assert result["mae"] == 0.2
    assert result["rmse"] == 0.3
    assert result["pred_file"] == "old.csv"
    assert update_data == {"mae": 0.2, "rmse": 0.3}

    updated_yaml = update_yaml_file(test_file, "metrics.new.nested.path", {"value": 42})
    assert updated_yaml["metrics"]["new"]["nested"]["path"] == {"value": 42}

    yaml_with_comments = """
metrics:
  discovery:  # Discovery metrics
    mae: 0.1  # Mean absolute error
    rmse: 0.2  # Root mean squared error
"""
    with open(test_file, mode="w") as file:
        file.write(yaml_with_comments)

    updated_yaml = update_yaml_file(
        test_file, "metrics.discovery", {"mae": 0.3, "rmse": 0.4}
    )
    assert updated_yaml["metrics"]["discovery"] == {"mae": 0.3, "rmse": 0.4}

    with open(test_file) as file:
        content = file.read()
    assert "discovery:  # Discovery metrics\n    mae: 0.3" in content, f"{content=}"

    commented_data = CommentedMap({"value": 1})
    commented_data.yaml_add_eol_comment("A comment", "value")
    updated_yaml = update_yaml_file(test_file, "new.path", commented_data)
    assert updated_yaml["new"]["path"]["value"] == 1
    with open(test_file) as file:
        content = file.read()
    assert "discovery:  # Discovery metrics\n    mae: 0.3" in content, f"{content=}"
    assert "value: 1  # A comment" in content, f"{content=}"

    # Callables ignore preserve_existing and may drop unspecified prior keys.
    with open(test_file, mode="w") as file:
        round_trip_yaml.dump(
            {"metrics": {"discovery": {"mae": 0.1, "pred_file": "old.csv"}}}, file
        )

    def replace_mae_only(section: dict[str, object]) -> dict[str, object]:
        """Return only the updated field (drops unspecified prior keys)."""
        assert section["pred_file"] == "old.csv"
        return {"mae": 0.9}

    updated = update_yaml_file(
        test_file, "metrics.discovery", replace_mae_only, preserve_existing=True
    )
    assert updated["metrics"]["discovery"] == {"mae": 0.9}

    # Dict updates honor preserve_existing.
    with open(test_file, mode="w") as file:
        round_trip_yaml.dump(
            {"metrics": {"discovery": {"mae": 0.1, "pred_file": "old.csv"}}}, file
        )
    updated = update_yaml_file(
        test_file, "metrics.discovery", {"mae": 0.5}, preserve_existing=True
    )
    assert updated["metrics"]["discovery"] == {"mae": 0.5, "pred_file": "old.csv"}
    updated = update_yaml_file(
        test_file, "metrics.discovery", {"mae": 0.6}, preserve_existing=False
    )
    assert updated["metrics"]["discovery"] == {"mae": 0.6}

    with pytest.raises(FileNotFoundError):
        update_yaml_file("non-existent.yml", "path", {"data": 1})

    for path in ("metrics..discovery", "metrics..", "metrics.discovery..", "."):
        with pytest.raises(ValueError, match="Invalid dotted_path="):
            update_yaml_file(test_file, path, {"data": 1})


# --- model artifact filenames / FileRef helpers ---


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
        canonical_scientific_notation(cast("Any", value))


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
    path = "models/mace/mace-mp-0/2026-07-01-discovery.csv.gz"
    assert make_file_ref(path) == {"name": path}
    assert make_file_ref(
        path, url="https://figshare.com/files/1", size=10, md5="a" * 32
    ) == {
        "name": path,
        "url": "https://figshare.com/files/1",
        "size": 10,
        "md5": "a" * 32,
    }
    for size, md5 in ((10, None), (None, "a" * 32)):
        with pytest.raises(ValueError, match="size and md5"):
            make_file_ref(path, size=size, md5=md5)


def test_file_ref_accessors() -> None:
    """Nested file refs expose names/URLs and can be traversed."""
    nested = {
        "name": "models/mace/mace-mp-0/2026-07-01-discovery.csv.gz",
        "url": "https://figshare.com/files/1",
    }
    assert file_ref_name(nested) == nested["name"]
    assert file_ref_url(nested) == nested["url"]
    assert file_ref_name("legacy/path.csv.gz") is None
    assert list(iter_file_refs({"metrics": {"pred_file": nested}})) == [
        (("metrics", "pred_file"), nested["name"])
    ]
    with pytest.raises(ValueError, match=r"Invalid FileRef at metrics\.pred_file"):
        list(iter_file_refs({"metrics": {"pred_file": "legacy/path.csv.gz"}}))


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
