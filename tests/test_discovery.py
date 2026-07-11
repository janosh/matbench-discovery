"""Tests for the unified calculator-backed WBM discovery pipeline."""

from __future__ import annotations

import copy
import os
from functools import partial
from types import SimpleNamespace
from typing import TYPE_CHECKING, Any, cast

import pandas as pd
import pytest
import yaml
from ase.build import bulk
from ase.calculators.emt import EMT
from pymatgen.core import Lattice, Structure
from pymatgen.entries.compatibility import MaterialsProject2020Compatibility
from pymatgen.entries.computed_entries import ComputedStructureEntry
from pymatviz.enums import Key

import matbench_discovery.discovery as discovery_core
import models.run_discovery as discovery_runner
from matbench_discovery.calculators import CALCULATORS
from matbench_discovery.discovery import (
    DiscoveryArtifacts,
    MergedDiscoveryRun,
    RelaxationRecord,
    RelaxationSettings,
    discovery_pred_col,
    merge_discovery_shards,
    partition_material_ids,
    read_discovery_shard,
    relax_atoms,
    run_discovery_shard,
    write_discovery_artifacts,
)
from matbench_discovery.energy import calc_energy_from_e_refs
from matbench_discovery.enums import MbdKey, Model, TestSubset

if TYPE_CHECKING:
    from pathlib import Path

    from ase import Atoms


def _bulk_atoms(material_id: str, element: str = "Cu", repeat: int = 1) -> Atoms:
    """Build a tagged periodic structure for discovery runner tests."""
    atoms = bulk(element) * (repeat, 1, 1)
    atoms.info[str(Key.mat_id)] = material_id
    return atoms


def _run_test_shard(
    tmp_path: Path,
    *,
    shard_index: int,
    n_shards: int,
    material_id: str,
) -> str:
    """Write one complete EMT shard and return its path."""
    atoms_by_id = {material_id: _bulk_atoms(material_id)}
    shard_path = tmp_path / f"shard-{shard_index:03d}-of-{n_shards:03d}.jsonl"
    run_discovery_shard(
        calculator=EMT(),
        model_key="emt",
        atoms_by_id=atoms_by_id,
        assigned_material_ids=[material_id],
        shard_path=str(shard_path),
        shard_index=shard_index,
        n_shards=n_shards,
        settings=RelaxationSettings(max_steps=0, cell_filter="none"),
        dataset_material_ids=["wbm-1", "wbm-2"] if n_shards == 2 else [material_id],
    )
    return str(shard_path)


def test_relaxation_settings_load_yaml_and_overrides() -> None:
    """Model YAML hyperparameters load exactly and explicit flags take precedence."""
    settings = RelaxationSettings.from_model("mace_mp_0")
    assert settings == RelaxationSettings(
        max_force=0.05,
        max_steps=500,
        ase_optimizer="FIRE",
        cell_filter="FrechetCellFilter",
    )
    overridden = RelaxationSettings.from_model(
        "mace_mp_0",
        max_force=0.01,
        max_steps=7,
        ase_optimizer="GOQN",
        cell_filter="none",
    )
    assert overridden == RelaxationSettings(0.01, 7, "GOQN", None)
    assert RelaxationSettings(cell_filter="NONE") == RelaxationSettings(
        cell_filter=None
    )


@pytest.mark.parametrize("model_key", CALCULATORS)
def test_all_calculators_have_valid_relaxation_settings(model_key: str) -> None:
    """Every shared-runner model has parseable relaxation hyperparameters."""
    assert isinstance(RelaxationSettings.from_model(model_key), RelaxationSettings)


@pytest.mark.parametrize(
    ("kwargs", "match"),
    [
        ({"max_force": 0}, "max_force"),
        ({"max_steps": -1}, "max_steps"),
        ({"ase_optimizer": "unknown"}, "optimizer"),
        ({"cell_filter": "unknown"}, "cell filter"),
    ],
)
def test_relaxation_settings_reject_invalid_values(
    kwargs: dict[str, Any], match: str
) -> None:
    """Invalid optimizer, filter, force, and step settings fail early."""
    with pytest.raises(ValueError, match=match):
        RelaxationSettings(**kwargs)


def test_partition_material_ids_is_deterministic_and_balanced() -> None:
    """Atom balancing is independent of mapping insertion order and deterministic."""
    atoms_by_id = {
        "wbm-3": _bulk_atoms("wbm-3", repeat=3),
        "wbm-1": _bulk_atoms("wbm-1", repeat=1),
        "wbm-4": _bulk_atoms("wbm-4", repeat=4),
        "wbm-2": _bulk_atoms("wbm-2", repeat=2),
    }
    reversed_atoms = dict(reversed(atoms_by_id.items()))
    shards = partition_material_ids(atoms_by_id, 2)
    assert shards == partition_material_ids(reversed_atoms, 2)
    shard_sizes = [
        sum(len(atoms_by_id[material_id]) for material_id in shard) for shard in shards
    ]
    assert shard_sizes == [5, 5]


@pytest.mark.parametrize(("task_id", "should_skip"), [(5, False), (6, True)])
def test_implicit_slurm_dry_run_uses_only_first_task(
    monkeypatch: pytest.MonkeyPatch, task_id: int, should_skip: bool
) -> None:
    """Large Slurm arrays do not over-shard or race the four-structure dry run."""
    monkeypatch.setenv("SLURM_ARRAY_TASK_COUNT", "20")
    monkeypatch.setenv("SLURM_ARRAY_TASK_MIN", "5")
    monkeypatch.setenv("SLURM_ARRAY_TASK_ID", str(task_id))
    assert discovery_runner._effective_shard_args(  # noqa: SLF001
        None, None, dry_run=True
    ) == (1, 0, should_skip)


@pytest.mark.parametrize(
    ("n_shards", "task_min", "task_max", "task_id", "expected"),
    [
        (None, 5, 6, 6, (2, 1)),
        (10, 3, 7, 7, (10, 7)),
        (None, 3, 7, 3, "Sparse Slurm arrays"),
    ],
)
def test_slurm_shard_selection(
    monkeypatch: pytest.MonkeyPatch,
    n_shards: int | None,
    task_min: int,
    task_max: int,
    task_id: int,
    expected: tuple[int, int] | str,
) -> None:
    """Slurm arrays preserve contiguous and explicit sparse shard mappings."""
    monkeypatch.setenv("SLURM_ARRAY_TASK_COUNT", "2")
    monkeypatch.setenv("SLURM_ARRAY_TASK_MIN", str(task_min))
    monkeypatch.setenv("SLURM_ARRAY_TASK_MAX", str(task_max))
    monkeypatch.setenv("SLURM_ARRAY_TASK_ID", str(task_id))
    if isinstance(expected, str):
        with pytest.raises(ValueError, match=expected):
            discovery_runner._slurm_shard_selection(n_shards, None)  # noqa: SLF001
    else:
        assert (
            discovery_runner._slurm_shard_selection(  # noqa: SLF001
                n_shards, None
            )
            == expected
        )


def test_emt_relaxation_smoke() -> None:
    """EMT produces a finite energy, structure, and convergence state."""
    record = relax_atoms(
        _bulk_atoms("wbm-cu"),
        EMT(),
        RelaxationSettings(max_force=0.5, max_steps=2, cell_filter="none"),
    )
    assert record.material_id == "wbm-cu"
    assert record.error is None
    assert record.structure is not None
    assert record.energy is not None
    assert record.n_steps <= 2
    single_point_record = relax_atoms(
        _bulk_atoms("wbm-cu-single-point"),
        EMT(),
        RelaxationSettings(max_steps=0, cell_filter=None),
    )
    assert single_point_record.converged is False


def test_shard_resume_preserves_success_and_failure_records(tmp_path: Path) -> None:
    """Completed and failed attempts are both durable and skipped on resume."""
    atoms_by_id = {
        "wbm-cu": _bulk_atoms("wbm-cu"),
        "wbm-si": _bulk_atoms("wbm-si", element="Si"),
    }
    shard_path = tmp_path / "shard-000-of-001.jsonl"
    settings = RelaxationSettings(max_steps=0, cell_filter="none")
    run_shard = partial(
        run_discovery_shard,
        model_key="emt",
        atoms_by_id=atoms_by_id,
        assigned_material_ids=list(atoms_by_id),
        shard_path=str(shard_path),
        shard_index=0,
        n_shards=1,
        settings=settings,
        dataset_material_ids=list(atoms_by_id),
    )
    first_run = run_shard(calculator=EMT())
    assert [record.material_id for record in first_run.records] == [
        "wbm-cu",
        "wbm-si",
    ]
    assert first_run.records[0].energy is not None
    assert first_run.records[1].error is not None

    size_before_resume = os.path.getsize(shard_path)
    resumed = run_shard(calculator=EMT())
    assert resumed.records == first_run.records
    assert os.path.getsize(shard_path) == size_before_resume


def test_shard_resume_handles_incomplete_events_and_header_change(
    tmp_path: Path,
) -> None:
    """Resume repairs incomplete events but rejects changed settings."""
    material_id = "wbm-1"
    _run_test_shard(tmp_path, shard_index=0, n_shards=1, material_id=material_id)
    shard_path = tmp_path / "shard-000-of-001.jsonl"
    shard_path.write_bytes(shard_path.read_bytes().rstrip(b"\n"))
    discovery_core._append_shard_event(  # noqa: SLF001
        str(shard_path),
        {"type": "run_metadata", "run_metadata": {"run_time_sec": 1}},
    )
    assert len(read_discovery_shard(str(shard_path)).run_metadata_segments) == 2

    with shard_path.open(mode="a", encoding="utf-8") as file:
        file.write('{"type": "record", "material_id":')

    _run_test_shard(tmp_path, shard_index=0, n_shards=1, material_id=material_id)
    assert read_discovery_shard(str(shard_path)).ignored_trailing_line is False

    with pytest.raises(ValueError, match="does not match this run"):
        run_discovery_shard(
            calculator=EMT(),
            model_key="emt",
            atoms_by_id={material_id: _bulk_atoms(material_id)},
            assigned_material_ids=[material_id],
            shard_path=str(shard_path),
            shard_index=0,
            n_shards=1,
            settings=RelaxationSettings(max_steps=1, cell_filter=None),
            dataset_material_ids=[material_id],
        )


def test_merge_rejects_incomplete_duplicate_and_wrong_coverage(
    tmp_path: Path,
) -> None:
    """Strict merge rejects missing shards, duplicate shard IDs, and foreign IDs."""
    shard_0 = _run_test_shard(tmp_path, shard_index=0, n_shards=2, material_id="wbm-1")
    shard_1 = _run_test_shard(tmp_path, shard_index=1, n_shards=2, material_id="wbm-2")
    merge_shards = partial(
        merge_discovery_shards,
        model_key="emt",
        expected_material_ids=["wbm-1", "wbm-2"],
    )
    merged_run = merge_shards([shard_0, shard_1])
    assert [record.material_id for record in merged_run.records] == ["wbm-1", "wbm-2"]
    assert merged_run.n_shards == 2

    incomplete_path = tmp_path / "incomplete" / "shard-000-of-001.jsonl"
    incomplete_atoms = {
        material_id: _bulk_atoms(material_id) for material_id in ("wbm-1", "wbm-2")
    }
    run_discovery_shard(
        calculator=EMT(),
        model_key="emt",
        atoms_by_id=incomplete_atoms,
        assigned_material_ids=list(incomplete_atoms),
        shard_path=str(incomplete_path),
        shard_index=0,
        n_shards=1,
        settings=RelaxationSettings(max_steps=0, cell_filter=None),
        dataset_material_ids=list(incomplete_atoms),
    )
    shard_lines = incomplete_path.read_text().splitlines()
    incomplete_path.write_text("\n".join(shard_lines[:2]) + "\n")
    with pytest.raises(ValueError, match="records do not match assignment"):
        merge_shards([str(incomplete_path)])

    with pytest.raises(ValueError, match="Expected shard indices"):
        merge_shards([shard_0])
    with pytest.raises(ValueError, match="Duplicate shard indices"):
        merge_shards([shard_0, shard_0])
    with pytest.raises(ValueError, match="Inconsistent shard metadata"):
        merge_discovery_shards(
            [shard_0, shard_1],
            model_key="emt",
            expected_material_ids=["wbm-1", "wbm-foreign"],
        )


def test_artifacts_match_legacy_mp2020_join_semantics(tmp_path: Path) -> None:
    """Unified MP2020 formation energies equal the former CSE mutation workflow."""
    material_id = "wbm-1"
    failed_id = "wbm-2"
    structure = Structure(
        Lattice.cubic(3.6),
        ["Cu"] * 4,
        [[0, 0, 0], [0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]],
    )
    raw_energy = -14.0
    record = RelaxationRecord(
        material_id=material_id,
        structure=structure.as_dict(),
        energy=raw_energy,
        converged=True,
        n_steps=4,
    )
    merged_run = MergedDiscoveryRun(
        model_key="emt",
        settings=RelaxationSettings(max_steps=4),
        records=(
            record,
            RelaxationRecord(
                material_id=failed_id,
                structure=None,
                energy=None,
                converged=False,
                n_steps=0,
                error="RuntimeError: failed",
            ),
        ),
        run_metadata={},
        n_shards=1,
    )
    df_summary = pd.DataFrame(
        {
            Key.mat_id: [material_id, failed_id],
            Key.formula: ["Cu4", "Cu"],
        }
    ).set_index(Key.mat_id)
    original_entry = ComputedStructureEntry(
        structure, -8.0, parameters={"run_type": "GGA"}
    )
    df_entries = pd.DataFrame(
        {Key.computed_structure_entry: [original_entry]}, index=[material_id]
    )
    compatibility = MaterialsProject2020Compatibility(check_potcar=False)
    ref_energies = {"Cu": -3.0}

    legacy_entry = copy.deepcopy(original_entry)
    legacy_entry._energy = raw_energy  # noqa: SLF001
    legacy_entry._structure = structure  # noqa: SLF001
    processed_entry = compatibility.process_entries(
        [legacy_entry], clean=True, verbose=False
    )[0]
    expected_e_form = calc_energy_from_e_refs(
        "Cu4", ref_energies, total_energy=processed_entry.energy
    )

    pred_file = tmp_path / "preds.csv.gz"
    geo_file = tmp_path / "geo.jsonl.gz"
    artifacts = write_discovery_artifacts(
        merged_run,
        pred_file_path=str(pred_file),
        geo_opt_file_path=str(geo_file),
        df_wbm_summary=df_summary,
        df_wbm_cse=df_entries,
        compatibility=MaterialsProject2020Compatibility(check_potcar=False),
        elemental_ref_energies=ref_energies,
    )
    pred_col = discovery_pred_col("emt")
    assert artifacts.predictions.loc[material_id, pred_col] == pytest.approx(
        expected_e_form
    )
    assert pd.read_csv(pred_file).loc[0, pred_col] == pytest.approx(expected_e_form)
    assert pd.isna(artifacts.predictions.loc[failed_id, pred_col])
    assert artifacts.predictions.loc[failed_id, "error"] == "RuntimeError: failed"
    assert artifacts.n_failed == 1
    df_geo = pd.read_json(geo_file, lines=True)
    assert len(df_geo) == 1
    assert df_geo.loc[0, Key.mat_id] == material_id
    assert df_geo.loc[0, artifacts.struct_col]["@class"] == "Structure"


@pytest.mark.parametrize(
    "cli_args",
    [
        ["--write-yaml"],
        ["--merge-shards", "--write-yaml", "--dry-run"],
        ["--merge-shards", "--shard-index", "0"],
    ],
)
def test_cli_rejects_unsafe_flag_combinations(cli_args: list[str]) -> None:
    """CLI rejects YAML writes outside complete merges and sharded merge runs."""
    with pytest.raises(SystemExit, match="2"):
        discovery_runner.main(["--model", "emt", *cli_args])


@pytest.mark.parametrize(
    "model_ref",
    [*discovery_core.ARCHIVED_DISCOVERY_MODELS, "equflash-29M-oam"],
)
def test_cli_rejects_archived_discovery_models(
    capsys: pytest.CaptureFixture[str], model_ref: str
) -> None:
    """Archived model names and YAML keys fail with their canonical reason."""
    model_key = Model.from_ref(model_ref).name
    with pytest.raises(SystemExit, match="2"):
        discovery_runner.main(["--model", model_ref, "--print-cmd", "--dry-run"])
    assert f"{model_key} discovery is archived:" in capsys.readouterr().err


def test_list_models_matches_calculator_registry(
    capsys: pytest.CaptureFixture[str],
) -> None:
    """Model listing contains every runnable calculator exactly once."""
    assert discovery_runner.main(["--list-models"]) == 0
    listed_models = {
        line.partition(":")[0] for line in capsys.readouterr().out.splitlines()
    }
    assert listed_models == set(CALCULATORS)


@pytest.mark.parametrize(
    ("model_key", "extra_args"),
    [
        ("mace_mp_0", ["--max-steps", "3"]),
        ("emt", []),
    ],
)
def test_dependency_isolated_print_cmd(
    capsys: pytest.CaptureFixture[str], model_key: str, extra_args: list[str]
) -> None:
    """Print commands include model dependencies and forwarded run options."""
    assert (
        discovery_runner.main(
            ["--model", model_key, "--print-cmd", "--dry-run", *extra_args]
        )
        == 0
    )
    command = capsys.readouterr().out
    assert "uv run --no-project" in command
    assert all(dependency in command for dependency in CALCULATORS[model_key].deps)
    assert f"models/run_discovery.py --model {model_key}" in command
    assert all(arg in command for arg in extra_args)
    assert "--dry-run" in command
    assert "--print-cmd" not in command


def test_runner_preserves_existing_artifact_columns(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Unified artifacts retain established YAML columns for published models."""
    assert discovery_runner._artifact_columns(Model.mace_mp_0) == (  # noqa: SLF001
        "e_form_per_atom_mace",
        "mace_structure",
    )
    assert discovery_runner._artifact_columns(None) == (None, None)  # noqa: SLF001
    existing_path = Model.mace_mp_0.metrics["discovery"]["pred_file"]
    make_artifact_data = partial(
        discovery_runner._artifact_yaml_data,  # noqa: SLF001
        Model.mace_mp_0,
        "discovery",
        column_key="pred_col",
        column="e_form_per_atom_mace",
    )
    for path, clears_url in (
        (f"{discovery_runner.ROOT}/{existing_path}", False),
        (f"{discovery_runner.ROOT}/models/mace/new-preds.csv.gz", True),
    ):
        artifact_data = make_artifact_data(path)
        assert ("pred_file_url" in artifact_data) is clears_url
        if clears_url:
            assert artifact_data["pred_file_url"] is None

    # simulate the backslash-separated relative path os.path.relpath returns on Windows
    monkeypatch.setattr(
        discovery_runner.os.path,
        "relpath",
        lambda _path, _start: existing_path.replace("/", "\\"),
    )
    artifact_data = make_artifact_data(f"{discovery_runner.ROOT}/{existing_path}")
    assert artifact_data["pred_file"] == existing_path
    assert "pred_file_url" not in artifact_data


def test_write_yaml_results_masks_outliers_and_updates_yaml(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """--write-yaml masks unrealistic outliers and writes subset metrics + paths."""
    import matbench_discovery.data as data_module

    material_ids = ["wbm-1", "wbm-2", "wbm-3", "wbm-4"]
    df_fake_wbm = pd.DataFrame(
        {
            MbdKey.each_true: [-0.1, 0.2, 0.0, -0.2],
            MbdKey.e_form_dft: [-1.0, -0.5, -0.8, -0.9],
            MbdKey.uniq_proto: [True, True, False, True],
        },
        index=pd.Index(material_ids, name=str(Key.mat_id)),
    )
    monkeypatch.setattr(data_module, "df_wbm", df_fake_wbm)

    pred_col = discovery_pred_col("emt")
    artifacts = DiscoveryArtifacts(
        pred_file_path=f"{tmp_path}/preds.csv.gz",
        geo_opt_file_path=f"{tmp_path}/geo.jsonl.gz",
        pred_col=pred_col,
        struct_col="structure",
        # wbm-3 is a 7.8 eV/atom outlier, wbm-4 failed to relax
        predictions=pd.DataFrame(
            {pred_col: [-1.05, -0.45, 7.0, None]}, index=material_ids
        ),
        n_success=3,
        n_failed=1,
    )
    yaml_path = tmp_path / "model.yml"
    old_discovery = {
        "pred_file": "models/old/preds.csv.gz",
        "pred_file_url": "https://example.com/old",
    }
    yaml_path.write_text(yaml.safe_dump({"metrics": {"discovery": old_discovery}}))
    mock_model = cast(
        "Model",
        SimpleNamespace(yaml_path=str(yaml_path), metrics={"discovery": old_discovery}),
    )
    discovery_runner._write_yaml_results(mock_model, artifacts)  # noqa: SLF001

    written = yaml.safe_load(yaml_path.read_text())
    discovery_yaml = written["metrics"]["discovery"]
    assert discovery_yaml["pred_file"] == artifacts.pred_file_path
    assert discovery_yaml["pred_file_url"] is None, "stale URL must be invalidated"
    assert discovery_yaml["pred_col"] == pred_col
    assert written["metrics"]["geo_opt"]["struct_col"] == "structure"
    assert set(discovery_yaml) >= {str(subset) for subset in TestSubset}
    full_metrics = discovery_yaml[str(TestSubset.full_test_set)]
    assert full_metrics[str(MbdKey.missing_preds)] == 2
    assert full_metrics["TP"] == 1.0
    assert isinstance(full_metrics["TP"], float), "counts are written as floats"
    # DAF denominator is the uniq-proto stable prevalence 2/3, Precision is 1
    uniq_metrics = discovery_yaml[str(TestSubset.uniq_protos)]
    assert uniq_metrics["DAF"] == pytest.approx(1.5)


def test_merge_paths_reuse_shard_run_date(tmp_path: Path) -> None:
    """A next-day run or merge resumes one prior directory and its artifact prefix."""
    # clearly-past dates so today's default shard dir can never shadow the glob
    old_shard_dir = tmp_path / "2020-01-02-wbm-IS2RE-FIRE-shards"
    old_shard_dir.mkdir()
    resolve_paths = partial(
        discovery_runner._resolve_output_paths,  # noqa: SLF001
        out_dir=str(tmp_path),
        settings=RelaxationSettings(),
        dry_run=False,
        shard_dir=None,
        pred_file=None,
        geo_opt_file=None,
    )
    shard_dir, pred_file, geo_opt_file = resolve_paths()
    assert shard_dir == str(old_shard_dir)
    assert pred_file == str(tmp_path / "2020-01-02-wbm-IS2RE-FIRE.csv.gz")
    assert geo_opt_file == str(tmp_path / "2020-01-02-wbm-IS2RE-FIRE.jsonl.gz")

    (tmp_path / "2020-01-01-wbm-IS2RE-FIRE-shards").mkdir()
    with pytest.raises(ValueError, match="Multiple discovery shard directories"):
        resolve_paths()


def test_cli_returns_failure_when_any_relaxation_fails(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A shard cannot report success when even one material fails."""

    def load_partially_failing_atoms(
        model_key: str, *, dry_run: bool = False
    ) -> dict[str, Atoms]:
        """Return one EMT-compatible and one incompatible structure."""
        assert model_key == "emt"
        assert dry_run is True
        return {
            "wbm-cu": _bulk_atoms("wbm-cu"),
            "wbm-si": _bulk_atoms("wbm-si", "Si"),
        }

    monkeypatch.setattr(
        discovery_runner, "load_wbm_atoms", load_partially_failing_atoms
    )
    assert (
        discovery_runner.main(
            [
                "--model",
                "emt",
                "--dry-run",
                "--max-steps",
                "0",
                "--out-dir",
                str(tmp_path),
            ]
        )
        == 1
    )
