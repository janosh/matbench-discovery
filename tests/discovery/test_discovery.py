"""Tests for the discovery relaxation pipeline and command-line runner."""

from __future__ import annotations

import copy
import os
from functools import partial
from types import SimpleNamespace
from typing import TYPE_CHECKING, Any, cast

import numpy as np
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
from matbench_discovery.data import file_ref_name, file_ref_url, make_file_ref
from matbench_discovery.discovery import (
    DISCOVERY_PRED_COL,
    DiscoveryArtifacts,
    MergedDiscoveryRun,
    RelaxationRecord,
    RelaxationSettings,
    merge_discovery_shards,
    read_discovery_shard,
    relax_atoms,
    run_discovery_shard,
    write_discovery_artifacts,
)
from matbench_discovery.energy import calc_energy_from_e_refs
from matbench_discovery.enums import MbdKey, Model, TestSubset
from matbench_discovery.hpc import partition_material_ids

if TYPE_CHECKING:
    from pathlib import Path

    from ase import Atoms


def bulk_atoms(material_id: str, element: str = "Cu", repeat: int = 1) -> Atoms:
    """Build a tagged periodic structure for discovery runner tests."""
    atoms = bulk(element) * (repeat, 1, 1)
    atoms.info[str(Key.mat_id)] = material_id
    return atoms


def run_single_point_shard(
    atoms_by_id: dict[str, Atoms], shard_path: str, **overrides: object
) -> discovery_core.LoadedShard:
    """Run an EMT shard with fixed-cell settings, allowing overrides."""
    kwargs: dict[str, Any] = {
        "calculator": EMT(),
        "model_key": "emt",
        "atoms_by_id": atoms_by_id,
        "assigned_material_ids": list(atoms_by_id),
        "shard_path": shard_path,
        "shard_index": 0,
        "n_shards": 1,
        "settings": RelaxationSettings(max_steps=0, cell_filter="none"),
        "dataset_material_ids": list(atoms_by_id),
    }
    kwargs |= overrides
    return run_discovery_shard(**kwargs)


def test_artifacts_match_legacy_mp2020_join_semantics(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """MP2020 artifacts retain valid entries and record individual rejections."""
    material_id = "wbm-1"
    failed_id = "wbm-2"
    rejected_id = "wbm-3"
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
                material_id=rejected_id,
                structure=structure.as_dict(),
                energy=-13.0,
                converged=True,
                n_steps=3,
            ),
            RelaxationRecord(material_id=failed_id, error="RuntimeError: failed"),
        ),
        run_metadata={},
        n_shards=1,
    )
    df_summary = pd.DataFrame(
        {
            Key.mat_id: [material_id, failed_id, rejected_id],
            Key.formula: ["Cu4", "Cu", "Cu4"],
        }
    ).set_index(Key.mat_id)
    original_entry = ComputedStructureEntry(
        structure, -8.0, parameters={"run_type": "GGA"}
    )
    df_entries = pd.DataFrame(
        {Key.computed_structure_entry: [original_entry, copy.deepcopy(original_entry)]},
        index=[material_id, rejected_id],
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

    process_entry = compatibility.process_entry

    def reject_selected_entry(
        entry: ComputedStructureEntry, **kwargs: bool
    ) -> ComputedStructureEntry | None:
        """Reject one entry while processing all others normally."""
        if entry.energy == -13.0:
            return None
        return process_entry(entry, **kwargs)

    monkeypatch.setattr(compatibility, "process_entry", reject_selected_entry)
    pred_file = tmp_path / "preds.csv.gz"
    geo_file = tmp_path / "geo.jsonl.gz"
    artifacts = write_discovery_artifacts(
        merged_run,
        pred_file_path=str(pred_file),
        geo_opt_file_path=str(geo_file),
        df_wbm_summary=df_summary,
        df_wbm_cse=df_entries,
        compatibility=compatibility,
        elemental_ref_energies=ref_energies,
    )
    pred_col = DISCOVERY_PRED_COL
    assert artifacts.predictions.loc[material_id, pred_col] == pytest.approx(
        expected_e_form
    )
    df_pred = pd.read_csv(pred_file)
    assert set(df_pred) == {"material_id", pred_col}
    assert df_pred.loc[0, pred_col] == pytest.approx(expected_e_form)
    assert pd.isna(artifacts.predictions.loc[failed_id, pred_col])
    assert artifacts.predictions.loc[failed_id, "error"] == "RuntimeError: failed"
    assert pd.isna(artifacts.predictions.loc[rejected_id, pred_col])
    assert (
        artifacts.predictions.loc[rejected_id, "error"]
        == "MP2020 compatibility rejected entry"
    )
    assert artifacts.n_failed == 2
    df_geo = pd.read_json(geo_file, lines=True)
    assert len(df_geo) == 1
    assert set(df_geo) == {"material_id", "structure"}
    assert df_geo.loc[0, "material_id"] == material_id
    assert df_geo.loc[0, "structure"]["@class"] == "Structure"


def test_runner_writes_canonical_artifact_metadata(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Regenerated artifacts invalidate stale URLs and record logical identity."""
    existing_path = file_ref_name(Model.mace_mp_0.metrics["discovery"]["pred_file"])
    assert existing_path is not None
    make_artifact_data = discovery_runner._artifact_yaml_data  # noqa: SLF001
    for path in (
        f"{discovery_runner.ROOT}/{existing_path}",
        f"{discovery_runner.ROOT}/models/mace/mace-mp-0/2026-07-02-discovery.csv.gz",
    ):
        artifact_data = make_artifact_data(path)
        assert file_ref_url(artifact_data["pred_file"]) is None
        assert "pred_file_url" not in artifact_data
        assert "pred_file_artifact" not in artifact_data

    artifact_data = discovery_runner._artifact_yaml_data(  # noqa: SLF001
        f"{discovery_runner.ROOT}/{existing_path}",
    )
    assert "pred_file_url" not in artifact_data
    assert "pred_file_artifact" not in artifact_data

    def windows_relpath(_path: str, _start: str) -> str:
        """Return a Windows-style relative artifact path."""
        return existing_path.replace("/", "\\")

    monkeypatch.setattr(discovery_runner.os.path, "relpath", windows_relpath)
    artifact_data = make_artifact_data(f"{discovery_runner.ROOT}/{existing_path}")
    assert artifact_data["pred_file"] == make_file_ref(existing_path)

    def raise_cross_drive(_path: str, _start: str) -> str:
        """Mimic relpath failing across Windows drives."""
        raise ValueError("path is on mount 'C:', start on mount 'D:'")

    monkeypatch.setattr(discovery_runner.os.path, "relpath", raise_cross_drive)
    external_path = "/mnt/c/tmp/preds.csv.gz"
    artifact_data = make_artifact_data(external_path)
    assert artifact_data["pred_file"] == make_file_ref(os.path.abspath(external_path))


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
    pred_col = DISCOVERY_PRED_COL
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
        "pred_file": make_file_ref(
            "models/mace/mace-mp-0/2026-07-01-discovery.csv.gz",
            url="https://example.com/old",
        ),
    }
    yaml_path.write_text(yaml.safe_dump({"metrics": {"discovery": old_discovery}}))
    mock_model = cast(
        "Model",
        SimpleNamespace(yaml_path=str(yaml_path), metrics={"discovery": old_discovery}),
    )
    run_metadata = {
        "hardware": "NVIDIA H200",
        "run_time_sec": 123.45,
        "max_rss_gb": 4.2,
        "hostnames": ["node-1"],  # audit-only fields must not leak into the YAML
    }
    discovery_runner._write_yaml_results(mock_model, artifacts, run_metadata)  # noqa: SLF001
    written = yaml.safe_load(yaml_path.read_text())
    discovery_yaml = written["metrics"]["discovery"]
    # normpath makes the comparison robust to Windows CI, where tmp_path sits on a
    # different drive than the repo (relpath impossible -> absolute native path)
    pred_name = file_ref_name(discovery_yaml["pred_file"])
    assert pred_name is not None
    assert os.path.normpath(pred_name) == os.path.normpath(artifacts.pred_file_path)
    assert file_ref_url(discovery_yaml["pred_file"]) is None, (
        "stale URL must be invalidated"
    )
    assert "pred_file_url" not in discovery_yaml
    assert "pred_col" not in discovery_yaml
    assert discovery_yaml["hardware"] == "NVIDIA H200"
    assert discovery_yaml["run_time_sec"] == pytest.approx(123.45)
    assert discovery_yaml["max_rss_gb"] == pytest.approx(4.2)
    assert "hostnames" not in discovery_yaml
    assert "max_gpu_mem_gb" not in discovery_yaml, "absent cost fields stay absent"
    geo_opt_yaml = written["metrics"]["geo_opt"]
    assert "struct_col" not in geo_opt_yaml
    assert "hardware" not in geo_opt_yaml, "cost provenance lives under discovery"
    assert set(discovery_yaml) >= {str(subset) for subset in TestSubset}
    full_metrics = discovery_yaml[str(TestSubset.full_test_set)]
    assert full_metrics[str(MbdKey.missing_preds)] == 2
    assert full_metrics["TP"] == 1.0
    assert isinstance(full_metrics["TP"], float), "counts are written as floats"
    # DAF denominator is the uniq-proto stable prevalence 2/3, Precision is 1
    uniq_metrics = discovery_yaml[str(TestSubset.uniq_protos)]
    assert uniq_metrics["DAF"] == pytest.approx(1.5)


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
    [*discovery_core.ARCHIVED_DISCOVERY_MODELS, "equflash-29m-oam"],
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
    assert listed_models == set(CALCULATORS) - set(
        discovery_core.ARCHIVED_DISCOVERY_MODELS
    )


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
    assert pred_file == str(tmp_path / "2020-01-02-discovery.csv.gz")
    assert geo_opt_file == str(tmp_path / "2020-01-02-geo-opt.jsonl.gz")

    (tmp_path / "2020-01-01-wbm-IS2RE-FIRE-shards").mkdir()
    with pytest.raises(ValueError, match="Multiple discovery shard directories"):
        resolve_paths()


def test_cli_returns_failure_when_any_relaxation_fails(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A shard cannot report success when even one material fails."""

    def load_partially_failing_atoms(
        model_key: str, *, dry_run: bool = False
    ) -> dict[str, Atoms]:
        """Return one EMT-compatible and one incompatible structure."""
        assert model_key == "emt"
        assert dry_run is True
        return {
            "wbm-cu": bulk_atoms("wbm-cu"),
            "wbm-si": bulk_atoms("wbm-si", "Si"),
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


def test_emt_relaxation_smoke() -> None:
    """EMT produces a finite energy, structure, and convergence state."""
    record = relax_atoms(
        bulk_atoms("wbm-cu"),
        EMT(),
        RelaxationSettings(max_force=0.5, max_steps=2, cell_filter="none"),
    )
    assert record.material_id == "wbm-cu"
    assert record.error is None
    assert record.structure is not None
    assert record.energy is not None
    assert record.n_steps <= 2
    single_point_record = relax_atoms(
        bulk_atoms("wbm-cu-single-point"),
        EMT(),
        RelaxationSettings(max_steps=0, cell_filter=None),
    )
    assert single_point_record.converged is False


def _run_test_shard(
    tmp_path: Path, *, shard_index: int, n_shards: int, material_id: str
) -> str:
    """Write one complete EMT shard and return its path."""
    shard_path = tmp_path / f"shard-{shard_index:03d}-of-{n_shards:03d}.jsonl"
    run_single_point_shard(
        {material_id: bulk_atoms(material_id)},
        str(shard_path),
        shard_index=shard_index,
        n_shards=n_shards,
        dataset_material_ids=(["wbm-1", "wbm-2"] if n_shards == 2 else [material_id]),
    )
    return str(shard_path)


def _run_two_test_shards(tmp_path: Path) -> list[str]:
    """Write a complete two-shard EMT run."""
    return [
        _run_test_shard(
            tmp_path,
            shard_index=shard_index,
            n_shards=2,
            material_id=material_id,
        )
        for shard_index, material_id in enumerate(("wbm-1", "wbm-2"))
    ]


def _merge_test_shards(shard_paths: list[str]) -> MergedDiscoveryRun:
    """Merge shards from the standard two-material EMT test run."""
    return merge_discovery_shards(
        shard_paths,
        model_key="emt",
        expected_material_ids=["wbm-1", "wbm-2"],
    )


def test_partition_material_ids_is_deterministic_and_balanced() -> None:
    """Atom balancing is independent of mapping insertion order and deterministic."""
    atoms_by_id = {
        "wbm-3": bulk_atoms("wbm-3", repeat=3),
        "wbm-1": bulk_atoms("wbm-1", repeat=1),
        "wbm-4": bulk_atoms("wbm-4", repeat=4),
        "wbm-2": bulk_atoms("wbm-2", repeat=2),
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
        (None, 3, 7, 3, "Partial Slurm arrays"),
        # a 1-based full array with explicit --n-shards must fail on every task,
        # not run shards 1-9 and silently never produce shard 0
        (10, 1, 10, 1, "SLURM_ARRAY_TASK_MAX=10 exceeds the last shard index 9"),
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


def test_shard_resume_preserves_success_and_failure_records(tmp_path: Path) -> None:
    """Completed and failed attempts are both durable and skipped on resume."""
    atoms_by_id = {
        "wbm-cu": bulk_atoms("wbm-cu"),
        "wbm-si": bulk_atoms("wbm-si", element="Si"),
    }
    shard_path = tmp_path / "shard-000-of-001.jsonl"
    first_run = run_single_point_shard(atoms_by_id, str(shard_path))
    assert [record.material_id for record in first_run.records] == [
        "wbm-cu",
        "wbm-si",
    ]
    assert first_run.records[0].energy is not None
    assert first_run.records[1].error is not None

    size_before_resume = os.path.getsize(shard_path)
    resumed = run_single_point_shard(atoms_by_id, str(shard_path))
    assert resumed.records == first_run.records
    assert os.path.getsize(shard_path) == size_before_resume


def test_shard_resume_handles_incomplete_events_and_header_change(
    tmp_path: Path,
) -> None:
    """Resume repairs incomplete events but rejects changed settings."""
    material_id = "wbm-1"
    shard_path = tmp_path / "shard-000-of-001.jsonl"
    run_single_point_shard({material_id: bulk_atoms(material_id)}, str(shard_path))
    shard_path.write_bytes(shard_path.read_bytes().rstrip(b"\n"))
    discovery_core._append_shard_event(  # noqa: SLF001
        str(shard_path),
        {"type": "run_metadata", "run_metadata": {"run_time_sec": 1}},
    )
    assert len(read_discovery_shard(str(shard_path)).run_metadata_segments) == 2

    with shard_path.open(mode="a", encoding="utf-8") as file:
        file.write('{"type": "record", "material_id":')

    run_single_point_shard({material_id: bulk_atoms(material_id)}, str(shard_path))
    assert read_discovery_shard(str(shard_path)).ignored_trailing_line is False

    with pytest.raises(ValueError, match="does not match this run"):
        run_single_point_shard(
            {material_id: bulk_atoms(material_id)},
            str(shard_path),
            settings=RelaxationSettings(max_steps=1, cell_filter=None),
        )
    with shard_path.open(mode="a", encoding="utf-8") as file:
        file.write("{malformed}\n")
    with pytest.raises(ValueError, match="Malformed JSON event"):
        read_discovery_shard(str(shard_path))


def test_shard_survives_unserializable_calculator_info(tmp_path: Path) -> None:
    """Arrays in atoms.info serialize; exotic objects downgrade to failure records."""
    array_atoms = bulk_atoms("wbm-array")
    array_atoms.info["forces"] = np.arange(3.0)
    exotic_atoms = bulk_atoms("wbm-exotic")
    exotic_atoms.info["handle"] = object()
    atoms_by_id = {"wbm-array": array_atoms, "wbm-exotic": exotic_atoms}
    shard = run_single_point_shard(atoms_by_id, f"{tmp_path}/shard-000-of-001.jsonl")
    records_by_id = {record.material_id: record for record in shard.records}
    array_record = records_by_id["wbm-array"]
    assert array_record.error is None
    assert array_record.structure is not None
    assert array_record.structure["properties"]["forces"] == [0.0, 1.0, 2.0]
    exotic_record = records_by_id["wbm-exotic"]
    assert exotic_record.structure is None
    assert exotic_record.error is not None
    assert "unserializable relaxation record" in exotic_record.error


@pytest.mark.parametrize("initial_content", ["", '{"type": "header", "model_'])
def test_shard_self_heals_interrupted_header(
    tmp_path: Path,
    initial_content: str,
) -> None:
    """A crash during the very first header append must not require manual repair."""
    shard_path = tmp_path / "shard-000-of-001.jsonl"
    shard_path.write_text(initial_content)
    shard = run_single_point_shard({"wbm-1": bulk_atoms("wbm-1")}, str(shard_path))
    assert [record.material_id for record in shard.records] == ["wbm-1"]


def test_merge_drops_cost_fields_unless_every_shard_reported_them(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Cost provenance is all-or-nothing across shards, mirroring the MD eval."""

    def fake_run_metadata(_run_time_sec: float) -> dict[str, float | str]:
        """Return deterministic cost metadata on every platform."""
        return {
            "hardware": "test hardware",
            "run_time_sec": 1.0,
            "max_rss_gb": 2.0,
            "max_gpu_mem_gb": 3.0,
        }

    monkeypatch.setattr(discovery_core, "_run_segment_metadata", fake_run_metadata)
    shard_0, shard_1 = _run_two_test_shards(tmp_path)
    complete_run = _merge_test_shards([shard_0, shard_1])
    assert set(discovery_core.COST_PROVENANCE_KEYS) <= set(complete_run.run_metadata)

    # strip shard 1's run_metadata events to simulate a crash before provenance write
    with open(shard_1, encoding="utf-8") as file:
        shard_1_lines = [
            line for line in file.read().splitlines() if '"run_metadata"' not in line
        ]
    with open(shard_1, mode="w", encoding="utf-8") as file:
        file.write("\n".join(shard_1_lines) + "\n")
    partial_run = _merge_test_shards([shard_0, shard_1])
    for cost_key in discovery_core.COST_PROVENANCE_KEYS:
        assert cost_key not in partial_run.run_metadata, cost_key
    assert partial_run.run_metadata["n_shards"] == 2  # audit fields remain


def test_merge_rejects_incomplete_duplicate_and_wrong_coverage(tmp_path: Path) -> None:
    """Strict merge rejects missing shards, duplicate shard IDs, and foreign IDs."""
    shard_0, shard_1 = _run_two_test_shards(tmp_path)
    merged_run = _merge_test_shards([shard_0, shard_1])
    assert [record.material_id for record in merged_run.records] == ["wbm-1", "wbm-2"]
    assert merged_run.n_shards == 2

    incomplete_path = tmp_path / "incomplete" / "shard-000-of-001.jsonl"
    incomplete_atoms = {
        material_id: bulk_atoms(material_id) for material_id in ("wbm-1", "wbm-2")
    }
    run_single_point_shard(incomplete_atoms, str(incomplete_path))
    shard_lines = incomplete_path.read_text().splitlines()
    incomplete_path.write_text("\n".join(shard_lines[:2]) + "\n")
    with pytest.raises(ValueError, match="records do not match assignment"):
        _merge_test_shards([str(incomplete_path)])

    with pytest.raises(ValueError, match="Expected shard indices"):
        _merge_test_shards([shard_0])
    with pytest.raises(ValueError, match="Duplicate shard indices"):
        _merge_test_shards([shard_0, shard_0])
    with pytest.raises(ValueError, match="Inconsistent shard metadata"):
        merge_discovery_shards(
            [shard_0, shard_1],
            model_key="emt",
            expected_material_ids=["wbm-1", "wbm-foreign"],
        )
