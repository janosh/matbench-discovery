"""Tests for the unified resumable PhononDB kappa pipeline."""

from __future__ import annotations

import gzip
import json
import os
from dataclasses import asdict, replace
from types import SimpleNamespace
from typing import TYPE_CHECKING, Any, cast

import numpy as np
import pandas as pd
import pytest
from ase import Atoms
from pymatviz.enums import Key

from matbench_discovery import data as mbd_data
from matbench_discovery.enums import MbdKey, Model
from matbench_discovery.metrics.phonons import write_metrics_to_yaml
from matbench_discovery.phonons.adapters import (
    FairchemKappaAdapter,
    RelaxationResult,
    StandardKappaAdapter,
    get_kappa_adapter,
)
from matbench_discovery.phonons.adapters.equflash import EquFlashKappaAdapter
from matbench_discovery.phonons.adapters.pet import PetKappaAdapter
from matbench_discovery.phonons.pipeline import (
    DRY_RUN_MAX_FC3_EVALUATIONS,
    KAPPA_PROTOCOL,
    KappaComputation,
    KappaRecord,
    KappaSettings,
    atomic_write_gzip_json,
    calculate_kappa_for_structure,
    checkpoint_digest,
    file_sha256,
    merge_kappa_shards,
    read_kappa_record,
    record_path,
    run_kappa_shard,
    should_calculate_conductivity,
    write_kappa_artifacts,
)
from matbench_discovery.phonons.schema import (
    normalize_kappa_dataframe,
    normalize_kappa_result,
    voigt_6_to_full_3x3,
)

if TYPE_CHECKING:
    from pathlib import Path

    from ase.calculators.calculator import Calculator
    from phono3py.api_phono3py import Phono3py

    from matbench_discovery.phonons.pipeline import MergedKappaRun


class ComputationStub:
    """Return deterministic successful calculations while counting invocations."""

    def __init__(self) -> None:
        """Initialize an empty call log."""
        self.material_ids: list[str] = []
        self.settings: list[KappaSettings] = []
        self.max_fc3_evaluations: list[int | None] = []

    def __call__(
        self,
        *,
        atoms: Atoms,
        calculator: Calculator,
        settings: KappaSettings,
        adapter: StandardKappaAdapter,
        log_file: str | None = None,
        max_fc3_evaluations: int | None = None,
    ) -> KappaComputation:
        """Return one normalized result and optional force sets."""
        del calculator, adapter, log_file
        material_id = str(atoms.info[Key.mat_id])
        self.material_ids.append(material_id)
        self.settings.append(settings)
        self.max_fc3_evaluations.append(max_fc3_evaluations)
        result = {
            str(Key.mat_id): material_id,
            str(Key.formula): atoms.get_chemical_formula(),
            str(Key.has_imag_ph_modes): False,
            str(MbdKey.kappa_tot_rta): [[1, 2, 3, 4, 5, 6]],
            "errors": [],
            "error_traceback": [],
        }
        forces = (
            {
                "fc2_set": np.full((1, len(atoms), 3), 2),
                "fc3_set": np.full((1, len(atoms), 3), 3),
            }
            if settings.save_forces
            else None
        )
        return KappaComputation(result=result, forces=forces, error=None)


class Phono3pyStub:
    """Minimal mutable phono3py state used by the orchestration test."""

    def __init__(self) -> None:
        """Initialize without force constants."""
        self.forces: np.ndarray | None = None

    def produce_fc3(self, *, symmetrize_fc3r: bool) -> None:
        """Validate FC3 production after forces are assigned."""
        assert symmetrize_fc3r is True
        assert self.forces is not None


class PipelineAdapterStub(StandardKappaAdapter):
    """Exercise the shared orchestration without invoking external phonon packages."""

    def __init__(self) -> None:
        """Initialize deterministic phono state and call counters."""
        self.phono3py = Phono3pyStub()
        self.n_fc3_calls = 0

    def relax(
        self,
        atoms: Atoms,
        calculator: Calculator,
        settings: KappaSettings,
        *,
        log_file: str | None,
    ) -> RelaxationResult:
        """Return the structure unchanged with stable symmetry."""
        del calculator, settings, log_file
        return RelaxationResult(
            atoms=atoms,
            initial_spg_num=225,
            final_spg_num=225,
            max_stress=None,
            reached_max_steps=False,
            n_steps=0,
        )

    def init_phono3py(self, atoms: Atoms, settings: KappaSettings) -> Phono3py:
        """Return the deterministic phono state."""
        del atoms, settings
        return cast("Phono3py", self.phono3py)

    def calculate_fc2(
        self,
        phono3py: Phono3py,
        calculator: Calculator,
        settings: KappaSettings,
        *,
        progress: dict[str, Any] | None = None,
    ) -> tuple[Phono3py, np.ndarray, np.ndarray]:
        """Return one imaginary mode and a deterministic FC2 force set."""
        del calculator, settings, progress
        return phono3py, np.ones((1, 1, 3)), np.array([[-1.0, 0.0, 0.0, 1.0]])

    def calculate_fc3(
        self,
        phono3py: Phono3py,
        calculator: Calculator,
        settings: KappaSettings,
        *,
        progress: dict[str, Any] | None = None,
        max_evaluations: int | None = None,
    ) -> np.ndarray:
        """Record and return one deterministic FC3 force set."""
        del phono3py, calculator, settings, progress, max_evaluations
        self.n_fc3_calls += 1
        return np.full((1, 1, 3), 3.0)


def make_atoms_by_id() -> dict[str, Atoms]:
    """Build a small insertion-ordered structure mapping with varied costs."""
    atoms_by_id = {}
    for material_idx, n_atoms in enumerate((1, 4, 2, 3)):
        atoms = Atoms("H" * n_atoms, positions=np.zeros((n_atoms, 3)))
        material_id = f"mp-{material_idx}"
        atoms.info[Key.mat_id] = material_id
        atoms_by_id[material_id] = atoms
    return atoms_by_id


@pytest.fixture
def shard_env(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> SimpleNamespace:
    """Provide a deterministic dataset, calculator, computation, and shard path."""
    from matbench_discovery.phonons import pipeline

    dataset_path = f"{tmp_path}/phonondb.extxyz"
    with open(dataset_path, mode="w", encoding="utf-8") as file:
        file.write("test dataset identity")
    computation = ComputationStub()
    monkeypatch.setattr(pipeline, "calculate_kappa_for_structure", computation)
    atoms_by_id = make_atoms_by_id()
    calculator = cast("Calculator", object())
    shard_dir = f"{tmp_path}/shards"

    def run(
        *,
        settings: KappaSettings,
        shard_index: int = 0,
        n_shards: int = 1,
        retry_failures: bool = False,
        dry_run: bool = False,
        checkpoint: str | None = None,
        dtype: str = "float64",
    ) -> tuple[KappaRecord, ...]:
        """Run one shard with the shared deterministic test environment."""
        return run_kappa_shard(
            calculator=calculator,
            model_key="test_model",
            atoms_by_id=atoms_by_id,
            dataset_path=dataset_path,
            shard_dir=shard_dir,
            shard_index=shard_index,
            n_shards=n_shards,
            settings=settings,
            retry_failures=retry_failures,
            dry_run=dry_run,
            checkpoint=checkpoint,
            dtype=dtype,
        )

    return SimpleNamespace(
        dataset_path=dataset_path,
        atoms_by_id=atoms_by_id,
        computation=computation,
        calculator=calculator,
        shard_dir=shard_dir,
        run=run,
    )


def update_temp_kappa_yaml(
    tmp_path: Path,
    initial_yaml: str,
    *,
    pred_file_path: str = "models/test/new.json.gz",
    run_metadata: dict[str, object] | None = None,
    force_file_path: str | None = None,
    run_info_path: str | None = None,
    replace_pred_file: bool = False,
) -> dict[str, Any]:
    """Update a temporary model YAML and return its kappa metrics mapping."""
    yaml_path = f"{tmp_path}/model.yml"
    with open(yaml_path, mode="w", encoding="utf-8") as file:
        file.write(initial_yaml)
    model = cast("Model", SimpleNamespace(yaml_path=yaml_path))
    write_metrics_to_yaml(
        model,
        {"srme": 0.25, "sre": 0.125},
        pred_file_path,
        run_metadata=run_metadata,
        force_file_path=force_file_path,
        run_info_path=run_info_path,
        replace_pred_file=replace_pred_file,
    )
    with open(yaml_path, encoding="utf-8") as file:
        metadata = mbd_data.round_trip_yaml.load(file)
    return cast("dict[str, Any]", metadata["metrics"]["phonons"]["kappa_103"])


def test_kappa_settings_validate_and_canonicalize() -> None:
    """Sparse versioned settings use defaults and canonicalize ASE aliases."""
    raw_settings = {
        "protocol": KAPPA_PROTOCOL,
        "ase_optimizer": "goqn",
        "ase_filter": "exp",
    }
    settings = KappaSettings.from_mapping(raw_settings)
    assert settings.ase_optimizer == "GOQN"
    assert settings.ase_filter == "ExpCellFilter"
    assert settings.temperatures == (300.0,)
    assert settings.is_plusminus == "auto"
    assert settings.digest == KappaSettings.from_mapping(raw_settings).digest

    raw_settings.pop("protocol")
    with pytest.raises(ValueError, match="protocol must be"):
        KappaSettings.from_mapping(raw_settings)


@pytest.mark.parametrize(
    ("field_name", "invalid_value", "exception_type", "error_match"),
    [
        ("displacement_distance", 0, ValueError, "positive and finite"),
        ("temperatures", (), ValueError, "positive finite"),
        ("max_steps", -1, ValueError, "non-negative integer"),
        ("batch_size", 0, ValueError, "positive integer"),
        ("max_atoms_per_batch", 0, ValueError, "positive integer or null"),
        ("ase_filter", "UnitCellFilter", ValueError, "Unsupported kappa ASE filter"),
        ("relaxation_mode", "three-stage", ValueError, "must be one of"),
        ("is_plusminus", "sometimes", TypeError, "boolean or 'auto'"),
    ],
)
def test_kappa_settings_reject_invalid_values(
    field_name: str,
    invalid_value: object,
    exception_type: type[Exception],
    error_match: str,
) -> None:
    """Invalid scientific and execution settings fail before a run starts."""
    with pytest.raises(exception_type, match=error_match):
        replace(KappaSettings(), **{field_name: invalid_value})


def test_checkpoint_digest_hashes_artifact_bytes(tmp_path: Path) -> None:
    """Checkpoint provenance follows file content rather than only its path."""
    checkpoint_path = f"{tmp_path}/model.pt"
    with open(checkpoint_path, mode="wb") as file:
        file.write(b"first checkpoint")
    first_digest = checkpoint_digest(checkpoint_path)
    checkpoint_stat = os.stat(checkpoint_path)
    with open(checkpoint_path, mode="wb") as file:
        file.write(b"other checkpoint")
    os.utime(
        checkpoint_path,
        ns=(checkpoint_stat.st_atime_ns, checkpoint_stat.st_mtime_ns),
    )
    assert checkpoint_digest(checkpoint_path) != first_digest


def test_generated_models_use_stable_registry_provenance() -> None:
    """Generated NequIP caches do not change merge identity across hosts."""
    from models import run_kappa

    model = Model.nequip_mp_l_0_1
    assert (
        run_kappa.resolved_artifact_identifier(model, None)
        == model.metadata["checkpoint_url"]
    )


def test_kappa_output_paths_reuse_exactly_one_prior_shard_tree(tmp_path: Path) -> None:
    """Implicit paths reuse one dated shard tree and reject ambiguous histories."""
    from models import run_kappa

    out_dir = f"{tmp_path}/outputs"
    first_shard_dir = f"{out_dir}/2026-07-08-phonondb-kappa-103-shards"
    os.makedirs(first_shard_dir)
    paths = run_kappa.resolve_output_paths(
        out_dir=out_dir,
        dry_run=False,
        shard_dir=None,
    )
    assert paths == (
        first_shard_dir,
        first_shard_dir.removesuffix("-shards") + ".json.gz",
    )
    custom_shard_dir = f"{out_dir}/custom-run"
    assert (
        run_kappa.resolve_output_paths(
            out_dir=out_dir, dry_run=False, shard_dir=custom_shard_dir
        )[1]
        == f"{custom_shard_dir}.json.gz"
    )

    os.makedirs(f"{out_dir}/2026-07-09-phonondb-kappa-103-shards")
    with pytest.raises(ValueError, match="Multiple kappa shard directories"):
        run_kappa.resolve_output_paths(
            out_dir=out_dir,
            dry_run=False,
            shard_dir=None,
        )


def test_kappa_cli_requires_declared_explicit_checkpoint(
    capsys: pytest.CaptureFixture[str],
) -> None:
    """Models without a downloadable artifact reject implicit calculator loading."""
    from models import run_kappa

    with pytest.raises(SystemExit, match="2"):
        run_kappa.main(["--model", "equiformer_v3_mp", "--print-cmd"])
    assert "requires --checkpoint" in capsys.readouterr().err
    with pytest.raises(SystemExit, match="2"):
        run_kappa.main(
            ["--model", "equiformer_v3_mp", "--merge-shards", "--write-yaml"]
        )
    assert "requires --checkpoint" in capsys.readouterr().err

    assert (
        run_kappa.main(
            [
                "--model",
                "equiformer_v3_mp",
                "--checkpoint",
                "/tmp/model.pt",
                "--print-cmd",
            ]
        )
        == 0
    )
    assert "--checkpoint /tmp/model.pt" in capsys.readouterr().out


def test_yaml_write_requires_verified_settings_and_dataset(tmp_path: Path) -> None:
    """Official metrics reject overridden protocols and noncanonical datasets."""
    from models import run_kappa

    dataset_path = f"{tmp_path}/phonondb.extxyz"
    with open(dataset_path, mode="w", encoding="utf-8") as file:
        file.write("canonical dataset")
    settings = KappaSettings.from_model("mace_mp_0")
    manifest = SimpleNamespace(
        settings=settings,
        dataset_hash=file_sha256(dataset_path),
    )
    merged_run = cast(
        "MergedKappaRun",
        SimpleNamespace(manifest=manifest, run_metadata={"n_failed": 0}),
    )
    run_kappa.validate_yaml_write_provenance(
        Model.mace_mp_0,
        merged_run,
        canonical_dataset_path=dataset_path,
    )

    manifest.settings = replace(settings, force_max=settings.force_max * 2)
    with pytest.raises(ValueError, match="overridden protocol"):
        run_kappa.validate_yaml_write_provenance(
            Model.mace_mp_0,
            merged_run,
            canonical_dataset_path=dataset_path,
        )
    manifest.settings = settings
    manifest.dataset_hash = "different"
    with pytest.raises(ValueError, match="noncanonical"):
        run_kappa.validate_yaml_write_provenance(
            Model.mace_mp_0,
            merged_run,
            canonical_dataset_path=dataset_path,
        )
    manifest.dataset_hash = file_sha256(dataset_path)
    merged_run.run_metadata["n_failed"] = 1
    with pytest.raises(ValueError, match="failed records"):
        run_kappa.validate_yaml_write_provenance(
            Model.mace_mp_0,
            merged_run,
            canonical_dataset_path=dataset_path,
        )


def test_dry_run_uses_and_merges_capped_protocol(shard_env: SimpleNamespace) -> None:
    """Dry runs execute and record the capped smoke-test protocol."""
    settings = KappaSettings(
        temperatures=(300, 600),
        max_steps=300,
        relaxation_mode="none",
        save_forces=True,
    )
    shard_env.run(settings=settings, dry_run=True)
    computation_stub = shard_env.computation
    assert computation_stub.settings
    assert all(
        run_settings.max_steps == 2 for run_settings in computation_stub.settings
    )
    assert all(
        run_settings.temperatures == (300.0,)
        for run_settings in computation_stub.settings
    )
    assert not any(
        run_settings.save_forces for run_settings in computation_stub.settings
    )
    assert set(computation_stub.max_fc3_evaluations) == {DRY_RUN_MAX_FC3_EVALUATIONS}

    merged_run = merge_kappa_shards(
        shard_env.shard_dir,
        model_key="test_model",
        expected_material_ids=tuple(shard_env.atoms_by_id),
    )
    assert merged_run.manifest.dry_run is True
    assert merged_run.manifest.settings.max_steps == 2
    assert merged_run.manifest.settings.save_forces is False


@pytest.mark.parametrize(
    (
        "has_imaginary_modes",
        "broken_symmetry",
        "ignore_imaginary_freqs",
        "conductivity_broken_symm",
        "expected",
    ),
    [
        (False, False, False, False, True),
        (True, False, False, False, False),
        (True, False, True, False, True),
        (False, True, False, False, False),
        (False, True, True, False, True),
        (True, True, True, False, True),
        (False, True, False, True, True),
    ],
)
def test_conductivity_gates_match_legacy_policy(
    has_imaginary_modes: bool,
    broken_symmetry: bool,
    ignore_imaginary_freqs: bool,
    conductivity_broken_symm: bool,
    expected: bool,
) -> None:
    """The ignore flag preserves the legacy policy of bypassing both gates."""
    settings = replace(
        KappaSettings(),
        ignore_imaginary_freqs=ignore_imaginary_freqs,
        conductivity_broken_symm=conductivity_broken_symm,
    )
    assert (
        should_calculate_conductivity(
            has_imaginary_modes=has_imaginary_modes,
            broken_symmetry=broken_symmetry,
            settings=settings,
        )
        is expected
    )


@pytest.mark.parametrize(
    ("ignore_imaginary_freqs", "expected_fc3_calls", "expected_skipped"),
    [(False, 0, True), (True, 1, False)],
)
def test_calculate_kappa_for_structure_runs_shared_physics_flow(
    monkeypatch: pytest.MonkeyPatch,
    ignore_imaginary_freqs: bool,
    expected_fc3_calls: int,
    expected_skipped: bool,
) -> None:
    """The shared flow executes FC2, gates FC3, and retains separate force sets."""
    from matbench_discovery.phonons import thermal_conductivity

    def calculate_conductivity_stub(
        phono3py: Phono3py, temperatures: tuple[float, ...]
    ) -> tuple[Phono3py, dict[str, Any], object]:
        """Return one finite conductivity result without solving the BTE."""
        assert temperatures == (300.0,)
        return (
            phono3py,
            {str(MbdKey.kappa_tot_rta): np.eye(3).tolist()},
            object(),
        )

    monkeypatch.setattr(
        thermal_conductivity, "calculate_conductivity", calculate_conductivity_stub
    )
    atoms = Atoms("H", cell=np.eye(3), positions=[[0, 0, 0]], pbc=True)
    atoms.info[Key.mat_id] = "mp-core"
    adapter = PipelineAdapterStub()
    computation = calculate_kappa_for_structure(
        atoms=atoms,
        calculator=cast("Calculator", object()),
        settings=KappaSettings(
            relaxation_mode="none",
            max_steps=0,
            ignore_imaginary_freqs=ignore_imaginary_freqs,
            save_forces=True,
        ),
        adapter=adapter,
    )

    assert computation.error is None
    assert adapter.n_fc3_calls == expected_fc3_calls
    assert computation.result.get("conductivity_skipped", False) is expected_skipped
    assert computation.forces is not None
    assert np.asarray(computation.forces["fc2_set"]).shape == (1, 1, 3)
    assert np.asarray(computation.forces["fc3_set"]).size == expected_fc3_calls * 3


def test_normalize_kappa_schema_aliases_and_voigt() -> None:
    """Legacy aliases and Voigt tensors become one canonical schema."""
    result = normalize_kappa_result(
        {
            "mp_id": "mp-1",
            str(Key.spg_num): 225,
            "relaxed_space_group_number": 221,
            str(MbdKey.kappa_tot_rta): [1, 2, 3, 4, 5, 6],
        }
    )
    assert result[str(Key.mat_id)] == "mp-1"
    assert result[str(Key.init_spg_num)] == 225
    assert result[str(Key.final_spg_num)] == 221
    assert result["broken_symmetry"] is True
    assert np.asarray(result[str(MbdKey.kappa_tot_rta)]).tolist() == [
        [1, 6, 5],
        [6, 2, 4],
        [5, 4, 3],
    ]
    for missing_spg in (np.nan, pd.NA):
        missing_result = normalize_kappa_result(
            {str(Key.init_spg_num): missing_spg, str(Key.final_spg_num): 225}
        )
        assert missing_result["broken_symmetry"] is False
    failed_result = normalize_kappa_result(
        {
            "errors": ["failed after partial conductivity output"],
            str(MbdKey.kappa_tot_rta): [1, 2, 3, 4, 5, 6],
        }
    )
    assert np.isnan(failed_result[str(MbdKey.kappa_tot_rta)])
    assert voigt_6_to_full_3x3([1, 2, 3]) == [1, 2, 3]


def test_normalize_kappa_dataframe_preserves_order() -> None:
    """Dataframe normalization retains row order while removing aliases."""
    import pandas as pd

    df_legacy = pd.DataFrame(
        {
            "mp_id": ["mp-2", "mp-1"],
            "initial_spg_num": [225, 221],
            str(MbdKey.kappa_tot_rta): [[1, 2, 3, 4, 5, 6]] * 2,
            "errors": [np.nan, []],
        }
    )
    df_normalized: pd.DataFrame = normalize_kappa_dataframe(df_legacy)
    assert list(df_normalized[str(Key.mat_id)]) == ["mp-2", "mp-1"]
    assert "mp_id" not in df_normalized
    assert np.asarray(df_normalized[str(MbdKey.kappa_tot_rta)].iloc[0]).shape == (
        3,
        3,
    )


@pytest.mark.parametrize(
    ("model_key", "adapter_type"),
    [
        ("mace_mp_0", StandardKappaAdapter),
        ("pet_oam_xl_1_0_0", PetKappaAdapter),
        ("equflash_29m_oam", EquFlashKappaAdapter),
        ("eqv2_s_dens_mp", FairchemKappaAdapter),
        ("dpa_4_0_pro_mptrj", StandardKappaAdapter),
        ("orb_v3", StandardKappaAdapter),
    ],
)
def test_get_kappa_adapter_dispatch(
    model_key: str, adapter_type: type[StandardKappaAdapter]
) -> None:
    """Backend families dispatch to their typed adapters."""
    assert isinstance(get_kappa_adapter(model_key), adapter_type)


def test_resumable_shards_strict_merge_and_artifacts(
    tmp_path: Path, shard_env: SimpleNamespace
) -> None:
    """Resume records, reject extras, and keep force sets separate."""
    shard_dir = shard_env.shard_dir
    atoms_by_id = shard_env.atoms_by_id
    settings = KappaSettings(relaxation_mode="none", max_steps=0, save_forces=True)
    computation_stub = shard_env.computation

    first_records = shard_env.run(settings=settings, n_shards=2)
    first_call_count = len(computation_stub.material_ids)
    assert first_call_count == len(first_records)
    shard_env.run(settings=settings, n_shards=2)
    assert len(computation_stub.material_ids) == first_call_count

    failed_record_path = record_path(shard_dir, first_records[0].material_id)
    failed_record = replace(
        read_kappa_record(failed_record_path), error="persisted failure"
    )
    atomic_write_gzip_json(failed_record_path, asdict(failed_record))
    shard_env.run(settings=settings, n_shards=2)
    assert len(computation_stub.material_ids) == first_call_count
    mismatched_retry = replace(
        failed_record,
        run_metadata={**failed_record.run_metadata, "hardware": "different hardware"},
    )
    atomic_write_gzip_json(failed_record_path, asdict(mismatched_retry))
    with pytest.raises(ValueError, match="different hardware or package versions"):
        shard_env.run(settings=settings, n_shards=2, retry_failures=True)
    atomic_write_gzip_json(failed_record_path, asdict(failed_record))
    shard_env.run(
        n_shards=2,
        settings=settings,
        retry_failures=True,
    )
    assert len(computation_stub.material_ids) == first_call_count + 1
    retried_record = read_kappa_record(failed_record_path)
    assert retried_record.run_metadata["retry_count"] == 1
    assert (
        retried_record.run_metadata["run_time_sec"]
        >= failed_record.run_metadata["run_time_sec"]
    )

    with pytest.raises(ValueError, match="missing="):
        merge_kappa_shards(shard_dir, model_key="test_model")
    shard_env.run(
        shard_index=1,
        n_shards=2,
        settings=settings,
    )
    mixed_record_path = record_path(shard_dir, first_records[0].material_id)
    original_record = read_kappa_record(mixed_record_path)
    mixed_metadata = dict(original_record.run_metadata)
    mixed_metadata["versions"] = {
        **cast("dict[str, str]", mixed_metadata["versions"]),
        "synthetic-backend": "different",
    }
    atomic_write_gzip_json(
        mixed_record_path,
        asdict(replace(original_record, run_metadata=mixed_metadata)),
    )
    with pytest.raises(ValueError, match="mixed or missing package versions"):
        merge_kappa_shards(shard_dir, model_key="test_model")
    atomic_write_gzip_json(mixed_record_path, asdict(original_record))

    merged_run = merge_kappa_shards(shard_dir, model_key="test_model")
    artifacts = write_kappa_artifacts(
        merged_run, pred_file_path=f"{tmp_path}/kappa.json.gz"
    )
    assert len(artifacts.predictions) == len(atoms_by_id)
    assert artifacts.force_file_path is not None
    assert artifacts.force_file_path == f"{tmp_path}/kappa-force-sets.json.gz"
    with gzip.open(artifacts.pred_file_path, mode="rt", encoding="utf-8") as file:
        prediction_rows = json.load(file)
    assert all("fc2_set" not in row and "fc3_set" not in row for row in prediction_rows)
    with gzip.open(artifacts.force_file_path, mode="rt", encoding="utf-8") as file:
        force_rows = json.load(file)
    assert len(force_rows) == len(atoms_by_id)
    with pytest.raises(ValueError, match="inside shard directory"):
        write_kappa_artifacts(
            merged_run,
            pred_file_path=f"{shard_dir}/records/final.json.gz",
        )

    atomic_write_gzip_json(f"{shard_dir}/records/extra.json.gz", {})
    with pytest.raises(ValueError, match="extra="):
        merge_kappa_shards(shard_dir, model_key="test_model")


def test_manifest_rejects_changed_settings(shard_env: SimpleNamespace) -> None:
    """A resumed run cannot silently change its scientific protocol."""
    settings = KappaSettings(relaxation_mode="none", max_steps=0)
    shard_env.run(settings=settings)
    with pytest.raises(ValueError, match=r"does not match.*settings"):
        shard_env.run(settings=replace(settings, force_max=2e-4))
    with pytest.raises(ValueError, match=r"does not match.*dtype"):
        shard_env.run(settings=settings, dtype="float32")


def test_merge_preserves_manifest_checkpoint_identity(
    shard_env: SimpleNamespace,
) -> None:
    """A merge validates records against the persisted checkpoint identity."""
    settings = KappaSettings(relaxation_mode="none", max_steps=0)
    checkpoint = "registry:model-v1"
    shard_env.run(settings=settings, checkpoint=checkpoint)

    merged_run = merge_kappa_shards(shard_env.shard_dir, model_key="test_model")
    assert len(merged_run.records) == len(shard_env.atoms_by_id)
    assert merged_run.manifest.checkpoint_hash == checkpoint_digest(checkpoint)


@pytest.mark.parametrize(
    "initial_yaml",
    ["metrics: {}\n", "metrics:\n  phonons: not available\n"],
)
def test_kappa_metric_yaml_replaces_missing_or_unavailable_phonons(
    tmp_path: Path, initial_yaml: str
) -> None:
    """Completed runs replace absent or unavailable phonon metadata."""
    kappa_metrics = update_temp_kappa_yaml(
        tmp_path, initial_yaml, replace_pred_file=True
    )
    assert kappa_metrics["κ_SRME"] == 0.25
    assert kappa_metrics["pred_file"] == "models/test/new.json.gz"


def test_kappa_metric_yaml_round_trip_updates_provenance(tmp_path: Path) -> None:
    """Complete-run metadata and sidecars round-trip through model YAML."""
    kappa_metrics = update_temp_kappa_yaml(
        tmp_path,
        (
            "metrics:\n"
            "  phonons:\n"
            "    kappa_103:\n"
            "      κ_SRME: 1.0\n"
            "      pred_file: models/test/new.json.gz\n"
            "      pred_file_url: https://example.com/old.json.gz\n"
            "      force_file_url: https://example.com/old-forces.json.gz\n"
            "      run_info_file_url: https://example.com/old-run-info.json\n"
        ),
        run_metadata={
            "hardware": "NVIDIA H200",
            "run_time_sec": 12.5,
            "max_rss_gb": 3.5,
            "max_gpu_mem_gb": 4.5,
        },
        force_file_path="models/test/forces.json.gz",
        run_info_path="models/test/run-info.json",
        replace_pred_file=True,
    )
    assert kappa_metrics["κ_SRME"] == 0.25
    assert kappa_metrics["κ_SRE"] == 0.125
    assert kappa_metrics["pred_file"] == "models/test/new.json.gz"
    assert kappa_metrics["pred_file_url"] is None
    assert kappa_metrics["force_file"] == "models/test/forces.json.gz"
    assert kappa_metrics["force_file_url"] is None
    assert kappa_metrics["run_info_file"] == "models/test/run-info.json"
    assert kappa_metrics["run_info_file_url"] is None
    assert kappa_metrics["hardware"] == "NVIDIA H200"
    assert kappa_metrics["run_time_sec"] == 12.5
    assert kappa_metrics["max_rss_gb"] == 3.5
    assert kappa_metrics["max_gpu_mem_gb"] == 4.5


def test_kappa_metric_yaml_clears_url_when_sidecar_path_changes(
    tmp_path: Path,
) -> None:
    """Changing a sidecar path invalidates its existing remote URL."""
    kappa_metrics = update_temp_kappa_yaml(
        tmp_path,
        (
            "metrics:\n"
            "  phonons:\n"
            "    kappa_103:\n"
            "      κ_SRME: 1.0\n"
            "      pred_file: models/test/pred.json.gz\n"
            "      pred_file_url: https://example.com/pred.json.gz\n"
            "      force_file: models/test/old-forces.json.gz\n"
            "      force_file_url: https://example.com/old-forces.json.gz\n"
        ),
        pred_file_path="models/test/pred.json.gz",
        force_file_path="models/test/new-forces.json.gz",
    )
    assert kappa_metrics["pred_file_url"] == "https://example.com/pred.json.gz"
    assert kappa_metrics["force_file"] == "models/test/new-forces.json.gz"
    assert kappa_metrics["force_file_url"] is None


def test_kappa_metric_yaml_replacement_clears_stale_sidecars(
    tmp_path: Path,
) -> None:
    """Replacing predictions removes provenance absent from the new run."""
    kappa_metrics = update_temp_kappa_yaml(
        tmp_path,
        (
            "metrics:\n"
            "  phonons:\n"
            "    kappa_103:\n"
            "      κ_SRME: 1.0\n"
            "      force_file: models/test/old-forces.json.gz\n"
            "      force_file_url: https://example.com/old-forces.json.gz\n"
            "      run_info_file: models/test/old-run-info.json\n"
            "      run_info_file_url: https://example.com/old-run-info.json\n"
            "      max_gpu_mem_gb: 9.0\n"
        ),
        replace_pred_file=True,
    )
    for stale_key in (
        "force_file",
        "force_file_url",
        "run_info_file",
        "run_info_file_url",
        "max_gpu_mem_gb",
    ):
        assert stale_key not in kappa_metrics
