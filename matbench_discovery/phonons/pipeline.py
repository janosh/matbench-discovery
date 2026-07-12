"""Resumable, sharded PhononDB thermal-conductivity benchmark pipeline."""

from __future__ import annotations

import copy
import glob
import gzip
import hashlib
import json
import math
import os
import socket
import sys
import tempfile
import time
import traceback
from collections.abc import Iterator, Mapping, Sequence
from contextlib import contextmanager
from dataclasses import asdict, dataclass, replace
from datetime import UTC, datetime
from typing import TYPE_CHECKING, Any, Literal
from urllib.parse import quote

import numpy as np
import pandas as pd
from filelock import FileLock
from pymatviz.enums import Key
from tqdm import tqdm

from matbench_discovery import ROOT
from matbench_discovery.ase_relax import canonical_filter_name, canonical_optimizer_name
from matbench_discovery.enums import DataFiles, Model
from matbench_discovery.hpc import (
    detect_hardware,
    merge_audit_metadata,
    partition_material_ids,
    peak_memory_gb,
    reset_gpu_peak_memory,
)
from matbench_discovery.hpc import package_versions as installed_package_versions
from matbench_discovery.phonons import check_imaginary_freqs
from matbench_discovery.phonons.adapters import StandardKappaAdapter, get_kappa_adapter
from matbench_discovery.phonons.schema import normalize_kappa_result

if TYPE_CHECKING:
    from ase import Atoms
    from ase.calculators.calculator import Calculator

KAPPA_RECORD_SCHEMA_VERSION = 1
KAPPA_MANIFEST_SCHEMA_VERSION = 2
PHONONDB_N_STRUCTURES = 103
KAPPA_MANIFEST_FILE = "manifest.json"
KAPPA_RECORD_DIR = "records"
KAPPA_FORCE_DIR = "forces"
DRY_RUN_MAX_FC3_EVALUATIONS = 8
KAPPA_PROTOCOL = "phonondb-v1"

RelaxationMode = Literal["single-stage", "two-stage", "none"]
PlusMinusSetting = bool | Literal["auto"]


@dataclass(frozen=True)
class KappaSettings:
    """Validated model-specific protocol for the PhononDB kappa benchmark."""

    displacement_distance: float = 0.01
    temperatures: tuple[float, ...] = (300.0,)
    ase_optimizer: str = "FIRE"
    ase_filter: str | None = "FrechetCellFilter"
    max_steps: int = 300
    force_max: float = 1e-4
    symprec: float = 1e-5
    enforce_relax_symm: bool = True
    conductivity_broken_symm: bool = False
    ignore_imaginary_freqs: bool = False
    is_plusminus: PlusMinusSetting = "auto"
    batch_size: int = 1
    max_atoms_per_batch: int | None = None
    relaxation_mode: RelaxationMode = "single-stage"
    save_forces: bool = False

    def __post_init__(self) -> None:
        """Canonicalize names and reject physically or structurally invalid values."""
        for field_name in ("displacement_distance", "force_max", "symprec"):
            value = getattr(self, field_name)
            if (
                isinstance(value, bool)
                or not isinstance(value, int | float)
                or not math.isfinite(value)
                or value <= 0
            ):
                raise ValueError(
                    f"{field_name} must be positive and finite, got {value}"
                )
        for field_name, minimum, requirement in (
            ("max_steps", 0, "a non-negative integer"),
            ("batch_size", 1, "a positive integer"),
        ):
            value = getattr(self, field_name)
            if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
                raise ValueError(f"{field_name} must be {requirement}, got {value}")
        if self.max_atoms_per_batch is not None and (
            isinstance(self.max_atoms_per_batch, bool)
            or not isinstance(self.max_atoms_per_batch, int)
            or self.max_atoms_per_batch < 1
        ):
            raise ValueError(
                "max_atoms_per_batch must be a positive integer or null, got "
                f"{self.max_atoms_per_batch}"
            )
        if not isinstance(self.ase_optimizer, str):
            raise TypeError("ase_optimizer must be a string")
        if self.ase_filter is not None and not isinstance(self.ase_filter, str):
            raise TypeError("ase_filter must be a string or null")
        if any(isinstance(temperature, bool) for temperature in self.temperatures):
            raise TypeError("temperatures must contain numbers, not booleans")
        temperatures = tuple(float(temperature) for temperature in self.temperatures)
        if not temperatures or any(
            not math.isfinite(temperature) or temperature <= 0
            for temperature in temperatures
        ):
            raise ValueError(
                f"temperatures must contain positive finite values, got {temperatures}"
            )
        if self.relaxation_mode not in {"single-stage", "two-stage", "none"}:
            raise ValueError(
                "relaxation_mode must be one of single-stage, two-stage, none; "
                f"got {self.relaxation_mode!r}"
            )
        for field_name in (
            "enforce_relax_symm",
            "conductivity_broken_symm",
            "ignore_imaginary_freqs",
            "save_forces",
        ):
            if not isinstance(getattr(self, field_name), bool):
                raise TypeError(f"{field_name} must be a boolean")
        if not (isinstance(self.is_plusminus, bool) or self.is_plusminus == "auto"):
            raise TypeError("is_plusminus must be a boolean or 'auto'")

        filter_name = canonical_filter_name(self.ase_filter)
        if filter_name not in {"FrechetCellFilter", "ExpCellFilter", None}:
            raise ValueError(f"Unsupported kappa ASE filter {filter_name!r}")
        object.__setattr__(self, "temperatures", temperatures)
        object.__setattr__(
            self, "ase_optimizer", canonical_optimizer_name(self.ase_optimizer)
        )
        object.__setattr__(self, "ase_filter", filter_name)

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any]) -> KappaSettings:
        """Construct settings from a versioned protocol and sparse overrides."""
        overrides = dict(data)
        protocol = overrides.pop("protocol", None)
        if protocol != KAPPA_PROTOCOL:
            raise ValueError(
                f"hyperparams.kappa protocol must be {KAPPA_PROTOCOL!r}, "
                f"got {protocol!r}"
            )
        return cls(**overrides)

    @classmethod
    def from_model(cls, model_ref: str) -> KappaSettings:
        """Load a model's versioned ``hyperparams.kappa`` protocol."""
        model = Model.from_ref(model_ref)
        hyperparams = model.metadata.get("hyperparams")
        if not isinstance(hyperparams, dict):
            raise TypeError(f"{model.name} hyperparams must be a mapping")
        kappa_data = hyperparams.get("kappa")
        if not isinstance(kappa_data, dict):
            raise TypeError(f"{model.name} has no verified hyperparams.kappa settings")
        return cls.from_mapping(kappa_data)

    @property
    def digest(self) -> str:
        """Return a deterministic SHA-256 hash of the protocol settings."""
        payload = json.dumps(
            json_ready(asdict(self)),
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        )
        return hashlib.sha256(payload.encode()).hexdigest()


@dataclass(frozen=True)
class KappaRunManifest:
    """Immutable identity and deterministic assignments for one sharded run."""

    model_key: str
    dataset_path: str
    dataset_hash: str
    material_ids: tuple[str, ...]
    material_id_hash: str
    n_shards: int
    assignments: tuple[tuple[str, ...], ...]
    settings: KappaSettings
    settings_hash: str
    adapter_name: str
    source_hash: str
    checkpoint_hash: str | None
    dtype: str
    device: str
    dry_run: bool = False
    schema_version: int = KAPPA_MANIFEST_SCHEMA_VERSION

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> KappaRunManifest:
        """Hydrate and validate a manifest from JSON data."""
        schema_version = int(data.get("schema_version", -1))
        settings_data = data.get("settings")
        if not isinstance(settings_data, dict):
            raise TypeError("Kappa manifest settings must be a mapping")
        manifest = cls(
            model_key=str(data["model_key"]),
            dataset_path=str(data["dataset_path"]),
            dataset_hash=str(data["dataset_hash"]),
            material_ids=tuple(map(str, data["material_ids"])),
            material_id_hash=str(data["material_id_hash"]),
            n_shards=int(data["n_shards"]),
            assignments=tuple(
                tuple(map(str, assignment)) for assignment in data["assignments"]
            ),
            settings=KappaSettings(**settings_data),
            settings_hash=str(data["settings_hash"]),
            adapter_name=str(data["adapter_name"]),
            source_hash=str(data["source_hash"]),
            checkpoint_hash=(
                str(data["checkpoint_hash"])
                if data.get("checkpoint_hash") is not None
                else None
            ),
            dtype=str(data["dtype"]),
            device=str(data["device"]),
            dry_run=bool(data.get("dry_run", False)),
            schema_version=schema_version,
        )
        manifest.validate()
        return manifest

    def validate(self) -> None:
        """Reject malformed assignments, hashes, and unsupported schema versions."""
        if self.schema_version != KAPPA_MANIFEST_SCHEMA_VERSION:
            raise ValueError(
                f"Unsupported kappa manifest schema {self.schema_version}, "
                f"expected {KAPPA_MANIFEST_SCHEMA_VERSION}"
            )
        if self.n_shards < 1 or len(self.assignments) != self.n_shards:
            raise ValueError(
                f"Manifest has {len(self.assignments)} assignments for "
                f"{self.n_shards} shards"
            )
        assigned_ids = [
            material_id for assignment in self.assignments for material_id in assignment
        ]
        if len(assigned_ids) != len(set(assigned_ids)):
            raise ValueError("Kappa manifest assignments contain duplicate IDs")
        if set(assigned_ids) != set(self.material_ids):
            raise ValueError("Kappa manifest assignments do not cover its material IDs")
        if self.material_id_hash != material_id_digest(self.material_ids):
            raise ValueError("Kappa manifest material_id_hash is invalid")
        if self.settings_hash != self.settings.digest:
            raise ValueError("Kappa manifest settings_hash is invalid")
        if self.dtype not in {"float32", "float64"}:
            raise ValueError(f"Kappa manifest has invalid dtype {self.dtype!r}")
        if self.device not in {"cpu", "cuda"}:
            raise ValueError(f"Kappa manifest has invalid device {self.device!r}")


@dataclass(frozen=True)
class KappaComputation:
    """One in-memory material result before its atomic record is written."""

    result: dict[str, Any]
    forces: dict[str, Any] | None
    error: str | None


@dataclass(frozen=True)
class KappaRecord:
    """One durable per-material result and its execution provenance."""

    material_id: str
    result: dict[str, Any]
    error: str | None
    run_metadata: dict[str, Any]
    force_file: str | None = None
    schema_version: int = KAPPA_RECORD_SCHEMA_VERSION

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> KappaRecord:
        """Hydrate and validate one gzip record."""
        result = data.get("result")
        metadata = data.get("run_metadata")
        if not isinstance(result, dict) or not isinstance(metadata, dict):
            raise TypeError("Kappa record result and run_metadata must be mappings")
        record = cls(
            material_id=str(data["material_id"]),
            result=normalize_kappa_result(result),
            error=str(data["error"]) if data.get("error") is not None else None,
            run_metadata=dict(metadata),
            force_file=(
                str(data["force_file"]) if data.get("force_file") is not None else None
            ),
            schema_version=int(data["schema_version"]),
        )
        if record.schema_version != KAPPA_RECORD_SCHEMA_VERSION:
            raise ValueError(
                f"Unsupported kappa record schema {record.schema_version}, "
                f"expected {KAPPA_RECORD_SCHEMA_VERSION}"
            )
        result_id = record.result.get(str(Key.mat_id))
        if result_id != record.material_id:
            raise ValueError(
                f"Kappa record ID mismatch: {record.material_id!r} != {result_id!r}"
            )
        return record


@dataclass(frozen=True)
class MergedKappaRun:
    """A complete, strictly validated set of PhononDB material records."""

    shard_dir: str
    manifest: KappaRunManifest
    records: tuple[KappaRecord, ...]
    run_metadata: dict[str, Any]


@dataclass(frozen=True)
class KappaArtifacts:
    """Final prediction, force-set, and provenance artifact paths."""

    pred_file_path: str
    run_info_path: str
    force_file_path: str | None
    predictions: pd.DataFrame
    n_failed: int


def material_id_digest(material_ids: Sequence[str]) -> str:
    """Return an order-sensitive digest for the canonical dataset ID sequence."""
    return hashlib.sha256("\n".join(map(str, material_ids)).encode()).hexdigest()


def file_sha256(path: str) -> str:
    """Stream a file into a SHA-256 digest."""
    with open(path, mode="rb") as file:
        return hashlib.file_digest(file, "sha256").hexdigest()


def checkpoint_digest(checkpoint: str | None) -> str | None:
    """Hash checkpoint bytes when local, otherwise hash its supplied identifier."""
    if checkpoint is None:
        return None
    if os.path.isfile(checkpoint):
        if sys.platform == "win32":
            return file_sha256(checkpoint)
        file_stats = os.stat(checkpoint)
        cache_path = f"{checkpoint}.sha256.json"
        try:
            with FileLock(f"{cache_path}.lock"):
                if os.path.isfile(cache_path):
                    with open(cache_path, encoding="utf-8") as file:
                        cache_data = json.load(file)
                    if isinstance(cache_data, dict):
                        cached_hash = cache_data.get("sha256")
                        if (
                            cache_data.get("size") == file_stats.st_size
                            and cache_data.get("mtime_ns") == file_stats.st_mtime_ns
                            and cache_data.get("ctime_ns") == file_stats.st_ctime_ns
                            and isinstance(cached_hash, str)
                            and len(cached_hash) == 64
                            and all(char in "0123456789abcdef" for char in cached_hash)
                        ):
                            return cached_hash
                digest = file_sha256(checkpoint)
                atomic_write_json(
                    cache_path,
                    {
                        "size": file_stats.st_size,
                        "mtime_ns": file_stats.st_mtime_ns,
                        "ctime_ns": file_stats.st_ctime_ns,
                        "sha256": digest,
                    },
                )
                return digest
        except (OSError, TypeError, ValueError):
            return file_sha256(checkpoint)
    return hashlib.sha256(checkpoint.encode()).hexdigest()


def kappa_source_digest(adapter: StandardKappaAdapter, model_key: str) -> str:
    """Hash shared runner, pipeline, calculator, schema, and adapter source files."""
    from matbench_discovery import ase_relax, calculators, hpc, phonons
    from matbench_discovery.phonons import schema, thermal_conductivity

    source_paths = {
        os.path.realpath(ase_relax.__file__),
        os.path.realpath(__file__),
        os.path.realpath(calculators.__file__),
        os.path.realpath(hpc.__file__),
        os.path.realpath(phonons.__file__),
        os.path.realpath(schema.__file__),
        os.path.realpath(thermal_conductivity.__file__),
        os.path.realpath(f"{ROOT}/models/run_kappa.py"),
    }
    if (calc_spec := calculators.CALCULATORS.get(model_key)) and calc_spec.project:
        source_paths.add(os.path.realpath(f"{ROOT}/{calc_spec.project}/pyproject.toml"))
    for adapter_class in type(adapter).mro():
        if (module := sys.modules.get(adapter_class.__module__)) and (
            module_file := getattr(module, "__file__", None)
        ):
            source_paths.add(os.path.realpath(module_file))

    source_hasher = hashlib.sha256()
    for source_path in sorted(source_paths):
        if not os.path.isfile(source_path):
            raise FileNotFoundError(f"Kappa source file not found: {source_path}")
        source_hasher.update(os.path.basename(source_path).encode())
        with open(source_path, mode="rb") as file:
            source_hasher.update(file.read())
    return source_hasher.hexdigest()


def material_id_from_atoms(atoms: Atoms) -> str:
    """Extract and canonicalize a PhononDB material ID from ASE metadata."""
    for key in (str(Key.mat_id), "mp_id"):
        material_id = atoms.info.get(key)
        if isinstance(material_id, str) and material_id:
            atoms.info[str(Key.mat_id)] = material_id
            return material_id
    raise ValueError(
        f"Structure has no material_id/mp_id in info keys {sorted(atoms.info)}"
    )


def load_phonondb_atoms(
    dataset_path: str | None = None,
    *,
    dry_run: bool = False,
    dry_run_size: int = 1,
) -> dict[str, Atoms]:
    """Load the canonical 103 structures and optionally retain a smoke-test subset."""
    from ase.io import read

    path = dataset_path or DataFiles.phonondb_pbe_103_structures.path
    atoms_list = read(path, format="extxyz", index=":")
    atoms_by_id: dict[str, Atoms] = {}
    for atoms in atoms_list:
        material_id = material_id_from_atoms(atoms)
        if material_id in atoms_by_id:
            raise ValueError(f"Duplicate PhononDB material ID {material_id!r}")
        atoms_by_id[material_id] = atoms
    if len(atoms_by_id) != PHONONDB_N_STRUCTURES:
        raise ValueError(
            f"Expected {PHONONDB_N_STRUCTURES} unique PhononDB structures, "
            f"got {len(atoms_by_id)}"
        )
    if not dry_run:
        return atoms_by_id
    if not 1 <= dry_run_size <= len(atoms_by_id):
        raise ValueError(
            f"dry_run_size must be in [1, {len(atoms_by_id)}], got {dry_run_size}"
        )
    smallest_structures = sorted(
        atoms_by_id.items(), key=lambda item: (len(item[1]), item[0])
    )[:dry_run_size]
    return dict(smallest_structures)


def build_kappa_manifest(
    *,
    model_key: str,
    atoms_by_id: Mapping[str, Atoms],
    dataset_path: str,
    settings: KappaSettings,
    adapter: StandardKappaAdapter,
    n_shards: int,
    dtype: str,
    device: str,
    checkpoint: str | None = None,
    dry_run: bool = False,
) -> KappaRunManifest:
    """Build the deterministic immutable manifest expected by every shard."""
    material_ids = tuple(atoms_by_id)
    assignments = tuple(map(tuple, partition_material_ids(atoms_by_id, n_shards)))
    effective_settings = dry_run_settings(settings) if dry_run else settings
    manifest = KappaRunManifest(
        model_key=model_key,
        dataset_path=os.path.abspath(dataset_path),
        dataset_hash=file_sha256(dataset_path),
        material_ids=material_ids,
        material_id_hash=material_id_digest(material_ids),
        n_shards=n_shards,
        assignments=assignments,
        settings=effective_settings,
        settings_hash=effective_settings.digest,
        adapter_name=adapter.name,
        source_hash=kappa_source_digest(adapter, model_key),
        checkpoint_hash=checkpoint_digest(checkpoint),
        dtype=dtype,
        device=device,
        dry_run=dry_run,
    )
    manifest.validate()
    return manifest


def ensure_kappa_manifest(shard_dir: str, expected_manifest: KappaRunManifest) -> None:
    """Atomically create a run manifest or require an exact existing match."""
    os.makedirs(shard_dir, exist_ok=True)
    manifest_path = f"{shard_dir}/{KAPPA_MANIFEST_FILE}"
    with FileLock(f"{manifest_path}.lock"):
        if os.path.isfile(manifest_path):
            existing_manifest = read_kappa_manifest(manifest_path)
            if existing_manifest != expected_manifest:
                differing_fields = [
                    field_name
                    for field_name in expected_manifest.__dataclass_fields__
                    if getattr(existing_manifest, field_name)
                    != getattr(expected_manifest, field_name)
                ]
                raise ValueError(
                    "Existing kappa run manifest does not match this run; "
                    f"different fields: {differing_fields}"
                )
        else:
            atomic_write_json(manifest_path, asdict(expected_manifest))


def read_kappa_manifest(path: str) -> KappaRunManifest:
    """Read and validate one immutable kappa run manifest."""
    with open(path, encoding="utf-8") as file:
        data = json.load(file)
    if not isinstance(data, dict):
        raise TypeError(f"Kappa manifest {path} must contain a mapping")
    return KappaRunManifest.from_dict(data)


def record_filename(material_id: str) -> str:
    """Return a reversible filesystem-safe filename for one material ID."""
    return f"{quote(material_id, safe='-_.')}.json.gz"


def record_path(shard_dir: str, material_id: str) -> str:
    """Return the durable record path for one material."""
    return f"{shard_dir}/{KAPPA_RECORD_DIR}/{record_filename(material_id)}"


def read_kappa_record(path: str) -> KappaRecord:
    """Read and validate one atomic gzip material record."""
    with gzip.open(path, mode="rt", encoding="utf-8") as file:
        data = json.load(file)
    if not isinstance(data, dict):
        raise TypeError(f"Kappa record {path} must contain a mapping")
    return KappaRecord.from_dict(data)


def should_calculate_conductivity(
    *,
    has_imaginary_modes: bool,
    broken_symmetry: bool,
    settings: KappaSettings,
) -> bool:
    """Apply the historical imaginary-mode and broken-symmetry gate policy."""
    if settings.ignore_imaginary_freqs:
        return True
    return not has_imaginary_modes and (
        settings.conductivity_broken_symm or not broken_symmetry
    )


def calculate_kappa_for_structure(
    *,
    atoms: Atoms,
    calculator: Calculator,
    settings: KappaSettings,
    adapter: StandardKappaAdapter,
    log_file: str | None = None,
    max_fc3_evaluations: int | None = None,
) -> KappaComputation:
    """Run the common relax → FC2 → gate → FC3 → conductivity flow."""
    material_id = material_id_from_atoms(atoms)
    atoms = copy.deepcopy(atoms)
    formula = atoms.info.get("name")
    result: dict[str, Any] = {
        str(Key.mat_id): material_id,
        str(Key.formula): str(formula) if formula else atoms.get_chemical_formula(),
        "errors": [],
        "error_traceback": [],
    }
    forces: dict[str, Any] | None = None
    error: str | None = None
    error_stage = "Relax"
    try:
        relaxation = adapter.relax(atoms, calculator, settings, log_file=log_file)
        atoms = relaxation.atoms
        result |= {
            str(Key.init_spg_num): relaxation.initial_spg_num,
            str(Key.final_spg_num): relaxation.final_spg_num,
            "max_stress": relaxation.max_stress,
            "reached_max_steps": relaxation.reached_max_steps,
            "n_steps": relaxation.n_steps,
            "broken_symmetry": (relaxation.initial_spg_num != relaxation.final_spg_num),
            "redirected_to_symmetry": relaxation.redirected_to_symmetry,
        }
        error_stage = "ForceConstant"
        phono3py = adapter.init_phono3py(atoms, settings)
        phono3py, fc2_set, frequencies = adapter.calculate_fc2(
            phono3py,
            calculator,
            settings,
        )
        has_imaginary_modes = check_imaginary_freqs(frequencies)
        result |= {
            str(Key.has_imag_ph_modes): has_imaginary_modes,
            str(Key.ph_freqs): frequencies,
        }
        should_calculate_conductivity_value = should_calculate_conductivity(
            has_imaginary_modes=has_imaginary_modes,
            broken_symmetry=bool(result["broken_symmetry"]),
            settings=settings,
        )
        fc3_set: np.ndarray | list[Any] = []
        if should_calculate_conductivity_value:
            fc3_set = adapter.calculate_fc3(
                phono3py,
                calculator,
                settings,
                max_evaluations=max_fc3_evaluations,
            )
            phono3py.forces = np.asarray(fc3_set)
            phono3py.produce_fc3(symmetrize_fc3r=True)
        if settings.save_forces:
            forces = {"fc2_set": fc2_set, "fc3_set": fc3_set}
        if should_calculate_conductivity_value:
            error_stage = "Conductivity"
            from matbench_discovery.phonons import thermal_conductivity as ltc

            _phono3py, kappa_data, _conductivity = ltc.calculate_conductivity(
                phono3py, temperatures=settings.temperatures
            )
            result.update(kappa_data)
        else:
            result["conductivity_skipped"] = True
    except Exception as exc:  # noqa: BLE001 - every material failure is persisted
        error = f"{error_stage}Error: {type(exc).__name__}: {exc}"
        result["errors"].append(error)
        result["error_traceback"].append(traceback.format_exc())
    return KappaComputation(normalize_kappa_result(result), forces, error)


def run_kappa_shard(
    *,
    calculator: Calculator,
    model_key: str,
    atoms_by_id: Mapping[str, Atoms],
    dataset_path: str,
    shard_dir: str,
    shard_index: int,
    n_shards: int,
    settings: KappaSettings,
    adapter: StandardKappaAdapter | None = None,
    checkpoint: str | None = None,
    dtype: str = "float64",
    device: str = "cpu",
    retry_failures: bool = False,
    dry_run: bool = False,
) -> tuple[KappaRecord, ...]:
    """Run or resume one atom-balanced shard with atomic records per material."""
    adapter = adapter or get_kappa_adapter(model_key)
    manifest = build_kappa_manifest(
        model_key=model_key,
        atoms_by_id=atoms_by_id,
        dataset_path=dataset_path,
        settings=settings,
        adapter=adapter,
        n_shards=n_shards,
        dtype=dtype,
        device=device,
        checkpoint=checkpoint,
        dry_run=dry_run,
    )
    settings = manifest.settings
    ensure_kappa_manifest(shard_dir, manifest)
    if not 0 <= shard_index < n_shards:
        raise ValueError(
            f"shard_index must be in [0, {n_shards - 1}], got {shard_index}"
        )

    os.makedirs(f"{shard_dir}/{KAPPA_RECORD_DIR}", exist_ok=True)
    if settings.save_forces:
        os.makedirs(f"{shard_dir}/{KAPPA_FORCE_DIR}", exist_ok=True)
    calculator = adapter.prepare_calculator(calculator, settings)
    hardware = detect_hardware()
    versions = package_versions(model_key)
    assigned_ids = manifest.assignments[shard_index]

    for material_id in tqdm(
        assigned_ids, desc=f"{model_key} kappa shard {shard_index}"
    ):
        output_path = record_path(shard_dir, material_id)
        with FileLock(f"{output_path}.lock"):
            existing_record = (
                read_kappa_record(output_path) if os.path.isfile(output_path) else None
            )
            if existing_record is not None:
                validate_record_provenance(existing_record, manifest)
                force_exists = (
                    existing_record.force_file is not None
                    and os.path.isfile(f"{shard_dir}/{existing_record.force_file}")
                )
                if existing_record.error is None and (
                    not settings.save_forces or force_exists
                ):
                    continue
                if existing_record.error is not None and not retry_failures:
                    continue
                if existing_record.error is not None and (
                    existing_record.run_metadata.get("hardware") != hardware
                    or existing_record.run_metadata.get("versions") != versions
                ):
                    raise ValueError(
                        "Cannot combine a retry with different hardware or "
                        "package versions"
                    )

            reset_gpu_peak_memory()
            start_time = time.perf_counter()
            computation = calculate_kappa_for_structure(
                atoms=atoms_by_id[material_id],
                calculator=calculator,
                settings=settings,
                adapter=adapter,
                log_file=(
                    f"{shard_dir}/relaxations/{quote(material_id, safe='-_.')}.log"
                ),
                max_fc3_evaluations=(
                    DRY_RUN_MAX_FC3_EVALUATIONS if manifest.dry_run else None
                ),
            )
            elapsed = time.perf_counter() - start_time
            metadata = material_run_metadata(
                elapsed,
                hardware=hardware,
                versions=versions,
                manifest=manifest,
                shard_index=shard_index,
                previous_metadata=(
                    existing_record.run_metadata
                    if existing_record is not None and existing_record.error is not None
                    else None
                ),
            )
            relative_force_path = None
            if computation.forces is not None:
                absolute_force_path = (
                    f"{shard_dir}/{KAPPA_FORCE_DIR}/{record_filename(material_id)}"
                )
                atomic_write_gzip_json(absolute_force_path, computation.forces)
                relative_force_path = os.path.relpath(absolute_force_path, shard_dir)
            record = KappaRecord(
                material_id=material_id,
                result=computation.result,
                error=computation.error,
                run_metadata=metadata,
                force_file=relative_force_path,
            )
            atomic_write_gzip_json(output_path, asdict(record))

    records = tuple(
        read_kappa_record(record_path(shard_dir, material_id))
        for material_id in assigned_ids
    )
    for record in records:
        validate_record_provenance(record, manifest)
    return records


def _manifest_provenance(manifest: KappaRunManifest) -> dict[str, Any]:
    """Return immutable manifest fields copied into record and run metadata."""
    return {
        "model_key": manifest.model_key,
        "dataset_hash": manifest.dataset_hash,
        "settings_hash": manifest.settings_hash,
        "source_hash": manifest.source_hash,
        "checkpoint_hash": manifest.checkpoint_hash,
        "dtype": manifest.dtype,
        "device": manifest.device,
    }


def validate_record_provenance(record: KappaRecord, manifest: KappaRunManifest) -> None:
    """Require a record's immutable run hashes to match its manifest."""
    expected = _manifest_provenance(manifest)
    mismatches = {
        key: (record.run_metadata.get(key), value)
        for key, value in expected.items()
        if record.run_metadata.get(key) != value
    }
    if mismatches:
        raise ValueError(
            f"Record {record.material_id!r} provenance does not match manifest: "
            f"{mismatches}"
        )


def merge_kappa_shards(
    shard_dir: str,
    *,
    model_key: str,
    expected_material_ids: Sequence[str] | None = None,
    require_103: bool = False,
) -> MergedKappaRun:
    """Strictly merge complete records, rejecting missing, extra, or mixed runs."""
    manifest = read_kappa_manifest(f"{shard_dir}/{KAPPA_MANIFEST_FILE}")
    if manifest.model_key != model_key:
        raise ValueError(
            f"Manifest model {manifest.model_key!r} does not match {model_key!r}"
        )

    expected_ids_source = (
        expected_material_ids
        if expected_material_ids is not None
        else manifest.material_ids
    )
    expected_ids = tuple(map(str, expected_ids_source))
    if len(expected_ids) != len(set(expected_ids)):
        raise ValueError("Expected kappa material IDs contain duplicates")
    if expected_ids != manifest.material_ids:
        raise ValueError(
            "Expected kappa material IDs do not exactly match the manifest"
        )
    if require_103 and len(expected_ids) != PHONONDB_N_STRUCTURES:
        raise ValueError(
            f"A complete kappa merge requires {PHONONDB_N_STRUCTURES} IDs, "
            f"got {len(expected_ids)}"
        )

    expected_paths = {
        os.path.normpath(record_path(shard_dir, material_id))
        for material_id in expected_ids
    }
    actual_paths = {
        os.path.normpath(path)
        for path in glob.glob(f"{shard_dir}/{KAPPA_RECORD_DIR}/*.json.gz")
    }
    if actual_paths != expected_paths:
        missing_paths = sorted(expected_paths - actual_paths)
        extra_paths = sorted(actual_paths - expected_paths)
        raise ValueError(
            "Kappa records do not match complete manifest coverage: "
            f"missing={missing_paths}, extra={extra_paths}"
        )

    records = tuple(
        read_kappa_record(record_path(shard_dir, material_id))
        for material_id in expected_ids
    )
    for expected_id, record in zip(expected_ids, records, strict=True):
        if record.material_id != expected_id:
            raise ValueError(
                f"Kappa record order/ID mismatch: {record.material_id!r} != "
                f"{expected_id!r}"
            )
        validate_record_provenance(record, manifest)

    if manifest.settings.save_forces:
        force_root = os.path.realpath(f"{shard_dir}/{KAPPA_FORCE_DIR}")
        for record in records:
            if record.error is None and record.force_file is None:
                raise ValueError(
                    f"Successful record {record.material_id!r} is missing force data"
                )
            if record.force_file is None:
                continue
            source_path = os.path.realpath(f"{shard_dir}/{record.force_file}")
            if os.path.commonpath(
                (force_root, source_path)
            ) != force_root or not os.path.isfile(source_path):
                raise ValueError(
                    f"Invalid force sidecar for {record.material_id}: {source_path}"
                )

    run_metadata: dict[str, Any] = {
        **merge_audit_metadata(
            [record.run_metadata for record in records], strict=True
        ),
        **_manifest_provenance(manifest),
        "n_shards": manifest.n_shards,
        "n_structures": len(records),
        "n_failed": sum(record.error is not None for record in records),
        "dataset_path": manifest.dataset_path,
        "settings": asdict(manifest.settings),
        "adapter": manifest.adapter_name,
    }
    return MergedKappaRun(
        shard_dir=os.path.abspath(shard_dir),
        manifest=manifest,
        records=records,
        run_metadata=run_metadata,
    )


def write_kappa_artifacts(
    merged_run: MergedKappaRun,
    *,
    pred_file_path: str,
    force_file_path: str | None = None,
    run_info_path: str | None = None,
) -> KappaArtifacts:
    """Write normalized predictions plus separate force and provenance sidecars."""
    artifact_stem = (
        pred_file_path.removesuffix(".json.gz")
        if pred_file_path.endswith(".json.gz")
        else os.path.splitext(pred_file_path)[0]
    )
    run_info_path = run_info_path or f"{artifact_stem}-run-info.json"
    candidate_force_path = (
        force_file_path or f"{artifact_stem}-force-sets.json.gz"
        if merged_run.manifest.settings.save_forces
        else None
    )
    artifact_paths = [
        os.path.realpath(path)
        for path in (pred_file_path, run_info_path, candidate_force_path)
        if path is not None
    ]
    if len(set(artifact_paths)) != len(artifact_paths):
        raise ValueError("Kappa artifact output paths must be distinct")
    shard_dir = os.path.realpath(merged_run.shard_dir)
    for artifact_path in artifact_paths:
        if os.path.commonpath((shard_dir, artifact_path)) == shard_dir:
            raise ValueError(
                f"Refusing to write final artifact inside shard directory: "
                f"{artifact_path}"
            )

    rows = [record.result for record in merged_run.records]
    atomic_write_gzip_json(pred_file_path, rows)

    atomic_write_json(
        run_info_path,
        {
            "run_metadata": merged_run.run_metadata,
            "manifest": asdict(merged_run.manifest),
        },
    )

    if candidate_force_path is not None:
        stream_force_artifact(
            candidate_force_path,
            records=merged_run.records,
            shard_dir=merged_run.shard_dir,
        )

    df_predictions = pd.DataFrame(rows)
    df_predictions = df_predictions.set_index(str(Key.mat_id), drop=False)
    return KappaArtifacts(
        pred_file_path=pred_file_path,
        run_info_path=run_info_path,
        force_file_path=candidate_force_path,
        predictions=df_predictions,
        n_failed=sum(record.error is not None for record in merged_run.records),
    )


@contextmanager
def _atomic_output_path(path: str) -> Iterator[str]:
    """Yield a temporary sibling path and atomically promote it on success."""
    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)
    file_descriptor, temporary_path = tempfile.mkstemp(
        prefix=f".{os.path.basename(path)}.", dir=parent or "."
    )
    os.close(file_descriptor)
    try:
        yield temporary_path
        os.replace(temporary_path, path)
    finally:
        if os.path.isfile(temporary_path):
            os.remove(temporary_path)


def stream_force_artifact(
    path: str,
    *,
    records: Sequence[KappaRecord],
    shard_dir: str,
) -> None:
    """Atomically stream per-material force sets into one gzip JSON array."""
    with _atomic_output_path(path) as temporary_path:
        n_written = 0
        with gzip.open(temporary_path, mode="wt", encoding="utf-8") as output_file:
            output_file.write("[")
            for record in records:
                if record.force_file is None:
                    continue
                source_path = f"{shard_dir}/{record.force_file}"
                with gzip.open(source_path, mode="rt", encoding="utf-8") as source_file:
                    force_data = json.load(source_file)
                if not isinstance(force_data, dict):
                    raise TypeError(
                        f"Force sidecar for {record.material_id} must be a mapping"
                    )
                if n_written:
                    output_file.write(",")
                json.dump(
                    json_ready({str(Key.mat_id): record.material_id, **force_data}),
                    output_file,
                    sort_keys=True,
                    allow_nan=False,
                )
                n_written += 1
            output_file.write("]")


def package_versions(model_key: str | None = None) -> dict[str, str]:
    """Collect core and model-backend package versions used by an execution."""
    from packaging.requirements import InvalidRequirement, Requirement

    from matbench_discovery.calculators import CALCULATORS

    package_names = {
        "ase",
        "jax",
        "matbench-discovery",
        "moyopy",
        "numpy",
        "phono3py",
        "phonopy",
        "scipy",
        "spglib",
        "tensorflow",
        "torch",
    }
    if model_key in CALCULATORS:
        for dependency in CALCULATORS[model_key].deps:
            try:
                package_names.add(Requirement(dependency).name)
            except InvalidRequirement:
                continue
    return installed_package_versions(package_names)


def material_run_metadata(
    run_time_sec: float,
    *,
    hardware: str,
    versions: Mapping[str, str],
    manifest: KappaRunManifest,
    shard_index: int,
    previous_metadata: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Collect timing, memory, run hashes, and audit provenance for one material."""
    current_memory = peak_memory_gb()
    if previous_metadata is not None:
        run_time_sec += float(previous_metadata.get("run_time_sec", 0))
        for memory_key in ("max_rss_gb", "max_gpu_mem_gb"):
            memory_values = [
                float(value)
                for value in (
                    current_memory.get(memory_key),
                    previous_metadata.get(memory_key),
                )
                if isinstance(value, int | float)
            ]
            if memory_values:
                current_memory[memory_key] = max(memory_values)
    metadata: dict[str, Any] = {
        **_manifest_provenance(manifest),
        "hardware": hardware,
        "run_time_sec": round(run_time_sec, 2),
        "shard_index": shard_index,
        "completed_at": datetime.now(UTC).isoformat(),
        "hostname": socket.gethostname(),
        "versions": dict(versions),
        **current_memory,
    }
    if previous_metadata is not None:
        metadata["retry_count"] = int(previous_metadata.get("retry_count", 0)) + 1
    for environment_key in ("SLURM_JOB_ID", "SLURM_ARRAY_TASK_ID"):
        if environment_value := os.getenv(environment_key):
            metadata[environment_key.casefold()] = environment_value
    return metadata


def dry_run_settings(settings: KappaSettings) -> KappaSettings:
    """Cap relaxation and force-sidecar work for a one-structure smoke test."""
    return replace(
        settings,
        max_steps=min(settings.max_steps, 2),
        temperatures=settings.temperatures[:1],
        save_forces=False,
    )


def json_ready(value: object) -> object:
    """Recursively convert scientific Python values to strict JSON data."""
    if isinstance(value, Mapping):
        return {str(key): json_ready(item) for key, item in value.items()}
    if isinstance(value, np.ndarray):
        return json_ready(value.tolist())  # ty: ignore[no-matching-overload]
    if isinstance(value, np.generic):
        return json_ready(value.item())
    if isinstance(value, tuple | list):
        return [json_ready(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, str | int | float | bool) or value is None:
        return value
    raise TypeError(f"Object of type {type(value).__name__} is not JSON serializable")


def atomic_write_json(path: str, data: object) -> None:
    """Atomically write strict JSON to a plain-text file."""
    with (
        _atomic_output_path(path) as temporary_path,
        open(temporary_path, mode="w", encoding="utf-8") as file,
    ):
        json.dump(json_ready(data), file, sort_keys=True, allow_nan=False)
        file.flush()
        os.fsync(file.fileno())


def atomic_write_gzip_json(path: str, data: object) -> None:
    """Atomically write strict JSON to a gzip file."""
    with (
        _atomic_output_path(path) as temporary_path,
        gzip.open(temporary_path, mode="wt", encoding="utf-8") as file,
    ):
        json.dump(json_ready(data), file, sort_keys=True, allow_nan=False)
        file.flush()
        os.fsync(file.fileno())
