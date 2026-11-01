"""Shared WBM discovery relaxation, sharding, merge, and artifact pipeline."""

from __future__ import annotations

import copy
import json
import os
import time
from dataclasses import asdict, dataclass, field, replace
from typing import TYPE_CHECKING, Any

import numpy as np
import pandas as pd
from filelock import FileLock
from pymatgen.core import Structure
from pymatgen.entries.compatibility import MaterialsProject2020Compatibility
from pymatgen.entries.computed_entries import ComputedStructureEntry
from pymatgen.io.ase import AseAtomsAdaptor
from pymatviz.enums import Key
from tqdm import tqdm

from matbench_discovery import atomic_output_path
from matbench_discovery.ase_relax import (
    canonical_filter_name,
    canonical_optimizer_name,
    resolve_cell_filter,
    resolve_optimizer,
)
from matbench_discovery.enums import DataFiles, Model, TestSubset
from matbench_discovery.hpc import (
    check_resume_labels,
    material_id_from_atoms,
    material_id_hash,
    memory_peaks,
    merge_audit_metadata,
    reset_gpu_peak_memory,
)

if TYPE_CHECKING:
    from collections.abc import Callable, Mapping, Sequence

    from ase import Atoms
    from ase.calculators.calculator import Calculator

# schema 2: segment_start events (process runtime), records with run time + memory
DISCOVERY_SHARD_SCHEMA_VERSION = 2
DISCOVERY_PRED_COL = "e_form_per_atom"
DISCOVERY_STRUCT_COL = "structure"
DISCOVERY_ID_COL = "material_id"
DEFAULT_DRY_RUN_STRUCTURES = 4
ARCHIVED_DISCOVERY_MODELS: dict[str, str] = {
    **dict.fromkeys(
        (
            "alchembert",
            "alignn",
            "bowsr",
            "cgcnn",
            "cgcnn_p",
            "ema_gnn",
            "esnet",
            "megnet",
            "voronoi_rf",
            "wrenformer",
        ),
        "archived direct predictor without an ASE calculator",
    ),
    **dict.fromkeys(
        ("equflash_29m_oam", "equflashv2_45m_oam"),
        "batched runtime requires an unsupported --no-deps install",
    ),
    **dict.fromkeys(
        ("equiformer_v3_mp", "equiformer_v3_oam"),
        "no stable checkpoint artifact or isolated runtime",
    ),
    "gnome": "model weights were not released",
}


@dataclass(frozen=True)
class RelaxationSettings:
    """ASE geometry-relaxation settings shared by all calculator-backed models."""

    max_force: float = 0.05
    max_steps: int = 500
    ase_optimizer: str = "FIRE"
    cell_filter: str | None = "FrechetCellFilter"

    def __post_init__(self) -> None:
        """Validate and canonicalize optimizer and cell-filter names."""
        if self.max_force <= 0 or not np.isfinite(self.max_force):
            raise ValueError(
                f"max_force must be positive and finite, got {self.max_force}"
            )
        if self.max_steps < 0:
            raise ValueError(f"max_steps must be non-negative, got {self.max_steps}")
        object.__setattr__(
            self, "ase_optimizer", canonical_optimizer_name(self.ase_optimizer)
        )
        object.__setattr__(self, "cell_filter", canonical_filter_name(self.cell_filter))

    @classmethod
    def from_model(
        cls,
        model_key: str,
        *,
        max_force: float | None = None,
        max_steps: int | None = None,
        ase_optimizer: str | None = None,
        cell_filter: str | None = None,
    ) -> RelaxationSettings:
        """Load settings from a model YAML and apply explicit overrides.

        Raises:
            ValueError: If model_key is neither a registered model nor the "emt"
                debug calculator (which has no model YAML).
        """
        hyperparams = (
            {}
            if model_key == "emt"
            else Model.from_ref(model_key).metadata.get("hyperparams", {})
        )
        if not isinstance(hyperparams, dict):
            raise TypeError(f"{model_key} hyperparams must be a mapping")
        evaluation = hyperparams.get("evaluation", {})
        if not isinstance(evaluation, dict):
            raise TypeError(f"{model_key} hyperparams.evaluation must be a mapping")

        overrides = {
            "max_force": max_force,
            "max_steps": max_steps,
            "ase_optimizer": ase_optimizer,
            "cell_filter": cell_filter,
        }
        casts: dict[str, Callable[[Any], Any]] = {
            "max_force": float,
            "max_steps": int,
            "ase_optimizer": str,
        }
        # explicit flag > YAML value (incl. a YAML null cell_filter) > field default
        kwargs = {
            name: casts.get(name, lambda value: value)(
                override if override is not None else evaluation[name]
            )
            for name, override in overrides.items()
            if override is not None or name in evaluation
        }
        return cls(**kwargs)


@dataclass(frozen=True)
class RelaxationRecord:
    """One attempted WBM relaxation, including failures for complete coverage.

    The defaults describe a failed attempt; successful relaxations set them all.
    ``run_time_sec`` is the attempt's wall time and ``max_rss_gb``/``max_gpu_mem_gb``
    the segment's running memory peaks after it (RSS includes loading WBM), so a
    later-killed segment's cost still counts after resuming.
    """

    material_id: str
    structure: dict[str, Any] | None = None
    energy: float | None = None
    converged: bool = False
    n_steps: int = 0
    error: str | None = None
    run_time_sec: float = 0.0
    max_rss_gb: float | None = None
    max_gpu_mem_gb: float | None = None

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> RelaxationRecord:
        """Hydrate a relaxation record from its JSON representation."""
        return cls(
            material_id=str(data["material_id"]),
            structure=data.get("structure"),
            energy=(
                float(energy)
                if isinstance(energy := data.get("energy"), int | float)
                else None
            ),
            converged=bool(data.get("converged", False)),
            n_steps=int(data.get("n_steps", 0)),
            error=str(error) if (error := data.get("error")) is not None else None,
            run_time_sec=float(data["run_time_sec"]),
            max_rss_gb=data.get("max_rss_gb"),
            max_gpu_mem_gb=data.get("max_gpu_mem_gb"),
        )


@dataclass(frozen=True)
class ShardHeader:
    """Immutable metadata identifying one deterministic discovery shard."""

    model_key: str
    shard_index: int
    n_shards: int
    assigned_material_ids: tuple[str, ...]
    dataset_size: int
    dataset_id_hash: str
    settings: RelaxationSettings
    schema_version: int = DISCOVERY_SHARD_SCHEMA_VERSION

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> ShardHeader:
        """Hydrate and validate a shard header from a JSON event."""
        settings_data = data.get("settings")
        if not isinstance(settings_data, dict):
            raise TypeError("Shard header settings must be a mapping")
        return cls(
            model_key=str(data["model_key"]),
            shard_index=int(data["shard_index"]),
            n_shards=int(data["n_shards"]),
            assigned_material_ids=tuple(map(str, data["assigned_material_ids"])),
            dataset_size=int(data["dataset_size"]),
            dataset_id_hash=str(data["dataset_id_hash"]),
            settings=RelaxationSettings(**settings_data),
            schema_version=int(data["schema_version"]),
        )


@dataclass
class ShardSegment:
    """One execution (OS process) of a resumable shard: its segment_start runtime
    (calculators.load_calculator_with_runtime) and the records it wrote.
    """

    runtime: dict[str, Any]
    material_ids: list[str] = field(default_factory=list)

    def cost_unit(
        self, records_by_id: Mapping[str, RelaxationRecord]
    ) -> dict[str, Any]:
        """Merge unit: runtime, summed record run time, final memory peaks."""
        records = [records_by_id[material_id] for material_id in self.material_ids]
        return {
            **self.runtime,
            "run_time_sec": sum(rec.run_time_sec for rec in records),
            "max_rss_gb": records[-1].max_rss_gb,
            "max_gpu_mem_gb": records[-1].max_gpu_mem_gb,
        }


@dataclass(frozen=True)
class LoadedShard:
    """Parsed shard header, records and execution segments."""

    header: ShardHeader
    records: tuple[RelaxationRecord, ...]
    segments: tuple[ShardSegment, ...]
    ignored_trailing_line: bool = False


@dataclass(frozen=True)
class MergedDiscoveryRun:
    """Strictly validated records and provenance from a complete shard set."""

    model_key: str
    settings: RelaxationSettings
    records: tuple[RelaxationRecord, ...]
    run_metadata: dict[str, Any]
    n_shards: int


@dataclass(frozen=True)
class DiscoveryArtifacts:
    """Final discovery prediction, geometry-optimization and run-info artifacts."""

    pred_file_path: str
    geo_opt_file_path: str
    run_info_path: str
    predictions: pd.DataFrame
    n_success: int
    n_failed: int


def load_wbm_atoms(
    model_key: str,
    *,
    dry_run: bool = False,
) -> dict[str, Atoms]:
    """Load WBM initial structures and optionally select a small smoke-test subset."""
    from matbench_discovery.data import ase_atoms_from_zip, df_wbm

    read_limit = 2048 if dry_run else None
    atoms_list = ase_atoms_from_zip(DataFiles.wbm_initial_atoms.path, limit=read_limit)
    atoms_by_id: dict[str, Atoms] = {}
    for atoms in atoms_list:
        material_id = material_id_from_atoms(atoms)
        if material_id in atoms_by_id:
            raise ValueError(f"Duplicate WBM structure {material_id!r}")
        atoms_by_id[material_id] = atoms

    expected_ids = list(map(str, df_wbm.index))
    if not dry_run:
        missing_ids = set(expected_ids) - set(atoms_by_id)
        extra_ids = set(atoms_by_id) - set(expected_ids)
        if missing_ids or extra_ids:
            raise ValueError(
                f"WBM atoms/summary coverage mismatch: {len(missing_ids)} missing, "
                f"{len(extra_ids)} extra"
            )
        return {material_id: atoms_by_id[material_id] for material_id in expected_ids}

    candidate_ids = [
        material_id for material_id in expected_ids if material_id in atoms_by_id
    ]
    if model_key == "emt":
        emt_elements = {"H", "Al", "Ni", "Cu", "Pd", "Ag", "Pt", "Au"}
        candidate_ids = [
            material_id
            for material_id in candidate_ids
            if set(atoms_by_id[material_id].get_chemical_symbols()) <= emt_elements
        ]
    selected_ids = candidate_ids[:DEFAULT_DRY_RUN_STRUCTURES]
    if not selected_ids:
        raise ValueError(f"No dry-run-compatible WBM structures found for {model_key}")
    return {material_id: atoms_by_id[material_id] for material_id in selected_ids}


def relax_atoms(
    atoms: Atoms,
    calculator: Calculator,
    settings: RelaxationSettings,
) -> RelaxationRecord:
    """Relax one structure with ASE and return a serializable success record."""
    material_id = material_id_from_atoms(atoms)
    atoms = atoms.copy()
    atoms.calc = calculator
    optimizer_steps = 0
    converged = False
    if settings.max_steps > 0:
        filter_cls = resolve_cell_filter(settings.cell_filter)
        relax_target = filter_cls(atoms) if filter_cls is not None else atoms
        optimizer_cls = resolve_optimizer(settings.ase_optimizer)
        optimizer = optimizer_cls(relax_target, logfile=None)  # ty: ignore[invalid-argument-type]
        converged = bool(
            optimizer.run(fmax=settings.max_force, steps=settings.max_steps)
        )
        optimizer_steps = optimizer.get_number_of_steps()

    energy = float(atoms.get_potential_energy())
    if not np.isfinite(energy):
        raise ValueError(f"{material_id} produced non-finite energy {energy}")
    structure = AseAtomsAdaptor.get_structure(atoms).as_dict()
    return RelaxationRecord(
        material_id=material_id,
        structure=structure,
        energy=energy,
        converged=converged,
        n_steps=optimizer_steps,
    )


def _json_default(value: object) -> object:
    """Convert NumPy values while rejecting unsupported JSON objects.

    Calculators may stuff arrays into ``atoms.info``, which
    ``AseAtomsAdaptor.get_structure`` copies into ``Structure.properties``.
    """
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()  # ty: ignore[no-matching-overload]  # shape is unknown
    raise TypeError(f"Object of type {type(value).__name__} is not JSON serializable")


def _append_shard_event(shard_path: str, event: Mapping[str, Any]) -> None:
    """Append and fsync one JSONL event."""
    payload = json.dumps(event, default=_json_default, allow_nan=False)
    with open(shard_path, mode="ab+") as file:
        file.seek(0, os.SEEK_END)
        if file.tell():
            file.seek(-1, os.SEEK_END)
            if file.read(1) != b"\n":
                file.write(b"\n")
        file.write(payload.encode() + b"\n")
        file.flush()
        os.fsync(file.fileno())


def _drop_truncated_line(shard_path: str) -> None:
    """Atomically drop a shard's truncated final line (an interrupted append)."""
    with open(shard_path, encoding="utf-8") as file:
        lines = file.readlines()
    with (
        atomic_output_path(shard_path) as tmp_path,
        open(tmp_path, mode="w", encoding="utf-8") as file,
    ):
        file.writelines(lines[:-1])
        file.flush()
        # fsync before the rename so a power loss cannot replace hours of durably
        # appended records with a partially written file
        os.fsync(file.fileno())


def read_discovery_shard(shard_path: str) -> LoadedShard:
    """Read one discovery JSONL shard and tolerate a truncated trailing event."""
    with open(shard_path, encoding="utf-8") as file:
        lines = file.readlines()
    if not lines:
        raise ValueError(f"Empty discovery shard {shard_path}")

    events: list[dict[str, Any]] = []
    ignored_trailing_line = False
    for line_index, line in enumerate(lines):
        if not line.strip():
            continue
        try:
            event = json.loads(line)
        except json.JSONDecodeError:
            if line_index != len(lines) - 1 or line.endswith("\n"):
                raise ValueError(
                    f"Malformed JSON event at line {line_index + 1} in {shard_path}"
                ) from None
            ignored_trailing_line = True
            break
        if not isinstance(event, dict):
            raise TypeError(f"Shard event at line {line_index + 1} must be a mapping")
        events.append(event)

    if not events or events[0].get("type") != "header":
        raise ValueError(f"First event in {shard_path} must be a header")
    header = ShardHeader.from_dict(events[0])
    if header.schema_version != DISCOVERY_SHARD_SCHEMA_VERSION:
        raise ValueError(
            f"Unsupported discovery shard schema {header.schema_version}, "
            f"expected {DISCOVERY_SHARD_SCHEMA_VERSION}"
        )

    records: list[RelaxationRecord] = []
    segments: list[ShardSegment] = []
    for line_number, event in enumerate(events[1:], start=2):
        event_type = event.get("type")
        where = f"{event_type} event {line_number} in {shard_path}"
        if event_type == "segment_start":
            if not isinstance(runtime := event.get("runtime"), dict):
                raise TypeError(f"{where} needs a runtime mapping")
            segments.append(ShardSegment(runtime=runtime))
        elif event_type != "record":
            raise ValueError(f"Unknown shard event type {event_type!r} in {shard_path}")
        elif not segments:
            raise ValueError(f"{where} precedes any segment_start")
        else:
            record = RelaxationRecord.from_dict(event)
            records.append(record)
            segments[-1].material_ids.append(record.material_id)

    if len(records) != len({record.material_id for record in records}):
        raise ValueError(f"Duplicate material records in {shard_path}")
    return LoadedShard(
        header=header,
        records=tuple(records),
        segments=tuple(segments),
        ignored_trailing_line=ignored_trailing_line,
    )


def _discard_interrupted_header(shard_path: str) -> None:
    """Delete an empty shard or one with an interrupted initial header append."""
    if not os.path.isfile(shard_path):
        return
    with open(shard_path, encoding="utf-8") as file:
        content_lines = [line for line in file if line.strip()]
    if not content_lines:
        os.remove(shard_path)
        return
    if len(content_lines) != 1:
        return
    try:
        json.loads(content_lines[0])
    except json.JSONDecodeError:
        if not content_lines[0].endswith("\n"):
            os.remove(shard_path)


def run_discovery_shard(
    *,
    calculator: Calculator,
    model_key: str,
    atoms_by_id: Mapping[str, Atoms],
    assigned_material_ids: Sequence[str],
    shard_path: str,
    shard_index: int,
    n_shards: int,
    settings: RelaxationSettings,
    dataset_material_ids: Sequence[str],
    runtime: Mapping[str, Any],
) -> LoadedShard:
    """Run or resume one shard, appending a durable record after each structure.

    Each execution appends a ``segment_start`` event holding ``runtime``
    (calculators.load_calculator_with_runtime), then records with their run time
    and memory peaks. A resume under other run labels than an earlier segment's
    raises (hpc.check_resume_labels).

    A lock file serializes concurrent writers (e.g. duplicate sbatch submissions or
    a Slurm requeue overlapping a lingering process), which would otherwise
    interleave duplicate records and poison the shard for every later read.
    """
    header = ShardHeader(
        model_key=model_key,
        shard_index=shard_index,
        n_shards=n_shards,
        assigned_material_ids=tuple(assigned_material_ids),
        dataset_size=len(dataset_material_ids),
        dataset_id_hash=material_id_hash(dataset_material_ids),
        settings=settings,
    )
    if shard_parent := os.path.dirname(shard_path):
        os.makedirs(shard_parent, exist_ok=True)
    with FileLock(f"{shard_path}.lock"):
        _discard_interrupted_header(shard_path)
        if os.path.isfile(shard_path):
            existing_shard = read_discovery_shard(shard_path)
            if existing_shard.header != header:
                raise ValueError(
                    f"Existing shard header in {shard_path} does not match this run"
                )
            if existing_shard.ignored_trailing_line:
                _drop_truncated_line(shard_path)
        else:
            _append_shard_event(shard_path, {"type": "header", **asdict(header)})
            existing_shard = read_discovery_shard(shard_path)

        completed_ids = {record.material_id for record in existing_shard.records}
        if unexpected_ids := completed_ids - set(assigned_material_ids):
            raise ValueError(
                f"Shard contains unassigned records: {sorted(unexpected_ids)}"
            )

        pending_ids = [
            material_id
            for material_id in assigned_material_ids
            if material_id not in completed_ids
        ]
        if not pending_ids:
            return existing_shard
        # every segment that wrote records must share the run labels (versions too,
        # so a resume after the env re-resolved fails now, not at the final merge)
        for idx, segment in enumerate(existing_shard.segments):
            if segment.material_ids:
                check_resume_labels(segment.runtime, runtime, f"{shard_path} {idx=}")

        reset_gpu_peak_memory(runtime)
        _append_shard_event(
            shard_path, {"type": "segment_start", "runtime": dict(runtime)}
        )
        for material_id in tqdm(pending_ids, desc=f"{model_key} discovery shard"):
            start_time = time.perf_counter()
            try:
                record = relax_atoms(atoms_by_id[material_id], calculator, settings)
            except Exception as exc:  # noqa: BLE001 - every failure needs a record
                record = RelaxationRecord(
                    material_id=material_id, error=f"{type(exc).__name__}: {exc}"
                )
                print(f"Failed to relax {material_id}: {record.error}")
            cost = {
                "run_time_sec": round(time.perf_counter() - start_time, 3),
                # failures may legitimately never reach the GPU
                **memory_peaks(runtime, succeeded=record.error is None),
            }
            record = replace(record, **cost)
            try:
                _append_shard_event(shard_path, {"type": "record", **asdict(record)})
            except (TypeError, ValueError) as exc:
                # a non-JSON-serializable relaxed structure must not become a poison
                # pill that crashes this run and every resume at the same material
                error = f"unserializable relaxation record: {exc}"
                record = replace(RelaxationRecord(material_id, error=error), **cost)
                print(f"Failed to serialize {material_id}: {record.error}")
                _append_shard_event(shard_path, {"type": "record", **asdict(record)})
        return read_discovery_shard(shard_path)


def merge_discovery_shards(
    shard_paths: Sequence[str],
    *,
    model_key: str,
    expected_material_ids: Sequence[str],
) -> MergedDiscoveryRun:
    """Strictly merge shards, rejecting duplicate, missing, or foreign materials."""
    if not shard_paths:
        raise ValueError("No discovery shard files supplied")
    shards = [read_discovery_shard(path) for path in sorted(shard_paths)]
    first_header = shards[0].header
    expected_ids = tuple(map(str, expected_material_ids))
    expected_id_set = set(expected_ids)
    expected_hash = material_id_hash(expected_ids)

    if len(expected_ids) != len(expected_id_set):
        raise ValueError("Expected material IDs contain duplicates")
    if first_header.model_key != model_key:
        raise ValueError(
            f"Shard model {first_header.model_key!r} does not match {model_key!r}"
        )

    shard_indices = [shard.header.shard_index for shard in shards]
    if len(shard_indices) != len(set(shard_indices)):
        raise ValueError(f"Duplicate shard indices: {shard_indices}")
    expected_indices = list(range(first_header.n_shards))
    if sorted(shard_indices) != expected_indices:
        raise ValueError(
            f"Expected shard indices {expected_indices}, got {sorted(shard_indices)}"
        )

    records_by_id: dict[str, RelaxationRecord] = {}
    for shard in shards:
        header = shard.header
        if (
            header.model_key != model_key
            or header.n_shards != first_header.n_shards
            or header.settings != first_header.settings
            or header.dataset_size != len(expected_ids)
            or header.dataset_id_hash != expected_hash
        ):
            raise ValueError(
                f"Inconsistent shard metadata for shard {header.shard_index}"
            )
        assigned_set = set(header.assigned_material_ids)
        record_ids = {record.material_id for record in shard.records}
        if record_ids != assigned_set:
            raise ValueError(
                f"Shard {header.shard_index} records do not match assignment: "
                f"missing={sorted(assigned_set - record_ids)}, "
                f"extra={sorted(record_ids - assigned_set)}"
            )
        for record in shard.records:
            if record.material_id in records_by_id:
                raise ValueError(f"Duplicate material record {record.material_id!r}")
            records_by_id[record.material_id] = record

    merged_ids = set(records_by_id)
    if merged_ids != expected_id_set:
        raise ValueError(
            "Shard assignment does not match expected WBM coverage: "
            f"missing={sorted(expected_id_set - merged_ids)}, "
            f"extra={sorted(merged_ids - expected_id_set)}"
        )

    # one merge unit per execution segment that wrote records
    cost_units = [
        segment.cost_unit(records_by_id)
        for shard in shards
        for segment in shard.segments
        if segment.material_ids
    ]
    run_metadata = merge_audit_metadata(cost_units, strict=True)
    run_metadata |= {
        "model_key": model_key,
        "n_shards": first_header.n_shards,
        "n_structures": len(expected_ids),
        "settings": asdict(first_header.settings),
    }
    return MergedDiscoveryRun(
        model_key=model_key,
        settings=first_header.settings,
        records=tuple(records_by_id[material_id] for material_id in expected_ids),
        run_metadata=run_metadata,
        n_shards=first_header.n_shards,
    )


def write_discovery_metrics(
    model: Model,
    df_wbm: pd.DataFrame,
    preds: pd.Series,
    *,
    uniq_proto_prevalence: float,
) -> dict[TestSubset, dict[str, str | float]]:
    """Score one model's WBM formation-energy predictions and write every test
    subset's metrics to its YAML (shared by the runner and scripts/evals/discovery.py).
    """
    from matbench_discovery.metrics import discovery as discovery_metrics

    metric_reference, model_preds = discovery_metrics.prepare_model_predictions(
        df_wbm, preds
    )
    metrics_by_subset = discovery_metrics.calc_discovery_metrics(
        metric_reference, model_preds, uniq_proto_prevalence=uniq_proto_prevalence
    )
    return discovery_metrics.write_all_metrics_to_yaml(
        model, metrics_by_subset, metric_reference, model_preds
    )


def _load_wbm_cse_frame(material_ids: Sequence[str] | None = None) -> pd.DataFrame:
    """Load canonical WBM entries, subsetting before expensive CSE hydration."""
    df_wbm_cse = pd.read_json(
        DataFiles.wbm_computed_structure_entries.path, lines=True
    ).set_index(Key.mat_id)
    if material_ids is not None:
        df_wbm_cse = df_wbm_cse.loc[df_wbm_cse.index.isin(material_ids)]
    df_wbm_cse[Key.computed_structure_entry] = [
        ComputedStructureEntry.from_dict(entry)
        for entry in tqdm(
            df_wbm_cse[Key.computed_structure_entry], desc="Hydrating WBM entries"
        )
    ]
    return df_wbm_cse


def write_discovery_artifacts(
    merged_run: MergedDiscoveryRun,
    *,
    pred_file_path: str,
    geo_opt_file_path: str,
    df_wbm_summary: pd.DataFrame | None = None,
    df_wbm_cse: pd.DataFrame | None = None,
    compatibility: MaterialsProject2020Compatibility | None = None,
    elemental_ref_energies: Mapping[str, float] | None = None,
) -> DiscoveryArtifacts:
    """Apply MP2020 corrections and write canonical prediction CSV and JSONL."""
    from matbench_discovery.data import df_wbm
    from matbench_discovery.energy import (
        calc_energy_from_e_refs,
        mp_elemental_ref_energies,
    )

    df_summary = df_wbm if df_wbm_summary is None else df_wbm_summary
    successful_ids = [
        record.material_id
        for record in merged_run.records
        if record.structure is not None and record.energy is not None
    ]
    df_entries = (
        _load_wbm_cse_frame(material_ids=successful_ids)
        if df_wbm_cse is None
        else df_wbm_cse.copy()
    )
    if Key.mat_id in df_entries:
        df_entries = df_entries.set_index(Key.mat_id)
    ref_energies = (
        mp_elemental_ref_energies
        if elemental_ref_energies is None
        else dict(elemental_ref_energies)
    )
    compatibility = compatibility or MaterialsProject2020Compatibility()
    output_rows: dict[str, dict[str, Any]] = {}
    geo_opt_rows: list[dict[str, Any]] = []
    for record in merged_run.records:
        row = output_rows[record.material_id] = {
            "energy": record.energy,
            "corrected_energy": np.nan,
            DISCOVERY_PRED_COL: np.nan,
            "converged": record.converged,
            "n_steps": record.n_steps,
            "error": record.error,
        }
        if record.structure is None or record.energy is None:
            continue
        if record.material_id not in df_entries.index:
            raise ValueError(
                f"No WBM computed structure entry for {record.material_id}"
            )
        entry = df_entries.loc[record.material_id, Key.computed_structure_entry]
        if df_wbm_cse is not None:
            entry = copy.deepcopy(entry)
        if isinstance(entry, dict):
            entry = ComputedStructureEntry.from_dict(entry)
        entry._energy = record.energy  # noqa: SLF001 - preserves legacy join semantics
        entry._structure = Structure.from_dict(record.structure)  # noqa: SLF001
        processed_entry = compatibility.process_entry(entry, clean=True)
        if processed_entry is None:
            row["error"] = "MP2020 compatibility rejected entry"
            continue

        if record.material_id not in df_summary.index:
            raise ValueError(f"No WBM summary row for {record.material_id}")
        formula = str(df_summary.loc[record.material_id, Key.formula])
        row["corrected_energy"] = float(processed_entry.energy)
        row[DISCOVERY_PRED_COL] = calc_energy_from_e_refs(
            formula,
            ref_energies=ref_energies,
            total_energy=float(processed_entry.energy),
        )
        geo_opt_rows.append(
            {
                DISCOVERY_ID_COL: record.material_id,
                DISCOVERY_STRUCT_COL: record.structure,
                "energy": record.energy,
                "converged": record.converged,
                "n_steps": record.n_steps,
            }
        )

    df_predictions = pd.DataFrame.from_dict(output_rows, orient="index")
    df_predictions.index.name = DISCOVERY_ID_COL
    numeric_cols = df_predictions.select_dtypes(include="number").columns
    df_predictions[numeric_cols] = df_predictions[numeric_cols].round(4)
    run_info_path = (
        f"{pred_file_path.removesuffix('.gz').removesuffix('.csv')}-run-info.json"
    )
    for output_path in (pred_file_path, geo_opt_file_path, run_info_path):
        if output_parent := os.path.dirname(output_path):
            os.makedirs(output_parent, exist_ok=True)
    # to_csv infers gzip from the .gz suffix; index.name is already DISCOVERY_ID_COL
    df_predictions[[DISCOVERY_PRED_COL]].to_csv(pred_file_path)
    pd.DataFrame(geo_opt_rows).to_json(
        geo_opt_file_path,
        orient="records",
        lines=True,
        compression="gzip" if geo_opt_file_path.endswith(".gz") else None,
    )

    # merged cost + audit provenance, so the YAML cost block is artifact-backed
    with open(run_info_path, mode="w", encoding="utf-8") as file:
        json.dump({"run_metadata": merged_run.run_metadata}, file, indent=2)
        file.write("\n")

    return DiscoveryArtifacts(
        pred_file_path=pred_file_path,
        geo_opt_file_path=geo_opt_file_path,
        run_info_path=run_info_path,
        predictions=df_predictions,
        n_success=len(geo_opt_rows),
        n_failed=len(merged_run.records) - len(geo_opt_rows),
    )


def dry_run_settings(settings: RelaxationSettings) -> RelaxationSettings:
    """Cap relaxation steps for a seconds-long end-to-end smoke test."""
    return replace(settings, max_steps=min(settings.max_steps, 2))
