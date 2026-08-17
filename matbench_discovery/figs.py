"""Export deterministic analysis payloads for the Svelte site.

Payloads are *data-only*: series arrays plus data-derived stats (MAE, AUC, F1, ...).
All presentation (axes, ref lines, legends, per-model colors, render order, default
visibility) lives inline in the Svelte pages that import these files.
site/src/figs/payloads.d.ts documents the expected payload shapes.

Static payloads use deterministic gzip. Multi-model payloads use canonical JSONL with
one shared ``_base`` record followed by one line per immutable ``model_key``.
"""

from __future__ import annotations

import contextlib
import gzip
import hashlib
import importlib.metadata
import json
import os
import posixpath
import re
import sys
import zlib
from enum import StrEnum
from typing import TYPE_CHECKING, Any, Final, TypeGuard

from matbench_discovery.enums import Model

if TYPE_CHECKING:
    from collections.abc import Collection, Mapping, Sequence

PAYLOAD_SCHEMA_VERSION: Final = 2


class PayloadMode(StrEnum):
    """Allowed update modes for multi-model JSONL payloads."""

    targeted = "targeted"
    full_roster = "full-roster"
    migrate_provenance = "migrate-provenance"
    migrate_model_key = "migrate-model-key"


def canonical_json(data: object) -> str:
    """Serialize JSON data canonically for hashing and exact comparisons."""
    return json.dumps(data, allow_nan=False, separators=(",", ":"), sort_keys=True)


def canonical_sha256(data: object) -> str:
    """Return SHA-256 of the UTF-8 canonical JSON representation of ``data``."""
    return hashlib.sha256(canonical_json(data).encode()).hexdigest()


def artifact_manifest(role: str, path: str) -> dict[str, str | int]:
    """Build a path-free content manifest for one benchmark or model artifact."""
    if not role or "/" in role or "\\" in role:
        raise ValueError(f"Artifact role must be a non-empty logical name: {role!r}")
    path = os.path.abspath(path)
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Payload input file not found: {path!r}")
    with open(path, "rb") as file:
        sha256 = hashlib.file_digest(file, "sha256").hexdigest()
        size = os.fstat(file.fileno()).st_size
    return {"role": role, "sha256": sha256, "size": size}


def model_payload_identity(
    model: Model, artifact_role: str, artifact_path: str
) -> dict[str, Any]:
    """Return mandatory stable identity and input provenance for one model record."""
    return {
        "model_key": model.key,
        "label": model.label,
        "input_artifacts": [artifact_manifest(artifact_role, artifact_path)],
    }


def discovery_model_identity(model: Model) -> dict[str, Any]:
    """Return payload identity for a model's discovery predictions."""
    return model_payload_identity(model, "discovery_predictions", model.discovery_path)


def _source_manifest(role: str, path: str) -> dict[str, str | int]:
    """Build a content manifest for one repo-owned recipe source file."""
    from matbench_discovery import ROOT, repo_relative_path

    relative_path = repo_relative_path(path)
    return artifact_manifest(role, f"{ROOT}/{relative_path}") | {"path": relative_path}


def build_payload_provenance(
    *,
    generator: str,
    benchmark_inputs: Mapping[str, str],
    source_files: Mapping[str, str],
    parameters: Mapping[str, Any],
    packages: Sequence[str],
) -> dict[str, Any]:
    """Describe exact shared inputs, recipe bytes, parameters, and runtime."""
    if not benchmark_inputs:
        raise ValueError("Payload provenance requires at least one benchmark input")
    if "generator" in source_files:
        raise ValueError("Recipe source role 'generator' is reserved")
    benchmark_manifest = sorted(
        (artifact_manifest(role, path) for role, path in benchmark_inputs.items()),
        key=lambda entry: str(entry["role"]),
    )
    source_manifest = sorted(
        (
            _source_manifest(role, path)
            for role, path in {"generator": generator, **source_files}.items()
        ),
        key=lambda entry: str(entry["role"]),
    )
    return {
        "benchmark_inputs": benchmark_manifest,
        "recipe": {
            "sources": source_manifest,
            "sha256": canonical_sha256(source_manifest),
        },
        "parameters": dict(parameters),
        "runtime": {
            "python": sys.version.split()[0],
            "packages": {
                package: importlib.metadata.version(package)
                for package in sorted(set(packages))
            },
        },
    }


def build_discovery_payload_provenance(
    *,
    generator: str,
    test_subset: str,
    parameters: Mapping[str, Any],
    benchmark_inputs: Mapping[str, str] | None = None,
    source_files: Mapping[str, str] | None = None,
    packages: Sequence[str] = (),
    prediction_round_decimals: int | None = 3,
) -> dict[str, Any]:
    """Build provenance with the shared discovery benchmark and loader inputs."""
    from matbench_discovery import ROOT
    from matbench_discovery.data import MAX_E_FORM_ERROR_THRESHOLD
    from matbench_discovery.enums import DataFiles

    return build_payload_provenance(
        generator=generator,
        benchmark_inputs={"wbm_summary": DataFiles.wbm_summary.path}
        | dict(benchmark_inputs or {}),
        source_files={"prediction_loader": f"{ROOT}/matbench_discovery/data.py"}
        | dict(source_files or {}),
        parameters={
            "max_formation_energy_error": MAX_E_FORM_ERROR_THRESHOLD,
            "prediction_round_decimals": prediction_round_decimals,
            "test_subset": test_subset,
            **parameters,
        },
        packages=("numpy", "pandas", "pymatviz", *packages),
    )


def computation_fingerprint(
    base: Mapping[str, Any], model_record: Mapping[str, Any]
) -> str:
    """Hash the complete shared and per-model computation identity."""
    provenance = base["provenance"]
    return canonical_sha256(
        {
            "benchmark_inputs": provenance["benchmark_inputs"],
            "recipe_sha256": provenance["recipe"]["sha256"],
            "parameters": provenance["parameters"],
            "runtime": provenance["runtime"],
            "input_artifacts": model_record["input_artifacts"],
        }
    )


# === IO ===
def write_json_gz(path: str, data: dict[str, Any]) -> int:
    """Write deterministic gzipped JSON; return compressed byte size.

    Skips the write when the existing file already decompresses to the same JSON:
    different zlib builds encode identical input to different (valid) gzip streams,
    so unconditional rewrites would churn committed payloads across CI runners.
    """
    payload = json.dumps(data, allow_nan=False, separators=(",", ":")).encode()
    # a missing/corrupt existing file falls through to the (re)write below. Per the
    # gzip docs, invalid files raise OSError (BadGzipFile), EOFError (truncation)
    # or zlib.error (corrupt deflate stream)
    with contextlib.suppress(OSError, EOFError, zlib.error), gzip.open(path) as file:
        if file.read() == payload:
            return os.path.getsize(path)
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)  # "." for bare filenames
    compressed = gzip.compress(payload, compresslevel=9, mtime=0)
    with open(path, "wb") as file:
        return file.write(compressed)


def _read_jsonl_records(path: str) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Read a schema-v2 base record followed by model records."""
    with open(path, encoding="utf-8") as file:
        records = [json.loads(line) for line in file if line.strip()]
    if not records or not isinstance(records[0], dict) or set(records[0]) != {"_base"}:
        raise ValueError(f"{path}: _base must be the first JSONL record")
    base, models = records[0]["_base"], records[1:]
    if not isinstance(base, dict) or not all(
        isinstance(model, dict) for model in models
    ):
        raise TypeError(f"{path}: JSONL records must be objects")
    model_keys = [model.get("model_key") for model in models]
    if model_keys != sorted(model_keys, key=str):
        raise ValueError(f"{path}: model records are not sorted by model_key")
    return base, models


def _validate_artifact_manifest(
    entries: object, *, source_files: bool = False
) -> list[dict[str, Any]]:
    """Validate and return a canonical artifact or recipe-source manifest."""
    if not isinstance(entries, list) or not entries:
        raise ValueError("Artifact manifests must be non-empty lists")
    expected_keys = {"role", "sha256", "size"} | ({"path"} if source_files else set())
    manifest_entries = [
        entry for entry in entries if _is_manifest_entry(entry, expected_keys)
    ]
    if len(manifest_entries) != len(entries):
        raise ValueError(f"Manifest entries must contain {sorted(expected_keys)}")
    roles = [entry["role"] for entry in manifest_entries]
    if roles != sorted(set(roles)):
        raise ValueError("Manifest roles must be unique and sorted")
    if any(
        not isinstance(role, str) or re.fullmatch(r"[A-Za-z0-9_.-]+", role) is None
        for role in roles
    ):
        raise ValueError(f"Invalid manifest roles: {roles}")
    if any(
        not isinstance(entry["sha256"], str)
        or re.fullmatch(r"[0-9a-f]{64}", entry["sha256"]) is None
        or type(entry["size"]) is not int
        or entry["size"] < 0
        for entry in manifest_entries
    ):
        raise ValueError("Manifest SHA-256 or byte size is invalid")
    if source_files and any(
        not isinstance(path := entry["path"], str)
        or not path
        or os.path.isabs(path)
        or "\\" in path
        or path != posixpath.normpath(path)
        or path in (".", "..")
        or path.startswith("../")
        for entry in manifest_entries
    ):
        raise ValueError("Recipe source paths must be repo-relative POSIX")
    return manifest_entries


def _is_manifest_entry(
    entry: object, expected_keys: set[str]
) -> TypeGuard[dict[str, Any]]:
    """Return whether ``entry`` is a manifest object with exactly expected keys."""
    return isinstance(entry, dict) and set(entry) == expected_keys


def _validate_payload(base: dict[str, Any], models: list[dict[str, Any]]) -> None:
    """Validate the strict schema-v2 payload contract."""
    if base.get("schema_version") != PAYLOAD_SCHEMA_VERSION:
        raise ValueError(f"Expected payload schema_version={PAYLOAD_SCHEMA_VERSION}")
    provenance = base.get("provenance")
    if not isinstance(provenance, dict):
        raise TypeError("Payload _base.provenance must be an object")
    expected_fields = {"benchmark_inputs", "recipe", "parameters", "runtime"}
    if set(provenance) - {"source_commit"} != expected_fields:
        raise ValueError(f"Payload provenance must contain {sorted(expected_fields)}")
    _validate_artifact_manifest(provenance["benchmark_inputs"])
    recipe = provenance["recipe"]
    if not isinstance(recipe, dict) or set(recipe) != {"sources", "sha256"}:
        raise ValueError("Payload recipe must contain sources and sha256")
    _validate_artifact_manifest(recipe["sources"], source_files=True)
    if recipe["sha256"] != canonical_sha256(recipe["sources"]):
        raise ValueError("Payload recipe SHA-256 does not match its source manifest")
    if not isinstance(provenance["parameters"], dict):
        raise TypeError("Payload provenance parameters must be an object")
    runtime = provenance["runtime"]
    if (
        not isinstance(runtime, dict)
        or set(runtime) != {"python", "packages"}
        or not isinstance(runtime["python"], str)
        or not isinstance(runtime["packages"], dict)
        or not all(
            isinstance(package, str) and isinstance(version, str)
            for package, version in runtime["packages"].items()
        )
    ):
        raise ValueError("Payload runtime must contain Python and package versions")
    source_commit = provenance.get("source_commit")
    if source_commit is not None and (
        not isinstance(source_commit, str)
        or re.fullmatch(r"[0-9a-f]{40}", source_commit) is None
    ):
        raise ValueError(f"Invalid informational source_commit: {source_commit!r}")
    for model in models:
        if "key" in model:
            raise ValueError("Legacy payload field 'key' is forbidden; use model_key")
        model_key = model.get("model_key")
        if (
            not isinstance(model_key, str)
            or re.fullmatch(r"[a-z0-9]+(?:[.-][a-z0-9]+)*", model_key) is None
        ):
            raise ValueError(f"Invalid canonical model_key: {model_key!r}")
        if not isinstance(model.get("label"), str) or not model["label"]:
            raise ValueError(f"Payload record {model_key!r} has no display label")
        _validate_artifact_manifest(model.get("input_artifacts"))
    if len(models) != len({model["model_key"] for model in models}):
        raise ValueError("Payload model_key values must be unique")


def read_jsonl_payload(path: str) -> dict[str, Any]:
    """Read and validate a schema-v2 JSONL figure payload."""
    base, models = _read_jsonl_records(path)
    _validate_payload(base, models)
    return {**base, "models": models}


def _base_without_source_commit(base: Mapping[str, Any]) -> dict[str, Any]:
    """Copy base data while excluding informational Git metadata."""
    provenance = dict(base["provenance"])
    provenance.pop("source_commit", None)
    return {**base, "provenance": provenance}


def _records_equal_except_label(
    left: Mapping[str, Any], right: Mapping[str, Any]
) -> bool:
    """Return whether two records have identical canonical non-label content."""
    return canonical_json(
        {key: value for key, value in left.items() if key != "label"}
    ) == canonical_json({key: value for key, value in right.items() if key != "label"})


def _assert_unchanged_record(
    path: str,
    old_base: Mapping[str, Any],
    new_base: Mapping[str, Any],
    old: Mapping[str, Any],
    new: Mapping[str, Any],
) -> None:
    """Require byte-stable derived data when the computation fingerprint is fixed."""
    if computation_fingerprint(old_base, old) == computation_fingerprint(
        new_base, new
    ) and not _records_equal_except_label(old, new):
        model_key = new["model_key"]
        raise ValueError(
            f"{path}: {model_key!r} changed derived data with an unchanged complete "
            "fingerprint. Only label may change."
        )


def write_jsonl_payload(
    path: str,
    payload: dict[str, Any],
    *,
    mode: PayloadMode,
    key_migration: tuple[str, str] | None = None,
    target_keys: Collection[str] | None = None,
) -> int:
    """Write a strict, deterministic, key-addressed multi-model JSONL payload."""
    if not isinstance(mode, PayloadMode):
        raise TypeError(f"mode must be a PayloadMode, got {mode!r}")
    if mode != PayloadMode.migrate_model_key and key_migration is not None:
        raise ValueError("key_migration is only valid in migrate-model-key mode")
    if mode == PayloadMode.migrate_model_key and key_migration is None:
        raise ValueError("migrate-model-key mode requires an OLD=NEW mapping")
    if "provenance" not in payload:
        raise ValueError("Payload generation must provide shared provenance")

    fresh_models = [
        {key: value for key, value in model.items() if key not in ("color", "visible")}
        for model in payload["models"]
    ]
    fresh_base = {key: value for key, value in payload.items() if key != "models"} | {
        "schema_version": PAYLOAD_SCHEMA_VERSION
    }
    _validate_payload(fresh_base, fresh_models)
    fresh_by_key = {model["model_key"]: model for model in fresh_models}
    if mode == PayloadMode.targeted:
        if target_keys is None:
            raise ValueError("targeted payload writes require selected model keys")
        unexpected = set(fresh_by_key) - set(target_keys)
        if unexpected:
            raise ValueError(
                f"{path}: targeted payload contains unselected model keys {unexpected}"
            )

    committed_base: dict[str, Any] | None = None
    committed_models: list[dict[str, Any]] = []
    if os.path.isfile(path):
        committed_base, committed_models = _read_jsonl_records(path)
    if (
        mode in (PayloadMode.targeted, PayloadMode.migrate_model_key)
        and committed_base is None
    ):
        raise FileNotFoundError(
            f"{path} not found: {mode.value} runs require an existing payload"
        )

    if committed_base is None:
        output_base = fresh_base
        output_models = fresh_models
    else:
        _validate_payload(committed_base, committed_models)
        committed_by_key = {model["model_key"]: model for model in committed_models}
        committed_keys = set(committed_by_key)
        fresh_keys = set(fresh_by_key)
        old_base = _base_without_source_commit(committed_base)
        new_base = _base_without_source_commit(fresh_base)
        base_changed = canonical_json(old_base) != canonical_json(new_base)
        if mode == PayloadMode.migrate_provenance:
            if committed_keys != fresh_keys:
                raise ValueError(
                    f"{path}: provenance migration cannot change the model roster; "
                    "run lifecycle changes separately with --full-roster"
                )
            for model_key, committed_model in committed_by_key.items():
                _assert_unchanged_record(
                    path,
                    committed_base,
                    fresh_base,
                    committed_model,
                    fresh_by_key[model_key],
                )
            if (
                canonical_json(old_base["provenance"])
                == canonical_json(new_base["provenance"])
                and base_changed
            ):
                raise ValueError(
                    f"{path}: shared derived data changed with unchanged provenance"
                )
            output_base = fresh_base
            output_models = fresh_models
        elif base_changed:
            raise ValueError(
                f"{path}: shared payload provenance or derived data changed. "
                "Regenerate explicitly with --migrate-provenance."
            )
        else:
            output_base = committed_base
            if mode != PayloadMode.migrate_model_key:
                for model_key in committed_keys & fresh_keys:
                    _assert_unchanged_record(
                        path,
                        committed_base,
                        fresh_base,
                        committed_by_key[model_key],
                        fresh_by_key[model_key],
                    )
            if mode == PayloadMode.targeted:
                output_models = list((committed_by_key | fresh_by_key).values())
            elif mode == PayloadMode.migrate_model_key:
                if key_migration is None:
                    raise ValueError(
                        "migrate-model-key mode requires an OLD=NEW mapping"
                    )
                old_key, new_key = key_migration
                new_model = Model.from_ref(new_key)
                if new_model.key != new_key or old_key not in new_model.key_aliases:
                    raise ValueError(
                        f"Key migration {old_key!r} -> {new_key!r} requires "
                        "model_key_aliases to contain the old key"
                    )
                if old_key in committed_by_key and new_key in committed_by_key:
                    raise ValueError(
                        f"Invalid key migration state for {old_key!r} -> {new_key!r}"
                    )
                expected = dict(committed_by_key)
                if old_key in expected:
                    expected[new_key] = dict(expected.pop(old_key), model_key=new_key)
                if fresh_keys != set(expected):
                    raise ValueError(f"{path}: key migration changed the model roster")
                for model_key, old in expected.items():
                    if not _records_equal_except_label(old, fresh_by_key[model_key]):
                        kind = "migrated" if model_key == new_key else "peer"
                        raise ValueError(f"{path}: key migration changed {kind} record")
                output_models = fresh_models
            else:
                for model_key in fresh_keys - committed_keys:
                    try:
                        aliases = set(Model.from_ref(model_key).key_aliases)
                    except ValueError:
                        aliases = set()
                    for old_key in (committed_keys - fresh_keys) & aliases:
                        raise ValueError(
                            f"{path}: model key replacement {old_key!r} -> "
                            f"{model_key!r} requires --migrate-model-key"
                        )
                output_models = fresh_models

    output_models.sort(key=lambda model: model["model_key"])
    records = [{"_base": output_base}, *output_models]
    body = "".join(canonical_json(record) + "\n" for record in records).encode()
    n_bytes = len(body)
    if os.path.isfile(path):
        with open(path, "rb") as file:
            if file.read() == body:
                return n_bytes
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "wb") as file:
        file.write(body)
    print(
        f"Wrote {os.path.basename(path)} "
        f"({n_bytes:,} bytes, {len(output_models)} models)"
    )
    return n_bytes


def write_site_payload(name: str, payload: dict[str, Any]) -> int:
    """Write one strict multi-model payload to ``site/src/figs/<name>.jsonl``."""
    from matbench_discovery import SITE_FIG_DATA
    from matbench_discovery.cli import cli_args, payload_mode

    mode = payload_mode()
    path = f"{SITE_FIG_DATA}/{name}.jsonl"
    return write_jsonl_payload(
        path,
        payload,
        mode=mode,
        key_migration=cli_args.migrate_model_key,
        target_keys={model.key for model in cli_args.models}
        if mode == PayloadMode.targeted
        else None,
    )
