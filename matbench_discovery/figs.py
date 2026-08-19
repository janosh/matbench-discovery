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
from typing import TYPE_CHECKING, Any, Final

from matbench_discovery.enums import Model

if TYPE_CHECKING:
    from collections.abc import Mapping, Sequence

PAYLOAD_SCHEMA_VERSION: Final = 2


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
    """Build separate computation identity and informational runtime audit data."""
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
        "identity": {
            "benchmark_inputs": benchmark_manifest,
            "recipe": {
                "sources": source_manifest,
                "sha256": canonical_sha256(source_manifest),
            },
            "parameters": dict(parameters),
        },
        "audit": {
            "runtime": {
                "python": sys.version.split()[0],
                "packages": {
                    package: importlib.metadata.version(package)
                    for package in sorted(set(packages))
                },
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
    """Build discovery identity from benchmark, generator, and caller recipe sources."""
    from matbench_discovery.data import MAX_E_FORM_ERROR_THRESHOLD
    from matbench_discovery.enums import DataFiles

    return build_payload_provenance(
        generator=generator,
        benchmark_inputs={"wbm_summary": DataFiles.wbm_summary.path}
        | dict(benchmark_inputs or {}),
        source_files=source_files or {},
        parameters={
            "max_formation_energy_error": MAX_E_FORM_ERROR_THRESHOLD,
            "prediction_round_decimals": prediction_round_decimals,
            "test_subset": test_subset,
            **parameters,
        },
        packages=("numpy", "pandas", "pymatviz", *packages),
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


def _validate_artifact_manifest(entries: object, *, source_files: bool = False) -> None:
    """Validate a canonical artifact or recipe-source manifest."""
    if not isinstance(entries, list):
        raise TypeError("Artifact manifests must be lists")
    if not entries:
        raise ValueError("Artifact manifests must be non-empty lists")
    expected_keys = {"role", "sha256", "size"} | ({"path"} if source_files else set())
    roles: list[str] = []
    for entry in entries:
        if not isinstance(entry, dict):
            raise TypeError("Manifest entries must be objects")
        if set(entry) != expected_keys:
            raise ValueError(f"Manifest entries must contain {sorted(expected_keys)}")
        manifest = {str(key): value for key, value in entry.items()}
        role = manifest["role"]
        if not isinstance(role, str) or re.fullmatch(r"[A-Za-z0-9_.-]+", role) is None:
            raise ValueError(f"Invalid manifest role: {role!r}")
        sha256 = manifest["sha256"]
        if not isinstance(sha256, str) or re.fullmatch(r"[0-9a-f]{64}", sha256) is None:
            raise ValueError(f"Invalid manifest SHA-256 for role {role!r}")
        if type(manifest["size"]) is not int or manifest["size"] < 0:
            raise ValueError(f"Invalid manifest byte size for role {role!r}")
        if source_files:
            path = manifest["path"]
            if (
                not isinstance(path, str)
                or not path
                or os.path.isabs(path)
                or "\\" in path
                or path != posixpath.normpath(path)
                or path in (".", "..")
                or path.startswith("../")
            ):
                raise ValueError(
                    f"Recipe source path must be repo-relative POSIX: {path!r}"
                )
        roles.append(role)
    if roles != sorted(set(roles)):
        raise ValueError("Manifest roles must be unique and sorted")


def _validate_payload(base: dict[str, Any], models: list[dict[str, Any]]) -> None:
    """Validate the strict schema-v2 payload contract."""
    if set(base) != {"schema_version", "identity", "audit", "derived"}:
        raise ValueError(
            "Payload _base must contain schema_version, identity, audit, and derived"
        )
    if base.get("schema_version") != PAYLOAD_SCHEMA_VERSION:
        raise ValueError(f"Expected payload schema_version={PAYLOAD_SCHEMA_VERSION}")
    identity = base["identity"]
    if not isinstance(identity, dict):
        raise TypeError("Payload _base.identity must be an object")
    expected_identity = {"benchmark_inputs", "recipe", "parameters"}
    if set(identity) != expected_identity:
        raise ValueError(f"Payload identity must contain {sorted(expected_identity)}")
    _validate_artifact_manifest(identity["benchmark_inputs"])
    recipe = identity["recipe"]
    if not isinstance(recipe, dict) or set(recipe) != {"sources", "sha256"}:
        raise ValueError("Payload recipe must contain sources and sha256")
    _validate_artifact_manifest(recipe["sources"], source_files=True)
    if recipe["sha256"] != canonical_sha256(recipe["sources"]):
        raise ValueError("Payload recipe SHA-256 does not match its source manifest")
    if not isinstance(identity["parameters"], dict):
        raise TypeError("Payload identity parameters must be an object")
    audit = base["audit"]
    if not isinstance(audit, dict) or set(audit) != {"runtime"}:
        raise ValueError("Payload audit must contain runtime")
    runtime = audit["runtime"]
    if not isinstance(runtime, dict) or set(runtime) != {"python", "packages"}:
        raise ValueError("Payload runtime must contain Python and package versions")
    if not isinstance(runtime["python"], str) or not runtime["python"]:
        raise ValueError("Payload runtime Python version must be a non-empty string")
    packages = runtime["packages"]
    if not isinstance(packages, dict):
        raise TypeError("Payload runtime packages must be an object")
    if not all(
        isinstance(package, str) and isinstance(version, str)
        for package, version in packages.items()
    ):
        raise TypeError("Payload runtime package names and versions must be strings")
    if not isinstance(base["derived"], dict):
        raise TypeError("Payload _base.derived must be an object")
    for model in models:
        if "key" in model:
            raise ValueError("Legacy payload field 'key' is forbidden; use model_key")
        model_key = model.get("model_key")
        if not isinstance(model_key, str):
            raise TypeError(f"Invalid canonical model_key: {model_key!r}")
        if re.fullmatch(r"[a-z0-9]+(?:[.-][a-z0-9]+)*", model_key) is None:
            raise ValueError(f"Invalid canonical model_key: {model_key!r}")
        label = model.get("label")
        if not isinstance(label, str) or not label:
            raise ValueError(f"Payload record {model_key!r} has no display label")
        _validate_artifact_manifest(model.get("input_artifacts"))
    if len(models) != len({model["model_key"] for model in models}):
        raise ValueError("Payload model_key values must be unique")


def read_jsonl_payload(path: str) -> dict[str, Any]:
    """Read and validate a schema-v2 JSONL figure payload."""
    base, models = _read_jsonl_records(path)
    _validate_payload(base, models)
    derived = base.pop("derived")
    return base | derived | {"models": models}


def write_jsonl_payload(
    path: str,
    payload: dict[str, Any],
    *,
    target_keys: set[str] | None,
) -> int:
    """Replace a full payload or splice records selected by immutable model key."""
    missing_sections = {"identity", "audit", "models"} - set(payload)
    if missing_sections:
        raise ValueError(f"Payload generation is missing {sorted(missing_sections)}")
    if "provenance" in payload:
        raise ValueError("Payload provenance must be split into identity and audit")

    fresh_models = [
        {key: value for key, value in model.items() if key not in ("color", "visible")}
        for model in payload["models"]
    ]
    fresh_base = {
        "schema_version": PAYLOAD_SCHEMA_VERSION,
        "identity": payload["identity"],
        "audit": payload["audit"],
        "derived": {
            key: value
            for key, value in payload.items()
            if key not in {"schema_version", "identity", "audit", "models"}
        },
    }
    _validate_payload(fresh_base, fresh_models)
    fresh_by_key = {model["model_key"]: model for model in fresh_models}
    if target_keys is not None:
        unexpected = set(fresh_by_key) - target_keys
        if unexpected:
            raise ValueError(
                f"{path}: targeted payload contains unselected model keys "
                f"{sorted(unexpected)}"
            )

    committed_base: dict[str, Any] | None = None
    committed_by_key: dict[str, dict[str, Any]] = {}
    if os.path.isfile(path):
        committed_base, committed_models = _read_jsonl_records(path)
        _validate_payload(committed_base, committed_models)
        committed_by_key = {model["model_key"]: model for model in committed_models}
    elif target_keys is not None:
        raise FileNotFoundError(
            f"{path} not found: targeted runs require an existing payload"
        )

    if committed_base is not None:
        identity_changed = canonical_json(committed_base["identity"]) != canonical_json(
            fresh_base["identity"]
        )
        if target_keys is not None and identity_changed:
            raise ValueError(
                f"{path}: targeted payload identity changed; use --full-roster"
            )
        if (
            target_keys is None
            and identity_changed
            and committed_by_key.keys() != fresh_by_key.keys()
        ):
            raise ValueError(
                f"{path}: payload identity and roster changed together; "
                "commit them separately"
            )
        if not identity_changed:
            if canonical_json(committed_base["derived"]) != canonical_json(
                fresh_base["derived"]
            ):
                raise ValueError(f"{path}: shared data changed with unchanged identity")
            for model_key in committed_by_key.keys() & fresh_by_key.keys():
                old_model = committed_by_key[model_key]
                new_model = fresh_by_key[model_key]
                expected = old_model | {"label": new_model["label"]}
                if old_model["input_artifacts"] == new_model[
                    "input_artifacts"
                ] and canonical_json(expected) != canonical_json(new_model):
                    raise ValueError(
                        f"{path}: {model_key!r} changed data with unchanged inputs"
                    )
        if target_keys is not None:
            for model_key in target_keys:
                committed_by_key.pop(model_key, None)
            fresh_models = list((committed_by_key | fresh_by_key).values())
    fresh_models.sort(key=lambda model: model["model_key"])
    records = [{"_base": fresh_base}, *fresh_models]
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
        f"({n_bytes:,} bytes, {len(fresh_models)} models)"
    )
    return n_bytes


def write_site_payload(
    name: str, payload: dict[str, Any], *, directory: str | None = None
) -> int:
    """Write one strict multi-model payload to ``<directory>/<name>.jsonl``."""
    from matbench_discovery import SITE_FIG_DATA
    from matbench_discovery.cli import cli_args, is_full_model_run

    target_keys = None
    if not is_full_model_run():
        target_keys = {model.key for model in cli_args.models}
    return write_jsonl_payload(
        f"{directory or SITE_FIG_DATA}/{name}.jsonl",
        payload,
        target_keys=target_keys,
    )
