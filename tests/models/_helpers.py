"""Helpers for model metadata schema, identity, and artifact-path checks."""

from __future__ import annotations

import functools
import os
import re
from glob import glob
from typing import Any

import yaml
from jsonschema import Draft7Validator

from matbench_discovery import PKG_DIR, ROOT
from matbench_discovery.data import (
    FILE_REF_KEYS,
    file_ref_name,
    parse_artifact_filename,
)

MODEL_KEY_PATTERN = re.compile(r"[a-z0-9]+(?:[.-][a-z0-9]+)*")
FAMILY_DIR_PATTERN = re.compile(r"[a-z0-9]+(?:_[a-z0-9]+)*")


@functools.cache
def load_model_schema() -> dict[str, Any]:
    """Load the model metadata JSON Schema."""
    with open(f"{ROOT}/tests/model-schema.yml", encoding="utf-8") as file:
        return yaml.safe_load(file)


@functools.cache
def load_modeling_tasks_schema() -> dict[str, Any]:
    """Load the modeling-tasks JSON Schema."""
    with open(f"{ROOT}/tests/modeling-tasks-schema.yml", encoding="utf-8") as file:
        return yaml.safe_load(file)


def model_yaml_paths() -> list[str]:
    """Return every top-level model metadata YAML path."""
    return sorted(glob(f"{ROOT}/models/[!_]*/[!_]*.yml"))


def non_aborted_model_yaml_paths() -> list[str]:
    """Return model YAML paths included in the generated Model enum."""
    yaml_paths: list[str] = []
    for yaml_path in model_yaml_paths():
        with open(yaml_path, encoding="utf-8") as file:
            metadata = yaml.safe_load(file)
        if metadata.get("lifecycle") == "aborted":
            continue
        yaml_paths.append(yaml_path)
    return yaml_paths


def validate_against_schema(
    instance: object, schema: dict[str, Any], label: str
) -> None:
    """Raise AssertionError when instance fails schema validation."""
    validator = Draft7Validator(schema)
    errors = sorted(validator.iter_errors(instance), key=lambda err: list(err.path))
    if not errors:
        return
    messages = "\n".join(
        f"  - {label} {'.'.join(str(part) for part in error.path)}: {error.message}"
        for error in errors
    )
    raise AssertionError(f"Schema validation failed for {label}:\n{messages}")


def validate_model_yaml(yaml_path: str) -> None:
    """Validate one model YAML file against tests/model-schema.yml."""
    with open(yaml_path, encoding="utf-8") as file:
        metadata = yaml.safe_load(file)
    validate_against_schema(metadata, load_model_schema(), yaml_path)


def validate_modeling_tasks_yaml() -> None:
    """Validate matbench_discovery/modeling-tasks.yml against its schema."""
    tasks_path = f"{PKG_DIR}/modeling-tasks.yml"
    with open(tasks_path, encoding="utf-8") as file:
        tasks = yaml.safe_load(file)
    validate_against_schema(tasks, load_modeling_tasks_schema(), tasks_path)


def enum_values_from_schema(schema: dict[str, Any], property_name: str) -> set[str]:
    """Return enum values for a top-level property, `$ref`, or array-of-`$ref`."""
    property_schema = schema["properties"][property_name]
    if "enum" in property_schema:
        return set(property_schema["enum"])
    ref = property_schema.get("$ref")
    if ref is None and property_schema.get("type") == "array":
        ref = property_schema["items"]["$ref"]
    ref_name = ref.rsplit("/", maxsplit=1)[-1]
    return set(schema["definitions"][ref_name]["enum"])


def declared_artifact_paths(value: object) -> list[str]:
    """Collect declared local artifact paths from nested model metadata."""
    if isinstance(value, list):
        return [path for item in value for path in declared_artifact_paths(item)]
    if not isinstance(value, dict):
        return []
    paths = [
        name
        for key, nested_value in value.items()
        if key in FILE_REF_KEYS and (name := file_ref_name(nested_value)) is not None
    ]
    return paths + [
        path
        for nested_value in value.values()
        for path in declared_artifact_paths(nested_value)
    ]


def assert_canonical_artifact_path(
    artifact_path: str, *, family_dir: str, model_key: str, yaml_path: str
) -> None:
    """Assert one artifact path uses canonical grammar and model ownership."""
    expected_artifact_dir = f"models/{family_dir}/{model_key}"
    assert os.path.dirname(artifact_path) == expected_artifact_dir, (
        f"Artifact {artifact_path!r} must live under {expected_artifact_dir!r} "
        f"({yaml_path})"
    )
    parse_artifact_filename(os.path.basename(artifact_path))
