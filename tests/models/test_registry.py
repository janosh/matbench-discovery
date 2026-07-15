"""Model metadata schema, identity, registry, and runner guardrails."""

from __future__ import annotations

import os
import re
from functools import cache
from glob import glob
from typing import Any

import pytest
import yaml
from jsonschema import Draft7Validator

from matbench_discovery import DATA_DIR, PKG_DIR, ROOT
from matbench_discovery.calculators import CALCULATORS
from matbench_discovery.data import DATASETS, iter_file_refs, parse_artifact_filename
from matbench_discovery.discovery import ARCHIVED_DISCOVERY_MODELS
from matbench_discovery.enums import ArchitectureType, Model, Open, Targets, Task

with open(f"{ROOT}/tests/model-schema.yml", encoding="utf-8") as file:
    MODEL_SCHEMA: dict[str, Any] = yaml.safe_load(file)

MODEL_YAML_PATHS = sorted(glob(f"{ROOT}/models/[!_]*/[!_]*.yml"))
OPEN_DATASETS = {
    dataset["name"]
    for dataset in DATASETS.values()
    if isinstance(dataset, dict) and dataset.get("open")
}
TEST_TASK_NAMES = {"IP2E", "IS2E", "IS2RE", "IS2RE_SR"}


def validate_against_schema(
    instance: object, schema: dict[str, Any], label: str
) -> None:
    """Raise AssertionError when instance fails schema validation."""
    errors = sorted(
        Draft7Validator(schema).iter_errors(instance), key=lambda err: list(err.path)
    )
    assert not errors, f"Schema validation failed for {label}:\n" + "\n".join(
        f"  - {'.'.join(str(part) for part in error.path)}: {error.message}"
        for error in errors
    )


@cache
def validate_model_yaml(yaml_path: str) -> dict[str, Any]:
    """Validate a model YAML against tests/model-schema.yml, returning its metadata."""
    with open(yaml_path, encoding="utf-8") as file:
        metadata = yaml.safe_load(file)
    validate_against_schema(metadata, MODEL_SCHEMA, yaml_path)
    return metadata


@pytest.mark.parametrize("yaml_path", MODEL_YAML_PATHS)
def test_model_yaml_schema_and_identity(yaml_path: str) -> None:
    """Each model YAML passes schema checks and uses canonical identity/paths."""
    # the schema enforces model_key's kebab-case pattern, so only check the family dir
    metadata = validate_model_yaml(yaml_path)
    model_key = metadata["model_key"]
    family_dir = os.path.basename(os.path.dirname(yaml_path))
    assert re.fullmatch(r"[a-z0-9]+(?:_[a-z0-9]+)*", family_dir)
    assert os.path.normpath(yaml_path) == os.path.normpath(
        f"{ROOT}/models/{family_dir}/{model_key}.yml"
    )
    expected_artifact_dir = f"models/{family_dir}/{model_key}"
    for _key_path, artifact_path in iter_file_refs(metadata):
        assert os.path.dirname(artifact_path) == expected_artifact_dir, (
            f"Artifact {artifact_path!r} must live under {expected_artifact_dir!r} "
            f"({yaml_path})"
        )
        parse_artifact_filename(os.path.basename(artifact_path))


def test_modeling_tasks_align_with_schema() -> None:
    """modeling-tasks.yml validates and its keys match model-schema metrics."""
    tasks_path = f"{PKG_DIR}/modeling-tasks.yml"
    schema_path = f"{ROOT}/tests/modeling-tasks-schema.yml"
    with open(tasks_path, encoding="utf-8") as file:
        modeling_tasks = yaml.safe_load(file)
    with open(schema_path, encoding="utf-8") as file:
        validate_against_schema(modeling_tasks, yaml.safe_load(file), tasks_path)
    metrics_schema = MODEL_SCHEMA["properties"]["metrics"]
    assert set(modeling_tasks) - {"cps"} == set(metrics_schema["properties"])


def test_datasets_yaml_matches_python_registry() -> None:
    """Training-set keys cannot drift between data/datasets.yml and Python."""
    with open(f"{DATA_DIR}/datasets.yml", encoding="utf-8") as file:
        datasets_from_file = yaml.safe_load(file)
    assert set(datasets_from_file) == set(DATASETS)


@pytest.mark.parametrize(
    ("py_values", "schema_property"),
    [
        ({member.value for member in ArchitectureType}, "architecture_types"),
        ({member.value for member in Open}, "openness"),
        ({member.value for member in Targets}, "targets"),
        (
            {task.value for task in Task if task.name not in TEST_TASK_NAMES},
            "train_task",
        ),
        (
            {task.value for task in Task if task.name in TEST_TASK_NAMES},
            "test_task",
        ),
    ],
)
def test_enums_match_schema(py_values: set[str], schema_property: str) -> None:
    """Python taxonomy/task enums stay aligned with model-schema vocabulary."""
    property_schema = MODEL_SCHEMA["properties"][schema_property]
    schema_values = property_schema.get("enum")
    if schema_values is None:
        ref = property_schema.get("$ref")
        if ref is None:
            ref = property_schema["items"]["$ref"]
        schema_values = MODEL_SCHEMA["definitions"][ref.rsplit("/", maxsplit=1)[-1]][
            "enum"
        ]
    assert py_values == set(schema_values)


def test_model_registry_identity() -> None:
    """Model YAMLs, enum members, registries, and family dirs stay aligned."""
    yaml_paths = [
        yaml_path
        for yaml_path in MODEL_YAML_PATHS
        if validate_model_yaml(yaml_path).get("lifecycle") != "aborted"
    ]
    enum_paths = sorted(model.yaml_path for model in Model)
    assert sorted(map(os.path.normpath, yaml_paths)) == sorted(
        map(os.path.normpath, enum_paths)
    )

    for model in Model:
        metadata = model.metadata
        model_key = metadata["model_key"]
        family_dir = os.path.basename(os.path.dirname(model.yaml_path))
        assert model.name == model_key.replace("-", "_").replace(".", "_")
        assert model.key == model_key
        assert "family" not in metadata
        assert model.family == family_dir
        # yaml_path is base_dir + rel_path, so rel_path pins the full path
        assert model.rel_path == f"{family_dir}/{model_key}.yml"

    registry_keys = set(CALCULATORS) | set(ARCHIVED_DISCOVERY_MODELS)
    assert registry_keys - {"emt"} <= Model.__members__.keys()
    active_families = {model.family for model in Model.active()}
    # normpath: Windows globs may end in `\`; basename of a trailing-sep path is "".
    model_dir_families = {
        os.path.basename(os.path.normpath(path))
        for path in glob(f"{ROOT}/models/[!_]*/")
    }
    # Families intentionally without an active Model enum entry.
    inactive_only_families = {"alignn_ff"}
    assert active_families | inactive_only_families == model_dir_families
    assert active_families.isdisjoint(inactive_only_families)


@pytest.mark.parametrize(
    ("model_key", "is_valid"),
    [
        ("model", True),
        ("model-1", True),
        ("model.1", True),
        ("Uppercase-model", False),
        ("model+plus", False),
        ("model..1", False),
        ("double--hyphen", False),
        ("equiformer_v3_mp", False),
        ("foo_bar", False),
    ],
)
def test_model_key_schema_requires_canonical_keys(
    model_key: str, is_valid: bool
) -> None:
    """Model key schema accepts only canonical lowercase kebab-case."""
    pattern = MODEL_SCHEMA["properties"]["model_key"]["pattern"]
    assert bool(re.search(pattern, model_key)) is is_valid


@pytest.mark.parametrize("task", ["diatomics", "discovery", "kappa"])
def test_shared_runner_is_only_executable(task: str) -> None:
    """Shared benchmark tasks have one runner and no per-model forks."""
    assert os.path.isfile(f"{ROOT}/models/run_{task}.py")
    assert not glob(f"{ROOT}/models/**/test_*_{task}.py", recursive=True)
    if task == "discovery":
        with open(f"{ROOT}/.github/pull_request_template.md") as file:
            assert "test_<arch_name>_discovery.py" not in file.read()


def test_runnable_kappa_models_have_complete_shared_contract() -> None:
    """Every calculator-backed phonon model has settings and backend dispatch."""
    from matbench_discovery.phonons.pipeline import KappaSettings

    configured_models = {
        model.name
        for model in Model
        if model.metadata.get("hyperparams", {}).get("evaluation", {}).get("kappa")
    }
    assert configured_models == set(CALCULATORS) - {"emt"}
    for model_key in configured_models:
        assert isinstance(KappaSettings.from_model(model_key), KappaSettings)

    prediction_models = {
        model.name
        for model in Model
        if model.metrics.get("phonons", {}).get("kappa_103")
    }
    assert prediction_models - configured_models == {"matris_v050_mptrj"}


@pytest.mark.parametrize("model", Model.active())
def test_active_model_metadata_policy(model: Model) -> None:
    """Active model metadata follows ordering, openness, and length policies."""
    training_sets = model.metadata["training_sets"]
    training_set_keys = set(training_sets)
    assert training_set_keys <= set(DATASETS), f"Invalid training set: {training_sets}"
    expected_order = [key for key in DATASETS if key in training_set_keys]
    assert training_sets == expected_order
    if training_set_keys <= OPEN_DATASETS:
        assert model.metadata["openness"].endswith("OD"), (
            f"{model.label} was only trained on open datasets but is "
            f"marked as {model.metadata['openness']}. Should be marked as OD."
        )

    assert model.label == model.metadata["model_name"]
    assert 3 <= len(model.label) <= 50
    model_version = model.metadata["model_version"]
    assert model_version is None or 1 <= len(model_version) <= 40
    assert len(model.metadata["authors"]) <= 30


@pytest.mark.parametrize("model", list(Model))
def test_runnable_models_have_reproducible_runners(model: Model) -> None:
    """Runnable models expose shared calculators or retained task scripts."""
    if model.is_active and model.metrics.get("discovery", {}).get("pred_file"):
        assert model.name in CALCULATORS or model.name in ARCHIVED_DISCOVERY_MODELS, (
            f"{model.name} has discovery results but no shared or archived runner state"
        )
    if (
        model.name in CALCULATORS
        or model.metadata.get("targets") == "E"
        or model.metadata.get("checkpoint_url") is None
        or model.metadata.get("lifecycle") in {"aborted", "superseded"}
    ):
        return
    model_dir = os.path.dirname(model.yaml_path)
    assert glob(f"{model_dir}/test_*.py") or glob(f"{model_dir}/test_*.ipynb"), (
        f"Missing test file in {model_dir}"
    )


def test_environments_drive_calculator_dependencies() -> None:
    """Registered calculators expose exactly their YAML environment blocks."""
    for calculator_key, calc_spec in CALCULATORS.items():
        if calculator_key == "emt":
            continue
        model = Model[calculator_key]
        environment = model.metadata["environment"]
        assert calc_spec.deps == tuple(environment["dependencies"])
        assert calc_spec.python_version == environment.get("python_version")
        assert calc_spec.find_links == tuple(environment.get("find_links", []))
        assert calc_spec.extra_index_url == tuple(
            environment.get("extra_index_urls", [])
        )
        assert calc_spec.project == environment.get("project")
