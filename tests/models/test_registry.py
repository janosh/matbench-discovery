"""Model metadata schema, identity, registry, and runner guardrails."""

from __future__ import annotations

import os
import re
import shutil
from glob import glob
from typing import TYPE_CHECKING

import pytest
import yaml

import scripts.apply_pr_models_overlay as overlay
from matbench_discovery import DATA_DIR, PKG_DIR, ROOT
from matbench_discovery.calculators import CALCULATORS
from matbench_discovery.data import DATASETS
from matbench_discovery.discovery import ARCHIVED_DISCOVERY_MODELS
from matbench_discovery.enums import ArchitectureType, Model, Open, Targets, Task
from tests.models._helpers import (
    FAMILY_DIR_PATTERN,
    MODEL_KEY_PATTERN,
    assert_canonical_artifact_path,
    declared_artifact_paths,
    enum_values_from_schema,
    load_model_schema,
    model_yaml_paths,
    non_aborted_model_yaml_paths,
    validate_model_yaml,
    validate_modeling_tasks_yaml,
)

if TYPE_CHECKING:
    from pathlib import Path

OPEN_DATASETS = {
    dataset["name"]
    for dataset in DATASETS.values()
    if isinstance(dataset, dict) and dataset.get("open")
}
MODEL_DIRS = sorted(glob(f"{ROOT}/models/[!_]*/"))
TEST_TASK_NAMES = {"IP2E", "IS2E", "IS2RE", "IS2RE_SR"}


@pytest.mark.parametrize("yaml_path", model_yaml_paths())
def test_model_yaml_schema_and_identity(yaml_path: str) -> None:
    """Each model YAML passes schema checks and uses canonical identity/paths."""
    validate_model_yaml(yaml_path)
    with open(yaml_path, encoding="utf-8") as file:
        metadata = yaml.safe_load(file)
    assert isinstance(metadata, dict)
    model_key = metadata["model_key"]
    family_dir = os.path.basename(os.path.dirname(yaml_path))
    assert MODEL_KEY_PATTERN.fullmatch(model_key)
    assert FAMILY_DIR_PATTERN.fullmatch(family_dir)
    assert os.path.normpath(yaml_path) == os.path.normpath(
        f"{ROOT}/models/{family_dir}/{model_key}.yml"
    )
    for artifact_path in declared_artifact_paths(metadata):
        assert_canonical_artifact_path(
            artifact_path,
            family_dir=family_dir,
            model_key=model_key,
            yaml_path=yaml_path,
        )


def test_modeling_tasks_align_with_schema() -> None:
    """modeling-tasks.yml validates and its keys match model-schema metrics."""
    validate_modeling_tasks_yaml()
    with open(f"{PKG_DIR}/modeling-tasks.yml", encoding="utf-8") as file:
        modeling_tasks = yaml.safe_load(file)
    metrics_schema = load_model_schema()["properties"]["metrics"]
    assert set(modeling_tasks) - {"cps"} == set(metrics_schema["properties"])


def test_overlaid_model_yaml_revalidates_schema(tmp_path: Path) -> None:
    """The ingest overlay path revalidates schema before trusted code proceeds."""
    trusted_root = os.path.join(tmp_path, "trusted")
    pr_root = os.path.join(tmp_path, "pr")
    os.makedirs(trusted_root)
    source_yaml = model_yaml_paths()[0]
    relative_path = os.path.relpath(source_yaml, ROOT).replace("\\", "/")
    destination = os.path.join(trusted_root, relative_path)
    os.makedirs(os.path.dirname(destination), exist_ok=True)
    shutil.copy2(source_yaml, destination)
    invalid_yaml = os.path.join(pr_root, relative_path)
    os.makedirs(os.path.dirname(invalid_yaml), exist_ok=True)
    with open(invalid_yaml, "w", encoding="utf-8") as file:
        file.write("model_key: INVALID\n")
    original_cwd = os.getcwd()
    try:
        os.chdir(trusted_root)
        assert overlay.main(str(pr_root), [relative_path]) == 0
        with pytest.raises(AssertionError, match="Schema validation failed"):
            validate_model_yaml(relative_path)
    finally:
        os.chdir(original_cwd)


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
    assert py_values == enum_values_from_schema(load_model_schema(), schema_property)


def test_model_yaml_enum_bijection() -> None:
    """Model YAMLs and the generated enum form an exact identity/path bijection."""
    yaml_paths = non_aborted_model_yaml_paths()
    enum_paths = sorted(f"{model.yaml_path}" for model in Model)
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
        assert model.rel_path == f"{family_dir}/{model_key}.yml"
        assert os.path.normpath(model.yaml_path) == os.path.normpath(
            f"{ROOT}/models/{family_dir}/{model_key}.yml"
        )


@pytest.mark.parametrize(
    "registry_key",
    sorted(key for key in (*CALCULATORS, *ARCHIVED_DISCOVERY_MODELS) if key != "emt"),
)
def test_registry_keys_resolve_to_model_enum(registry_key: str) -> None:
    """CALCULATORS and archived discovery keys resolve to Model enum members."""
    assert registry_key in Model.__members__


@pytest.mark.parametrize(
    ("model_key", "is_valid"),
    [
        ("equiformer-v3-mp", True),
        ("mace-mp-0", True),
        ("esen-30m-mp", True),
        ("cgcnn-p", True),
        ("dpa-4.0-pro-mptrj", True),
        ("chgnet-0.3.0", True),
        ("pet-oam-xl-1.0.0", True),
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
    pattern = load_model_schema()["properties"]["model_key"]["pattern"]
    assert bool(re.search(pattern, model_key)) is is_valid


@pytest.mark.parametrize("task", ["diatomics", "discovery", "kappa"])
def test_shared_runner_is_only_executable(task: str) -> None:
    """Shared benchmark tasks have one runner and no per-model forks."""
    assert os.path.isfile(f"{ROOT}/models/run_{task}.py")
    assert not glob(f"{ROOT}/models/**/test_*_{task}.py", recursive=True)


def test_runnable_kappa_models_have_complete_shared_contract() -> None:
    """Every calculator-backed phonon model has settings and backend dispatch."""
    from matbench_discovery.phonons.pipeline import KappaSettings

    configured_models = {
        model.name
        for model in Model
        if isinstance(
            model.metadata.get("hyperparams", {}).get("evaluation", {}).get("kappa"),
            dict,
        )
    }
    assert configured_models == set(CALCULATORS) - {"emt"}
    for model_key in configured_models:
        assert isinstance(KappaSettings.from_model(model_key), KappaSettings)

    prediction_models = {
        model.name
        for model in Model
        if isinstance(model.metrics.get("phonons"), dict)
        and isinstance(model.metrics["phonons"].get("kappa_103"), dict)
    }
    assert prediction_models - configured_models == {"matris_v050_mptrj"}


@pytest.mark.parametrize("model", Model.active())
def test_active_model_training_set_policy(model: Model) -> None:
    """Training sets follow DATASETS order; open-only models must be marked OD."""
    training_sets = model.metadata["training_sets"]
    training_set_keys = set(training_sets)
    assert training_set_keys <= set(DATASETS), f"Invalid training set: {training_sets}"
    expected_order = [key for key in DATASETS if key in training_set_keys]
    assert training_sets == expected_order
    if training_set_keys <= OPEN_DATASETS:
        openness = model.metadata.get("openness")
        assert openness is not None, (
            f"{model.label} was only trained on open datasets but has no "
            "openness metadata. Should be marked as OD."
        )
        assert openness.endswith("OD"), (
            f"{model.label} was only trained on open datasets but is "
            f"marked as {openness}. Should be marked as OD."
        )

    assert model.label == model.metadata["model_name"]
    assert 3 <= len(model.label) <= 50
    model_version = model.metadata["model_version"]
    assert model_version is None or 1 <= len(model_version) <= 40
    assert 1 <= len(model.metadata["authors"]) <= 30
    repo = model.metadata["repo"]
    assert repo is None or repo.startswith(("http://", "https://"))


def test_active_model_count_matches_family_dirs() -> None:
    """Active families match model dirs; inactive-only families are allowlisted."""
    active_families = {
        os.path.basename(os.path.dirname(model.yaml_path))
        for model in Model
        if model.is_active
    }
    model_dir_families = {os.path.basename(path.rstrip("/")) for path in MODEL_DIRS}
    # Families intentionally without an active Model enum entry.
    inactive_only_families = {"alignn_ff"}
    assert active_families | inactive_only_families == model_dir_families
    assert active_families.isdisjoint(inactive_only_families)


@pytest.mark.parametrize("model", list(Model))
def test_runnable_models_have_reproducible_runners(model: Model) -> None:
    """Runnable models expose shared calculators or retained task scripts."""
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


def test_discovery_contributor_policy_forbids_runner_forks() -> None:
    """Contributor policy forbids per-model discovery runners."""
    with open(f"{ROOT}/.github/pull_request_template.md") as file:
        pr_template = file.read()
    assert "test_<arch_name>_discovery.py" not in pr_template


@pytest.mark.parametrize("model", Model.active())
def test_active_discovery_models_have_shared_or_archived_runner(model: Model) -> None:
    """Active discovery models are shared-runner-backed or explicitly archived."""
    discovery_metrics = model.metrics.get("discovery")
    if not isinstance(discovery_metrics, dict):
        return
    if not discovery_metrics.get("pred_file"):
        return
    assert model.name in CALCULATORS or model.name in ARCHIVED_DISCOVERY_MODELS, (
        f"{model.name} has discovery results but no shared or archived runner state"
    )
