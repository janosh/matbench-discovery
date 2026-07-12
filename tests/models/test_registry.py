"""Tests for model metadata YAML files and registry."""

import os
import re
from glob import glob

import pytest
import yaml

from matbench_discovery import ROOT
from matbench_discovery.calculators import CALCULATORS
from matbench_discovery.data import DATASETS
from matbench_discovery.discovery import ARCHIVED_DISCOVERY_MODELS
from matbench_discovery.enums import Model

OPEN_DATASETS = {
    dataset["name"]
    for dataset in DATASETS.values()
    if isinstance(dataset, dict) and dataset.get("open")
}

# Get model directories for testing
MODEL_DIRS = sorted(glob(f"{ROOT}/models/[!_]*/"))


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
        if isinstance(model.metadata.get("hyperparams", {}).get("kappa"), dict)
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


def test_model_dirs_have_metadata() -> None:
    """Test that all model directories have required metadata."""
    required_types = {
        "authors": list,  # dict with name, affiliation, orcid?, email?
        "date_added": str,
        "model_name": str,
        "model_version": str,
        "repo": str,
        "training_set": list,
    }

    # Count completed models
    completed_models = [model for model in Model if model.is_complete]
    assert len(completed_models) >= len(MODEL_DIRS) - 5, (
        "Missing metadata for some models"
    )

    for model in completed_models:
        model_dir = f"{ROOT}/models/{model.key}/"
        for key, expected_type in required_types.items():
            assert key in model.metadata, f"Required {key=} missing in {model_dir}"
            actual_value = model.metadata[key]
            err_msg = f"Invalid {key=}, expected {expected_type} in {model_dir}"
            assert isinstance(actual_value, expected_type), err_msg

            if key != "training_set":
                continue
            training_sets = actual_value
            # allow either string key or dict
            training_set_keys = set(training_sets)
            assert training_set_keys <= set(DATASETS), (
                f"Invalid training set: {training_sets}"
            )
            # Check if model was trained only on open datasets
            if training_set_keys <= OPEN_DATASETS and not model.metadata[
                "openness"
            ].endswith("OD"):
                # if so, check that the model is marked as OD (open data)
                raise ValueError(
                    f"{model.label} was only trained on open datasets but is "
                    f"marked as {model.metadata['openness']}. Should be marked as OD."
                )

        authors = model.metadata["authors"]
        metadata_model_name = model.metadata["model_name"]
        model_version = model.metadata["model_version"]
        repo = model.metadata["repo"]
        assert model.label == metadata_model_name, (
            f"{model.label=} != {metadata_model_name=}"
        )

        # make sure all keys are valid
        assert 3 <= len(model.label) < 50, (
            f"Invalid name={model.label!r} not between 3 and 50 characters"
        )
        assert 1 <= len(model_version) < 30, (
            f"Invalid {model_version=} not between 1 and 30 characters"
        )
        assert 1 < len(authors) < 30, f"{len(authors)=} not between 1 and 30"
        assert repo == "missing" or repo.startswith("https://"), (
            f"Invalid {repo=} not starting with https://"
        )


def test_model_dirs_have_reproducible_runners() -> None:
    """Require runnable models to have shared calculator or retained task coverage."""
    for model in Model:
        if (
            model.name in CALCULATORS
            or model.metadata.get("targets") == "E"
            or model.metadata.get("checkpoint_url") == "missing"
            or model.metadata.get("status") in {"aborted", "superseded"}
        ):
            continue
        model_dir = os.path.dirname(model.yaml_path)
        assert glob(f"{model_dir}/test_*.py") or glob(f"{model_dir}/test_*.ipynb"), (
            f"Missing test file in {model_dir}"
        )


def test_discovery_contributor_policy_forbids_runner_forks() -> None:
    """Contributor policy forbids per-model discovery runners."""
    with open(f"{ROOT}/.github/pull_request_template.md") as file:
        pr_template = file.read()
    assert "test_<arch_name>_discovery.py" not in pr_template


def test_active_discovery_models_have_reproducible_runner() -> None:
    """Active discovery models are shared-runner-backed or explicitly archived."""
    assert set(ARCHIVED_DISCOVERY_MODELS) <= {model.name for model in Model}
    for model in Model.active():
        discovery_metrics = model.metrics.get("discovery")
        if not isinstance(discovery_metrics, dict) or not discovery_metrics.get(
            "pred_file"
        ):
            continue
        assert model.name in CALCULATORS or model.name in ARCHIVED_DISCOVERY_MODELS, (
            f"{model.name} has discovery results but no shared or archived runner state"
        )


@pytest.mark.parametrize(
    "model_key, is_valid",
    [
        ("equiformer-v3-mp", True),  # kebab-cased
        ("mace-mp-0", True),
        ("eSEN-30m-mp", True),  # uppercase and digits allowed
        ("cgcnn+p", True),  # + allowed
        ("dpa-4.0-pro-mptrj", True),  # dots allowed
        ("equiformer_v3_mp", False),  # underscores rejected
        ("foo_bar", False),
    ],
)
def test_model_key_schema_forbids_underscores(model_key: str, is_valid: bool) -> None:
    """model_key schema pattern rejects underscores to keep model URLs param-cased."""
    with open(f"{ROOT}/tests/model-schema.yml") as file:
        pattern = yaml.safe_load(file)["properties"]["model_key"]["pattern"]
    # JSON Schema applies `pattern` as an (anchored here) regex, same as re.search
    assert bool(re.search(pattern, model_key)) is is_valid
