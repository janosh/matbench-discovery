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


def test_model_dirs_have_metadata() -> None:
    """Test that all model directories have required metadata."""
    required = {
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
        for key, expected in required.items():
            assert key in model.metadata, f"Required {key=} missing in {model_dir}"

            if key == "training_set":
                training_sets = model.metadata[key]
                # allow either string key or dict
                assert isinstance(training_sets, list)
                assert set(training_sets) <= {*DATASETS}, (
                    f"Invalid training set: {training_sets}"
                )
                # Check if model was trained only on open datasets
                if set(training_sets) <= OPEN_DATASETS and not model.metadata[
                    "openness"
                ].endswith("OD"):
                    # if so, check that the model is marked as OD (open data)
                    raise ValueError(
                        f"{model.label} was only trained on open datasets but is "
                        f"marked as {model.metadata['openness']}. Should be marked as "
                        "OD."
                    )
                continue

            actual_val = model.metadata[key]
            if type(expected) is type:
                err_msg = f"Invalid {key=}, expected {expected} in {model_dir}"
                assert isinstance(actual_val, expected), err_msg

        authors, date_added, yml_model_name, model_version, repo = (
            model.metadata[key] for key in list(required)[:-1]
        )
        assert model.label == yml_model_name, f"{model.label=} != {yml_model_name=}"

        # make sure all keys are valid
        for name in model.label if isinstance(model.label, list) else [model.label]:
            assert 3 <= len(name) < 50, (
                f"Invalid {name=} not between 3 and 50 characters"
            )
        assert 1 <= len(model_version) < 30, (
            f"Invalid {model_version=} not between 1 and 30 characters"
        )
        assert isinstance(date_added, str), f"Invalid {date_added=} not a string"
        assert isinstance(authors, list)
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
            or model.metadata.get("status") == "aborted"
        ):
            continue
        model_dir = os.path.dirname(model.yaml_path)
        assert glob(f"{model_dir}/test_*.py") or glob(f"{model_dir}/test_*.ipynb"), (
            f"Missing test file in {model_dir}"
        )


def test_discovery_uses_only_shared_runner() -> None:
    """Runnable calculators use one runner and contributor policy forbids forks."""
    assert os.path.isfile(f"{ROOT}/models/run_discovery.py")
    assert os.path.isfile(f"{ROOT}/models/run_diatomics.py")
    assert not glob(f"{ROOT}/models/**/test_*_discovery.py", recursive=True)
    assert set(CALCULATORS).isdisjoint(ARCHIVED_DISCOVERY_MODELS)
    for model_key in set(CALCULATORS) - {"emt"}:
        assert Model.from_ref(model_key).name == model_key
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


def test_model_enum() -> None:
    """Test Model enum functionality."""
    # Skip file existence checks in CI environment
    for model in Model:
        assert os.path.isfile(model.yaml_path)
    for model in Model.active():
        assert "/models/" in model.discovery_path

    # Test model properties that don't depend on file existence
    assert Model.mace_mp_0.label == "MACE-MP-0"
    assert Model.mace_mp_0.name == Model.mace_mp_0.value == "mace_mp_0"
    assert Model.active() == tuple(model for model in Model if model.is_complete)
    assert not Model.alphanet_mptrj.is_complete
    assert not Model.dpa_3_1_mptrj.is_complete
    model_keys = {model.key for model in Model}
    for model in Model:
        if model.metadata.get("status") == "superseded":
            assert model.metadata["superseded_by"] in model_keys


@pytest.mark.parametrize(
    "input_value, expected_model",
    [
        # Exact matches
        ("mace_mp_0", Model.mace_mp_0),
        ("eqv2_s_dens_mp", Model.eqv2_s_dens_mp),
        # Dash conversion
        ("mace-mp-0", Model.mace_mp_0),
        ("eqV2-s-dens-mp", Model.eqv2_s_dens_mp),
        # Case insensitive
        ("MACE-MP-0", Model.mace_mp_0),
        ("EQV2-S-DENS-MP", Model.eqv2_s_dens_mp),
        # Mixed separators
        ("mace-mp_0", Model.mace_mp_0),
        ("mace_mp-0", Model.mace_mp_0),
    ],
)
def test_model_missing_valid_inputs(input_value: str, expected_model: Model) -> None:
    """Test that _missing_ method correctly handles valid inputs."""
    assert Model._missing_(input_value) is expected_model


@pytest.mark.parametrize(
    "input_value",
    [123, None, [], {}, "nonexistent", "mace-mp-1", "eqv2-s-dens", "", "   "],
)
def test_model_missing_invalid_inputs(
    input_value: str | int | None | list | dict,
) -> None:
    """Test that _missing_ method returns None for invalid inputs."""
    assert Model._missing_(input_value) is None


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
