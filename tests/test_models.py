import os
from glob import glob

import pytest
import yaml

from matbench_discovery import DATA_DIR, ROOT
from matbench_discovery.enums import Model
from matbench_discovery.models import model_is_compliant, validate_model_metadata

with open(f"{DATA_DIR}/datasets.yml", encoding="utf-8") as file:
    DATASETS = yaml.safe_load(file)

OPEN_DATASETS = {dataset["name"] for dataset in DATASETS.values() if dataset["open"]}

# Get model directories for testing
MODEL_DIRS = sorted(glob(f"{ROOT}/models/[!_]*/"))


def parse_version(version: str) -> tuple[int, ...]:
    """Parse version string into tuple of integers."""
    return tuple(map(int, version.split(".")))


def test_model_dirs_have_metadata() -> None:
    """Test that all model directories have required metadata."""
    required = {
        "authors": list,  # dict with name, affiliation, orcid?, email?
        "date_added": str,
        "model_name": str,
        "model_version": str,
        "repo": str,
        "training_set": {
            "title": str,  # name of training set
            "url": str,  # url to e.g. figshare
            "n_structures": int,  # number of structures
        },
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
                openness = model.metadata["openness"].endswith("OD")
                if set(training_sets) <= OPEN_DATASETS and not openness:
                    # if so, check that the model is marked as OD (open data)
                    raise ValueError(
                        f"{model.label} was only trained on open datasets but is "
                        f"marked as {model.metadata['openness']}. Should be marked as "
                        "OD."
                    )

            actual_val = model.metadata[key]
            if isinstance(expected, dict) and key != "training_set":
                missing_keys = {*expected} - {*actual_val}
                assert not missing_keys, f"{missing_keys=} under {key=} in {model_dir}"
                continue

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
        assert 1 < len(model_version) < 30, (
            f"Invalid {model_version=} not between 1 and 30 characters"
        )
        assert isinstance(date_added, str), f"Invalid {date_added=} not a string"
        assert isinstance(authors, list)
        assert 1 < len(authors) < 30, f"{len(authors)=} not between 1 and 30"
        assert repo == "missing" or repo.startswith("https://"), (
            f"Invalid {repo=} not starting with https://"
        )


def test_model_dirs_have_test_scripts() -> None:
    for model_dir in MODEL_DIRS:
        test_scripts = glob(f"{model_dir}*test_*.py")
        test_nbs = glob(f"{model_dir}*test_*.ipynb")
        assert len(test_scripts + test_nbs) > 0, f"Missing test file in {model_dir}"


def test_model_enum() -> None:
    """Test Model enum functionality."""
    # Skip file existence checks in CI environment
    for model in Model:
        assert os.path.isfile(model.yaml_path)
        assert "/models/" in model.discovery_path

    # Test model properties that don't depend on file existence
    assert Model.mace_mp_0.label == "MACE-MP-0"
    assert Model.mace_mp_0.name == Model.mace_mp_0.value == "mace_mp_0"


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
    "model, is_compliant",
    [
        (Model.megnet, True),
        (Model.eqv2_m_omat_salex_mp, False),
        (Model.eqv2_s_dens_mp, True),
        (Model.orb_v2, False),
        (Model.wrenformer, True),
        (Model.voronoi_rf, True),
        (Model.gnome, False),
        (Model.mattersim_v1_5m, False),
    ],
)
def test_model_is_compliant(model: Model, is_compliant: bool) -> None:
    """Test model compliance checking."""
    assert model.is_compliant is is_compliant
    # Also test the function directly for consistency
    assert model_is_compliant(model.metadata) is is_compliant


def test_model_is_compliant_invalid_training_set() -> None:
    """Test that model_is_compliant raises TypeError for non-list training_set."""
    metadata = {"openness": "OSOD", "training_set": "MPtrj", "model_name": "test_model"}
    with pytest.raises(TypeError, match="expected list of training sets"):
        model_is_compliant(metadata)


@pytest.mark.parametrize(
    "metadata, error_match",
    [
        ({"status": "incomplete"}, "has status != 'complete'"),
        ({"model_type": "InvalidType"}, "is not a valid ModelType"),
    ],
)
def test_validate_model_metadata_errors(
    metadata: dict[str, str], error_match: str
) -> None:
    """Test validate_model_metadata raises ValueError for invalid metadata."""
    with pytest.raises(ValueError, match=error_match):
        validate_model_metadata(metadata, "test_model.yml")
