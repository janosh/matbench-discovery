import os
from glob import glob

import pytest
import yaml

from matbench_discovery import DATA_DIR, __version__
from matbench_discovery.enums import Model
from matbench_discovery.models import MODEL_DIRS, MODEL_METADATA, model_is_compliant

with open(f"{DATA_DIR}/datasets.yml", encoding="utf-8") as file:
    DATASETS = yaml.safe_load(file)

OPEN_DATASETS = {dataset["title"] for dataset in DATASETS.values() if dataset["open"]}


def parse_version(version: str) -> tuple[int, ...]:
    """Parse version string into tuple of integers."""
    return tuple(map(int, version.split(".")))


def test_model_dirs_have_metadata() -> None:
    """Test that all model directories have required metadata."""
    required = {
        "authors": list,  # dict with name, affiliation, orcid?, email?
        "date_added": str,
        "matbench_discovery_version": str,
        "model_name": str,
        "model_version": str,
        "repo": str,
        "training_set": {
            "title": str,  # name of training set
            "url": str,  # url to e.g. figshare
            "n_structures": int,  # number of structures
        },
    }

    assert len(MODEL_METADATA) >= len(MODEL_DIRS), "Missing metadata for some models"

    for model_name, metadata in MODEL_METADATA.items():
        model_dir = metadata["model_dir"]
        for key, expected in required.items():
            assert key in metadata, f"Required {key=} missing in {model_dir}"

            if key == "training_set":
                training_sets = metadata[key]
                # allow either string key or dict
                assert isinstance(training_sets, list)
                assert set(training_sets) <= {*DATASETS}, (
                    f"Invalid training set: {training_sets}"
                )
                # Check if model was trained only on open datasets
                openness = metadata["openness"].endswith("OD")
                if set(training_sets) <= OPEN_DATASETS and not openness:
                    # if so, check that the model is marked as OD (open data)
                    raise ValueError(
                        f"{model_name} was only trained on open datasets but is "
                        f"marked as {metadata['openness']}. Should be marked as "
                        "OD."
                    )

            actual_val = metadata[key]
            if isinstance(expected, dict) and key != "training_set":
                missing_keys = {*expected} - {*actual_val}
                assert not missing_keys, f"{missing_keys=} under {key=} in {model_dir}"
                continue

            if type(expected) is type:
                err_msg = f"Invalid {key=}, expected {expected} in {model_dir}"
                assert isinstance(actual_val, expected), err_msg

        authors, date_added, mbd_version, yml_model_name, model_version, repo = (
            metadata[key] for key in list(required)[:-1]
        )
        assert model_name == yml_model_name, f"{model_name=} != {yml_model_name=}"

        # make sure all keys are valid
        for name in model_name if isinstance(model_name, list) else [model_name]:
            assert 3 <= len(name) < 50, (
                f"Invalid {name=} not between 3 and 50 characters"
            )
        assert 1 < len(model_version) < 30, (
            f"Invalid {model_version=} not between 1 and 30 characters"
        )
        # TODO increase max allowed version when updating package
        assert (
            parse_version("1.0.0")
            <= parse_version(mbd_version)
            <= parse_version(__version__)
        ), f"Invalid matbench-discovery version: {mbd_version}"
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
    "model_key, is_compliant",
    [
        (Model.megnet, True),
        (Model.eqv2_m, False),
        (Model.eqv2_s_dens, True),
        (Model.orb_v2, False),
        (Model.wrenformer, True),
        (Model.voronoi_rf, True),
        (Model.gnome, False),
        (Model.mattersim_v1_5m, False),
    ],
)
def test_model_is_compliant(model_key: Model, is_compliant: bool) -> None:
    assert model_is_compliant(MODEL_METADATA[model_key.label]) is is_compliant
