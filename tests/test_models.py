import os
from glob import glob

from matbench_discovery import ROOT, __version__
from matbench_discovery.models import MODEL_DIRS, MODEL_METADATA
from matbench_discovery.preds import Model


def parse_version(v: str) -> tuple[int, ...]:
    return tuple(map(int, v.split(".")))


def test_model_dirs_have_metadata() -> None:
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
            actual_val = metadata[key]
            if key == "training_set":
                # allow either string key or dict
                assert isinstance(actual_val, dict | str | list)
            if (isinstance(expected, dict) and key != "training_set") or (
                key == "training_set" and isinstance(actual_val, dict)
            ):
                missing_keys = {*expected} - {*actual_val}  # type: ignore[misc]
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
            assert (
                3 <= len(name) < 50
            ), f"Invalid {name=} not between 3 and 50 characters"
        assert (
            1 < len(model_version) < 30
        ), f"Invalid {model_version=} not between 1 and 30 characters"
        # TODO increase max allowed version when updating package
        assert (
            parse_version("1.0.0")
            <= parse_version(mbd_version)
            <= parse_version(__version__)
        ), f"Invalid matbench-discovery version: {mbd_version}"
        assert isinstance(date_added, str), f"Invalid {date_added=} not a string"
        assert isinstance(authors, list)
        assert 1 < len(authors) < 30, f"{len(authors)=} not between 1 and 30"
        assert repo.startswith(
            "https://"
        ), f"Invalid {repo=} not starting with https://"


def test_model_dirs_have_test_scripts() -> None:
    for model_dir in MODEL_DIRS:
        test_scripts = glob(f"{model_dir}*test_*.py")
        test_nbs = glob(f"{model_dir}*test_*.ipynb")
        assert len(test_scripts + test_nbs) > 0, f"Missing test file in {model_dir}"


def test_model_enum() -> None:
    for model_key in Model:
        model_yaml_path = f"{ROOT}/models/{model_key.url}"
        assert os.path.isfile(model_key.path)
        assert os.path.isfile(model_yaml_path) or model_key.url is None
