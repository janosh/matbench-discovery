from glob import glob

from matbench_discovery.models import MODEL_DIRS, MODEL_METADATA


def test_model_dirs_have_metadata() -> None:
    required = {
        "authors": list,  # dict with name, affiliation, orcid?, email?
        "date_added": str,
        "matbench_discovery_version": float,
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
        for key in required:
            assert key in metadata, f"Required {key=} missing in {model_dir}"
            if isinstance(required[key], dict):
                missing_keys = {*required[key]} - {*metadata[key]}  # type: ignore[misc]
                assert (
                    not missing_keys
                ), f"Missing sub-keys {missing_keys} of {key=} in {model_dir}"
                continue

            err_msg = f"Invalid {key=}, expected {required[key]} in {model_dir}"
            assert isinstance(metadata[key], required[key]), err_msg  # type: ignore[arg-type]

        authors, date_added, mbd_version, yml_model_name, model_version, repo = (
            metadata[key] for key in list(required)[:-1]
        )
        assert model_name == yml_model_name, f"{model_name=} != {yml_model_name=}"

        # make sure all keys are valid
        for name in model_name if isinstance(model_name, list) else [model_name]:
            assert (
                3 < len(name) < 50
            ), f"Invalid {name=} not between 3 and 50 characters"
        assert (
            1 < len(model_version) < 15
        ), f"Invalid {model_version=} not between 1 and 15 characters"
        # TODO increase max allowed version when updating package
        assert (
            1 <= mbd_version <= 1
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
