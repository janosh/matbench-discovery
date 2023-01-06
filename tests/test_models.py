import os
from datetime import date
from glob import glob

import yaml

from matbench_discovery import ROOT

MODEL_DIRS = glob(f"{ROOT}/models/*/")


def test_model_dirs_have_metadata() -> None:
    required = (
        "authors",
        "date_added",
        "matbench_discovery_version",
        "model_name",
        "model_version",
        "repo",
    )
    for model_dir in MODEL_DIRS:
        md_file = f"{model_dir}metadata.yml"
        assert os.path.isfile(md_file), f"Missing metadata file: {md_file}"

        # make sure all required keys are non-empty
        with open(md_file) as yml_file:
            metadata = yaml.full_load(yml_file)

        for key in required:
            assert metadata.get(key), f"Empty {key=} in {md_file}"

        authors, date_added, mbd_version, model_name, model_version, repo = (
            metadata[key] for key in required
        )

        # make sure all keys are valid
        assert (
            3 < len(model_name) < 50
        ), f"Invalid {model_name=} not between 3 and 50 characters"
        assert (
            1 < len(model_version) < 15
        ), f"Invalid {model_version=} not between 1 and 15 characters"
        # TODO increase max version when releasing new versions
        assert (
            1 <= mbd_version <= 1
        ), f"Invalid matbench-discovery version: {mbd_version}"
        assert isinstance(date_added, date), f"Invalid {date_added=} not a string"
        assert (
            isinstance(authors, list) and 1 < len(authors) < 30
        ), "authors not list or not between 1 and 30 authors"
        assert repo.startswith(
            "https://"
        ), f"Invalid {repo=} not starting with https://"


def test_model_dirs_have_test_scripts() -> None:
    for model_dir in MODEL_DIRS:
        test_scripts = glob(f"{model_dir}*test_*.py")
        test_nbs = glob(f"{model_dir}*test_*.ipynb")
        assert len(test_scripts + test_nbs) > 0, f"Missing test file in {model_dir}"
