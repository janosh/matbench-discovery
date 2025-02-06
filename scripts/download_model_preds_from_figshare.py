"""Upload model prediction files for different modeling tasks and different models to
Figshare articles via API.

This script creates/updates a Figshare article containing model predictions from all
models in the Matbench Discovery benchmark. This includes both energy predictions,
ML-relaxed structures, and symmetry analysis files.
"""

import argparse
import os
import tomllib
from collections.abc import Sequence
from typing import Any, Final

import yaml
from tqdm import tqdm

import matbench_discovery.figshare as figshare
from matbench_discovery import PKG_DIR, ROOT
from matbench_discovery.data import Model, round_trip_yaml
from matbench_discovery.models import MODEL_METADATA

with open(f"{PKG_DIR}/modeling-tasks.yml") as file:
    MODELING_TASKS: Final = yaml.safe_load(file)

with open(f"{ROOT}/pyproject.toml", mode="rb") as toml_file:
    pyproject = tomllib.load(toml_file)["project"]


def parse_args(args: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--models",
        nargs="+",
        default=[m.name for m in Model],
        choices=[m.name for m in Model],
        help=(
            "Space-separated list of model names to download. If not provided, all "
            "models will be downloaded."
        ),
    )
    parser.add_argument(
        "--tasks",
        nargs="+",
        choices=list(MODELING_TASKS),
        default=list(MODELING_TASKS),
        help=(
            "Space-separated list of modeling tasks to download. Defaults to all tasks."
        ),
    )
    parser.add_argument(
        "-n",
        "--dry-run",
        action="store_true",
        help="Print what would be downloaded without actually downloading",
    )

    return parser.parse_args(args)


def download_one_modeling_task_article(
    task: str, models: list[str], *, dry_run: bool = False
) -> None:
    """Update or create a Figshare article for a modeling task."""
    article_id = figshare.ARTICLE_IDS[f"model_preds_{task}"]

    dry_run = True

    if article_id is not None:
        # Check if article exists and is accessible
        if figshare.article_exists(article_id):
            print(f"\nFound existing article for {task=} with ID {article_id}")
        else:
            raise ValueError(f"\nArticle {article_id} for {task=} not found")

    article_url = f"{figshare.ARTICLE_URL_PREFIX}/{article_id}"
    print(f"Now downloading article at {article_url}")

    if dry_run:
        print("\nDry run mode - no files will be downloaded")

    existing_files = figshare.get_existing_files(article_id)
    print(f"Found {len(existing_files)} existing files:")
    for idx, (file_name, file_data) in enumerate(existing_files.items(), start=1):
        print(f"{idx}. {file_name}: {file_data.get('id')}")

    downloaded_files: dict[str, str] = {}

    for model_name in tqdm(models):
        model: Model = getattr(Model, model_name)
        if not os.path.isfile(model.yaml_path):
            print(
                f"Warning: missing model metadata file {model.yaml_path}, skipping..."
            )
            continue

        with open(model.yaml_path) as file:
            model_data = round_trip_yaml.load(file)

        metrics = model_data.get("metrics", {})
        metric_data = metrics.get(task, {})
        if not isinstance(metric_data, dict):
            continue

        # Recursively find all keys ending in _file in the metric_data dictionary
        def find_file_keys(data: dict[str, Any], prefix: str = "") -> dict[str, str]:
            """Find all keys ending in _file and their values in a nested dictionary."""
            result: dict[str, str] = {}
            for key, value in data.items():
                full_key = f"{prefix}.{key}" if prefix else key
                if isinstance(value, dict):
                    result |= find_file_keys(value, full_key)
                elif isinstance(value, str) and key.endswith("_file"):
                    result[full_key] = value
            return result

        for key_path, rel_file_path in find_file_keys(metric_data).items():
            file_path = f"{ROOT}/{rel_file_path}"
            if not os.path.isfile(file_path):
                print(
                    f"Warning: {task} file for {model_name} not found, "
                    f"expected at {file_path}"
                )
                continue

            # Check if file already exists and has same hash
            file_hash, _ = figshare.get_file_hash_and_size(file_path)
            filename = file_path.removeprefix(f"{ROOT}/")

    print(f"Downloaded {len(downloaded_files)} files:")
    if downloaded_files:
        for idx, (filename, url) in enumerate(downloaded_files.items(), start=1):
            print(f"{idx}. {filename}: {url}")
    else:
        print("\nNo files were downloaded.")


def main(args: Sequence[str] | None = None) -> int:
    """Main function to upload model prediction files to Figshare."""
    parsed_args = parse_args(args)
    models_to_update = parsed_args.models or sorted(MODEL_METADATA)
    tasks_to_update = parsed_args.tasks
    if dry_run := parsed_args.dry_run:
        print("\nDry run mode - no files will be uploaded")
    print(f"Updating {len(models_to_update)} models: {', '.join(models_to_update)}")
    print(f"Updating {len(tasks_to_update)} tasks: {', '.join(tasks_to_update)}")

    for task in tasks_to_update:
        try:
            download_one_modeling_task_article(task, models_to_update, dry_run=dry_run)
        except Exception as exc:  # prompt to delete article if something went wrong
            state = {
                key: locals().get(key)
                for key in ("task", "model_name", "models_to_update", "tasks_to_update")
            }
            exc.add_note(f"Upload failed with {state=}")
            raise

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
