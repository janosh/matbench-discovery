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

import matbench_discovery.remote.figshare as figshare
from matbench_discovery import PKG_DIR, ROOT
from matbench_discovery.data import round_trip_yaml
from matbench_discovery.enums import Model

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
        nargs="*",
        type=Model,  # type: ignore[arg-type]
        choices=Model,
        default=list(Model),
        help="Models to analyze. If none specified, analyzes all models.",
    )
    parser.add_argument(
        "--tasks",
        nargs="+",
        choices=list(MODELING_TASKS),
        default=list(MODELING_TASKS),
        help=(
            "Space-separated list of modeling tasks to update. Defaults to all tasks."
        ),
    )
    parser.add_argument(
        "-n",
        "--dry-run",
        action="store_true",
        help="Print what would be uploaded without actually uploading",
    )

    return parser.parse_args(args)


def get_article_metadata(task: str) -> dict[str, Sequence[object]]:
    """Get metadata for creating a new Figshare article for a modeling task."""
    task_info = MODELING_TASKS[task]
    return {
        "title": f"Matbench Discovery - Model Predictions for {task_info['label']}",
        "description": f"""
        This dataset contains model predictions from various models evaluated on the
        Matbench Discovery benchmark for the {task_info["label"].lower()} task.

        Task description: {task_info["description"]}

        For more information about the benchmark and models, visit:
        https://github.com/janosh/matbench-discovery.
        """.strip(),
        "defined_type": "dataset",
        "tags": [*pyproject["keywords"], f"task-{task}"],
        "categories": list(figshare.CATEGORIES),
    }


def update_one_modeling_task_article(
    task: str, models: list[Model], *, dry_run: bool = False
) -> None:
    """Update or create a Figshare article for a modeling task."""
    article_id = figshare.ARTICLE_IDS[f"model_preds_{task}"]

    if article_id is not None:
        # Check if article exists and is accessible
        if figshare.article_exists(article_id):
            print(f"\nFound existing article for {task=} with ID {article_id}")
        else:
            print(f"\nArticle {article_id} for {task=} not found")
            article_id = None

    if article_id is None:
        if dry_run:
            print(f"\nWould create new article for {task=}")
            article_id = 0
        else:
            metadata = get_article_metadata(task)
            article_id = figshare.create_article(metadata)
            print(
                f"\n⚠️ Created new Figshare article for {task=} with {article_id=}"
                f"\nUpdate FIGSHARE_ARTICLE_IDS in {__file__} with this ID!"
            )

    article_url = f"{figshare.ARTICLE_URL_PREFIX}/{article_id}"
    print(f"Now updating article at {article_url}")

    if dry_run:
        print("\nDry run mode - no files will be uploaded")

    existing_files = figshare.get_existing_files(article_id)
    print(f"Found {len(existing_files)} existing files:")
    for idx, (file_name, file_data) in enumerate(existing_files.items(), start=1):
        print(f"{idx}. {file_name}: {file_data.get('id')}")

    updated_files: dict[str, str] = {}  # files that were re-uploaded
    new_files: dict[str, str] = {}  # files that didn't exist before

    for model in tqdm(models):
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
                    f"Warning: {task} file for {model.name} not found, "
                    f"expected at {file_path}"
                )
                continue

            # Check if file already exists and has same hash
            file_hash, _ = figshare.get_file_hash_and_size(file_path)
            filename = file_path.removeprefix(f"{ROOT}/")

            if filename in existing_files:
                file_id, md5_hash = (
                    existing_files[filename][key] for key in ("id", "computed_md5")
                )
                if file_hash == md5_hash:
                    file_url = f"{figshare.DOWNLOAD_URL_PREFIX}/{file_id}"
                    # Update model metadata if URL not present
                    url_key = f"{key_path}_url"  # append _url to YAML key
                    if url_key not in metric_data:
                        *parts, last = url_key.split(".")
                        target = metric_data
                        for part in parts:
                            target = target[part]
                        target[last] = file_url
                else:
                    # Upload modified file
                    if not dry_run:
                        file_id = figshare.upload_file(
                            article_id,
                            file_path,
                            file_name=filename,
                        )
                        file_url = f"{figshare.DOWNLOAD_URL_PREFIX}/{file_id}"
                        updated_files[filename] = file_url
                        *parts, last = key_path.split(".")
                        target = metric_data
                        for part in parts:
                            target = target[part]
                        target[f"{last}_url"] = file_url
            else:
                # Upload new file
                if not dry_run:
                    file_id = figshare.upload_file(
                        article_id,
                        file_path,
                        file_name=filename,
                    )
                    file_url = f"{figshare.DOWNLOAD_URL_PREFIX}/{file_id}"
                    new_files[filename] = file_url
                    *parts, last = key_path.split(".")
                    target = metric_data
                    for part in parts:
                        target = target[part]
                    target[f"{last}_url"] = file_url

        # Save updated model metadata if changed
        if not dry_run:
            with open(model.yaml_path, mode="w") as file:
                round_trip_yaml.dump(model_data, file)

    print(f"Newly added: {len(new_files)}")
    print(f"Updated: {len(updated_files)}")

    if new_files or updated_files:
        if new_files:
            print("\nNewly added files:")
            for idx, (filename, url) in enumerate(new_files.items(), start=1):
                print(f"{idx}. {filename}: {url}")

        if updated_files:
            print("\nUpdated files:")
            for idx, (filename, url) in enumerate(updated_files.items(), start=1):
                print(f"{idx}. {filename}: {url}")
    else:
        print("\nNo files were added or updated.")


def main(args: Sequence[str] | None = None) -> int:
    """Main function to upload model prediction files to Figshare."""
    parsed_args = parse_args(args)
    models_to_update = parsed_args.models
    tasks_to_update = parsed_args.tasks
    if dry_run := parsed_args.dry_run:
        print("\nDry run mode - no files will be uploaded")
    print(f"Updating {len(models_to_update)} models: {', '.join(models_to_update)}")
    print(f"Updating {len(tasks_to_update)} tasks: {', '.join(tasks_to_update)}")

    for task in tasks_to_update:
        try:
            update_one_modeling_task_article(task, models_to_update, dry_run=dry_run)
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
