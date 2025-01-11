"""Upload model prediction files from all models to a Figshare article.

This script creates/updates a Figshare article containing model predictions from all
models in the Matbench Discovery benchmark. This includes both energy predictions and
ML-relaxed structures.
"""

import argparse
import hashlib
import os
from collections.abc import Sequence
from typing import Final

import requests
import yaml
from tqdm import tqdm

from matbench_discovery import PKG_DIR, ROOT
from matbench_discovery.data import Model, round_trip_yaml
from matbench_discovery.figshare import (
    BASE_URL,
    create_article,
    make_request,
    upload_file_to_figshare,
)
from matbench_discovery.models import MODEL_METADATA

with open(f"{PKG_DIR}/modeling-tasks.yml") as file:
    MODELING_TASKS: Final = yaml.safe_load(file)

# Maps modeling tasks to their Figshare article IDs. New figshare articles will be
# created if the ID is None. Be sure to paste the new article ID into the
# FIGSHARE_ARTICLE_IDS dict below! It'll be printed by this script.
ARTICLE_URL_PREFIX = "https://figshare.com/articles/dataset"
FIGSHARE_ARTICLE_IDS = {
    "discovery": 28187990,
    "geo_opt": 28187999,
    "phonons": None,
}


def get_file_hash_and_size(
    file_name: str, chunk_size: int = 10_000_000
) -> tuple[str, int]:
    """Get the md5 hash and size of a file."""
    md5 = hashlib.md5()  # noqa: S324
    size = 0
    with open(file_name, mode="rb") as file:
        while data := file.read(chunk_size):
            size += len(data)
            md5.update(data)
    return md5.hexdigest(), size


def get_existing_files(article_id: int) -> dict[str, tuple[int, str]]:
    """Get a mapping of filenames to (file_id, md5) for files already in the article."""
    try:
        files = make_request("GET", f"{BASE_URL}/account/articles/{article_id}/files")
        return {file["name"]: (file["id"], file["computed_md5"]) for file in files}
    except requests.HTTPError as exc:
        if exc.response.status_code == 404:
            return {}
        raise


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
            "Space-separated list of model names to update. If not provided, all "
            "models will be updated."
        ),
    )
    parser.add_argument(
        "--tasks",
        nargs="+",
        choices=list(MODELING_TASKS),
        default=["discovery"],
        help=(
            "Space-separated list of modeling tasks to update. Defaults to 'discovery'."
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
        "tags": [
            "materials-science",
            "machine-learning",
            "crystal-structure-prediction",
            f"task-{task}",
        ],
        "categories": [
            25162,  # Structure and dynamics of materials
            25144,  # Inorganic materials (incl. nanomaterials)
            25186,  # Cheminformatics and Quantitative Structure-Activity Relationships
        ],
    }


def update_one_modeling_task_article(
    task: str, models: list[str], *, dry_run: bool = False
) -> None:
    """Update or create a Figshare article for a modeling task."""
    article_id = FIGSHARE_ARTICLE_IDS[task]

    if article_id is not None:
        # Check if article exists and is accessible
        try:
            make_request("GET", f"{BASE_URL}/account/articles/{article_id}")
            print(f"\nFound existing article for {task} task with ID {article_id}")
        except requests.HTTPError as exc:
            if exc.response.status_code == 404:
                print(f"\nArticle {article_id} for {task} task not found")
                article_id = None
            else:
                raise

    if article_id is None:
        if dry_run:
            print(f"\nWould create new article for {task} task")
            article_id = 0
        else:
            metadata = get_article_metadata(task)
            article_id = create_article(metadata)
            print(
                f"\n⚠️ Created new Figshare article for {task=} with {article_id=}"
                f"\nUpdate FIGSHARE_ARTICLE_IDS in {__file__} with this ID!"
            )

    article_url = f"{ARTICLE_URL_PREFIX}/{article_id}"
    print(f"Updating article with ID {article_id} at {article_url}")

    if dry_run:
        print("\nDry run mode - no files will be uploaded")

    existing_files = get_existing_files(article_id)
    print(f"Found {len(existing_files)} existing files")

    uploaded_files: dict[str, str] = {}
    new_files: dict[str, str] = {}

    pbar = tqdm(models)
    for model_name in pbar:
        model = getattr(Model, model_name)
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

        pred_file = metric_data.get("pred_file")
        if not pred_file:
            continue

        if not os.path.isfile(f"{ROOT}/{pred_file}"):
            print(
                f"Warning: {task} predictions for {model_name} not "
                f"found, expected at {pred_file}"
            )
            continue

        pbar.set_description(f"Processing {model_name}")
        model_updated = False

        filename = os.path.basename(pred_file)
        file_path = f"{ROOT}/{pred_file}"

        # Check if file already exists and has same hash
        file_hash, _ = get_file_hash_and_size(file_path)
        if filename in existing_files:
            file_id, stored_hash = existing_files[filename]
            if file_hash == stored_hash:
                file_url = f"https://figshare.com/ndownloader/files/{file_id}"
                uploaded_files[filename] = file_url
                # Update model metadata if URL not present
                if "pred_file_url" not in metric_data:
                    metric_data["pred_file_url"] = file_url
                    model_updated = True
                continue

        # Upload new or modified file
        if dry_run:
            file_url = "DRY_RUN_URL"
        else:
            file_id = upload_file_to_figshare(article_id, file_path)
            file_url = f"https://figshare.com/ndownloader/files/{file_id}"

        uploaded_files[filename] = file_url
        new_files[filename] = file_url

        # Update model metadata
        metric_data["pred_file_url"] = file_url
        model_updated = True

        # Save updated model metadata if changed
        if model_updated and not dry_run:
            with open(model.yaml_path, "w") as file:
                round_trip_yaml.dump(model_data, file)

    print(f"\nTotal files: {len(uploaded_files)}")
    print(f"Newly added: {len(new_files)}")

    print("\nAll uploaded files:")
    for filename, url in uploaded_files.items():
        print(f"{filename}: {url}")

    if new_files:
        print("\nNewly added files:")
        for filename, url in new_files.items():
            print(f"{filename}: {url}")
    else:
        print("\nNo new files were added.")


def main(args: Sequence[str] | None = None) -> int:
    """Main function to upload model prediction files to Figshare."""
    parsed_args = parse_args(args)
    models_to_update = parsed_args.models or sorted(MODEL_METADATA)
    tasks_to_update = parsed_args.tasks
    print(f"Updating {len(models_to_update)} models: {', '.join(models_to_update)}")
    print(f"Updating {len(tasks_to_update)} tasks: {', '.join(tasks_to_update)}")

    for task in tasks_to_update:
        try:
            update_one_modeling_task_article(
                task, models_to_update, dry_run=parsed_args.dry_run
            )
        except Exception as exc:  # prompt to delete article if something went wrong
            msg = f"Upload failed for {task=}"
            if model_name := locals().get("model_name", ""):
                msg += f" and {model_name=}"
            exc.add_note(msg)
            raise

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
