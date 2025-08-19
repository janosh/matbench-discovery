"""Upload model prediction files for different modeling tasks and different models to
Figshare articles via API.

This script creates/updates a Figshare article containing model predictions from all
models in the Matbench Discovery benchmark. This includes both energy predictions,
ML-relaxed structures, and symmetry analysis files.
"""

import os
import tomllib
from collections.abc import Sequence
from typing import Any, Literal

import yaml
from tqdm import tqdm

import matbench_discovery.remote.figshare as figshare
from matbench_discovery import PKG_DIR, ROOT
from matbench_discovery.cli import cli_parser
from matbench_discovery.data import round_trip_yaml
from matbench_discovery.enums import Model

with open(f"{PKG_DIR}/modeling-tasks.yml", encoding="utf-8") as file:
    MODELING_TASKS = yaml.safe_load(file)
# remove 'cps' task as it's a dynamic metric with changing weights
# no point in uploading to figshare
MODELING_TASKS.pop("cps", None)

with open(f"{ROOT}/pyproject.toml", mode="rb") as toml_file:
    pyproject = tomllib.load(toml_file)["project"]


def process_exclusion_prefixes(items: list[str], all_items: list[str]) -> list[str]:
    """Process items with exclusion prefixes (!) and return the final list.

    Args:
        items: List of items, some of which may be prefixed with '!' for exclusion
        all_items: Complete list of all possible items

    Returns:
        List of items after processing exclusions
    """
    include_items, exclude_items = [], []
    for item in items:
        if isinstance(item, str) and item.startswith("!"):
            exclude_items.append(item.removeprefix("!"))
        else:
            include_items.append(item)

    # If there are explicit inclusions, use those
    # Otherwise, start with all items and remove exclusions
    result = include_items or all_items.copy()
    for item in exclude_items:
        if item in result:
            result.remove(item)

    return result


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


def should_process_file(
    key: str, file_type: Literal["all", "analysis", "pred"]
) -> bool:
    """Filter files by type."""
    return file_type == "all" or key.endswith(f"{file_type}_file")


def update_one_modeling_task_article(
    task: str,
    models: list[Model],
    *,
    dry_run: bool = False,
    file_type: Literal["all", "analysis", "pred"] = "all",
    force_reupload: bool = False,
    interactive: bool = True,
) -> None:
    """Update or create a Figshare article for a modeling task."""
    article_id = figshare.ARTICLE_IDS[f"model_preds_{task}"]
    article_is_new = False

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
            article_is_new = True
            print(
                f"\n⚠️ Created new Figshare article for {task=} with {article_id=}"
                f"\nUpdate figshare.ARTICLE_IDS with this ID!"
            )

    article_url = f"{figshare.ARTICLE_URL_PREFIX}/{article_id}"
    print(f"Now updating article at {article_url}")

    if dry_run:
        print("\nDry run mode - no files will be uploaded")

    existing_files = figshare.get_existing_files(article_id)
    print(f"Found {len(existing_files)} existing files:")
    for idx, (file_name, file_data) in enumerate(existing_files.items(), start=1):
        print(f"{idx}. {file_name}: {file_data.get('id')}")

    # files that were skipped because they already exist
    skipped_files: dict[str, tuple[str, Model]] = {}  # filename -> (url, model)
    updated_files: dict[str, tuple[str, Model]] = {}  # files that were re-uploaded
    new_files: dict[str, tuple[str, Model]] = {}  # files that didn't exist before
    deleted_files: dict[str, int] = {}  # filename -> id (files that were deleted)

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
                elif (
                    isinstance(value, str)
                    and key.endswith("_file")
                    and should_process_file(key, file_type)
                ):
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

            filename = file_path.removeprefix(f"{ROOT}{os.sep}")

            # First check if the exact same file already exists
            if not force_reupload and not dry_run:
                file_hash, _ = figshare.get_file_hash_and_size(file_path)
                exists, file_id = figshare.file_exists_with_same_hash(
                    article_id, filename, file_hash
                )

                if exists and file_id is not None:
                    file_url = f"{figshare.DOWNLOAD_URL_PREFIX}/{file_id}"
                    skipped_files[filename] = (file_url, model)

                    # Update model metadata if URL not present
                    url_key = f"{key_path}_url"  # append _url to YAML key
                    if url_key not in metric_data:
                        *parts, last = url_key.split(".")
                        target = metric_data
                        for part in parts:
                            target = target[part]
                        target[last] = file_url

                    continue

            # Check for similar files that should be deleted
            similar_files = figshare.find_similar_files(filename, existing_files)
            if similar_files and not dry_run:
                print(f"\nFound similar files for {filename}:")
                for idx, (similar_name, similar_id) in enumerate(
                    similar_files, start=1
                ):
                    print(f"{idx}. {similar_name} (ID: {similar_id})")

                # Ask for user confirmation to delete similar files
                if interactive:
                    confirm = input("Delete these files before uploading? [y/N] ")
                    if confirm.lower() == "y":
                        for similar_name, similar_id in similar_files:
                            if similar_id is not None and figshare.delete_file(
                                article_id, similar_id
                            ):
                                deleted_files[similar_name] = similar_id
                                # keep local view in sync to classify uploads correctly
                                existing_files.pop(similar_name, None)
                                print(f"Deleted similar file: {similar_name}")
                else:
                    print("Skipping deletion of similar files (non-interactive mode)")

            # Upload file if it doesn't exist or force_reupload is True
            if not dry_run:
                file_id, was_uploaded = figshare.upload_file_if_needed(
                    article_id,
                    file_path,
                    file_name=filename,
                    force_reupload=force_reupload,
                )
                file_url = f"{figshare.DOWNLOAD_URL_PREFIX}/{file_id}"

                if filename in existing_files:
                    updated_files[filename] = (file_url, model)
                else:
                    new_files[filename] = (file_url, model)

                # Update model metadata with URL
                *parts, last = key_path.split(".")
                target = metric_data
                for part in parts:
                    target = target[part]
                target[f"{last}_url"] = file_url

        # Save updated model metadata if changed
        if not dry_run:
            with open(model.yaml_path, mode="w") as file:
                round_trip_yaml.dump(model_data, file)

    # Extract unique models from the file dictionaries
    new_models = {model.name for _, model in new_files.values()}
    updated_models = {model.name for _, model in updated_files.values()}

    print(f"Newly added: {len(new_files)} [{', '.join(new_models)}]")
    print(f"Updated: {len(updated_files)} [{', '.join(updated_models)}]")
    print(f"Skipped (already exists with same hash): {len(skipped_files)}")
    print(f"Deleted: {len(deleted_files)}")

    if new_files or updated_files or skipped_files or deleted_files:
        if deleted_files:
            print("\nDeleted files:")
            for idx, (filename, file_id) in enumerate(deleted_files.items(), start=1):
                print(f"{idx}. {filename}: {file_id}")

        if new_files:
            print("\nNewly added files:")
            for idx, (filename, (url, model)) in enumerate(new_files.items(), start=1):
                print(f"{idx}. {model.name} {filename}: {url}")

        if updated_files:
            print("\nUpdated files:")
            for idx, (filename, (url, model)) in enumerate(
                updated_files.items(), start=1
            ):
                print(f"{idx}. {model.name} {filename}: {url}")

        if skipped_files:
            print("\nSkipped files (already exist with same hash):")
            for idx, (filename, (url, model)) in enumerate(
                skipped_files.items(), start=1
            ):
                print(f"{idx}. {model.name} {filename}: {url}")

        # Publish the article if any new files were added, it's not a dry run,
        # and the article wasn't newly created
        if (new_files or updated_files) and not dry_run and not article_is_new:
            print(f"\nFiles were added or updated. Publishing article {article_url}")
            figshare.publish_article(article_id)
        elif article_is_new:
            print(
                "\n⚠️ Article was newly created. Please review it at "
                f"{article_url} before publishing manually."
            )
    else:
        print("\nNo files were added or updated.")


def main(raw_args: Sequence[str] | None = None) -> int:
    """Main function to upload model prediction files to Figshare.

    Args:
        raw_args: Command line arguments. If None, sys.argv[1:] will be used.

    Returns:
        int: Exit code (0 for success).
    """
    # Add figshare-specific arguments to the central CLI parser
    figshare_group = cli_parser.add_argument_group(
        "figshare", "Arguments for Figshare upload functionality"
    )
    figshare_group.add_argument(
        "--tasks",
        nargs="+",
        type=str,
        default=list(MODELING_TASKS),
        help="Space-separated list of modeling tasks to update. Defaults to all tasks. "
        "Prefix with '!' to exclude. Note: exclamation mark needs to be "
        "backslash-escaped in shell.",
    )
    figshare_group.add_argument(
        "--file-type",
        choices=["all", "analysis", "pred"],
        default="all",
        help="Type of files to upload: analysis, pred or all (default)",
    )
    figshare_group.add_argument(
        "--force-reupload",
        action="store_true",
        help="Force reupload of files even if they already exist with the same hash",
    )
    figshare_group.add_argument(
        "--no-interactive",
        action="store_true",
        default=False,
        help="Disable interactive prompts for file deletion (files will be skipped)",
    )

    args, _unknown = cli_parser.parse_known_args(raw_args)

    # Process exclusion prefixes for tasks
    all_tasks = list(MODELING_TASKS)
    args.tasks = process_exclusion_prefixes(args.tasks, all_tasks)

    # Process exclusion prefixes for models
    if hasattr(args, "models") and args.models is not None:
        # Convert Model enum instances to strings for exclusion processing
        model_names = [model.name for model in args.models]
        all_models = [model.name for model in Model]
        processed_models = process_exclusion_prefixes(model_names, all_models)
        # Convert back to Model enum instances
        args.models = [Model[model] for model in processed_models]

    models_to_update = args.models
    tasks_to_update = args.tasks
    if dry_run := args.dry_run:
        print("\nDry run mode - no files will be uploaded")
    print(
        f"Updating {len(models_to_update)} models: "
        f"{', '.join(model.name for model in models_to_update)}"
    )
    print(f"Updating {len(tasks_to_update)} tasks: {', '.join(tasks_to_update)}")
    print(f"File type filter: {args.file_type}")
    if args.force_reupload:
        print("Force reupload: True - will reupload files even if they already exist")

    for task in tasks_to_update:
        try:
            update_one_modeling_task_article(
                task,
                models_to_update,
                dry_run=dry_run,
                file_type=args.file_type,
                force_reupload=args.force_reupload,
                interactive=not args.no_interactive,
            )
        except Exception as exc:  # prompt to delete article if something went wrong
            state = {
                key: locals().get(key)
                for key in ("task", "models_to_update", "tasks_to_update")
            }
            exc.add_note(f"Upload failed with {state=}")
            raise

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
