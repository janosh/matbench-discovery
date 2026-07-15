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

from matbench_discovery import PKG_DIR, ROOT, repo_relative_path
from matbench_discovery.cli import cli_parser, models_arg
from matbench_discovery.data import iter_file_refs, round_trip_yaml
from matbench_discovery.enums import Model
from matbench_discovery.remote import figshare

with open(f"{ROOT}/pyproject.toml", mode="rb") as toml_file:
    pyproject = tomllib.load(toml_file)["project"]


def process_exclusion_prefixes(items: list[str], all_items: list[str]) -> list[str]:
    """Apply ``!`` exclusions to explicit items or the complete item list."""
    include_items = [item for item in items if not item.startswith("!")]
    excluded_items = {item[1:] for item in items if item.startswith("!")}
    # If there are explicit inclusions, use those
    # Otherwise, start with all items and remove exclusions
    return [item for item in include_items or all_items if item not in excluded_items]


def get_article_metadata(task_info: dict[str, Any]) -> dict[str, Sequence[object]]:
    """Get metadata for creating a new Figshare article for a modeling task."""
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
        "tags": [*pyproject["keywords"], f"task-{task_info['label']}"],
        "categories": list(figshare.CATEGORIES),
    }


def resolve_artifact_path(
    model_yaml_path: str,
    relative_file_path: str,
    root_dir: str = ROOT,
) -> str:
    """Resolve an artifact file confined to its model architecture directory."""
    if os.path.isabs(relative_file_path):
        raise ValueError(f"Artifact path must be relative: {relative_file_path!r}")

    file_path = os.path.abspath(f"{root_dir}/{relative_file_path}")
    model_dir = os.path.realpath(os.path.dirname(model_yaml_path))
    if os.path.commonpath((model_dir, os.path.realpath(file_path))) != model_dir:
        raise ValueError(
            f"Artifact path escapes model directory {model_dir!r}: "
            f"{relative_file_path!r}"
        )
    if not os.path.isfile(file_path):
        raise FileNotFoundError(file_path)
    return file_path


def set_file_ref_url(
    data: dict[str, Any],
    key_path: tuple[str, ...],
    file_url: str,
    *,
    size: int,
    md5: str,
) -> None:
    """Set ``url``, ``size`` and ``md5`` on a nested FileRef."""
    *parts, last = key_path
    target = data
    for part in parts:
        target = target[part]
    ref = target.get(last)
    if not isinstance(ref, dict):
        raise TypeError(f"Expected FileRef object at {'.'.join(key_path)}, got {ref!r}")
    target[last] = {**ref, "url": file_url, "size": size, "md5": md5}


def update_one_modeling_task_article(
    task: str,
    models: list[Model],
    *,
    modeling_tasks: dict[str, dict[str, str]],
    dry_run: bool = False,
    file_type: Literal["all", "analysis", "pred"] = "all",
    force_reupload: bool = False,
    interactive: bool = True,
) -> None:
    """Update or create a Figshare article for a modeling task."""
    if (task_info := modeling_tasks.get(task)) is None:
        raise KeyError(
            f"Missing task metadata for {task!r}. "
            f"Available tasks: {list(modeling_tasks)}"
        )
    article_id = figshare.ARTICLE_IDS[f"model_preds_{task}"]
    article_is_new = False

    # Check if article exists and is accessible
    if article_id is not None and figshare.article_exists(article_id):
        print(f"\nFound existing article for {task=} with ID {article_id}")
    elif article_id is not None:
        print(f"\nArticle {article_id} for {task=} not found")
        article_id = None

    if article_id is None:
        if dry_run:
            print(f"\nWould create new article for {task=}")
            article_id = 0
        else:
            article_id = figshare.create_article(get_article_metadata(task_info))
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

        # Ingestion schema-validates model YAMLs before archival, so malformed FileRefs
        # intentionally fail fast here rather than being silently skipped.
        for key_parts, rel_file_path in iter_file_refs(metric_data):
            if file_type != "all" and key_parts[-1] != f"{file_type}_file":
                continue
            try:
                file_path = resolve_artifact_path(model.yaml_path, rel_file_path)
            except FileNotFoundError:
                print(
                    f"Warning: {task} file for {model.name} not found, "
                    f"expected at {ROOT}/{rel_file_path}"
                )
                continue

            filename = repo_relative_path(file_path)

            # Hash feeds the exact-match check (skipped when forcing) and the upload
            # (skipped on dry runs); a forced dry run needs neither. Reused below.
            file_hash: str | None = None
            file_size: int | None = None
            if not (dry_run and force_reupload):
                file_hash, file_size = figshare.get_file_hash_and_size(file_path)

            # First check if the exact same file already exists
            if not force_reupload and file_hash is not None and file_size is not None:
                exists, file_id = figshare.file_exists_with_same_hash(
                    article_id, filename, file_hash
                )

                if exists and file_id is not None:
                    file_url = f"{figshare.DOWNLOAD_URL_PREFIX}/{file_id}"
                    skipped_files[filename] = (file_url, model)
                    set_file_ref_url(
                        metric_data,
                        key_parts,
                        file_url,
                        size=file_size,
                        md5=file_hash,
                    )
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
                if not interactive:
                    print("Skipping deletion of similar files (non-interactive mode)")
                elif (
                    input("Delete these files before uploading? [y/N] ").lower() == "y"
                ):
                    for similar_name, similar_id in similar_files:
                        if similar_id is not None and figshare.delete_file(
                            article_id, similar_id
                        ):
                            deleted_files[similar_name] = similar_id
                            # keep local view in sync to classify uploads correctly
                            existing_files.pop(similar_name, None)
                            print(f"Deleted similar file: {similar_name}")

            # Upload file if it doesn't exist or force_reupload is True
            if not dry_run:
                if file_hash is None or file_size is None:
                    raise RuntimeError(f"Missing hash/size for {file_path}")
                file_id, _was_uploaded = figshare.upload_file_if_needed(
                    article_id,
                    file_path,
                    file_name=filename,
                    force_reupload=force_reupload,
                )
                file_url = f"{figshare.DOWNLOAD_URL_PREFIX}/{file_id}"

                target_files = (
                    updated_files if filename in existing_files else new_files
                )
                target_files[filename] = (file_url, model)

                # Update model metadata with URL
                set_file_ref_url(
                    metric_data,
                    key_parts,
                    file_url,
                    size=file_size,
                    md5=file_hash,
                )

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

        for heading, files in (
            ("Newly added files", new_files),
            ("Updated files", updated_files),
            ("Skipped files (already exist with same hash)", skipped_files),
        ):
            if files:
                print(f"\n{heading}:")
                for idx, (filename, (url, model)) in enumerate(files.items(), start=1):
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
    """Upload selected model prediction artifacts to Figshare."""
    with open(f"{PKG_DIR}/modeling-tasks.yml", encoding="utf-8") as file:
        modeling_tasks = yaml.safe_load(file)
    # remove 'cps' task as it's a dynamic metric with changing weights
    # no point in uploading to figshare
    modeling_tasks.pop("cps", None)

    # Add figshare-specific arguments to the central CLI parser
    figshare_group = cli_parser.add_argument_group(
        "figshare", "Arguments for Figshare upload functionality"
    )
    figshare_group.add_argument(
        "--tasks",
        nargs="+",
        type=str,
        default=list(modeling_tasks),
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

    models_arg.type = str
    models_arg.choices = None
    models_arg.default = [model.name for model in Model.active()]
    args, _unknown = cli_parser.parse_known_args(raw_args)

    # Process exclusion prefixes for tasks
    args.tasks = process_exclusion_prefixes(args.tasks, list(modeling_tasks))

    # Process exclusion prefixes for models
    processed_models = process_exclusion_prefixes(
        args.models, [model.name for model in Model]
    )
    args.models = [Model(model_name) for model_name in processed_models]

    if dry_run := args.dry_run:
        print("\nDry run mode - no files will be uploaded")
    print(
        f"Updating {len(args.models)} models: "
        f"{', '.join(model.name for model in args.models)}"
    )
    print(f"Updating {len(args.tasks)} tasks: {', '.join(args.tasks)}")
    print(f"File type filter: {args.file_type}")
    if args.force_reupload:
        print("Force reupload: True - will reupload files even if they already exist")

    for task in args.tasks:
        try:
            update_one_modeling_task_article(
                task,
                args.models,
                modeling_tasks=modeling_tasks,
                dry_run=dry_run,
                file_type=args.file_type,
                force_reupload=args.force_reupload,
                interactive=not args.no_interactive,
            )
        except Exception as exc:  # prompt to delete article if something went wrong
            exc.add_note(f"Upload failed for {task=}")
            raise

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
