"""Update Matbench Discovery Data Files Figshare article via API."""

import argparse
import os
import tomllib
from collections.abc import Sequence
from typing import Any

from tqdm import tqdm

from matbench_discovery import DATA_DIR, PKG_DIR, ROOT
from matbench_discovery.data import round_trip_yaml
from matbench_discovery.enums import DataFiles
from matbench_discovery.remote import figshare


def main(
    yaml_path: str,
    article_id: int | None,
    pyproject: dict[str, Any],
    max_file_size: int = 500 * 1024**2,  # in bytes, files larger than this should
    # be uploaded manually to avoid timeouts
    files: Sequence[DataFiles] | None = None,
    *,
    dry_run: bool = False,
) -> int:
    """Keep Data Files Figshare article in sync with data-files.yml file.

    Args:
        yaml_path (str): Path to the YAML file containing the data files.
        article_id (int | None): ID of the Figshare article to update. If None, a new
            article will be created.
        pyproject (dict[str, Any]): Dictionary containing the project metadata.
        max_file_size (int): Maximum file size in bytes to be uploaded automatically.
        files (list[DataFiles] | None): Which DataFiles to upload (must exist
            locally). Defaults to all local ones.
        dry_run (bool): Only print planned creates, uploads, publishes and YAML
            writes. Reads the remote file listing but changes nothing.
    """
    existing_yaml: dict[str, dict[str, str]] = {}
    if os.path.isfile(yaml_path):
        with open(yaml_path, encoding="utf-8") as file:
            existing_yaml = round_trip_yaml.load(file)

    # validate every file before touching Figshare so a bad request changes nothing
    def is_local(data_file: DataFiles) -> bool:
        """Whether the data file exists locally."""
        return os.path.isfile(f"{DATA_DIR}/{data_file.rel_path}")

    if files is None:
        files = [data_file for data_file in DataFiles if is_local(data_file)]
    elif missing := [data_file.name for data_file in files if not is_local(data_file)]:
        raise FileNotFoundError(f"Requested --files missing locally: {missing}")
    if undescribed := [
        data_file.name
        for data_file in files
        if not (existing_yaml.get(data_file.name) or {}).get("description")
    ]:
        raise ValueError(
            f"Add a description for {undescribed} in {yaml_path} before uploading"
        )

    article_is_new = False
    if article_id is not None:
        if figshare.article_exists(article_id):
            article_url = f"{figshare.ARTICLE_URL_PREFIX}/{article_id}"
            print(f"\nFound existing article at {article_url}")
        else:
            print(f"\n{article_id=} not found")
            article_id = None

    description = """
    Matbench Discovery evaluates machine learning interatomic potentials across
    materials-science tasks, from crystal stability prediction to finite-temperature
    dynamics.

    These files provide training-set structures, benchmark inputs, and reference
    data for task-specific evaluations.

    For a description of each file, see
    https://matbench-discovery.materialsproject.org/data.

    The original MPtrj training set, with energies, forces, stresses and (partial)
    magmoms, is available at https://figshare.com/articles/dataset/23713842.
    """.replace("\n", " ").strip()
    metadata = {
        "title": "Matbench Discovery - Data Files",
        "description": description,
        "defined_type": "dataset",  # seems to be default anyway
        "tags": pyproject["keywords"],
        "categories": list(figshare.CATEGORIES),
        "references": list(pyproject["urls"].values()),
    }

    if article_id is None:
        article_is_new = True
        if dry_run:
            print("\n[dry-run] Would create a new Figshare article")
        else:
            article_id = figshare.create_article(metadata)
            print(
                f"\n⚠️ Created new Figshare article with {article_id=}"
                f"\nSet ARTICLE_IDS['data_files'] in {figshare.__file__} to this ID!"
            )
    article_url = f"{figshare.ARTICLE_URL_PREFIX}/{article_id}"

    try:
        # remote listing lets us skip re-uploading unchanged files
        existing_files = (
            [] if article_id is None else figshare.list_article_files(article_id)
        )
        print(f"Found {len(existing_files)} existing files on Figshare")
        files_by_name = {file["name"]: file for file in existing_files}

        # copy so entries for files not touched by this run survive
        files_in_article: dict[str, dict[str, str]] = existing_yaml.copy()
        updated_files: dict[str, str] = {}  # re-uploaded
        new_files: dict[str, str] = {}  # first upload of a registry entry

        for data_file in (pbar := tqdm(files)):
            pbar.set_description(f"Processing {data_file.name}")
            file_path = f"{DATA_DIR}/{data_file.rel_path}"

            # copy so we don't mutate existing_yaml
            file_data = existing_yaml[data_file.name].copy()
            file_data.setdefault("path", data_file.rel_path)

            # Always record the freshly computed hash: reusing a stored md5 here once
            # let re-uploads leave stale registry checksums (#357).
            # Files upload under their DATA_DIR-relative path (e.g. wbm/x.csv.gz), so
            # the remote listing must be matched by that name, not the basename.
            file_name = data_file.rel_path
            file_hash, file_size = figshare.get_file_hash_and_size(file_path)
            file_data["md5"] = file_hash
            files_in_article[data_file.name] = file_data

            if (
                existing_file := files_by_name.get(file_name)
            ) and file_hash == existing_file["computed_md5"]:
                file_id = existing_file["id"]  # local file matches remote, skip upload
            elif file_size > max_file_size:
                # can't upload -> remote artifact unchanged -> keep the previous
                # consistent (md5, url) pair (fresh md5 + old/no url would fail every
                # checksummed download). The run after the manual upload hits the
                # hash-match branch above and refreshes md5+url together.
                print(
                    f"\n⚠️  Skipping {file_name} ({file_size / 1024**2:.1f} MB)"
                    f"\nFile exceeds {max_file_size / 1024**2:.0f} MB limit. "
                    "Please upload manually at "
                    f"https://figshare.com/account/articles/{article_id}"
                )
                files_in_article[data_file.name] = existing_yaml[data_file.name]
                continue
            else:  # upload new or modified file, tracked for the summary print below
                # registry entries without a URL were never uploaded before
                is_update = "url" in existing_yaml[data_file.name]
                dest = updated_files if is_update else new_files
                if dry_run or article_id is None:
                    print(f"\n[dry-run] Would upload {file_path} as {file_name!r}")
                    dest[data_file.name] = "(not uploaded: dry run)"
                    continue
                file_id = figshare.upload_file(
                    article_id, file_path, file_name=file_name
                )
                dest[data_file.name] = f"{figshare.DOWNLOAD_URL_PREFIX}/{file_id}"
            file_data["url"] = f"{figshare.DOWNLOAD_URL_PREFIX}/{file_id}"

        if new_files or updated_files:
            if new_files:
                print("\nNewly added files:")
                for idx, (data_file, url) in enumerate(new_files.items(), start=1):
                    print(f"{idx}. {data_file}: {url}")

            if updated_files:
                print("\nUpdated files:")
                for idx, (data_file, url) in enumerate(updated_files.items(), start=1):
                    print(f"{idx}. {data_file}: {url}")

            # newly created articles need manual review before their first publish
            if dry_run or article_id is None:  # no article ID only in dry runs
                action = (
                    "leave new article unpublished" if article_is_new else "publish"
                )
                print(f"\n[dry-run] Would {action} {article_url}")
            elif article_is_new:
                print(
                    "\n⚠️ Article was newly created. Please review it at "
                    f"{article_url} before publishing manually."
                )
            else:
                print(
                    f"\nFiles were added or updated. Publishing article {article_url}"
                )
                figshare.publish_article(article_id)
        else:
            print("\nNo files were added or updated.")

        if dry_run:
            print(f"[dry-run] Would update {yaml_path}")
        else:
            with open(yaml_path, mode="w", encoding="utf-8") as file:
                round_trip_yaml.dump(files_in_article, file)

    except Exception as exc:  # add context to exception for better debugging
        state = {
            key: locals().get(key) for key in ("article_id", "file_path", "data_file")
        }
        exc.add_note(f"Upload failed with {state=}")
        raise

    return 0


def build_parser() -> argparse.ArgumentParser:
    """Build the CLI parser (strict, so typos like --dryrun fail instead of
    silently uploading and publishing).
    """
    parser = argparse.ArgumentParser(
        description="Upload data files to Figshare", allow_abbrev=False
    )
    parser.add_argument(
        "--files",
        type=DataFiles,
        nargs="+",
        help="DataFiles to upload (must exist locally). Defaults to all local files.",
    )
    parser.add_argument(
        "-n",
        "--dry-run",
        action="store_true",
        help="Print planned uploads/publishes without changing Figshare or the YAML",
    )
    return parser


if __name__ == "__main__":
    args = build_parser().parse_args()

    with open(f"{ROOT}/pyproject.toml", mode="rb") as toml_file:
        pyproject = tomllib.load(toml_file)["project"]

    raise SystemExit(
        main(
            yaml_path=f"{PKG_DIR}/data-files.yml",
            article_id=figshare.ARTICLE_IDS["data_files"],
            pyproject=pyproject,
            files=args.files,
            dry_run=args.dry_run,
        )
    )
