"""Update Matbench Discovery Data Files Figshare article via API."""

import os
import tomllib
from typing import Any

from tqdm import tqdm

import matbench_discovery.figshare as figshare
from matbench_discovery import DATA_DIR, PKG_DIR, ROOT
from matbench_discovery.data import DataFiles, round_trip_yaml

__author__ = "Janosh Riebesell"
__date__ = "2023-04-27"


def main(
    yaml_path: str,
    article_id: int | None,
    pyproject: dict[str, Any],
    max_file_size: int = 500 * 1024**2,  # in bytes, files larger than this should
    # be uploaded manually to avoid timeouts
) -> int:
    """Keep Data Files Figshare article in sync with data-files.yml file.

    Args:
        yaml_path (str): Path to the YAML file containing the data files.
        article_id (int | None): ID of the Figshare article to update. If None, a new
            article will be created.
        pyproject (dict[str, Any]): Dictionary containing the project metadata.
        max_file_size (int): Maximum file size in bytes to be uploaded automatically.
    """
    if article_id is not None:
        # Check if article exists and is accessible
        if figshare.article_exists(article_id):
            article_url = f"{figshare.ARTICLE_URL_PREFIX}/{article_id}"
            print(f"\nFound existing article at {article_url}")
        else:
            print(f"\n{article_id=} not found")
            article_id = None

    DESCRIPTION = """
    These are the Matbench Discovery data files, a benchmark for machine learning
    interatomic potentials (MLIP) on inorganic crystals. We evaluate models on
    multiple tasks including crystal stability prediction from unrelaxed structures,
    geometry optimization, modeling of harmonic and anharmonic phonons with more
    soon to come.

    The data files include relaxed structures of the MP training set, initial+relaxed
    structures of the WBM test set, both in pymatgen.Structure and ase.Atoms format.

    For a description of each file, see
    https://matbench-discovery.materialsproject.org/data#--direct-download.

    The original MPtrj training set containing 1.3M structures with their energies,
    forces, stresses and (partial) magmoms is available at https://figshare.com/articles/dataset/23713842.
    """.replace("\n", " ").strip()
    metadata = {
        "title": "Matbench Discovery - Data Files",
        "description": DESCRIPTION,
        "defined_type": "dataset",  # seems to be default anyway
        "tags": pyproject["keywords"],
        "categories": list(figshare.CATEGORIES),
        "references": list(pyproject["urls"].values()),
    }

    if article_id is None:
        article_id = figshare.create_article(metadata)
        print(
            f"\n⚠️ Created new Figshare article with {article_id=}"
            f"\nUpdate FIGSHARE_ARTICLE_ID in {__file__} with this ID!"
        )

    try:
        # Load existing YAML data if available
        existing_yaml: dict[str, dict[str, str]] = {}
        if os.path.isfile(yaml_path):
            with open(yaml_path) as file:
                existing_yaml = round_trip_yaml.load(file)

        # Get existing files from Figshare to avoid re-uploading unchanged files
        existing_files = figshare.list_article_files(article_id)
        print(f"Found {len(existing_files)} existing files on Figshare")

        # Create lookup dict for faster file checks
        files_by_name = {file["name"]: file for file in existing_files}

        # copy existing_data to preserve all existing entries
        files_in_article: dict[str, dict[str, str]] = existing_yaml.copy()
        updated_files: dict[str, str] = {}  # files that were re-uploaded
        new_files: dict[str, str] = {}  # files that didn't exist before

        pbar = tqdm(DataFiles)
        for data_file in pbar:
            pbar.set_description(f"Processing {data_file.name}")
            file_path = f"{DATA_DIR}/{data_file.rel_path}"

            if not os.path.isfile(file_path):
                print(f"Warning: {file_path} does not exist, skipping...")
                continue

            # Get existing data or create new entry (make copy to not modify original)
            file_data = existing_yaml.get(data_file.name, {}).copy()

            # Only set defaults for new entries
            if not file_data:
                file_data["path"] = data_file.rel_path
                file_data["description"] = "Description needed"

            # Check if file needs to be uploaded
            filename = os.path.basename(file_path)
            file_size = os.path.getsize(file_path)

            if file_size > max_file_size:
                print(
                    f"\n⚠️  Skipping {filename} ({file_size / 1024**2:.1f} MB)"
                    f"\nFile exceeds {max_file_size / 1024**2:.0f} MB limit. "
                    "Please upload manually at "
                    f"https://figshare.com/account/articles/{article_id}"
                )
                continue

            # Use stored MD5 if available, otherwise compute it
            if "md5" not in file_data:
                file_hash, _ = figshare.get_file_hash_and_size(file_path)
                file_data["md5"] = file_hash
            file_hash = file_data["md5"]

            if (
                existing_file := files_by_name.get(filename)
            ) and file_hash == existing_file["computed_md5"]:
                file_url = f"{figshare.DOWNLOAD_URL_PREFIX}/{existing_file['id']}"
                file_data["url"] = file_url
                files_in_article[data_file.name] = file_data
                continue

            # Upload new or modified file
            file_id = figshare.upload_file(article_id, file_path)
            file_url = f"{figshare.DOWNLOAD_URL_PREFIX}/{file_id}"

            file_data["url"] = file_url
            files_in_article[data_file.name] = file_data

            # Track whether file is new or updated
            if data_file.name in existing_yaml:
                updated_files[data_file.name] = file_url
            else:
                new_files[data_file.name] = file_url

        if new_files or updated_files:
            if new_files:
                print("\nNewly added files:")
                for idx, (data_file, url) in enumerate(new_files.items(), start=1):
                    print(f"{idx}. {data_file}: {url}")

            if updated_files:
                print("\nUpdated files:")
                for idx, (data_file, url) in enumerate(updated_files.items(), start=1):
                    print(f"{idx}. {data_file}: {url}")
        else:
            print("\nNo files were added or updated.")

        # Write updated YAML file
        with open(yaml_path, mode="w") as file:
            round_trip_yaml.dump(files_in_article, file)

    except Exception as exc:  # add context to exception for better debugging
        state = {
            key: locals().get(key) for key in ("article_id", "file_path", "data_file")
        }
        exc.add_note(f"Upload failed with {state=}")
        raise

    return 0


if __name__ == "__main__":
    with open(f"{ROOT}/pyproject.toml", mode="rb") as toml_file:
        pyproject = tomllib.load(toml_file)["project"]

    # upload/update data files to figshare
    main(
        yaml_path=f"{PKG_DIR}/data-files.yml",
        article_id=figshare.ARTICLE_IDS["data_files"],
        pyproject=pyproject,
    )
