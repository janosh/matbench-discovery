"""Create a new Figshare article via API.

Adapted from this example notebook:
https://colab.research.google.com/drive/13CAM8mL1u7ZsqNhfZLv7bNb1rdhMI64d?usp=sharing
Found notebook in docs: https://help.figshare.com/article/how-to-use-the-figshare-api
"""

import os
import tomllib  # needs python 3.11
from typing import Any

from tqdm import tqdm

from matbench_discovery import DATA_DIR, PKG_DIR, ROOT
from matbench_discovery.data import DataFiles, round_trip_yaml
from matbench_discovery.figshare import (
    BASE_URL,
    create_article,
    make_request,
    upload_file_to_figshare,
)

__author__ = "Janosh Riebesell"
__date__ = "2023-04-27"


def main(pyproject: dict[str, Any], yaml_path: str) -> int:
    """Main function to upload all files in FILE_PATHS to the same Figshare article."""
    pkg_name, version = pyproject["name"], pyproject["version"]
    # category IDs can be found at https://api.figshare.com/v2/categories
    categories = {
        25162: "Structure and dynamics of materials",
        25144: "Inorganic materials (incl. nanomaterials)",
        25186: "Cheminformatics and Quantitative Structure-Activity Relationships",
    }

    DESCRIPTION = f"""
    These are the v{version} data files for Matbench Discovery,
    {pyproject["description"].lower()}. It contains relaxed structures of the MP
    training set, initial+relaxed structures of the WBM test set, plus several
    checkpoints  for models trained on this data specifically for this benchmark.

    For a description of each file, see
    https://matbench-discovery.materialsproject.org/contribute#--direct-download.

    The original MPtrj training set containing 1.3M structures with their energies,
    forces, stresses and (partial) magmoms is available at https://figshare.com/articles/dataset/23713842.
    """.replace("\n", " ").strip()
    metadata = {
        "title": f"{pkg_name.title()} v{version}",
        "description": DESCRIPTION,
        "defined_type": "dataset",  # seems to be default anyway
        "tags": pyproject["keywords"],
        "categories": list(categories),
        "references": list(pyproject["urls"].values()),
    }
    try:
        article_id = create_article(metadata)
        uploaded_files: dict[str, dict[str, str]] = {}

        existing_data: dict[str, dict[str, str]] = {}
        if os.path.isfile(yaml_path):  # Load existing YAML data if available
            with open(yaml_path) as file:
                existing_data = round_trip_yaml.load(file)

        pbar = tqdm(DataFiles, desc="Uploading to Figshare")
        for key in pbar:
            pbar.set_postfix(file=key)
            file_path = f"{DATA_DIR}/{key.rel_path}"
            file_id = upload_file_to_figshare(article_id, file_path)
            file_url = f"https://figshare.com/ndownloader/files/{file_id}"

            # Create new entry while preserving existing description if available
            uploaded_files[key] = {
                "url": file_url,
                "path": key.rel_path,
                "description": existing_data.get(key, {}).get(
                    "description", "Description needed"
                ),
            }

        print("\nUploaded files:")
        for key, file_info in uploaded_files.items():
            print(f"{key}: {file_info['url']}")

        # Manually copy MPtrj (special case since separate Figshare article from
        # Matbench Discovery)
        if mp_trj := existing_data.get("mp_trj"):
            uploaded_files["mp_trj"] = mp_trj

        # Preserve _links section from existing data or use default
        if _links := existing_data.get("_links"):
            uploaded_files["_links"] = _links

        with open(yaml_path, mode="w") as file:
            round_trip_yaml.dump(uploaded_files, file)

    except Exception as exc:  # prompt to delete article if something went wrong
        if file_path := str(locals().get("file_path", "")):
            exc.add_note(f"{file_path=}")
        answer, article_id = "", int(locals().get("article_id", 0))
        while article_id and answer.lower() not in ("y", "n"):
            answer = input("Delete article? [y/n] ")
        if answer.lower() == "y":
            make_request("DELETE", f"{BASE_URL}/account/articles/{article_id}")

    return 0


if __name__ == "__main__":
    with open(f"{ROOT}/pyproject.toml", mode="rb") as toml_file:
        pyproject = tomllib.load(toml_file)["project"]

    figshare_yaml_path = f"{PKG_DIR}/data-files.yml"
    # upload all data files to figshare with current pyproject.toml version
    main(pyproject, figshare_yaml_path)
