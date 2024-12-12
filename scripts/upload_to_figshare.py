"""Create a new Figshare article via API.

Adapted from this example notebook:
https://colab.research.google.com/drive/13CAM8mL1u7ZsqNhfZLv7bNb1rdhMI64d?usp=sharing
Found notebook in docs: https://help.figshare.com/article/how-to-use-the-figshare-api
"""

import hashlib
import json
import os
import tomllib  # needs python 3.11
from typing import Any

import requests
from tqdm import tqdm

from matbench_discovery import DATA_DIR, PKG_DIR, ROOT
from matbench_discovery.data import DataFiles, round_trip_yaml

__author__ = "Janosh Riebesell"
__date__ = "2023-04-27"

with open(f"{ROOT}/site/.env") as file:
    # TOKEN: length 128, alphanumeric (e.g. 271431c6a94ff7...)
    TOKEN = file.read().split("figshare_token=")[1].split("\n")[0]

BASE_URL = "https://api.figshare.com/v2"


def make_request(
    method: str, url: str, *, data: Any = None, binary: bool = False
) -> Any:
    """Make a token-authorized HTTP request to the Figshare API.

    Args:
        method (str): HTTP method (GET, POST, PUT, DELETE).
        url (str): URL to send the request to.
        data (Any, optional): Data to send in the request body. Defaults to None.
        binary (bool, optional): Whether the data is binary. Defaults to False.

    Returns:
        Any: JSON response data or binary data.

    Raises:
        HTTPError: If the request fails. Error will contain the response body.
    """
    headers = {"Authorization": f"token {TOKEN}"}
    if data is not None and not binary:
        data = json.dumps(data)
    response = requests.request(method, url, headers=headers, data=data, timeout=10)
    try:
        response.raise_for_status()
        try:
            data = json.loads(response.content)
        except ValueError:
            data = response.content
    except requests.HTTPError as exc:
        exc.add_note(f"body={response.content.decode()}")
        raise

    return data


def create_article(metadata: dict[str, str | int | float]) -> int:
    """Create a new Figshare article with given metadata and return the article ID."""
    result = make_request("POST", f"{BASE_URL}/account/articles", data=metadata)
    print(f"Created article: {result['location']}\n")
    result = make_request("GET", result["location"])
    return result["id"]


def get_file_hash_and_size(
    file_name: str, chunk_size: int = 10_000_000
) -> tuple[str, int]:
    """Get the md5 hash and size of a file. File is read in chunks of chunk_size bytes.
    Default chunk size is 10_000_000 ~= 10MB.
    """
    md5 = hashlib.md5()  # noqa: S324
    size = 0
    with open(file_name, mode="rb") as file:
        while data := file.read(chunk_size):
            size += len(data)
            md5.update(data)
    return md5.hexdigest(), size


def upload_file_to_figshare(article_id: int, file_path: str) -> int:
    """Upload a file to Figshare and return the file ID."""
    # Initiate new upload
    md5, size = get_file_hash_and_size(file_path)
    data = dict(name=os.path.basename(file_path), md5=md5, size=size)
    endpoint = f"{BASE_URL}/account/articles/{article_id}/files"
    result = make_request("POST", endpoint, data=data)
    file_info = make_request("GET", result["location"])

    # Upload parts
    url = file_info["upload_url"]
    result = make_request("GET", url)
    with open(file_path, mode="rb") as file:
        for part in tqdm(result["parts"], desc=file_path):
            # Upload part
            u_data = file_info.copy()
            u_data.update(part)
            url = f"{u_data['upload_url']}/{part['partNo']}"
            file.seek(part["startOffset"])
            chunk = file.read(part["endOffset"] - part["startOffset"] + 1)
            make_request("PUT", url, data=chunk, binary=True)

    # Complete upload
    make_request("POST", f"{endpoint}/{file_info['id']}")
    return file_info["id"]


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
    {pyproject['description'].lower()}. It contains relaxed structures of the MP
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
