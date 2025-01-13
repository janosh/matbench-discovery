"""Helper functions for uploading files to Figshare via their API."""

import hashlib
import json
import os
from collections.abc import Mapping, Sequence
from typing import Any

import requests
from tqdm import tqdm

from matbench_discovery import ROOT

ENV_PATH = f"{ROOT}/site/.env"
BASE_URL = "https://api.figshare.com/v2"

FIGSHARE_TOKEN = os.getenv("FIGSHARE_TOKEN")
if not FIGSHARE_TOKEN and os.path.isfile(ENV_PATH):
    with open(ENV_PATH) as file:
        # TOKEN: length 128, alphanumeric (e.g. 271431c6a94ff7...)
        FIGSHARE_TOKEN = file.read().split("figshare_token=")[1].split("\n")[0]


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
    headers = {"Authorization": f"token {FIGSHARE_TOKEN}"}
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


def create_article(
    metadata: Mapping[str, Sequence[object]], *, verbose: bool = True
) -> int:
    """Create a new Figshare article with given metadata and return the article ID.

    Args:
        metadata (dict): Article metadata including title, description, etc.
        verbose (bool, optional): Whether to print the article URL and title.
            Defaults to True.

    Returns:
        int: The ID of the created article.
    """
    result = make_request("POST", f"{BASE_URL}/account/articles", data=metadata)
    if verbose:
        print(f"Created article: {result['location']} with title {metadata['title']}\n")
    result = make_request("GET", result["location"])
    return result["id"]


def get_file_hash_and_size(
    file_name: str, chunk_size: int = 10_000_000
) -> tuple[str, int]:
    """Get the md5 hash and size of a file.

    Args:
        file_name (str): Path to the file.
        chunk_size (int, optional): Size of chunks to read in bytes. Defaults to 10MB.

    Returns:
        tuple[str, int]: MD5 hash and file size in bytes.
    """
    md5 = hashlib.md5()  # noqa: S324
    size = 0
    with open(file_name, mode="rb") as file:
        while data := file.read(chunk_size):
            size += len(data)
            md5.update(data)
    return md5.hexdigest(), size


def upload_file_to_figshare(article_id: int, file_path: str) -> int:
    """Upload a file to Figshare and return the file ID.

    Args:
        article_id (int): ID of the article to upload to.
        file_path (str): Path to the file to upload.

    Returns:
        int: The ID of the uploaded file.
    """
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
