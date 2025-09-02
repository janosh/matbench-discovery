"""Helper functions for uploading files to Figshare via their API."""

import difflib
import hashlib
import json
import os
from collections.abc import Mapping, Sequence
from typing import Any, Final, cast

import requests
from tqdm import tqdm

from matbench_discovery import ROOT

ENV_PATH: Final[str] = f"{ROOT}/site/.env"
BASE_URL: Final[str] = "https://api.figshare.com/v2"
CONFIG: Final[dict[str, float]] = {"timeout": 30.0}

# Maps modeling tasks to their Figshare article IDs. New figshare articles will be
# created if the ID is None. Be sure to paste the new article ID into the
# ARTICLE_IDS dict below! It'll be printed by this script.
ARTICLE_URL_PREFIX: Final = "https://figshare.com/articles/dataset"
DOWNLOAD_URL_PREFIX: Final = "https://figshare.com/files"
ARTICLE_IDS: Final[dict[str, int | None]] = {
    "model_preds_discovery": 28187990,
    "model_preds_geo_opt": 28642406,
    "model_preds_phonons": 28347251,
    "model_preds_diatomics": 28437344,
    "data_files": 22715158,
}

# category IDs can be found at https://api.figshare.com/v2/categories
CATEGORIES: Final[dict[int, str]] = {
    25162: "Structure and dynamics of materials",
    25144: "Inorganic materials (incl. nanomaterials)",
    25186: "Cheminformatics and Quantitative Structure-Activity Relationships",
}


FIGSHARE_TOKEN = os.getenv("FIGSHARE_TOKEN")
if not FIGSHARE_TOKEN and os.path.isfile(ENV_PATH):
    with open(ENV_PATH) as file:
        # TOKEN: length 128, alphanumeric (e.g. 271431c6a94ff7...)
        FIGSHARE_TOKEN = file.read().split("figshare_token=")[1].split("\n")[0]


def set_timeout(timeout_seconds: float) -> None:
    """Set the HTTP request timeout used for all Figshare API calls.

    Args:
        timeout_seconds (float): Timeout in seconds. Must be greater than 0.
    """
    if timeout_seconds <= 0:
        raise ValueError("timeout_seconds must be greater than 0")
    CONFIG["timeout"] = timeout_seconds


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
    response = requests.request(
        method,
        url,
        headers=headers,
        data=data,
        timeout=CONFIG["timeout"],
    )
    try:
        response.raise_for_status()
        try:
            return json.loads(response.content)
        except ValueError:
            return response.content
    except requests.HTTPError as exc:
        exc.add_note(f"body={response.content.decode()}")
        raise


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


def upload_file(article_id: int, file_path: str, file_name: str = "") -> int:
    """Upload a file to Figshare and return the file ID.

    Args:
        article_id (int): ID of the article to upload to.
        file_path (str): Path to the file to upload.
        file_name (str, optional): Name as it will appear in Figshare. Defaults to the
            file path relative to repo's root dir: file_path.removeprefix(ROOT).

    Returns:
        int: The ID of the uploaded file.
    """
    # Initiate new upload
    md5, size = get_file_hash_and_size(file_path)
    file_name = file_name or file_path.removeprefix(f"{ROOT}/")
    data = dict(name=file_name, md5=md5, size=size)
    endpoint = f"{BASE_URL}/account/articles/{article_id}/files"
    result = make_request("POST", endpoint, data=data)
    file_info = make_request("GET", result["location"])

    # Upload parts with nested progress bar showing bytes and percent
    url = file_info["upload_url"]
    parts_info = make_request("GET", url)
    with (
        open(file_path, mode="rb") as file,
        tqdm(
            total=size,
            unit="B",
            unit_scale=True,
            unit_divisor=1024,
            position=1,
            desc=f"Uploading {file_name}",
            leave=False,
        ) as pbar,
    ):
        for part in parts_info["parts"]:
            # Upload part
            part_url = f"{file_info['upload_url']}/{part['partNo']}"
            file.seek(part["startOffset"])
            chunk_len = part["endOffset"] - part["startOffset"] + 1
            chunk = file.read(chunk_len)
            make_request("PUT", part_url, data=chunk, binary=True)
            pbar.update(len(chunk))
            pbar.set_postfix_str(f"{pbar.n / 1024**2:.2f}/{size / 1024**2:.2f} MB")

    # Complete upload
    make_request("POST", f"{endpoint}/{file_info['id']}")
    return file_info["id"]


def article_exists(article_id: int | str) -> bool:
    """Check if a Figshare article exists and is accessible.

    Args:
        article_id (int | str): The ID or URL of the article to check.

    Returns:
        bool: True if the article exists and is accessible, False otherwise.
    """
    article_url = (
        f"{BASE_URL}/account/articles/{article_id}"
        if isinstance(article_id, int)
        else article_id
    )
    try:
        make_request("GET", article_url)
    except requests.HTTPError as exc:
        if exc.response.status_code == 404:
            return False
        exc.add_note(f"{article_url=}")
        raise
    else:
        return True


def list_article_files(article_id: int) -> list[dict[str, Any]]:
    """Get a list of files in a Figshare article.

    Args:
        article_id (int): ID of the article to list files from.

    Returns:
        list[dict[str, Any]]: List of file information dictionaries. Each dictionary
            contains keys like 'name', 'id', 'size', 'computed_md5', etc.
            Empty list if article doesn't exist.

    Raises:
        requests.HTTPError: If the request fails for any reason other than 404.
    """
    try:
        return make_request("GET", f"{BASE_URL}/account/articles/{article_id}/files")
    except requests.HTTPError as exc:
        if exc.response.status_code == 404:
            return []
        raise


def get_existing_files(article_id: int) -> dict[str, dict[str, Any]]:
    """Get a mapping of filenames to dict with file details (usually id and md5 hash)
    for files already in the article.
    """
    try:
        files = make_request("GET", f"{BASE_URL}/account/articles/{article_id}/files")
        return {file.pop("name"): file for file in files}
    except requests.HTTPError as exc:
        if exc.response.status_code == 404:
            return {}
        raise


def file_exists_with_same_hash(
    article_id: int, file_name: str, file_hash: str
) -> tuple[bool, int | None]:
    """Check if a file with the same name and hash already exists in the article.

    Args:
        article_id (int): ID of the article to check.
        file_name (str): Name of the file to check.
        file_hash (str): MD5 hash of the file to check.

    Returns:
        tuple[bool, int | None]: A tuple containing:
            - bool: True if a file with the same name and hash exists, False otherwise.
            - int | None: The file ID if it exists, None otherwise.
    """
    existing_files = get_existing_files(article_id)
    if file_name in existing_files:
        existing_file = existing_files[file_name]
        if existing_file.get("computed_md5") == file_hash:
            return True, existing_file.get("id")
    return False, None


def find_similar_files(
    filename: str,
    existing_files: dict[str, dict[str, Any]],
    similarity_threshold: float = 0.7,
) -> list[tuple[str, int]]:
    """Find similar files using string similarity.

    Files must have same model family/subfolder, exceed similarity threshold,
    and have matching task types (kappa/phonon/discovery/geo) if present.

    Args:
        filename: File being uploaded
        existing_files: Existing files dictionary
        similarity_threshold: Min similarity ratio (0-1), default 0.7

    Returns:
        List of (name, id) tuples for similar files
    """
    parts = filename.split("/")
    if len(parts) < 3:
        return []

    model_family, model_subfolder, base_filename = parts[1], parts[2], parts[-1]
    task_type = _extract_task_type(base_filename)

    similar_files: list[tuple[str, int]] = []
    for existing_name, file_data in existing_files.items():
        existing_parts = existing_name.split("/")

        # Skip if not enough parts or different model family/subfolder
        if (
            len(existing_parts) < 3
            or existing_parts[1] != model_family
            or existing_parts[2] != model_subfolder
        ):
            continue

        existing_basename = existing_parts[-1]
        existing_task_type = _extract_task_type(existing_basename)

        # Skip if different task types
        if task_type and existing_task_type and task_type != existing_task_type:
            continue

        # Check similarity
        similarity = difflib.SequenceMatcher(
            None, base_filename, existing_basename
        ).ratio()
        if similarity >= similarity_threshold:
            file_id = cast("int", file_data.get("id"))
            similar_files.append((existing_name, file_id))

    return similar_files


def _extract_task_type(filename: str) -> str:
    """Extract task type (kappa, phonon, discovery, geo) from filename."""
    task_types = ["kappa", "phonon", "discovery", "geo"]
    parts = filename.split("-")
    return next(
        (task for part in parts for task in task_types if part.startswith(task)), ""
    )


def delete_file(article_id: int, file_id: int) -> bool:
    """Delete a file from a Figshare article.

    Args:
        article_id (int): ID of the article containing the file.
        file_id (int): ID of the file to delete.

    Returns:
        bool: True if the file was successfully deleted, False otherwise.
    """
    url = f"{BASE_URL}/account/articles/{article_id}/files/{file_id}"
    try:
        make_request("DELETE", url)  # should return None
        print(f"Successfully deleted file with ID {file_id} from article {article_id}")
        return True  # noqa: TRY300
    except Exception as exc:
        print(f"Failed to delete file with ID {file_id}: {exc}")
        return False


def upload_file_if_needed(
    article_id: int,
    file_path: str,
    file_name: str = "",
    *,
    force_reupload: bool = False,
) -> tuple[int, bool]:
    """Upload a file to Figshare if it doesn't already exist with the same hash.

    Args:
        article_id (int): ID of the article to upload to.
        file_path (str): Path to the file to upload.
        file_name (str, optional): Name as it will appear in Figshare. Defaults to the
            file path relative to repo's root dir: file_path.removeprefix(ROOT).
        force_reupload (bool, optional): If True, delete and reupload the file even if
            it already exists with the same hash. Defaults to False.

    Returns:
        tuple[int, bool]: A tuple containing:
            - int: The file ID.
            - bool: True if the file was uploaded, False if it already existed.
    """
    file_name = file_name or file_path.removeprefix(f"{ROOT}/")
    file_hash, _ = get_file_hash_and_size(file_path)

    # Check if file already exists with same hash
    exists, file_id = file_exists_with_same_hash(article_id, file_name, file_hash)

    if exists and file_id is not None:
        if force_reupload:
            print(
                f"{file_name=} exists but force_reupload=True, deleting and reuploading"
            )
            if delete_file(article_id, file_id):
                # Upload the file after successful deletion
                file_id = upload_file(article_id, file_path, file_name)
                return file_id, True
            print(f"Failed to delete existing file {file_name}, skipping upload")
            return file_id, False
        print(f"File {file_name} already exists with same hash, skipping upload")
        return file_id, False

    # Upload the file
    file_id = upload_file(article_id, file_path, file_name)
    return file_id, True


def publish_article(article_id: int, *, verbose: bool = True) -> bool:
    """Publish a Figshare article, making it publicly available.

    Args:
        article_id (int): ID of the article to publish.
        verbose (bool, optional): Whether to print status messages. Defaults to True.

    Returns:
        bool: True if the article was successfully published.
    """
    post_url = f"{BASE_URL}/account/articles/{article_id}/publish"
    try:
        make_request("POST", post_url)
        if verbose:
            article_url = f"{ARTICLE_URL_PREFIX}/{article_id}"
            print(f"Successfully published article {article_id} at {article_url}")
    except Exception as exc:
        if verbose:
            print(f"Failed to publish article {article_id}: {exc}")
        return False
    return True
