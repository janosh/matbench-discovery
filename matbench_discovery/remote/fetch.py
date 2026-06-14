"""Files download functions."""

import builtins
import os
import shutil
import sys
import tarfile
import traceback

import requests


def download_file(
    file_path: str, url: str, headers: dict[str, str] | None = None
) -> None:
    """Download the file from the given URL to the given file path.
    Prints rather than raises if the file cannot be downloaded.

    Args:
        file_path: Local path to write the downloaded file to.
        url: URL to download from.
        headers: Optional HTTP headers, e.g. an Authorization bearer token for
            gated HuggingFace checkpoints.
    """
    file_dir = os.path.dirname(file_path)
    if file_dir:
        os.makedirs(file_dir, exist_ok=True)
    tmp_file_path = f"{file_path}.part"
    download_finished = False

    # Convert any Figshare URL variant to the API download endpoint to avoid WAF
    # Handles: figshare.com/files/ID, figshare.com/ndownloader/files/ID,
    # and ndownloader.figshare.com/files/ID
    if "figshare.com" in url and "/files/" in url:
        file_id = url.rsplit("/files/", maxsplit=1)[-1].split("?", maxsplit=1)[0]
        url = f"https://api.figshare.com/v2/file/download/{file_id}"

    try:
        # Stream large files to avoid loading entire file into memory
        response = requests.get(url, timeout=600, stream=True, headers=headers)
        response.raise_for_status()

        with open(tmp_file_path, mode="wb") as file:
            file.writelines(response.iter_content(chunk_size=8192))

        download_finished = True
        os.replace(tmp_file_path, file_path)
    except (OSError, requests.RequestException):
        error_msg = traceback.format_exc()
        if not download_finished:
            try:
                os.remove(tmp_file_path)
            except FileNotFoundError:
                pass
            except OSError:
                error_msg += (
                    f"\nFailed to remove partial download at {tmp_file_path=}.\n"
                    f"{traceback.format_exc()}"
                )
        print(f"Error downloading {url=}\nto {file_path=}.\n{error_msg}")


def extract_tar_if_needed(tar_path: str) -> str:
    """Extract a .tar archive into a sibling directory named after the archive
    (minus the .tar suffix) unless that directory already exists, and return that
    directory. Extraction goes to a temp directory moved into place only when
    complete, so interrupted extractions can't be mistaken for finished ones.
    """
    if not tar_path.endswith(".tar"):
        raise ValueError(f"Expected a .tar file, got {tar_path!r}")
    extract_dir = tar_path.removesuffix(".tar")
    if os.path.isdir(extract_dir):
        return extract_dir

    # sibling staging dir (not tempfile.mkdtemp: os.replace needs same filesystem);
    # clear any leftovers from a previously interrupted extraction
    tmp_dir = f"{extract_dir}.extracting"
    shutil.rmtree(tmp_dir, ignore_errors=True)
    with tarfile.open(tar_path) as archive:
        # the "data" extraction filter (PEP 706) only exists in Python 3.11.4+;
        # older 3.11 patches (allowed by requires-python >=3.11) lack the kwarg
        if hasattr(tarfile, "data_filter"):
            archive.extractall(tmp_dir, filter="data")
        else:
            archive.extractall(tmp_dir)  # noqa: S202
    os.replace(tmp_dir, extract_dir)
    print(f"Extracted {tar_path!r} to {extract_dir!r}")
    return extract_dir


def maybe_auto_download_file(url: str, abs_path: str, label: str | None = None) -> None:
    """Download file if not exist and user confirms or auto-download is enabled."""
    if os.path.isfile(abs_path):
        return

    # whether to auto-download model prediction files without prompting
    auto_download_files = os.getenv("MBD_AUTO_DOWNLOAD_FILES", "true").lower() == "true"

    is_ipython = hasattr(builtins, "__IPYTHON__")
    # default to 'y' if auto-download enabled or not in interactive session (TTY
    # or iPython)
    answer = (
        "y"
        if auto_download_files or not (is_ipython or sys.stdin.isatty())
        else input(
            f"{abs_path!r} associated with {label=} does not exist. Download it "
            "now? This will cache the file for future use. [y/n] "
        )
    )
    if answer.lower().strip() == "y":
        print(f"Downloading {label!r} from {url!r} to {abs_path!r}")
        download_file(abs_path, url)
