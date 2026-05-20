"""Files download functions."""

import builtins
import os
import sys
import traceback

import requests


def download_file(file_path: str, url: str) -> None:
    """Download the file from the given URL to the given file path.
    Prints rather than raises if the file cannot be downloaded.
    """
    file_dir = os.path.dirname(file_path)
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
        response = requests.get(url, timeout=600, stream=True)
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
