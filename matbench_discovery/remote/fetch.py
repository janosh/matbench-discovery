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
    try:
        response = requests.get(url, timeout=5)

        response.raise_for_status()

        with open(file_path, "wb") as file:
            file.write(response.content)
    except requests.exceptions.RequestException:
        print(f"Error downloading {url=}\nto {file_path=}.\n{traceback.format_exc()}")


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
