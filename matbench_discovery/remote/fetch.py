"""Files download functions."""

import builtins
import hashlib
import os
import sys
import traceback

import requests


def _headers_for_url(
    url: str, headers: dict[str, str] | None = None
) -> dict[str, str] | None:
    """Return request headers, adding HuggingFace bearer auth when available."""
    request_headers = dict(headers or {})
    if "huggingface.co" in url and "Authorization" not in request_headers:
        token = os.getenv("HF_TOKEN") or os.getenv("HUGGING_FACE_HUB_TOKEN")
        if token:
            request_headers["Authorization"] = f"Bearer {token}"
    return request_headers or None


def download_file(
    file_path: str,
    url: str,
    headers: dict[str, str] | None = None,
    md5: str | None = None,
) -> None:
    """Download the file from the given URL to the given file path.
    Prints rather than raises if the file cannot be downloaded.

    Args:
        file_path: Local path to write the downloaded file to.
        url: URL to download from.
        headers: Optional HTTP headers, e.g. an Authorization bearer token for
            gated HuggingFace checkpoints.
        md5: Optional expected MD5 checksum. On mismatch, the download is discarded
            (any previously cached file_path is left unchanged) and an error printed.
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

    headers = _headers_for_url(url, headers)
    try:
        # Stream large files to avoid loading entire file into memory
        response = requests.get(url, timeout=600, stream=True, headers=headers)
        response.raise_for_status()

        file_hash = hashlib.md5()  # noqa: S324
        with open(tmp_file_path, mode="wb") as file:
            for chunk in response.iter_content(chunk_size=8192):
                file.write(chunk)
                file_hash.update(chunk)

        if md5 and (actual_md5 := file_hash.hexdigest()) != md5:
            os.remove(tmp_file_path)
            print(
                f"MD5 mismatch for {url=}: expected {md5}, got {actual_md5}. "
                f"Discarded the download, {file_path=} left unchanged."
            )
            return

        # set flag before replace: if only the rename fails, the fully-downloaded
        # .part file is deliberately kept so it doesn't have to be re-fetched
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


def maybe_auto_download_file(
    url: str,
    abs_path: str,
    label: str | None = None,
    headers: dict[str, str] | None = None,
    md5: str | None = None,
) -> None:
    """Download a missing or checksum-stale file after confirmation."""
    if os.path.isfile(abs_path):
        if md5 is None:
            return
        with open(abs_path, mode="rb") as file:
            cached_md5 = hashlib.file_digest(file, "md5").hexdigest()
        if cached_md5 == md5:
            return
        print(
            f"Cached file {abs_path!r} has MD5 {cached_md5}, expected {md5}; "
            "re-downloading"
        )

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
        download_file(abs_path, url, headers=headers, md5=md5)
