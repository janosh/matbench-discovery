"""Files download functions."""

import builtins
import hashlib
import os
import sys
import tempfile
from urllib.parse import urlsplit

import requests

from matbench_discovery import file_digest


def download_file(
    file_path: str,
    url: str,
    headers: dict[str, str] | None = None,
    md5: str | None = None,
) -> None:
    """Atomically download a file, raising HTTP, checksum, and filesystem errors.

    Args:
        file_path: Local path to write the downloaded file to.
        url: URL to download from.
        headers: Optional HTTP headers, e.g. an Authorization bearer token for
            gated HuggingFace checkpoints.
        md5: Optional expected MD5 checksum. On mismatch, the download is discarded
            and any previously cached file_path is left unchanged.
    """
    file_dir = os.path.dirname(file_path)
    if file_dir:
        os.makedirs(file_dir, exist_ok=True)
    # Convert any Figshare URL variant to the API download endpoint to avoid WAF
    # Handles: figshare.com/files/ID, figshare.com/ndownloader/files/ID,
    # and ndownloader.figshare.com/files/ID
    parsed = urlsplit(url)
    hostname = parsed.hostname or ""
    if (
        hostname == "figshare.com" or hostname.endswith(".figshare.com")
    ) and "/files/" in parsed.path:
        file_id = parsed.path.rsplit("/files/", maxsplit=1)[-1]
        url = f"https://api.figshare.com/v2/file/download/{file_id}"

    is_huggingface = parsed.scheme == "https" and (
        hostname == "huggingface.co" or hostname.endswith(".huggingface.co")
    )
    request_headers = dict(headers or {})
    if (
        is_huggingface
        and not any(key.casefold() == "authorization" for key in request_headers)
        and (token := os.getenv("HF_TOKEN") or os.getenv("HUGGING_FACE_HUB_TOKEN"))
    ):
        request_headers["Authorization"] = f"Bearer {token}"
    file_descriptor, tmp_file_path = tempfile.mkstemp(
        dir=file_dir or ".", prefix=f".{os.path.basename(file_path)}.", suffix=".part"
    )
    os.close(file_descriptor)
    download_finished = False
    try:
        # Stream large files to avoid loading entire file into memory
        file_hash = hashlib.md5()  # noqa: S324
        with (
            requests.get(
                url, timeout=600, stream=True, headers=request_headers or None
            ) as response,
            open(tmp_file_path, mode="wb") as file,
        ):
            response.raise_for_status()
            for chunk in response.iter_content(chunk_size=8192):
                file.write(chunk)
                file_hash.update(chunk)

        if md5 and (actual_md5 := file_hash.hexdigest()) != md5:
            raise ValueError(
                f"MD5 mismatch for {url=}: expected {md5}, got {actual_md5}. "
                f"{file_path=} left unchanged."
            )

        # set flag before replace: if only the rename fails, the fully-downloaded
        # .part file is deliberately kept so it doesn't have to be re-fetched
        download_finished = True
        os.replace(tmp_file_path, file_path)
    except Exception as exc:
        exc.add_note(f"Downloading {url=} to {file_path=}")
        if (
            isinstance(exc, requests.HTTPError)
            and exc.response is not None
            and exc.response.status_code in (401, 403)
            and is_huggingface
        ):
            exc.add_note(
                "For gated HuggingFace repos, accept the model license and set "
                "HF_TOKEN in the environment. Check the token's access to this repo."
            )
        if download_finished:
            exc.add_note(f"Completed download retained at {tmp_file_path!r}")
        else:
            try:
                os.remove(tmp_file_path)
            except FileNotFoundError:
                pass
            except OSError as cleanup_error:
                exc.add_note(
                    f"Failed to remove partial download {tmp_file_path=}: "
                    f"{cleanup_error}"
                )
        raise


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
        if (cached_md5 := file_digest(abs_path, "md5")) == md5:
            return
        print(
            f"Cached file {abs_path!r} has MD5 {cached_md5}, expected {md5}; "
            "re-downloading"
        )

    auto_download_files = os.getenv("MBD_AUTO_DOWNLOAD_FILES", "true").lower() == "true"

    # default to 'y' if auto-download is on or no interactive session (TTY/iPython)
    is_ipython = hasattr(builtins, "__IPYTHON__")
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
