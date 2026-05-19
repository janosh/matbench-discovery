import os
from collections.abc import Iterator
from pathlib import Path
from queue import Queue
from threading import Barrier, Thread
from typing import Self
from unittest.mock import patch

import pytest
import requests

from matbench_discovery.remote.fetch import download_file, maybe_auto_download_file


def make_mock_response(content: bytes, status_code: int = 200) -> requests.Response:
    """Create a mock requests.Response with given content and status code."""
    response = requests.Response()
    response.status_code = status_code
    response._content = content  # noqa: SLF001
    response.iter_content = (
        lambda chunk_size=1, decode_unicode=False: [content]  # noqa: ARG005
    )
    return response


@pytest.mark.parametrize(
    "input_url, expected_url",
    [
        # Standard figshare.com/files/ format
        (
            "https://figshare.com/files/12345",
            "https://api.figshare.com/v2/file/download/12345",
        ),
        # ndownloader path variant
        (
            "https://figshare.com/ndownloader/files/99999",
            "https://api.figshare.com/v2/file/download/99999",
        ),
        # ndownloader subdomain variant
        (
            "https://ndownloader.figshare.com/files/55555",
            "https://api.figshare.com/v2/file/download/55555",
        ),
        # Query params stripped
        (
            "https://ndownloader.figshare.com/files/55555?access_token=abc",
            "https://api.figshare.com/v2/file/download/55555",
        ),
        # Non-figshare URL unchanged
        ("https://example.com/files/test.gz", "https://example.com/files/test.gz"),
    ],
)
def test_figshare_url_conversion(
    input_url: str,
    expected_url: str,
    tmp_path: Path,
) -> None:
    """Test that Figshare URL variants are converted to the API download endpoint."""
    dest = tmp_path / "out.gz"
    test_content = b"mock data"
    with patch(
        "requests.get", return_value=make_mock_response(test_content)
    ) as mock_get:
        download_file(str(dest), input_url)
        assert mock_get.call_args[0][0] == expected_url


def test_download_file(tmp_path: Path, capsys: pytest.CaptureFixture) -> None:
    """Test download_file function."""
    url = "https://example.com/test.txt"
    test_content = b"test content"
    dest_path = tmp_path / "test.txt"

    with patch("requests.get", return_value=make_mock_response(test_content)):
        download_file(str(dest_path), url)
        assert dest_path.read_bytes() == test_content

    # Mock failed request
    with patch("requests.get", return_value=make_mock_response(b"Not found", 404)):
        download_file(str(dest_path), url)  # Should print error but not raise

    stdout, stderr = capsys.readouterr()
    assert f"Error downloading {url=}" in stdout
    assert stderr == ""


def test_download_file_keeps_existing_file_on_stream_error(
    tmp_path: Path, capsys: pytest.CaptureFixture
) -> None:
    """Failed streamed downloads should not corrupt existing cached files."""
    url = "https://example.com/test.txt"
    dest_path = tmp_path / "test.txt"
    dest_path.write_bytes(b"old content")

    response = make_mock_response(b"")

    def broken_iter_content(
        chunk_size: int = 8192, decode_unicode: bool = False  # noqa: ARG001
    ) -> list[bytes]:
        raise requests.ConnectionError("stream failed")

    response.iter_content = broken_iter_content

    with patch("requests.get", return_value=response):
        download_file(str(dest_path), url)

    stdout, stderr = capsys.readouterr()
    assert f"Error downloading {url=}" in stdout
    assert stderr == ""
    assert dest_path.read_bytes() == b"old content"
    assert not list(tmp_path.glob("test.txt.*.part"))


def test_download_file_cleans_partial_file_after_truncated_stream(
    tmp_path: Path, capsys: pytest.CaptureFixture
) -> None:
    """Mid-stream truncation should not replace an existing cached file."""
    url = "https://example.com/test.txt"
    dest_path = tmp_path / "test.txt"
    dest_path.write_bytes(b"old content")

    response = make_mock_response(b"")

    def truncated_iter_content(
        chunk_size: int = 8192, decode_unicode: bool = False  # noqa: ARG001
    ) -> Iterator[bytes]:
        yield b"partial content"
        raise requests.exceptions.ChunkedEncodingError("stream ended early")

    response.iter_content = truncated_iter_content

    with patch("requests.get", return_value=response):
        download_file(str(dest_path), url)

    stdout, stderr = capsys.readouterr()
    assert f"Error downloading {url=}" in stdout
    assert "stream ended early" in stdout
    assert stderr == ""
    assert dest_path.read_bytes() == b"old content"
    assert not list(tmp_path.glob("test.txt.*.part"))


def test_download_file_cleans_partial_file_after_write_error(
    tmp_path: Path, capsys: pytest.CaptureFixture
) -> None:
    """Write failures such as ENOSPC should clean up temporary files."""
    url = "https://example.com/test.txt"
    dest_path = tmp_path / "test.txt"
    dest_path.write_bytes(b"old content")
    real_open = open

    class FailingWriteFile:
        def __init__(self, file_path: str) -> None:
            self.file = real_open(file_path, mode="wb")

        def __enter__(self) -> Self:
            return self

        def __exit__(self, *_exc_info: object) -> None:
            self.file.close()

        def write(self, data: bytes) -> int:
            self.file.write(data[:1])
            raise OSError("No space left on device")

    def failing_open(
        file: str, mode: str = "r", *args: object, **kwargs: object
    ) -> object:
        if mode == "wb" and str(file).startswith(str(dest_path)):
            return FailingWriteFile(str(file))
        return real_open(file, mode, *args, **kwargs)

    with (
        patch("requests.get", return_value=make_mock_response(b"test content")),
        patch("builtins.open", side_effect=failing_open),
    ):
        download_file(str(dest_path), url)

    stdout, stderr = capsys.readouterr()
    assert f"Error downloading {url=}" in stdout
    assert "No space left on device" in stdout
    assert stderr == ""
    assert dest_path.read_bytes() == b"old content"
    assert not list(tmp_path.glob("test.txt.*.part"))


def test_download_file_removes_part_file_after_empty_response(
    tmp_path: Path, capsys: pytest.CaptureFixture
) -> None:
    """HTTP 200 with an empty body should not create a cached file."""
    url = "https://example.com/test.txt"
    dest_path = tmp_path / "test.txt"

    with patch("requests.get", return_value=make_mock_response(b"")):
        download_file(str(dest_path), url)

    stdout, stderr = capsys.readouterr()
    assert f"Error downloading {url=}" in stdout
    assert "Downloaded empty file" in stdout
    assert stderr == ""
    assert not dest_path.exists()
    assert not list(tmp_path.glob("test.txt.*.part"))


def test_maybe_auto_download_file_raises_when_request_fails_before_streaming(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Request failures before streaming should leave no cache artifacts."""
    url = "https://example.com/file.txt"
    abs_path = tmp_path / "test" / "file.txt"
    abs_path.parent.mkdir()

    monkeypatch.setenv("MBD_AUTO_DOWNLOAD_FILES", "true")
    with (
        patch("requests.get", side_effect=requests.ConnectionError("connect failed")),
        pytest.raises(FileNotFoundError, match="Download failed for 'test'"),
    ):
        maybe_auto_download_file(
            url, str(abs_path), label="test", raise_on_failure=True
        )

    assert not abs_path.exists()
    assert not list((tmp_path / "test").glob("file.txt.*.part"))


def test_download_file_ignores_empty_keepalive_chunks(tmp_path: Path) -> None:
    """Empty streaming chunks should be ignored before real content arrives."""
    url = "https://example.com/test.txt"
    dest_path = tmp_path / "test.txt"
    response = make_mock_response(b"")

    def iter_content_with_keepalive(
        chunk_size: int = 8192, decode_unicode: bool = False  # noqa: ARG001
    ) -> list[bytes]:
        return [b"", b"test content"]

    response.iter_content = iter_content_with_keepalive

    with patch("requests.get", return_value=response):
        download_file(str(dest_path), url)

    assert dest_path.read_bytes() == b"test content"
    assert not list(tmp_path.glob("test.txt.*.part"))


def test_download_file_uses_distinct_temp_files_for_concurrent_downloads(
    tmp_path: Path,
) -> None:
    """Concurrent downloads should not share a sidecar .part file."""
    url = "https://example.com/test.txt"
    dest_path = tmp_path / "test.txt"
    barrier = Barrier(2)
    responses: Queue[requests.Response] = Queue()
    replace_sources: list[str] = []
    errors: list[BaseException] = []
    real_replace = os.replace

    for content in (b"first content", b"second content"):
        response = make_mock_response(b"")

        def iter_content(
            chunk_size: int = 8192,  # noqa: ARG001
            decode_unicode: bool = False,  # noqa: ARG001
            *,
            chunk: bytes = content,
        ) -> list[bytes]:
            barrier.wait(timeout=5)
            return [chunk]

        response.iter_content = iter_content
        responses.put(response)

    def mock_get(*_args: object, **_kwargs: object) -> requests.Response:
        return responses.get_nowait()

    def recording_replace(src: str, dst: str) -> None:
        replace_sources.append(src)
        real_replace(src, dst)

    def worker() -> None:
        try:
            download_file(str(dest_path), url)
        except BaseException as exc:  # pragma: no cover - re-raised below
            errors.append(exc)

    threads = [Thread(target=worker) for _ in range(2)]
    with patch("requests.get", side_effect=mock_get), patch(
        "os.replace", side_effect=recording_replace
    ):
        for thread in threads:
            thread.start()
        for thread in threads:
            thread.join()

    assert errors == []
    assert len(set(replace_sources)) == 2
    assert all(source.endswith(".part") for source in replace_sources)
    assert f"{dest_path}.part" not in replace_sources
    assert dest_path.read_bytes() in {b"first content", b"second content"}
    assert not list(tmp_path.glob("test.txt.*.part"))


def test_maybe_auto_download_file_replaces_empty_cache(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture
) -> None:
    """Empty cache files should be treated as invalid and downloaded again."""
    url = "https://example.com/file.txt"
    abs_path = tmp_path / "test" / "file.txt"
    abs_path.parent.mkdir()
    abs_path.touch()

    monkeypatch.setenv("MBD_AUTO_DOWNLOAD_FILES", "true")
    with patch("requests.get", return_value=make_mock_response(b"test content")):
        maybe_auto_download_file(url, str(abs_path), label="test")

    stdout, stderr = capsys.readouterr()
    assert f"Downloading 'test' from {url!r}" in stdout
    assert stderr == ""
    assert abs_path.read_bytes() == b"test content"


def test_maybe_auto_download_file_raises_after_failed_download(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Download helpers should fail fast if no non-empty file was written."""
    url = "https://example.com/file.txt"
    abs_path = tmp_path / "test" / "file.txt"
    abs_path.parent.mkdir()

    monkeypatch.setenv("MBD_AUTO_DOWNLOAD_FILES", "true")
    with (
        patch("requests.get", return_value=make_mock_response(b"Not found", 404)),
        pytest.raises(FileNotFoundError, match="Download failed for 'test'"),
    ):
        maybe_auto_download_file(
            url, str(abs_path), label="test", raise_on_failure=True
        )

    assert not abs_path.exists()


def test_maybe_auto_download_file_soft_fails_by_default(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Optional callers should retain the legacy try-download-then-check pattern."""
    url = "https://example.com/file.txt"
    abs_path = tmp_path / "test" / "file.txt"
    abs_path.parent.mkdir()

    monkeypatch.setenv("MBD_AUTO_DOWNLOAD_FILES", "true")
    with patch("requests.get", return_value=make_mock_response(b"Not found", 404)):
        maybe_auto_download_file(url, str(abs_path), label="test")

    assert not abs_path.exists()


def test_maybe_auto_download_file(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture
) -> None:
    """Test auto-download behavior of maybe_auto_download_file function."""
    url = "https://example.com/file.txt"
    abs_path = f"{tmp_path}/test/file.txt"
    os.makedirs(os.path.dirname(abs_path), exist_ok=True)

    mock_response = make_mock_response(b"test content")

    # Test 1: Auto-download enabled (default)
    monkeypatch.setenv("MBD_AUTO_DOWNLOAD_FILES", "true")
    with patch("requests.get", return_value=mock_response):
        maybe_auto_download_file(url, abs_path, label="test")
        stdout, _ = capsys.readouterr()
        assert f"Downloading 'test' from {url!r}" in stdout
        assert os.path.isfile(abs_path)

    # Test 2: Auto-download disabled
    os.remove(abs_path)
    monkeypatch.setenv("MBD_AUTO_DOWNLOAD_FILES", "false")
    assert not os.path.isfile(abs_path)

    # Mock user input 'n' to skip download
    with (
        patch("requests.get", return_value=mock_response),
        patch("builtins.input", return_value="n"),
        patch("sys.stdin.isatty", return_value=True),  # force interactive mode
    ):
        maybe_auto_download_file(url, abs_path, label="test")
        assert not os.path.isfile(abs_path)

    # Test 3: Auto-download disabled but user confirms
    with (
        patch("requests.get", return_value=mock_response),
        patch("builtins.input", return_value="y"),
        patch("sys.stdin.isatty", return_value=True),  # force interactive mode
    ):
        maybe_auto_download_file(url, abs_path, label="test")
        stdout, _ = capsys.readouterr()
        assert f"Downloading 'test' from {url!r}" in stdout
        assert os.path.isfile(abs_path)

    # Test 4: File already exists (no download attempt)
    with patch("requests.get") as mock_get:
        maybe_auto_download_file(url, abs_path, label="test")
        mock_get.assert_not_called()

    # Test 5: Non-interactive session (auto-download)
    os.remove(abs_path)
    with (
        patch("requests.get", return_value=mock_response),
        patch("sys.stdin.isatty", return_value=False),
        patch("builtins.input", return_value="y"),  # Fallback input mock
    ):
        maybe_auto_download_file(url, abs_path, label="test")
        stdout, _ = capsys.readouterr()
        assert f"Downloading 'test' from {url!r}" in stdout
        assert os.path.isfile(abs_path)

    # Test 6: IPython session with auto-download disabled
    os.remove(abs_path)
    with (
        patch("requests.get", return_value=mock_response),
        patch("builtins.input", return_value="n"),
        patch("sys.stdin.isatty", return_value=True),  # force interactive mode
    ):
        maybe_auto_download_file(url, abs_path, label="test")
        assert not os.path.isfile(abs_path)
