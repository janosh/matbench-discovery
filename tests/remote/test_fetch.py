"""Tests for remote file download helpers."""

import os
from collections.abc import Iterator
from contextlib import nullcontext
from pathlib import Path
from unittest.mock import patch

import pytest
import requests

from matbench_discovery.remote.fetch import download_file, maybe_auto_download_file


def make_mock_response(content: bytes, status_code: int = 200) -> requests.Response:
    """Create a mock requests.Response with given content and status code."""
    response = requests.Response()
    response.status_code = status_code
    response._content = content  # noqa: SLF001
    response.iter_content = (  # ty: ignore[invalid-assignment]
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


def test_download_file_current_directory(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Relative filenames in the current directory should not call makedirs('')."""
    url = "https://example.com/test.txt"
    dest_path = tmp_path / "test.txt"

    with (
        patch("requests.get", return_value=make_mock_response(b"test content")),
        patch("os.makedirs") as mock_makedirs,
    ):
        monkeypatch.chdir(tmp_path)
        download_file(dest_path.name, url)

    mock_makedirs.assert_not_called()
    assert dest_path.read_bytes() == b"test content"


@pytest.mark.parametrize("token_env", ["HF_TOKEN", "HUGGING_FACE_HUB_TOKEN"])
def test_download_file_adds_huggingface_token(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, token_env: str
) -> None:
    """HuggingFace downloads use bearer auth when a token env var is present."""
    url = "https://huggingface.co/org/repo/resolve/main/file.csv.gz"
    dest_path = tmp_path / "file.csv.gz"
    monkeypatch.delenv("HF_TOKEN", raising=False)
    monkeypatch.delenv("HUGGING_FACE_HUB_TOKEN", raising=False)
    monkeypatch.setenv(token_env, "hf_test")

    with patch("requests.get", return_value=make_mock_response(b"test")) as mock_get:
        download_file(str(dest_path), url)

    assert dest_path.read_bytes() == b"test"
    assert mock_get.call_args.kwargs["headers"] == {"Authorization": "Bearer hf_test"}


def test_download_file_keeps_completed_part_file_on_replace_error(
    tmp_path: Path, capsys: pytest.CaptureFixture
) -> None:
    """Completed downloads should survive final replace failures."""
    url = "https://example.com/test.txt"
    dest_path = tmp_path / "test.txt"
    part_path = Path(f"{dest_path}.part")
    dest_path.write_bytes(b"old content")

    with (
        patch("requests.get", return_value=make_mock_response(b"new content")),
        patch("os.replace", side_effect=PermissionError("replace denied")),
    ):
        download_file(str(dest_path), url)

    stdout, stderr = capsys.readouterr()
    assert f"Error downloading {url=}" in stdout
    assert "replace denied" in stdout
    assert stderr == ""
    assert dest_path.read_bytes() == b"old content"
    assert part_path.read_bytes() == b"new content"


@pytest.mark.parametrize(
    "stream_chunks, remove_error",
    [
        ((), None),
        ((b"partial content",), None),
        ((b"partial content",), PermissionError("cannot remove part file")),
    ],
)
def test_download_file_keeps_existing_file_on_stream_error(
    stream_chunks: tuple[bytes, ...],
    remove_error: OSError | None,
    tmp_path: Path,
    capsys: pytest.CaptureFixture,
) -> None:
    """Failed streamed downloads should not corrupt existing cached files."""
    url = "https://example.com/test.txt"
    dest_path = tmp_path / "test.txt"
    dest_path.write_bytes(b"old content")
    response = make_mock_response(b"")

    def broken_iter_content(**_kwargs: object) -> Iterator[bytes]:
        yield from stream_chunks
        raise requests.ConnectionError("stream failed")

    response.iter_content = broken_iter_content  # ty: ignore[invalid-assignment]
    remove_ctx = (
        patch("os.remove", side_effect=remove_error) if remove_error else nullcontext()
    )
    with patch("requests.get", return_value=response), remove_ctx:
        download_file(str(dest_path), url)

    stdout, stderr = capsys.readouterr()
    assert f"Error downloading {url=}" in stdout
    assert "stream failed" in stdout
    assert stderr == ""
    assert dest_path.read_bytes() == b"old content"
    if remove_error:
        assert "Failed to remove partial download" in stdout
        assert "cannot remove part file" in stdout
        assert os.path.isfile(f"{dest_path}.part")
    else:
        assert not os.path.isfile(f"{dest_path}.part")


@pytest.mark.parametrize(
    ("auto_download", "stdin_isatty", "answer", "should_download"),
    [
        (None, True, "n", True),
        ("true", True, "n", True),
        ("false", True, "n", False),
        ("false", True, "y", True),
        ("false", False, "n", True),
    ],
    ids=["auto_unset", "auto_enabled", "declined", "confirmed", "non_interactive"],
)
def test_maybe_auto_download_file_prompt_modes(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture,
    auto_download: str | None,
    stdin_isatty: bool,
    answer: str,
    *,
    should_download: bool,
) -> None:
    """maybe_auto_download_file honors env, prompt and non-interactive defaults."""
    url = "https://example.com/file.txt"
    abs_path = f"{tmp_path}/test/file.txt"
    os.makedirs(os.path.dirname(abs_path), exist_ok=True)
    mock_response = make_mock_response(b"test content")

    if auto_download is None:
        monkeypatch.delenv("MBD_AUTO_DOWNLOAD_FILES", raising=False)
    else:
        monkeypatch.setenv("MBD_AUTO_DOWNLOAD_FILES", auto_download)
    with (
        patch("requests.get", return_value=mock_response),
        patch("builtins.input", return_value=answer),
        patch("sys.stdin.isatty", return_value=stdin_isatty),
    ):
        maybe_auto_download_file(url, abs_path, label="test")

    stdout, stderr = capsys.readouterr()
    assert stderr == ""
    if should_download:
        assert f"Downloading 'test' from {url!r}" in stdout
        assert os.path.isfile(abs_path)
        assert Path(abs_path).read_bytes() == b"test content"
    else:
        assert stdout == ""
        assert not os.path.isfile(abs_path)


def test_maybe_auto_download_file_skips_existing_file(tmp_path: Path) -> None:
    """maybe_auto_download_file does not request an already cached file."""
    url = "https://example.com/file.txt"
    abs_path = f"{tmp_path}/test/file.txt"
    os.makedirs(os.path.dirname(abs_path), exist_ok=True)
    Path(abs_path).write_bytes(b"cached")

    with patch("requests.get") as mock_get:
        maybe_auto_download_file(url, abs_path, label="test")
    mock_get.assert_not_called()


def test_maybe_auto_download_file_forwards_huggingface_token(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Auto-download forwards HuggingFace bearer auth to the underlying request."""
    url = "https://huggingface.co/org/repo/resolve/main/file.csv.gz"
    abs_path = f"{tmp_path}/test/file.csv.gz"

    monkeypatch.setenv("MBD_AUTO_DOWNLOAD_FILES", "true")
    monkeypatch.delenv("HUGGING_FACE_HUB_TOKEN", raising=False)
    monkeypatch.setenv("HF_TOKEN", "hf_test")
    with patch("requests.get", return_value=make_mock_response(b"test")) as mock_get:
        maybe_auto_download_file(url, abs_path, label="test")

    assert os.path.isfile(abs_path)
    assert mock_get.call_args.kwargs["headers"] == {"Authorization": "Bearer hf_test"}
