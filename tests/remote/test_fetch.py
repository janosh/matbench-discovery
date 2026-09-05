"""Tests for remote file download helpers."""

import hashlib
import os
from collections.abc import Iterator
from contextlib import nullcontext
from pathlib import Path
from unittest.mock import patch

import pytest
import requests

from matbench_discovery.remote.fetch import download_file, maybe_auto_download_file
from tests.utils import make_mock_response


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
        (
            "https://example.com/figshare.com/files/123",
            "https://example.com/figshare.com/files/123",
        ),
        (
            "https://figshare.com.example.org/files/123",
            "https://figshare.com.example.org/files/123",
        ),
        (
            "https://figshare.com@other.example/files/123",
            "https://figshare.com@other.example/files/123",
        ),
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


def test_download_file(tmp_path: Path) -> None:
    """Downloads verify checksums and propagate HTTP failures, preserving cache."""
    url = "https://example.com/test.txt"
    test_content = b"test content"
    dest_path = tmp_path / "test.txt"

    with patch("requests.get", return_value=make_mock_response(test_content)):
        download_file(str(dest_path), url, md5=hashlib.md5(test_content).hexdigest())  # noqa: S324
        assert dest_path.read_bytes() == test_content
        assert not list(tmp_path.glob(".*.part"))

    with (
        patch("requests.get", return_value=make_mock_response(b"Not found", 404)),
        pytest.raises(requests.HTTPError, match="404") as error,
    ):
        download_file(str(dest_path), url)
    assert url in " ".join(error.value.__notes__)
    assert dest_path.read_bytes() == test_content
    with pytest.raises(ValueError, match="Invalid IPv6 URL"):
        download_file(str(dest_path), "https://[invalid")
    assert not list(tmp_path.glob(".*.part"))


def test_download_file_md5_mismatch_discards_download(tmp_path: Path) -> None:
    """A download failing md5 verification is discarded; cached file survives."""
    url = "https://example.com/test.txt"
    dest_path = tmp_path / "test.txt"
    dest_path.write_bytes(b"old content")
    expected_md5 = "0" * 32

    with (
        patch("requests.get", return_value=make_mock_response(b"new content")),
        pytest.raises(ValueError, match=f"MD5 mismatch.*expected {expected_md5}"),
    ):
        download_file(str(dest_path), url, md5=expected_md5)
    assert dest_path.read_bytes() == b"old content"
    assert not list(tmp_path.glob(".*.part"))


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
@pytest.mark.parametrize("status_code", [200, 401, 403, 404])
@pytest.mark.parametrize(
    ("url", "use_token"),
    [
        ("https://huggingface.co/org/repo/resolve/main/file.csv.gz", True),
        ("https://cdn-lfs.huggingface.co/file", True),
        ("https://huggingface.co.example.org/file", False),
        ("https://example.org/huggingface.co/file", False),
        ("https://example.org/file?source=huggingface.co", False),
        ("https://huggingface.co@other.example/file", False),
        ("http://huggingface.co/file", False),
    ],
)
def test_download_file_adds_huggingface_token(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    token_env: str,
    status_code: int,
    url: str,
    use_token: bool,
) -> None:
    """Only HTTPS HuggingFace hosts receive automatically supplied bearer tokens."""
    dest_path = tmp_path / "file.csv.gz"
    monkeypatch.delenv("HF_TOKEN", raising=False)
    monkeypatch.delenv("HUGGING_FACE_HUB_TOKEN", raising=False)
    monkeypatch.setenv(token_env, "hf_test")

    with (
        patch(
            "requests.get", return_value=make_mock_response(b"test", status_code)
        ) as mock_get,
        (
            pytest.raises(requests.HTTPError) if status_code != 200 else nullcontext()
        ) as error,
    ):
        download_file(str(dest_path), url)

    if error is None:
        assert dest_path.read_bytes() == b"test"
    else:
        notes = " ".join(error.value.__notes__)
        assert ("HF_TOKEN" in notes) is (use_token and status_code in (401, 403))
        assert not dest_path.is_file()
    assert mock_get.call_args.kwargs["headers"] == (
        {"Authorization": "Bearer hf_test"} if use_token else None
    )


def test_download_file_keeps_completed_part_file_on_replace_error(
    tmp_path: Path,
) -> None:
    """Completed downloads should survive final replace failures."""
    url = "https://example.com/test.txt"
    dest_path = tmp_path / "test.txt"
    dest_path.write_bytes(b"old content")

    with (
        patch("requests.get", return_value=make_mock_response(b"new content")),
        patch("os.replace", side_effect=PermissionError("replace denied")),
        pytest.raises(PermissionError, match="replace denied") as error,
    ):
        download_file(str(dest_path), url)
    assert dest_path.read_bytes() == b"old content"
    part_paths = list(tmp_path.glob(".*.part"))
    assert len(part_paths) == 1
    assert part_paths[0].read_bytes() == b"new content"
    assert repr(str(part_paths[0])) in " ".join(error.value.__notes__)


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
    with (
        patch("requests.get", return_value=response),
        remove_ctx,
        pytest.raises(requests.ConnectionError, match="stream failed") as error,
    ):
        download_file(str(dest_path), url)
    notes = " ".join(error.value.__notes__)
    assert url in notes
    assert dest_path.read_bytes() == b"old content"
    if remove_error:
        assert "Failed to remove partial download" in notes
        assert "cannot remove part file" in notes
        assert len(list(tmp_path.glob(".*.part"))) == 1
    else:
        assert not list(tmp_path.glob(".*.part"))


def test_overlapping_downloads_use_independent_temporary_files(tmp_path: Path) -> None:
    """Interleaved downloads cannot truncate or replace each other's temporary file."""
    destination = tmp_path / "model.ckpt"
    response = make_mock_response(b"")

    def interleaved_content(**_kwargs: object) -> Iterator[bytes]:
        """Complete a second download while the first is still streaming."""
        yield b"first "
        download_file(str(destination), "https://example.com/second")
        assert destination.read_bytes() == b"second download"
        yield b"download"

    response.iter_content = interleaved_content  # ty: ignore[invalid-assignment]
    with patch(
        "requests.get", side_effect=[response, make_mock_response(b"second download")]
    ):
        download_file(str(destination), "https://example.com/first")
    assert destination.read_bytes() == b"first download"
    assert list(tmp_path.iterdir()) == [destination]


@pytest.mark.parametrize("auth_header", ["Authorization", "authorization"])
def test_explicit_authorization_is_preserved(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    auth_header: str,
) -> None:
    """Explicit authentication takes precedence regardless of header casing."""
    monkeypatch.setenv("HF_TOKEN", "hf_environment")
    headers = {auth_header: "Bearer explicit"}
    with patch("requests.get", return_value=make_mock_response(b"data")) as mock_get:
        download_file(
            str(tmp_path / "model"), "https://huggingface.co/model", headers=headers
        )
    assert mock_get.call_args.kwargs["headers"] == headers


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


@pytest.mark.parametrize(
    ("md5", "expected_content", "should_download"),
    [
        (None, b"cached", False),
        (hashlib.md5(b"cached").hexdigest(), b"cached", False),  # noqa: S324
        (hashlib.md5(b"fresh").hexdigest(), b"fresh", True),  # noqa: S324
    ],
    ids=["unchecked", "checksum_match", "checksum_mismatch"],
)
def test_maybe_auto_download_file_validates_existing_file(
    tmp_path: Path,
    md5: str | None,
    expected_content: bytes,
    should_download: bool,
) -> None:
    """Existing files are reused only when their optional checksum matches."""
    url = "https://example.com/file.txt"
    abs_path = f"{tmp_path}/test/file.txt"
    os.makedirs(os.path.dirname(abs_path), exist_ok=True)
    Path(abs_path).write_bytes(b"cached")

    with patch("requests.get", return_value=make_mock_response(b"fresh")) as mock_get:
        maybe_auto_download_file(url, abs_path, label="test", md5=md5)
    assert Path(abs_path).read_bytes() == expected_content
    assert mock_get.called is should_download


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
