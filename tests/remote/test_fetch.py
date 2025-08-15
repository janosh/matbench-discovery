import os
from pathlib import Path
from unittest.mock import patch

import pytest
import requests

from matbench_discovery.remote.fetch import maybe_auto_download_file


def test_download_file(tmp_path: Path, capsys: pytest.CaptureFixture) -> None:
    """Test download_file function."""

    from matbench_discovery.remote.fetch import download_file

    url = "https://example.com/test.txt"
    test_content = b"test content"
    dest_path = tmp_path / "test.txt"

    # Mock successful request
    mock_response = requests.Response()
    mock_response.status_code = 200
    mock_response._content = test_content  # noqa: SLF001

    with patch("requests.get", return_value=mock_response):
        download_file(str(dest_path), url)
        assert dest_path.read_bytes() == test_content

    # Mock failed request
    mock_response = requests.Response()
    mock_response.status_code = 404
    mock_response._content = b"Not found"  # noqa: SLF001

    with patch("requests.get", return_value=mock_response):
        download_file(str(dest_path), url)  # Should print error but not raise

    stdout, stderr = capsys.readouterr()
    assert f"Error downloading {url=}" in stdout
    assert stderr == ""


def test_maybe_auto_download_file(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture
) -> None:
    """Test auto-download behavior of maybe_auto_download_file function."""
    url = "https://example.com/file.txt"
    abs_path = f"{tmp_path}/test/file.txt"
    os.makedirs(os.path.dirname(abs_path), exist_ok=True)

    # Mock successful request
    mock_response = requests.Response()
    mock_response.status_code = 200
    mock_response._content = b"test content"  # noqa: SLF001

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
