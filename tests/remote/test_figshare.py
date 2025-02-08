"""Unit tests for Figshare API helper functions."""

from pathlib import Path
from typing import Any
from unittest.mock import MagicMock, patch

import pytest
import requests

import matbench_discovery.remote.figshare as figshare


@pytest.mark.parametrize(
    "content,expected,binary",
    [
        (b'{"key": "value"}', {"key": "value"}, False),  # JSON
        (b"binary_data", b"binary_data", True),  # Binary
        (b"plain text", b"plain text", False),  # Text
        (b"", b"", False),  # Empty
        (b"{invalid-json}", b"{invalid-json}", False),  # Invalid JSON
    ],
)
def test_make_request(content: bytes, expected: Any, binary: bool) -> None:
    """Test make_request with various response types."""
    mock_response = MagicMock(content=content)

    with patch("requests.request", return_value=mock_response):
        assert figshare.make_request("GET", "test_url", binary=binary) == expected


@pytest.mark.parametrize(
    "error_content,status_code",
    [
        (b"Not found", 404),
        (b'{"error": "Invalid token"}', 401),
        (b"", 500),
        (b'{"message": "Rate limit exceeded"}', 429),
    ],
)
def test_make_request_errors(error_content: bytes, status_code: int) -> None:
    """Test make_request error handling with various HTTP error codes."""
    mock_response = MagicMock(content=error_content, status_code=status_code)
    mock_response.raise_for_status.side_effect = requests.HTTPError(
        response=mock_response
    )

    err_msg = f"body={error_content.decode()}"
    with (
        patch("requests.request", return_value=mock_response),
        pytest.raises(requests.HTTPError, match=err_msg),
    ):
        figshare.make_request("GET", "test_url")


@pytest.mark.parametrize("status_code,expected", [(200, True), (404, False)])
def test_article_exists(status_code: int, expected: bool) -> None:
    """Test true and false case for article_exists."""
    mock_response = MagicMock(content=b'{"id": 12345}')
    if status_code != 200:
        mock_response.raise_for_status.side_effect = requests.HTTPError(
            response=MagicMock(status_code=status_code)
        )

    with patch("requests.request", return_value=mock_response):
        assert figshare.article_exists(12345) == expected


# Various error codes that should raise exceptions
@pytest.mark.parametrize("status_code", [401, 403, 500, 503])
def test_article_exists_errors(status_code: int) -> None:
    """Test article_exists raises exceptions for non-404 errors."""
    mock_response = MagicMock()
    mock_response.raise_for_status.side_effect = requests.HTTPError(
        response=MagicMock(status_code=status_code)
    )

    err_msg = "article_url='https://api.figshare.com/v2/account/articles/12345'"
    with (
        patch("requests.request", return_value=mock_response),
        pytest.raises(requests.HTTPError, match=err_msg),
    ):
        figshare.article_exists(12345)


@pytest.mark.parametrize(
    "metadata,article_id",
    [
        ({"title": "Test Article"}, 12345),
        ({"title": "Test", "description": "Desc"}, 67890),
        ({"title": "Test", "tags": ["tag1", "tag2"]}, 11111),
        ({"title": "Test", "categories": [1, 2], "keywords": ["key"]}, 22222),
    ],
)
def test_create_article_variants(
    metadata: dict[str, Any], article_id: int, capsys: pytest.CaptureFixture
) -> None:
    """Test article creation with different metadata combinations."""
    with patch(
        "matbench_discovery.remote.figshare.make_request",
        side_effect=[{"location": "loc"}, {"id": article_id}],
    ):
        assert figshare.create_article(metadata, verbose=True) == article_id

        stdout, stderr = capsys.readouterr()
        assert stdout == f"Created article: loc with title {metadata['title']}\n\n"
        assert stderr == ""


@pytest.mark.parametrize(
    "test_data,expected_size,expected_md5",
    [
        (b"test data", 9, "eb733a00c0c9d336e65691a37ab54293"),
        (b"", 0, "d41d8cd98f00b204e9800998ecf8427e"),  # Empty
        (b"hello world", 11, "5eb63bbbe01eeed093cb22bb8f5acdc3"),  # Regular
        (b"a" * 1000, 1000, "cabe45dcc9ae5b66ba86600cca6b8ba8"),  # Large
    ],
)
def test_get_file_hash_and_size_variants(
    test_data: bytes,
    expected_size: int,
    expected_md5: str,
    tmp_path: Path,
) -> None:
    """Test file hash and size calculation with different file contents."""
    test_file = tmp_path / "test_file"
    test_file.write_bytes(test_data)

    md5, size = figshare.get_file_hash_and_size(str(test_file))
    assert size == expected_size
    assert md5 == expected_md5


def test_get_file_hash_and_size_large_file(tmp_path: Path) -> None:
    """Test hash and size calculation for large files using chunked reading."""
    chunks = [b"chunk1", b"chunk2", b"chunk3"]
    test_file = tmp_path / "large_file"
    test_file.write_bytes(b"".join(chunks))

    md5, size = figshare.get_file_hash_and_size(str(test_file), chunk_size=5)
    assert size == sum(len(chunk) for chunk in chunks)
    assert md5 == "2aca0a9378723b1bed59975523ed50cd"


@pytest.mark.parametrize(
    "file_parts,file_name",
    [
        (  # Default file name (from path)
            [{"partNo": 1, "startOffset": 0, "endOffset": 9}],
            "",
        ),
        (  # Custom file name
            [{"partNo": 1, "startOffset": 0, "endOffset": 9}],
            "custom_name.txt",
        ),
        (  # Multiple parts with default name
            [
                {"partNo": 1, "startOffset": 0, "endOffset": 4},
                {"partNo": 2, "startOffset": 5, "endOffset": 9},
            ],
            "",
        ),
        (  # Multiple parts with custom name
            [
                {"partNo": 1, "startOffset": 0, "endOffset": 4},
                {"partNo": 2, "startOffset": 5, "endOffset": 9},
            ],
            "renamed.dat",
        ),
    ],
)
def test_upload_file_to_figshare_variants(
    file_parts: list[dict[str, int]],
    file_name: str,
    tmp_path: Path,
) -> None:
    """Test file upload with different file parts configurations and file names."""
    test_file = tmp_path / "upload_test_file"
    test_file.write_bytes(b"test data")

    mock_responses = {
        "POST": {"location": "file_location"},
        "GET": {"id": 67890, "upload_url": "upload_url", "parts": file_parts},
        "PUT": None,  # PUT requests return None on success
    }

    def mock_make_request(method: str, url: str, **kwargs: Any) -> Any:
        """Mock request handler that returns appropriate response for each method."""
        if method == "GET" and url == "upload_url":
            return {"parts": file_parts}
        if method == "POST" and "data" in kwargs:
            # Verify the file name in the POST request
            data = kwargs["data"]
            assert data["name"] == file_name or test_file.name
        return mock_responses[method]

    with (
        patch("matbench_discovery.remote.figshare.ROOT", str(tmp_path)),
        patch(
            "matbench_discovery.remote.figshare.make_request",
            side_effect=mock_make_request,
        ),
    ):
        assert figshare.upload_file(12345, str(test_file), file_name=file_name) == 67890


DUMMY_FILES = [
    {"name": "file1.txt", "id": 1, "md5": "abc123", "size": 100, "status": "ok"},
    {"name": "file2.txt", "id": 2, "md5": "def456", "size": 200, "status": "ok"},
]


@pytest.mark.parametrize("files", [[], DUMMY_FILES])  # Empty and non-empty
def test_list_article_files(files: list[dict[str, Any]]) -> None:
    """Test list_article_files with various file configurations."""
    with patch("matbench_discovery.remote.figshare.make_request", return_value=files):
        assert figshare.list_article_files(12345) == files


def test_list_article_files_errors(capsys: pytest.CaptureFixture) -> None:
    """Test list_article_files HTTP error handling."""
    mock_response = MagicMock()
    mock_response.raise_for_status.side_effect = requests.HTTPError(
        response=MagicMock(status_code=404)
    )

    with patch("requests.request", return_value=mock_response):
        # should return empty list for 404 errors
        assert figshare.list_article_files(12345) == []

    stdout, stderr = capsys.readouterr()
    assert stdout == stderr == ""

    for status_code in (401, 403, 500, 503):
        mock_response = MagicMock()
        mock_response.raise_for_status.side_effect = requests.HTTPError(
            response=MagicMock(status_code=status_code)
        )

        with (
            patch("requests.request", return_value=mock_response),
            pytest.raises(requests.HTTPError, match="\nbody="),
        ):
            figshare.list_article_files(12345)


@pytest.mark.parametrize(
    "files,expected",
    [
        ([], {}),  # Empty list case
        (  # Single file case
            [{"name": "test.txt", "id": 1, "computed_md5": "abc123"}],
            {"test.txt": {"id": 1, "computed_md5": "abc123"}},
        ),
        (  # Multiple files case
            [
                {"name": "file1.txt", "id": 1, "computed_md5": "abc123"},
                {"name": "file2.txt", "id": 2, "computed_md5": "def456"},
            ],
            {
                "file1.txt": {"id": 1, "computed_md5": "abc123"},
                "file2.txt": {"id": 2, "computed_md5": "def456"},
            },
        ),
        (  # Files with same name (should use last one)
            [
                {"name": "test.txt", "id": 1, "computed_md5": "abc123"},
                {"name": "test.txt", "id": 2, "computed_md5": "def456"},
            ],
            {"test.txt": {"id": 2, "computed_md5": "def456"}},
        ),
    ],
)
def test_get_existing_files(
    files: list[dict[str, Any]], expected: dict[str, dict[str, Any]]
) -> None:
    """Test get_existing_files with various file configurations."""
    with patch("matbench_discovery.remote.figshare.make_request", return_value=files):
        assert figshare.get_existing_files(12345) == expected


def test_get_existing_files_404() -> None:
    """Test get_existing_files returns empty dict for 404 errors."""
    mock_response = MagicMock()
    mock_response.raise_for_status.side_effect = requests.HTTPError(
        response=MagicMock(status_code=404)
    )

    with patch("requests.request", return_value=mock_response):
        assert figshare.get_existing_files(12345) == {}
