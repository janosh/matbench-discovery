"""Unit tests for Figshare API helper functions."""

from pathlib import Path
from typing import Any
from unittest.mock import MagicMock, patch

import pytest
import requests

from matbench_discovery.figshare import (
    create_article,
    get_file_hash_and_size,
    make_request,
    upload_file_to_figshare,
)


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
        assert make_request("GET", "test_url", binary=binary) == expected


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

    with (
        patch("requests.request", return_value=mock_response),
        pytest.raises(requests.HTTPError) as exc_info,
    ):
        make_request("GET", "test_url")

    assert f"body={error_content.decode()}" in str(exc_info.value.__notes__)


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
        "matbench_discovery.figshare.make_request",
        side_effect=[{"location": "loc"}, {"id": article_id}],
    ):
        assert create_article(metadata, verbose=True) == article_id

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

    md5, size = get_file_hash_and_size(str(test_file))
    assert size == expected_size
    assert md5 == expected_md5


def test_get_file_hash_and_size_large_file(tmp_path: Path) -> None:
    """Test hash and size calculation for large files using chunked reading."""
    chunks = [b"chunk1", b"chunk2", b"chunk3"]
    test_file = tmp_path / "large_file"
    test_file.write_bytes(b"".join(chunks))

    md5, size = get_file_hash_and_size(str(test_file), chunk_size=5)
    assert size == sum(len(chunk) for chunk in chunks)
    assert md5 == "2aca0a9378723b1bed59975523ed50cd"


@pytest.mark.parametrize(
    "file_parts",
    [
        [{"partNo": 1, "startOffset": 0, "endOffset": 9}],  # Single part
        [  # Multiple parts
            {"partNo": 1, "startOffset": 0, "endOffset": 4},
            {"partNo": 2, "startOffset": 5, "endOffset": 9},
        ],
        [{"partNo": 1, "startOffset": 0, "endOffset": 0}],  # Empty file
        [  # Many small parts
            {"partNo": idx + 1, "startOffset": idx * 2, "endOffset": (idx + 1) * 2 - 1}
            for idx in range(50)
        ],
    ],
)
def test_upload_file_to_figshare_variants(
    file_parts: list[dict[str, int]], tmp_path: Path
) -> None:
    """Test file upload with different file parts configurations."""
    test_file = tmp_path / "upload_test_file"
    test_file.write_bytes(b"test data")

    mock_responses = {
        "POST": {"location": "file_location"},
        "GET": {"id": 67890, "upload_url": "upload_url", "parts": file_parts},
        "PUT": None,  # PUT requests return None on success
    }

    def mock_make_request(method: str, url: str, **_kwargs: Any) -> Any:
        """Mock request handler that returns appropriate response for each method."""
        if method == "GET" and url == "upload_url":
            return {"parts": file_parts}
        return mock_responses[method]

    with patch(
        "matbench_discovery.figshare.make_request", side_effect=mock_make_request
    ):
        assert upload_file_to_figshare(12345, str(test_file)) == 67890
