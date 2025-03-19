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


@pytest.mark.parametrize(
    "status_code,expected_result",
    [
        (204, True),  # Success case
        (200, True),  # Alternative success case
    ],
)
def test_delete_file_success(status_code: int, expected_result: bool) -> None:
    """Test delete_file function with successful deletion."""
    _mock_response = MagicMock(status_code=status_code)

    with patch(
        "matbench_discovery.remote.figshare.make_request", return_value=None
    ) as mock_make_request:
        result = figshare.delete_file(12345, 67890)
        assert result == expected_result
        mock_make_request.assert_called_once_with(
            "DELETE", f"{figshare.BASE_URL}/account/articles/12345/files/67890"
        )


def test_delete_file_error() -> None:
    """Test delete_file function with error during deletion."""
    with patch(
        "matbench_discovery.remote.figshare.make_request",
        side_effect=Exception("API Error"),
    ) as mock_request:
        assert figshare.delete_file(12345, 67890) is False
        mock_request.assert_called_once()


@pytest.mark.parametrize(
    "force_reupload,file_exists,expected_upload,expected_delete",
    [
        (False, False, True, False),  # File doesn't exist, should upload
        (False, True, False, False),  # File exists, shouldn't upload or delete
        (True, False, True, False),  # File doesn't exist, should upload, no delete
        (
            True,
            True,
            True,
            True,
        ),  # File exists, force reupload, should delete and upload
    ],
)
def test_upload_file_if_needed(
    force_reupload: bool,
    file_exists: bool,
    expected_upload: bool,
    expected_delete: bool,
    tmp_path: Path,
) -> None:
    """Test upload_file_if_needed with various combinations of parameters."""
    test_file = tmp_path / "test_file.txt"
    test_file.write_text("test content")
    file_id = 67890 if file_exists else None

    mock_delete = MagicMock(return_value=True)
    mock_upload = MagicMock(return_value=12345)

    with patch.multiple(
        "matbench_discovery.remote.figshare",
        get_file_hash_and_size=MagicMock(return_value=("test_hash", 12)),
        file_exists_with_same_hash=MagicMock(return_value=(file_exists, file_id)),
        delete_file=mock_delete,
        upload_file=mock_upload,
    ):
        result_id, was_uploaded = figshare.upload_file_if_needed(
            54321, str(test_file), force_reupload=force_reupload
        )

        # Verify expected behavior
        assert mock_delete.called == expected_delete
        assert mock_upload.called == expected_upload
        assert was_uploaded == expected_upload
        assert result_id == (12345 if expected_upload else file_id)


def test_upload_file_if_needed_delete_failure(tmp_path: Path) -> None:
    """Test upload_file_if_needed when delete fails but force_reupload is True."""
    test_file = tmp_path / "test_file.txt"
    test_file.write_text("test content")

    mock_delete = MagicMock(return_value=False)
    mock_upload = MagicMock()

    with patch.multiple(
        "matbench_discovery.remote.figshare",
        get_file_hash_and_size=MagicMock(return_value=("test_hash", 12)),
        file_exists_with_same_hash=MagicMock(return_value=(True, 67890)),
        delete_file=mock_delete,
        upload_file=mock_upload,
    ):
        result_id, was_uploaded = figshare.upload_file_if_needed(
            54321, str(test_file), force_reupload=True
        )

        # Verify expected behavior
        assert mock_delete.called
        assert not mock_upload.called
        assert not was_uploaded
        assert result_id == 67890


@pytest.mark.parametrize(
    "success,verbose",
    [
        (True, True),  # Successful publish with verbose output
        (True, False),  # Successful publish without verbose output
        (False, True),  # Failed publish with verbose output
        (False, False),  # Failed publish without verbose output
    ],
)
def test_publish_article(
    success: bool, verbose: bool, capsys: pytest.CaptureFixture
) -> None:
    """Test publish_article function with different combinations of parameters."""
    article_id = 12345
    err_msg = "Test error"

    def mock_make_request_side_effect(*_args: Any, **_kwargs: Any) -> Any:
        """Either return successful result or raise exception based on success param."""
        if success:
            return  # Successful response for POST
        raise Exception(err_msg)  # noqa: TRY002

    with patch(
        "matbench_discovery.remote.figshare.make_request",
        side_effect=mock_make_request_side_effect,
    ):
        assert figshare.publish_article(article_id, verbose=verbose) is success

        stdout, _ = capsys.readouterr()
        if success and verbose:
            expected_url = f"{figshare.ARTICLE_URL_PREFIX}/{article_id}"
            assert (
                f"Successfully published article {article_id} at {expected_url}"
                in stdout
            )
        elif not success and verbose:
            assert f"Failed to publish article {article_id}: {err_msg}" in stdout
        else:
            assert stdout == ""


@pytest.mark.parametrize(
    "filename,existing_files,expected_similar,threshold",
    [
        # Empty files case
        ("models/model1/ver1/file-kappa-103.json.gz", {}, [], 0.7),
        # Different model family - no match
        (
            "models/model1/ver1/file-kappa-103.json.gz",
            {"models/model2/ver1/file-kappa-103.json.gz": {"id": 123}},
            [],
            0.7,
        ),
        # Different task type - no match
        (
            "models/model1/ver1/file-kappa-103.json.gz",
            {"models/model1/ver1/file-phonon-50.json.gz": {"id": 123}},
            [],
            0.7,
        ),
        # Different subfolder - no match
        (
            "models/model1/ver1/file-kappa-103.json.gz",
            {"models/model1/ver2/file-kappa-103.json.gz": {"id": 123}},
            [],
            0.7,
        ),
        # High similarity, same task - match
        (
            "models/model1/ver1/file-kappa-103-v1.json.gz",
            {"models/model1/ver1/file-kappa-103-v2.json.gz": {"id": 123}},
            [("models/model1/ver1/file-kappa-103-v2.json.gz", 123)],
            0.7,
        ),
        # Multiple matches
        (
            "models/model1/ver1/file-kappa-103.json.gz",
            {
                "models/model1/ver1/file1-kappa-103.json.gz": {"id": 123},
                "models/model1/ver1/file2-kappa-103.json.gz": {"id": 456},
                "models/model2/ver1/file-kappa-103.json.gz": {"id": 789},
            },
            [
                ("models/model1/ver1/file1-kappa-103.json.gz", 123),
                ("models/model1/ver1/file2-kappa-103.json.gz", 456),
            ],
            0.7,
        ),
        # Higher threshold excludes matches
        (
            "models/model1/ver1/file-kappa-103.json.gz",
            {"models/model1/ver1/similar-kappa-103.json.gz": {"id": 123}},
            [],
            0.95,
        ),
    ],
)
def test_find_similar_files(
    filename: str,
    existing_files: dict[str, dict[str, Any]],
    expected_similar: list[tuple[str, int]],
    threshold: float,
) -> None:
    """Test find_similar_files with various scenarios."""
    with (
        patch("difflib.SequenceMatcher") as mock_matcher,
        patch(
            "matbench_discovery.remote.figshare._extract_task_type"
        ) as mock_extract_task_type,
    ):
        # Configure similarity values
        mock_seq_matcher = MagicMock()
        mock_seq_matcher.ratio.return_value = 0.8  # Default high similarity

        # For high threshold test
        if threshold > 0.9:
            mock_seq_matcher.ratio.return_value = 0.9  # Not enough for 0.95

        mock_matcher.return_value = mock_seq_matcher

        # Configure task type extraction
        def extract_task(filename: str) -> str:
            for task in ["kappa", "phonon", "discovery", "geo"]:
                if task in filename:
                    return task
            return ""

        mock_extract_task_type.side_effect = extract_task

        # Run test
        result = sorted(
            figshare.find_similar_files(filename, existing_files, threshold),
            key=lambda x: x[0],
        )
        assert result == sorted(expected_similar, key=lambda x: x[0])
