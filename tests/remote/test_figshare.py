"""Unit tests for Figshare API helper functions."""

from pathlib import Path
from typing import Any
from unittest.mock import MagicMock, patch

import pytest
import requests

from matbench_discovery.remote import figshare


@pytest.mark.parametrize(
    "content,expected,binary",
    [
        (b'{"key": "value"}', {"key": "value"}, False),  # JSON
        (b"binary_data", b"binary_data", True),  # Binary
        (b"{invalid-json}", b"{invalid-json}", False),  # Invalid JSON
    ],
)
def test_make_request(
    content: bytes, expected: dict[str, str] | bytes, binary: bool
) -> None:
    """Test make_request with various response types."""
    mock_response = MagicMock(content=content)

    with patch("requests.request", return_value=mock_response):
        assert figshare.make_request("GET", "test_url", binary=binary) == expected


@pytest.mark.parametrize(
    "error_content,status_code",
    [
        (b'{"error": "Invalid token"}', 401),
        (b"", 500),
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


def test_article_exists_errors() -> None:
    """Test article_exists raises exceptions for non-404 errors."""
    mock_response = MagicMock()
    mock_response.raise_for_status.side_effect = requests.HTTPError(
        response=MagicMock(status_code=500)
    )

    err_msg = "article_url='https://api.figshare.com/v2/account/articles/12345'"
    with (
        patch("requests.request", return_value=mock_response),
        pytest.raises(requests.HTTPError, match=err_msg),
    ):
        figshare.article_exists(12345)


def test_create_article(capsys: pytest.CaptureFixture) -> None:
    """Article creation follows the returned location and reports its title."""
    metadata = {"title": "Test", "description": "Desc", "tags": ["tag1", "tag2"]}
    article_id = 12345
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

    md5, size = figshare.get_file_hash_and_size(str(test_file), chunk_size=5)
    assert size == expected_size
    assert md5 == expected_md5


@pytest.mark.parametrize(
    "file_parts,file_name",
    [
        (  # Default file name (from path)
            [{"partNo": 1, "startOffset": 0, "endOffset": 9}],
            "",
        ),
        (  # Custom name and multiple parts
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

    def mock_make_request(
        method: str, url: str, **kwargs: dict[str, str | int] | bytes | bool
    ) -> dict[str, Any] | None:
        """Mock request handler that returns appropriate response for each method."""
        if method == "GET" and url == "upload_url":
            return {"parts": file_parts}
        if method == "POST" and isinstance(data := kwargs.get("data"), dict):
            # Verify the file name in the POST request
            assert data["name"] == (file_name or test_file.name)
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
    with patch(
        "matbench_discovery.remote.figshare.make_request", return_value=files
    ) as request:
        assert figshare.list_article_files(12345) == files
    request.assert_called_once_with(
        "GET",
        f"{figshare.BASE_URL}/account/articles/12345/files?page_size=1000&page=1",
    )


def test_list_article_files_paginates() -> None:
    """Fetch subsequent pages until Figshare returns a short page."""
    first_page = [
        {"name": f"file-{file_idx}.txt", "id": file_idx} for file_idx in range(1000)
    ]
    with patch(
        "matbench_discovery.remote.figshare.make_request",
        side_effect=[first_page, DUMMY_FILES],
    ) as request:
        assert figshare.list_article_files(12345) == [*first_page, *DUMMY_FILES]

    base_url = f"{figshare.BASE_URL}/account/articles/12345/files?page_size=1000"
    assert [mock_call.args for mock_call in request.call_args_list] == [
        ("GET", f"{base_url}&page=1"),
        ("GET", f"{base_url}&page=2"),
    ]


def test_list_article_files_errors(capsys: pytest.CaptureFixture) -> None:
    """Test list_article_files HTTP error handling."""
    mock_response = MagicMock()
    mock_response.raise_for_status.side_effect = requests.HTTPError(
        response=MagicMock(status_code=404)
    )

    with patch("requests.request", return_value=mock_response):
        # should return empty list for 404 errors
        assert figshare.list_article_files(12345) == []
        assert figshare.get_existing_files(12345) == {}

    stdout, stderr = capsys.readouterr()
    assert stdout == stderr == ""

    mock_response.raise_for_status.side_effect = requests.HTTPError(
        response=MagicMock(status_code=500)
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
        (  # Multiple files and duplicate names use the final entry
            [
                {"name": "file1.txt", "id": 1, "computed_md5": "abc123"},
                {"name": "test.txt", "id": 2, "computed_md5": "def456"},
                {"name": "test.txt", "id": 3, "computed_md5": "ghi789"},
            ],
            {
                "file1.txt": {"id": 1, "computed_md5": "abc123"},
                "test.txt": {"id": 3, "computed_md5": "ghi789"},
            },
        ),
    ],
)
def test_get_existing_files(
    files: list[dict[str, Any]], expected: dict[str, dict[str, Any]]
) -> None:
    """Test get_existing_files with various file configurations."""
    original_files = [file.copy() for file in files]
    with patch("matbench_discovery.remote.figshare.make_request", return_value=files):
        assert figshare.get_existing_files(12345) == expected
    assert files == original_files


def test_file_exists_with_same_hash_reuses_inventory() -> None:
    """A supplied article inventory avoids another Figshare listing request."""
    existing_files = {"file.txt": {"id": 123, "computed_md5": "abc"}}
    with patch.object(figshare, "get_existing_files") as get_existing_files:
        result = figshare.file_exists_with_same_hash(
            12345, "file.txt", "abc", existing_files=existing_files
        )
    assert result == (True, 123)
    get_existing_files.assert_not_called()


def test_delete_file_success() -> None:
    """Successful deletion calls the expected article endpoint."""
    with patch(
        "matbench_discovery.remote.figshare.make_request", return_value=None
    ) as mock_make_request:
        assert figshare.delete_file(12345, 67890) is True
        mock_make_request.assert_called_once_with(
            "DELETE", f"{figshare.BASE_URL}/account/articles/12345/files/67890"
        )


def test_delete_file_error() -> None:
    """Test delete_file function with error during deletion."""
    with patch(
        "matbench_discovery.remote.figshare.make_request",
        side_effect=requests.RequestException("API Error"),
    ) as mock_request:
        assert figshare.delete_file(12345, 67890) is False
        mock_request.assert_called_once()


@pytest.mark.parametrize(
    "force_reupload,file_exists,delete_success,expected_upload,expected_delete",
    [
        (False, False, True, True, False),
        (False, True, True, False, False),
        (True, False, True, True, False),
        (True, True, True, True, True),
        (True, True, False, False, True),
    ],
)
def test_upload_file_if_needed(
    force_reupload: bool,
    file_exists: bool,
    delete_success: bool,
    expected_upload: bool,
    expected_delete: bool,
    tmp_path: Path,
) -> None:
    """Test upload_file_if_needed with various combinations of parameters."""
    test_file = tmp_path / "test_file.txt"
    test_file.write_text("test content")
    file_id = 67890 if file_exists else None

    mock_delete = MagicMock(return_value=delete_success)
    mock_upload = MagicMock(return_value=12345)
    existing_files: dict[str, dict[str, Any]] = {}

    with patch.multiple(
        "matbench_discovery.remote.figshare",
        get_file_hash_and_size=MagicMock(return_value=("test_hash", 12)),
        file_exists_with_same_hash=MagicMock(return_value=(file_exists, file_id)),
        delete_file=mock_delete,
        upload_file=mock_upload,
    ):
        result_id, was_uploaded = figshare.upload_file_if_needed(
            54321,
            str(test_file),
            file_name="test_file.txt",
            force_reupload=force_reupload,
            existing_files=existing_files,
        )

        # Verify expected behavior
        assert mock_delete.called == expected_delete
        assert mock_upload.called == expected_upload
        assert was_uploaded == expected_upload
        assert result_id == (12345 if expected_upload else file_id)
        expected_files = (
            {"test_file.txt": {"id": 12345, "computed_md5": "test_hash"}}
            if expected_upload
            else {}
        )
        assert existing_files == expected_files


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

    with patch(
        "matbench_discovery.remote.figshare.make_request",
        side_effect=None if success else requests.RequestException(err_msg),
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
    assert (
        figshare.find_similar_files(filename, existing_files, threshold)
        == expected_similar
    )
