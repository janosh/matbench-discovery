"""Unit tests for Figshare API helper functions."""

import hashlib
from functools import partial
from pathlib import Path
from types import ModuleType
from typing import Any
from unittest.mock import MagicMock, patch

import pytest
import requests
from ruamel.yaml import YAML

import scripts.upload_data_files_to_figshare as upload_data
import scripts.upload_model_preds_to_figshare as upload_models
from matbench_discovery.enums import Model
from matbench_discovery.remote import figshare

ARTICLE_URL = f"{figshare.ARTICLE_URL_PREFIX}/12345"
KAPPA_FILE = "models/model1/ver1/file-kappa-103.json.gz"


@pytest.mark.parametrize(
    "uploader", [upload_data, upload_models], ids=["data", "models"]
)
def test_upload_preserves_unicode_yaml(
    uploader: ModuleType, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Archive updates preserve UTF-8 metadata under a Windows default encoding."""
    yaml_path = tmp_path / "metadata.yml"
    original = "description: κ — Å — 声子\nmetrics: {}\n"
    yaml_path.write_text(original, encoding="utf-8")
    monkeypatch.setattr(
        uploader, "open", partial(open, encoding="cp1252"), raising=False
    )
    monkeypatch.setattr(uploader, "round_trip_yaml", YAML())
    monkeypatch.setattr(figshare, "article_exists", lambda _article_id: True)
    monkeypatch.setattr(figshare, "get_existing_files", lambda _article_id: {})
    monkeypatch.setattr(figshare, "list_article_files", lambda _article_id: [])
    monkeypatch.setattr(figshare, "make_request", MagicMock(side_effect=AssertionError))
    if uploader is upload_data:
        assert (
            uploader.main(str(yaml_path), 123, {"keywords": [], "urls": {}}, files=[])
            == 0
        )
    else:
        monkeypatch.setattr(Model, "yaml_path", property(lambda _model: str(yaml_path)))
        uploader.update_one_modeling_task_article(
            "phonons", [Model.mace_mp_0], modeling_tasks={"phonons": {}}
        )
    assert yaml_path.read_text(encoding="utf-8") == original


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

    md5, size = figshare.get_file_hash_and_size(str(test_file))
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

    # computed_md5 satisfies the post-upload verification in upload_file
    file_md5 = hashlib.md5(test_file.read_bytes(), usedforsecurity=False).hexdigest()
    mock_responses = {
        "POST": {"location": "file_location"},
        "GET": {
            "id": 67890,
            "upload_url": "upload_url",
            "parts": file_parts,
            "computed_md5": file_md5,
        },
        "PUT": None,  # PUT requests return None on success
    }

    def mock_make_request(
        method: str, url: str, **kwargs: dict[str, str | int] | bytes | bool
    ) -> dict[str, Any] | None:
        """Mock request handler that returns appropriate response for each method."""
        if method == "GET" and url == "upload_url":
            return {"parts": file_parts}
        if method == "POST" and isinstance(data := kwargs.get("data"), dict):
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


@pytest.mark.parametrize(
    ("remote_responses", "expected_error"),
    [
        ([{"computed_md5": "abc123"}], None),  # matches on first poll
        ([{}, {"computed_md5": "abc123"}], None),  # md5 appears once hashing finishes
        ([{"computed_md5": "corrupt"}], "Figshare stored"),
        ([{}, {}, {}], "reported no checksum"),  # never hashed
    ],
)
def test_verify_upload(
    remote_responses: list[dict[str, str]], expected_error: str | None
) -> None:
    """Post-upload verification accepts a matching checksum and rejects the rest."""
    with (
        patch(
            "matbench_discovery.remote.figshare.make_request",
            side_effect=remote_responses,
        ),
        patch("matbench_discovery.remote.figshare.time.sleep"),  # skip the backoff
    ):
        if expected_error:
            with pytest.raises(ValueError, match=expected_error):
                figshare.verify_upload(1, 2, "abc123", attempts=len(remote_responses))
        else:
            figshare.verify_upload(1, 2, "abc123", attempts=len(remote_responses))


DUMMY_FILES = [
    {"name": "file1.txt", "id": 1, "md5": "abc123", "size": 100, "status": "ok"},
    {"name": "file2.txt", "id": 2, "md5": "def456", "size": 200, "status": "ok"},
]


@pytest.mark.parametrize(
    "pages",
    [
        [[]],
        [DUMMY_FILES],
        [[{"name": f"file-{idx}.txt", "id": idx} for idx in range(1000)], DUMMY_FILES],
    ],
    ids=["empty", "single-page", "paginated"],
)
def test_list_article_files(pages: list[list[dict[str, Any]]]) -> None:
    """Collect all files and request consecutive pages until a short page arrives."""
    with patch(
        "matbench_discovery.remote.figshare.make_request",
        side_effect=pages,
    ) as request:
        assert figshare.list_article_files(12345) == [
            file for page in pages for file in page
        ]

    base_url = f"{figshare.BASE_URL}/account/articles/12345/files?page_size=1000"
    assert [mock_call.args for mock_call in request.call_args_list] == [
        ("GET", f"{base_url}&page={page}") for page in range(1, len(pages) + 1)
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


@pytest.mark.parametrize(
    "side_effect,expected",
    [(None, True), (requests.RequestException("API Error"), False)],
)
def test_delete_file(side_effect: Exception | None, expected: bool) -> None:
    """Deletion hits the file endpoint and reports API errors as failure."""
    with patch(
        "matbench_discovery.remote.figshare.make_request", side_effect=side_effect
    ) as mock_request:
        assert figshare.delete_file(12345, 67890) is expected
    mock_request.assert_called_once_with(
        "DELETE", f"{figshare.BASE_URL}/account/articles/12345/files/67890"
    )


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
    "success,verbose,expected_stdout",
    [
        (True, True, f"Successfully published article 12345 at {ARTICLE_URL}"),
        (False, True, "Failed to publish article 12345: Test error"),
        (True, False, ""),  # verbose=False stays silent either way
        (False, False, ""),
    ],
)
def test_publish_article(
    success: bool, verbose: bool, expected_stdout: str, capsys: pytest.CaptureFixture
) -> None:
    """Publishing returns its outcome and reports it only when verbose."""
    with patch(
        "matbench_discovery.remote.figshare.make_request",
        side_effect=None if success else requests.RequestException("Test error"),
    ):
        assert figshare.publish_article(12345, verbose=verbose) is success

    stdout, _ = capsys.readouterr()
    assert expected_stdout in stdout
    assert verbose or stdout == ""


@pytest.mark.parametrize(
    "filename,existing_files,expected_similar,threshold",
    [
        (KAPPA_FILE, {}, [], 0.7),  # no candidates
        # a match needs same family, same subfolder, same task and high similarity
        (KAPPA_FILE, {"models/model2/ver1/f-kappa-103.json.gz": {"id": 1}}, [], 0.7),
        (KAPPA_FILE, {"models/model1/ver1/f-phonon-50.json.gz": {"id": 1}}, [], 0.7),
        (KAPPA_FILE, {"models/model1/ver2/f-kappa-103.json.gz": {"id": 1}}, [], 0.7),
        (
            KAPPA_FILE,
            {f"models/model1/ver1/nested/{KAPPA_FILE.rsplit('/', 1)[-1]}": {"id": 1}},
            [],
            0.7,
        ),
        (
            "models/model1/ver1/file-kappa-103-v1.json.gz",
            {"models/model1/ver1/file-kappa-103-v2.json.gz": {"id": 123}},
            [("models/model1/ver1/file-kappa-103-v2.json.gz", 123)],
            0.7,
        ),
        (  # every qualifying candidate is returned
            KAPPA_FILE,
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
        (  # same candidate as the match case above, excluded by a stricter threshold
            KAPPA_FILE,
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


@pytest.mark.parametrize("suffix_idx", range(6))
def test_similar_files_keep_artifact_roles_and_directories(suffix_idx: int) -> None:
    """Match dates without conflating complementary artifacts, settings or folders."""
    suffixes = [
        "phonons-kappa-103.json.gz",
        "phonons-kappa-103-forces.json.gz",
        "phonons-kappa-103-phonons.json.gz",
        "phonons-kappa-103-run-info.json",
        "geo-opt-symprec=1e-2-moyo=0.12.0.csv.gz",
        "geo-opt-symprec=1e-5-moyo=0.12.0.csv.gz",
    ]
    directory = "models/mace/mace-mp-0/harmonic"
    candidates = {
        f"{folder}/2026-09-05-{suffix}": {"id": idx}
        for folder in (directory, f"{directory}/other")
        for idx, suffix in enumerate(suffixes)
    }
    suffix = suffixes[suffix_idx]
    assert figshare.find_similar_files(
        f"{directory}/2026-09-06-{suffix}", candidates
    ) == [(f"{directory}/2026-09-05-{suffix}", suffix_idx)]
