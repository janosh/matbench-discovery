"""Tests for syncing the Data Files Figshare article without network access."""

import hashlib
from pathlib import Path

import pytest

import scripts.upload_data_files_to_figshare as upload_data
from matbench_discovery.data import round_trip_yaml
from matbench_discovery.enums import DataFiles

DATA_FILE = DataFiles.wbm_summary
PYPROJECT = {"keywords": ["mlip"], "urls": {"Repo": "https://example.com"}}


@pytest.mark.parametrize(
    ("remote_matches", "dry_run", "expect_upload"),
    [
        (True, False, False),  # remote file named by rel_path is recognized as same
        (False, False, True),
        (False, True, False),  # dry run plans but never uploads/publishes/writes
    ],
    ids=["unchanged", "upload", "dry-run"],
)
def test_sync_matches_remote_by_rel_path(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
    remote_matches: bool,
    dry_run: bool,
    expect_upload: bool,
) -> None:
    """Remote files match by upload name; uploads, publishes, writes follow."""
    local_file = tmp_path / DATA_FILE.rel_path
    local_file.parent.mkdir(parents=True)
    local_file.write_bytes(b"wbm summary")
    local_md5 = hashlib.md5(b"wbm summary").hexdigest()  # noqa: S324
    yaml_path = tmp_path / "data-files.yml"
    registry = {DATA_FILE.name: {"description": "WBM summary"}}
    write_registry(yaml_path, registry)

    calls: list[tuple[str, object]] = []
    remote_file = {
        "name": DATA_FILE.rel_path,  # uploads keep their DATA_DIR-relative name
        "computed_md5": local_md5 if remote_matches else "0" * 32,
        "id": 7,
    }
    monkeypatch.setattr(upload_data, "DATA_DIR", str(tmp_path))
    figshare = upload_data.figshare
    monkeypatch.setattr(figshare, "article_exists", lambda _id: True)
    monkeypatch.setattr(figshare, "list_article_files", lambda _id: [remote_file])

    def fake_upload(article_id: int, file_path: str, file_name: str = "") -> int:
        """Record the upload name and return a new file ID."""
        calls.append(("upload", (article_id, file_path, file_name)))
        return 8

    monkeypatch.setattr(figshare, "upload_file", fake_upload)
    monkeypatch.setattr(
        figshare, "publish_article", lambda _article_id: calls.append(("publish", 1))
    )

    exit_code = upload_data.main(
        str(yaml_path), 1, PYPROJECT, files=[DATA_FILE], dry_run=dry_run
    )

    assert exit_code == 0
    if expect_upload:
        assert calls == [
            ("upload", (1, f"{tmp_path}/{DATA_FILE.rel_path}", DATA_FILE.rel_path)),
            ("publish", 1),
        ]
    else:
        assert calls == []
    with open(yaml_path, encoding="utf-8") as file:
        written = round_trip_yaml.load(file)
    if dry_run:
        assert written == registry
        assert "[dry-run] Would upload" in capsys.readouterr().out
    else:
        file_id = 8 if expect_upload else 7
        assert written[DATA_FILE.name]["md5"] == local_md5
        assert written[DATA_FILE.name]["url"].endswith(f"/{file_id}")


def write_registry(yaml_path: Path, registry: dict[str, dict[str, str]]) -> None:
    """Write a minimal data-files.yml registry."""
    with open(yaml_path, mode="w", encoding="utf-8") as file:
        round_trip_yaml.dump(registry, file)


@pytest.mark.parametrize(
    ("registry", "write_file", "error", "match"),
    [
        # an explicitly requested file that is missing locally must not be skipped
        ({DATA_FILE.name: {"description": "WBM summary"}}, False, FileNotFoundError,
         "Requested --files missing locally"),
        # new entries need a real description, no "Description needed" placeholder
        ({DATA_FILE.name: {"description": ""}}, True, ValueError, "Add a description"),
    ],
)  # fmt: skip
def test_sync_rejects_missing_or_undescribed_files(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    registry: dict[str, dict[str, str]],
    write_file: bool,
    error: type[Exception],
    match: str,
) -> None:
    """Invalid requests fail before any Figshare call."""
    if write_file:
        local_file = tmp_path / DATA_FILE.rel_path
        local_file.parent.mkdir(parents=True)
        local_file.write_bytes(b"wbm summary")
    yaml_path = tmp_path / "data-files.yml"
    write_registry(yaml_path, registry)
    monkeypatch.setattr(upload_data, "DATA_DIR", str(tmp_path))
    # reaching Figshare would raise AttributeError instead of the validation error
    monkeypatch.delattr(upload_data.figshare, "article_exists")
    with pytest.raises(error, match=match):
        upload_data.main(str(yaml_path), 1, PYPROJECT, files=[DATA_FILE])


@pytest.mark.parametrize(
    ("argv", "expected"),
    [
        ([], {"files": None, "dry_run": False}),
        (["-n", "--files", DATA_FILE.name], {"files": [DATA_FILE], "dry_run": True}),
        (["--dry"], "unrecognized arguments: --dry"),  # no prefix abbreviations
        (["--files", "not_a_data_file"], "invalid DataFiles value"),
    ],
)
def test_cli_parses_or_rejects_args(
    capsys: pytest.CaptureFixture[str], argv: list[str], expected: dict | str
) -> None:
    """Valid args parse; typos exit with usage errors before any upload."""
    parser = upload_data.build_parser()
    if isinstance(expected, str):
        with pytest.raises(SystemExit, match="2"):
            parser.parse_args(argv)
        assert expected in capsys.readouterr().err
    else:
        assert vars(parser.parse_args(argv)) == expected
