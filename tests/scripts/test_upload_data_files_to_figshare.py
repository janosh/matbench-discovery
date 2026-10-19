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
    yaml_path.write_text("{}\n", encoding="utf-8")

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
            ("upload", (1, str(local_file), DATA_FILE.rel_path)),
            ("publish", 1),
        ]
    else:
        assert calls == []
    with open(yaml_path, encoding="utf-8") as file:
        written = round_trip_yaml.load(file)
    if dry_run:
        assert written == {}
        assert "[dry-run] Would upload" in capsys.readouterr().out
    else:
        file_id = 8 if expect_upload else 7
        assert written[DATA_FILE.name]["md5"] == local_md5
        assert written[DATA_FILE.name]["url"].endswith(f"/{file_id}")


@pytest.mark.parametrize(
    ("argv", "dry_run"),
    [([], False), (["-n"], True), (["--dry-run"], True), (["--dryrun"], None)],
)
def test_cli_rejects_unknown_flags(argv: list[str], dry_run: bool | None) -> None:
    """Typos in flags exit instead of falling through to a real upload."""
    parser = upload_data.build_parser()
    if dry_run is None:
        with pytest.raises(SystemExit, match="2"):
            parser.parse_args(argv)
    else:
        assert parser.parse_args(argv).dry_run is dry_run
