"""Tests for archiving model prediction files to Figshare without network access."""

from contextlib import nullcontext
from pathlib import Path
from typing import Literal
from unittest.mock import PropertyMock

import pytest

import scripts.upload_model_preds_to_figshare as upload
from matbench_discovery import repo_relative_path
from matbench_discovery.data import round_trip_yaml
from matbench_discovery.enums import Model

MODEL = Model.mace_mp_0
PRED_REL = "models/arch/model/geo-opt.jsonl.gz"
ANALYSIS_REL = "models/arch/model/geo-opt-symprec=1e-2.csv.gz"
REPRO_REL = "models/arch/model/old-run.jsonl.gz"  # reproducibility.runs[0].file
TASKS = {"geo_opt": {"label": "Geo-opt", "description": "Relaxed structures"}}


@pytest.mark.parametrize(
    ("items", "expected"),
    [
        (["a", "c"], ["a", "c"]),  # explicit inclusions
        (["!b"], ["a", "c"]),  # exclusions start from all items
        (["a", "b", "!b"], ["a"]),  # exclusions beat inclusions
        (["!d"], None),  # typos (also behind !) are usage errors
    ],
)
def test_process_exclusion_prefixes(
    items: list[str], expected: list[str] | None
) -> None:
    """Explicit items win over the full list; ``!`` prefixes always exclude."""
    with nullcontext() if expected else pytest.raises(SystemExit):
        assert upload.process_exclusion_prefixes(items, ["a", "b", "c"]) == expected


@pytest.mark.parametrize(
    ("file_type", "same_hash"),
    [("all", False), ("analysis", False), ("all", True)],
)
def test_upload_writes_file_refs_back_to_yaml(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    file_type: Literal["all", "analysis", "pred"],
    same_hash: bool,
) -> None:
    """Uploads (or same-hash matches) write url+size+md5 only for the filtered refs,
    including reproducibility run files nested in lists (file_type=all only).
    """
    for rel_path in (PRED_REL, ANALYSIS_REL, REPRO_REL):
        (tmp_path / rel_path).parent.mkdir(parents=True, exist_ok=True)
        (tmp_path / rel_path).write_bytes(b"artifact")
    metrics = {
        "pred_file": {"name": PRED_REL},
        "symprec=1e-2": {"analysis_file": {"name": ANALYSIS_REL}, "rmsd": 0.1},
        "reproducibility": {"runs": [{"label": "old", "file": {"name": REPRO_REL}}]},
    }
    yaml_path = f"{tmp_path}/models/arch/model.yml"
    with open(yaml_path, mode="w", encoding="utf-8") as file:
        round_trip_yaml.dump({"metrics": {"geo_opt": metrics}}, file)
    monkeypatch.setattr(Model, "yaml_path", PropertyMock(return_value=yaml_path))
    monkeypatch.setattr(
        upload, "resolve_artifact_path", lambda _yaml, rel: f"{tmp_path}/{rel}"
    )
    monkeypatch.setattr(
        upload,
        "repo_relative_path",
        lambda path: repo_relative_path(path, str(tmp_path)),
    )
    figshare = upload.figshare
    monkeypatch.setattr(figshare, "article_exists", lambda _id: True)
    monkeypatch.setattr(
        figshare, "get_file_hash_and_size", lambda path: (f"md5-{path[-12:]}", 8)
    )
    monkeypatch.setattr(figshare, "find_similar_files", lambda *_args: [])
    uploads: list[str] = []
    published: list[int] = []
    existing_ids = {PRED_REL: 5, ANALYSIS_REL: 6, REPRO_REL: 7}
    monkeypatch.setattr(figshare, "get_existing_files", lambda _id: {})
    monkeypatch.setattr(
        figshare,
        "file_exists_with_same_hash",
        lambda _id, name, *_args, **_kwargs: (
            (True, existing_ids[name]) if same_hash else (False, None)
        ),
    )

    def fake_upload(_article_id: int, _path: str, **kwargs: object) -> tuple[int, bool]:
        """Record the upload name and return a fresh file ID."""
        uploads.append(str(kwargs["file_name"]))
        return 10 + len(uploads), True

    monkeypatch.setattr(figshare, "upload_file_if_needed", fake_upload)
    monkeypatch.setattr(figshare, "publish_article", published.append)
    upload.update_one_modeling_task_article(
        "geo_opt", [MODEL], modeling_tasks=TASKS, file_type=file_type
    )

    with open(yaml_path, encoding="utf-8") as file:
        written = round_trip_yaml.load(file)["metrics"]["geo_opt"]
    refs = {
        PRED_REL: written["pred_file"],
        ANALYSIS_REL: written["symprec=1e-2"]["analysis_file"],
        REPRO_REL: written["reproducibility"]["runs"][0]["file"],
    }
    kinds = ((PRED_REL, "pred"), (ANALYSIS_REL, "analysis"), (REPRO_REL, "run"))
    selected = [rel_path for rel_path, kind in kinds if file_type in ("all", kind)]
    assert uploads == ([] if same_hash else selected)
    # unchanged files only get their existing URL recorded, so nothing to publish
    assert published == (
        [] if same_hash else [figshare.ARTICLE_IDS["model_preds_geo_opt"]]
    )
    assert written["symprec=1e-2"]["rmsd"] == 0.1  # sibling metrics survive
    for idx, rel_path in enumerate(selected, start=1):
        file_id = existing_ids[rel_path] if same_hash else 10 + idx
        assert refs[rel_path] == {
            "name": rel_path,
            "url": f"{figshare.DOWNLOAD_URL_PREFIX}/{file_id}",
            "size": 8,
            "md5": f"md5-{rel_path[-12:]}",
        }
    for rel_path in set(refs) - set(selected):
        assert refs[rel_path] == {"name": rel_path}  # filtered out: untouched
