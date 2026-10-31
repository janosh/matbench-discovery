"""Test enums module."""

import hashlib
import os
from enum import auto
from pathlib import Path
from typing import Any
from unittest.mock import patch

import pytest
import requests

import scripts.generate_model_enum as model_enum_generator
from matbench_discovery import DATA_DIR, ROOT, file_digest
from matbench_discovery.enums import (
    ArchitectureType,
    DataFiles,
    Files,
    LabelEnum,
    MbdKey,
    Model,
    Open,
    Task,
    TestSubset,
    resolve_verified_file,
)
from matbench_discovery.remote import figshare
from tests.utils import make_mock_response, patch_model_metadata


@pytest.mark.parametrize("enum_cls", [ArchitectureType, MbdKey, Open, Task, TestSubset])
def test_label_enum_values_and_labels(enum_cls: type[LabelEnum]) -> None:
    """Every labeled-enum member has a non-empty, unique string value and label."""
    values = [member.value for member in enum_cls]
    labels = [member.label for member in enum_cls]
    for value, label in zip(values, labels, strict=True):
        assert isinstance(value, str), f"{value=}"
        assert isinstance(label, str), f"{label=}"
        assert value != ""
        assert label != "", f"{value=}"
    assert len(values) == len(set(values)), "Values must be unique"
    assert len(labels) == len(set(labels)), "Labels must be unique"


# HTML tags are part of the label, so the exact string is checked
E_FORM_DFT_LABEL = (
    "DFT E<sub>form</sub> "
    "<span style='font-size: 0.8em; font-weight: lighter;'>(eV/atom)</span>"
)


@pytest.mark.parametrize(
    ("member", "value", "label"),
    [
        (MbdKey.e_form_dft, "e_form_per_atom_mp2020_corrected", E_FORM_DFT_LABEL),
        (
            MbdKey.each_true,
            "e_above_hull_mp2020_corrected_ppd_mp",
            "E<sub>MP hull dist</sub>",
        ),
        (Task.S2E, "S2E", "structure to energy"),
        (Task.S2EFS, "S2EFS", "structure to energy, force, stress"),
        (Task.S2EFSM, "S2EFSM", "structure to energy, force, stress, magmoms"),
        (
            Task.IS2RE_SR,
            "IS2RE-SR",
            "initial structure to relaxed energy with structure relaxation",
        ),
        (Open.OSOD, "OSOD", "open source, open data"),
        (Open.CSCD, "CSCD", "closed source, closed data"),
        (TestSubset.uniq_protos, "unique_prototypes", "Unique Structure Prototypes"),
        (TestSubset.full_test_set, "full_test_set", "Full Test Set"),
    ],
    ids=lambda arg: arg.name if isinstance(arg, LabelEnum) else None,
)
def test_label_enum_member_spot_checks(
    member: LabelEnum, value: str, label: str
) -> None:
    """Spot-check a few member values and labels the site and metrics code rely on."""
    assert member == value
    assert member.label == label


def test_files_enum() -> None:
    """Test error handling in Files enum."""
    assert Files.base_dir == DATA_DIR

    class SubFiles(Files, base_dir="foo"):
        test_file = auto(), "test/file.txt"

        @property
        def url(self) -> str:
            """URL associated with the file."""
            return "https://example.com/file.txt"

        @property
        def label(self) -> str:
            """Label associated with the file."""
            return "test"

    assert SubFiles.base_dir == "foo"

    test_file = SubFiles.test_file
    assert repr(test_file) == "SubFiles.test_file"
    assert str(test_file) == "test_file"


@pytest.mark.parametrize("enum_cls", [Model, DataFiles])
def test_files_members_are_distinct(enum_cls: type[Files]) -> None:
    """Files members must compare by value instead of collapsing into one element.

    Regression: Files.__new__ called str.__new__(cls) without the value, leaving every
    member's string content "" so the inherited str.__eq__/__hash__ made all members
    compare equal and hash to 0 (sets/dicts/`in` collapsed to one element).
    """
    members = list(enum_cls)
    # a set of members must keep every distinct member, not collapse to one
    assert len({*members}) == len(members)

    first, second = members[0], members[1]
    assert first != second
    # each member equals (and hashes like) its own value string
    assert first == first.value
    assert hash(first) == hash(first.value)
    assert repr(first) == f"{enum_cls.__name__}.{first.name}"
    assert str(first) == first.name == first.value


def test_diatomics_reference_md5_matches_bundled_copy() -> None:
    """The registered diatomics reference md5 matches the site's bundled copy."""
    diatomics = DataFiles.diatomics_dft_reference
    assert diatomics.yaml[diatomics.name]["md5"] == file_digest(
        f"{ROOT}/site/src/lib/diatomics-dft.json.gz", "md5"
    )


def test_data_files_members_match_yaml_registry() -> None:
    """Every public data-files YAML entry has exactly one enum member."""
    yaml_data = DataFiles.mp_energies.yaml
    public_keys = {key for key in yaml_data if not key.startswith("_")}
    assert set(DataFiles.__members__) == public_keys
    for data_file in DataFiles:
        assert data_file.rel_path == yaml_data[data_file.name]["path"]


def test_data_files_path_raises_when_md5_download_fails(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """DataFiles.path surfaces failed checksum verification to callers."""
    data_file = DataFiles.mp_energies
    expected_md5 = "0" * 32
    monkeypatch.setattr(DataFiles, "_base_dir", str(tmp_path))
    monkeypatch.setitem(
        data_file.__dict__,
        "yaml",
        {
            data_file.name: {
                "description": "test data file",
                "md5": expected_md5,
                "path": data_file.rel_path,
                "url": "https://example.com/file.csv.gz",
            }
        },
    )

    with (
        patch("requests.get", return_value=make_mock_response(b"bad data")),
        pytest.raises(ValueError, match=f"MD5 mismatch.*expected {expected_md5}"),
    ):
        _ = data_file.path

    assert not os.path.isfile(f"{tmp_path}/{data_file.rel_path}")


@pytest.mark.parametrize("repair_succeeds", [True, False])
def test_resolve_verified_file_rechecks_stale_cache(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, repair_succeeds: bool
) -> None:
    """A stale artifact is returned only after a successful verified repair."""
    from matbench_discovery import enums as enums_mod

    abs_path = f"{tmp_path}/artifact.csv"
    Path(abs_path).write_bytes(b"stale")
    expected_md5 = hashlib.md5(b"fresh", usedforsecurity=False).hexdigest()
    actual_md5 = hashlib.md5(b"stale", usedforsecurity=False).hexdigest()

    def repair_file(_url: str, repair_path: str, **_kwargs: object) -> None:
        """Replace stale bytes only in the successful repair case."""
        if repair_succeeds:
            Path(repair_path).write_bytes(b"fresh")

    monkeypatch.setattr(enums_mod, "maybe_auto_download_file", repair_file)

    def resolve() -> str:
        """Resolve the test artifact with its expected checksum."""
        return resolve_verified_file(
            abs_path,
            url="https://example.com/artifact.csv",
            label="artifact",
            md5=expected_md5,
        )

    if repair_succeeds:
        assert resolve() == abs_path
        assert Path(abs_path).read_bytes() == b"fresh"
    else:
        with pytest.raises(
            ValueError, match=f"expected MD5 {expected_md5}, got {actual_md5}"
        ):
            resolve()


@pytest.mark.parametrize("subdirectory", ["", "nested/", "../", "../mace-mp-0-extra/"])
def test_model_md_path_returns_path_for_dict_md(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, subdirectory: str
) -> None:
    """Resolve nested local predictions but reject paths outside the model directory."""
    from matbench_discovery import enums as enums_mod

    model = Model.mace_mp_0
    rel_path = f"models/mace/mace-mp-0/{subdirectory}2026-09-06-md-metrics.csv.gz"
    patch_model_metadata(
        monkeypatch,
        model,
        {
            "model_name": "MACE",
            "model_key": "mace-mp-0",
            "metrics": {"md": {"pred_file": {"name": rel_path}}},
        },
    )
    download_calls: list[tuple[object, ...]] = []
    monkeypatch.setattr(
        enums_mod,
        "maybe_auto_download_file",
        lambda *args, **_kwargs: download_calls.append(args),
    )
    monkeypatch.setattr(enums_mod, "ROOT", str(tmp_path))
    expected_path = f"{tmp_path}/{rel_path}"
    Path(expected_path).parent.mkdir(parents=True, exist_ok=True)
    Path(expected_path).write_bytes(b"md")
    if subdirectory.startswith("../"):
        with pytest.raises(ValueError, match="must live under"):
            _ = model.md_path
    else:
        assert model.md_path == expected_path
    assert download_calls == []


@pytest.mark.parametrize("harmonic_rerun", [True, False])
def test_kappa_103_phonon_path_layouts(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, *, harmonic_rerun: bool
) -> None:
    """The phonon sidecar resolves from kappa_103.harmonic or kappa_103 itself."""
    model = Model.mace_mp_0
    sub_dir = "harmonic/" if harmonic_rerun else ""
    rel_path = f"models/mace/mace-mp-0/{sub_dir}2026-09-05-phonons-kappa-103.json.gz"
    phonon_ref = {"phonon_file": {"name": rel_path}}
    kappa_103 = {"harmonic": phonon_ref} if harmonic_rerun else phonon_ref
    metrics = {"phonons": {"kappa_103": kappa_103}}
    metadata = {"model_name": "MACE", "model_key": "mace-mp-0", "metrics": metrics}
    patch_model_metadata(monkeypatch, model, metadata)
    monkeypatch.setattr("matbench_discovery.enums.ROOT", str(tmp_path))
    Path(f"{tmp_path}/{rel_path}").parent.mkdir(parents=True)
    Path(f"{tmp_path}/{rel_path}").write_bytes(b"phonons")
    assert model.kappa_103_phonon_path == f"{tmp_path}/{rel_path}"


@pytest.mark.parametrize("data_file", DataFiles)
@pytest.mark.network
def test_data_files_enum_urls(
    data_file: DataFiles, url_session: requests.Session
) -> None:
    """Test that each URL in data-files.yml is a reachable Figshare download URL."""
    name, url = data_file.name, data_file.url
    assert "figshare.com/files/" in url, (
        f"URL for {name} is not a Figshare download URL: {url}"
    )
    check_url(url_session, url)


@pytest.fixture(scope="session")
def figshare_data_file_md5s(url_session: requests.Session) -> dict[str, str]:
    """Map Figshare file id -> computed_md5 for the repo's data-files article."""
    article_id = figshare.ARTICLE_IDS["data_files"]
    response = url_session.get(
        f"https://api.figshare.com/v2/articles/{article_id}/files?page_size=500",
        timeout=TIMEOUT,
    )
    response.raise_for_status()
    return {str(file["id"]): file["computed_md5"] for file in response.json()}


@pytest.mark.parametrize("data_file", DataFiles)
@pytest.mark.network
def test_data_files_md5_matches_figshare(
    data_file: DataFiles, figshare_data_file_md5s: dict[str, str]
) -> None:
    """Each declared md5 in data-files.yml must match Figshare's computed_md5 for the
    current artifact, so registry drift (like the stale checksums behind #357) fails
    CI instead of making the file un-downloadable (download_file discards mismatches).
    """
    file_id = data_file.url.rsplit("/", maxsplit=1)[-1]
    if (computed_md5 := figshare_data_file_md5s.get(file_id)) is None:
        # file lives in an external Figshare article this repo doesn't control
        # (e.g. mp_trj_original in the original MPtrj article), can't enforce md5
        pytest.skip(f"{data_file.name} file id {file_id} not in data-files article")
    declared_md5 = data_file.yaml[data_file.name].get("md5")
    assert declared_md5 == computed_md5, (
        f"data-files.yml md5 for {data_file.name} ({declared_md5}) does not match "
        f"Figshare computed_md5 ({computed_md5}) for file id {file_id}. Update the "
        "registry entry (or re-run scripts/upload_data_files_to_figshare.py)."
    )


@pytest.mark.parametrize(
    ("metrics", "expected"),
    [(None, {}), ({"md": {}}, {"md": {}})],
)
def test_model_metadata_accessors(
    monkeypatch: pytest.MonkeyPatch,
    metrics: dict[str, Any] | None,
    expected: dict[str, Any],
) -> None:
    """label/key/metrics read YAML metadata; bare `metrics:` reads as no metrics."""
    model = Model.alignn
    metadata = {"model_name": "Fake Net", "model_key": "fake-net", "metrics": metrics}
    patch_model_metadata(monkeypatch, model, metadata)
    assert (model.label, model.key, model.metrics) == ("Fake Net", "fake-net", expected)
    assert model.md_path is None  # no md pred_file means no MD predictions


def test_generated_model_enum_is_current() -> None:
    """Committed Model members match the YAML-driven generator."""
    with open(f"{ROOT}/matbench_discovery/enums.py", encoding="utf-8") as file:
        source = file.read()
    assert model_enum_generator.generate_source(source) == source


def test_model_enum_generator_rejects_global_key_collisions(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Canonical keys are unique across active and aborted YAMLs."""
    for lifecycle in ("active", "aborted"):
        path = tmp_path / "models" / lifecycle / "model.yml"
        path.parent.mkdir(parents=True)
        path.write_text(f"model_key: model-a\nlifecycle: {lifecycle}\n")
    monkeypatch.setattr(model_enum_generator, "ROOT", str(tmp_path))
    source = (
        f"class Model:\n{model_enum_generator.BEGIN_MARKER}\n"
        f"{model_enum_generator.END_MARKER}\n"
    )
    with pytest.raises(ValueError, match="collides"):
        model_enum_generator.generate_source(source)


@pytest.mark.parametrize(
    "input_value, expected_model",
    [
        # Exact matches
        ("mace_mp_0", Model.mace_mp_0),
        ("eqv2_s_dens_mp", Model.eqv2_s_dens_mp),
        # Dash conversion
        ("mace-mp-0", Model.mace_mp_0),
        ("eqv2-s-dens-mp", Model.eqv2_s_dens_mp),
        ("chgnet-0.3.0", Model.chgnet_0_3_0),
        ("cgcnn-p", Model.cgcnn_p),
        # Case insensitive
        ("MACE-MP-0", Model.mace_mp_0),
        ("EQV2-S-DENS-MP", Model.eqv2_s_dens_mp),
        # Mixed separators
        ("mace-mp_0", Model.mace_mp_0),
        ("mace_mp-0", Model.mace_mp_0),
        (123, None),
        (None, None),
        ([], None),
        ({}, None),
        ("nonexistent", None),
        ("mace-mp-1", None),
        ("eqv2-s-dens", None),
        ("", None),
        ("   ", None),
    ],
)
def test_model_missing(input_value: object, expected_model: Model | None) -> None:
    """Model._missing_ normalizes valid references and rejects invalid ones."""
    assert Model._missing_(input_value) is expected_model


def test_model_md_path_passes_huggingface_token(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """MD prediction downloads from gated HuggingFace repos use bearer auth."""
    url = "https://huggingface.co/org/repo/resolve/main/md.csv.gz"
    rel_path = "models/alignn/alignn/nested/2026-09-06-md-metrics.csv.gz"
    abs_path = f"{tmp_path}/{rel_path}"
    monkeypatch.setenv("HF_TOKEN", "hf_secret")
    monkeypatch.setattr("matbench_discovery.enums.ROOT", str(tmp_path))
    patch_model_metadata(
        monkeypatch,
        Model.alignn,
        {
            "model_name": "Gated MD",
            "model_key": "alignn",
            "metrics": {"md": {"pred_file": {"name": rel_path, "url": url}}},
        },
    )

    with patch("requests.get", return_value=make_mock_response(b"md")) as mock_get:
        assert Model.alignn.md_path == abs_path

    assert os.path.isfile(abs_path)
    assert mock_get.call_args.args == (url,)
    assert mock_get.call_args.kwargs["headers"] == {"Authorization": "Bearer hf_secret"}


TIMEOUT = 30
URL_WORKERS = 16  # HTTP checks wait on the network, not the runner's CPU cores.
VALID_URL_STATUSES = {200, 202, 403, 429, 500, 502, 503, 504}


def check_url(session: requests.Session, url: str) -> None:
    """Assert a URL resolves - the adapter retries transients, so only 404s fail."""
    # 200: OK, 202: figshare async, 403: restricted, 429/5xx: rate limited or server
    # outage after retries, which says nothing about whether the URL is valid
    status = session.head(url, allow_redirects=True, timeout=TIMEOUT).status_code
    assert status in VALID_URL_STATUSES, f"unexpected {status=} for {url}"


@pytest.mark.network
def test_model_prediction_urls(url_session: requests.Session) -> None:
    """Test that all model prediction file URLs are valid."""
    import concurrent.futures

    tasks: dict[str, str] = {}
    for model in Model.active():
        if model.name == Model.mace_mpa_0.name:
            continue
        tasks[model.metadata["pr_url"]] = model.name
        metrics = model.metrics
        if not metrics:
            continue
        # walk nested metrics dicts to collect every FileRef URL with its dotted path
        stack: list[tuple[str, dict[str, Any]]] = [("", metrics)]
        while stack:
            parent_key, dct = stack.pop()
            for key, val in dct.items():
                key_path = f"{parent_key}.{key}" if parent_key else key
                if key == "url" and isinstance(val, str):
                    tasks[val] = f"{model.name}.{key_path}"
                elif isinstance(val, dict):
                    stack.append((key_path, val))

    errors: list[Exception] = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=URL_WORKERS) as executor:
        futures = {
            executor.submit(check_url, url_session, url): (url, desc)
            for url, desc in tasks.items()
        }
        for future in concurrent.futures.as_completed(futures):
            url, desc = futures[future]
            try:
                future.result()
            except (AssertionError, requests.RequestException) as exc:
                exc.add_note(f"Failed to validate\n{url}\n{desc}")
                errors.append(exc)

    if errors:
        msg = f"{len(errors)}/{len(tasks)} URLs failed validation"
        raise ExceptionGroup(msg, errors)
