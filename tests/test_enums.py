"""Test enums module."""

import os
from enum import auto
from pathlib import Path
from typing import Any
from unittest.mock import patch

import pytest
import requests
import requests.adapters

from matbench_discovery import DATA_DIR
from matbench_discovery.enums import (
    DataFiles,
    Files,
    MbdKey,
    Model,
    ModelType,
    Open,
    Task,
    TestSubset,
)
from matbench_discovery.remote import figshare
from matbench_discovery.remote.fetch import maybe_auto_download_file


def make_mock_response(content: bytes) -> requests.Response:
    """Create a successful streaming response with fixed byte content."""
    response = requests.Response()
    response.status_code = 200
    response._content = content  # noqa: SLF001
    response.iter_content = lambda chunk_size=8192: [content]  # ty: ignore[invalid-assignment] # noqa: ARG005
    return response


def test_mbd_key() -> None:
    """Test MbdKey enum."""
    # Test basic enum functionality
    assert MbdKey.e_form_dft == "e_form_per_atom_mp2020_corrected"
    # HTML tags are part of the label, so we need to test the exact string
    assert MbdKey.e_form_dft.label == (
        "DFT E<sub>form</sub> "
        "<span style='font-size: 0.8em; font-weight: lighter;'>(eV/atom)</span>"
    )
    assert MbdKey.each_true == "e_above_hull_mp2020_corrected_ppd_mp"
    assert MbdKey.each_true.label == "E<sub>MP hull dist</sub>"
    assert MbdKey.init_protostructure_spglib == (
        "protostructure_spglib_initial_structure"
    )
    assert MbdKey.protostructure_spglib == "protostructure_spglib"

    # Test that all keys have values and labels
    for key in MbdKey:
        assert isinstance(key.value, str), f"{key=}"
        assert isinstance(key.label, str), f"{key=}"
        assert key.value != ""
        assert key.label != ""

    # Test uniqueness of values and labels
    values = [key.value for key in MbdKey]
    labels = [key.label for key in MbdKey]
    assert len(values) == len(set(values)), "Values must be unique"
    assert len(labels) == len(set(labels)), "Labels must be unique"


def test_task() -> None:
    """Test Task enum."""
    # Test basic enum functionality
    assert Task.S2E == "S2E"
    assert Task.S2E.label == "structure to energy"
    assert Task.S2EFS == "S2EFS"
    assert Task.S2EFS.label == "structure to energy, force, stress"

    # Test that all tasks have values and labels
    for task in Task:
        assert isinstance(task.value, str)
        assert isinstance(task.label, str)
        assert task.value != ""
        assert task.label != ""

    # Test task descriptions make sense
    assert "energy" in Task.S2E.label
    assert "force" in Task.S2EF.label
    assert "stress" in Task.S2EFS.label
    assert "magmoms" in Task.S2EFSM.label


def test_model_type() -> None:
    """Test ModelType enum."""
    # Test basic enum functionality
    assert ModelType.GNN == "GNN"
    assert ModelType.GNN.label == "Graph Neural Network"
    assert ModelType.RF == "RF"
    assert ModelType.RF.label == "Random Forest"

    # Test that all model types have values and labels
    for model_type in ModelType:
        assert isinstance(model_type.value, str)
        assert isinstance(model_type.label, str)
        assert model_type.value != ""
        assert model_type.label != ""

    # Test model type descriptions make sense
    assert "Neural" in ModelType.GNN.label
    assert "Forest" in ModelType.RF.label
    assert "Transformer" in ModelType.Transformer.label
    assert "Fingerprint" in ModelType.Fingerprint.label


def test_open() -> None:
    """Test Open enum."""
    # Test basic enum functionality
    assert Open.OSOD == "OSOD"
    assert Open.OSOD.label == "open source, open data"
    assert Open.CSCD == "CSCD"
    assert Open.CSCD.label == "closed source, closed data"

    # Test that all openness types have values and labels
    for open_type in Open:
        assert isinstance(open_type.value, str)
        assert isinstance(open_type.label, str)
        assert open_type.value != ""
        assert open_type.label != ""

    # Test openness descriptions make sense
    assert "open source" in Open.OSOD.label
    assert "open data" in Open.OSOD.label
    assert "closed source" in Open.CSCD.label
    assert "closed data" in Open.CSCD.label


def test_test_subset() -> None:
    """Test TestSubset enum."""
    # Test basic enum functionality
    assert TestSubset.uniq_protos == "unique_prototypes"
    assert TestSubset.uniq_protos.label == "Unique Structure Prototypes"
    assert TestSubset.most_stable_10k == "most_stable_10k"
    assert TestSubset.most_stable_10k.label == "10k Most Stable Materials"

    # Test that all test subsets have values and labels
    for subset in TestSubset:
        assert isinstance(subset.value, str)
        assert isinstance(subset.label, str)
        assert subset.value != ""
        assert subset.label != ""

    # Test subset descriptions make sense
    assert "Unique" in TestSubset.uniq_protos.label
    assert "Stable" in TestSubset.most_stable_10k.label
    assert "Full" in TestSubset.full_test_set.label


def test_files_enum() -> None:
    """Test error handling in Files enum."""
    assert Files.base_dir == DATA_DIR

    # Test custom base_dir
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

    # Test __repr__ and __str__ methods
    test_file = SubFiles.test_file
    assert repr(test_file) == "SubFiles.test_file"
    assert str(test_file) == "test_file"

    # Test invalid label lookup
    label = "invalid-label"
    with pytest.raises(ValueError, match=f"{label=} not found in Files"):
        Files.from_label(label)


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


def test_data_files_enum() -> None:
    """Test DataFiles enum functionality."""
    # Test __repr__ and __str__ for DataFiles
    assert repr(DataFiles.mp_energies) == "DataFiles.mp_energies"
    assert str(DataFiles.mp_energies) == "mp_energies"

    # Test that paths are constructed correctly
    assert DataFiles.mp_energies.rel_path == "mp/2025-02-01-mp-energies.csv.gz"
    assert DataFiles.mp_patched_phase_diagram.rel_path == "mp/2023-02-07-ppd-mp.pkl.gz"
    assert DataFiles.mp_energies.name == "mp_energies"
    assert DataFiles.mp_energies.url.startswith("https://figshare.com/files/")

    # Test that multiple files exist and have correct attributes
    assert DataFiles.wbm_summary.rel_path == "wbm/2023-12-13-wbm-summary.csv.gz"
    assert DataFiles.wbm_summary.path == f"{DATA_DIR}/wbm/2023-12-13-wbm-summary.csv.gz"
    assert DataFiles.wbm_summary.url.startswith("https://figshare.com/files/")


def test_data_files_path_raises_when_md5_download_fails(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture
) -> None:
    """DataFiles.path surfaces failed checksum verification to callers."""
    data_file = DataFiles.mp_energies
    expected_md5 = "0" * 32
    monkeypatch.setattr(DataFiles, "_base_dir", str(tmp_path))
    monkeypatch.setenv("MBD_AUTO_DOWNLOAD_FILES", "true")
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
        pytest.raises(FileNotFoundError, match="Failed to download and verify"),
    ):
        _ = data_file.path

    stdout, stderr = capsys.readouterr()
    assert f"expected {expected_md5}" in stdout
    assert stderr == ""
    assert not os.path.isfile(f"{tmp_path}/{data_file.rel_path}")


@pytest.mark.parametrize("md_value", [None, "not available", 42])
def test_model_md_path_returns_none_for_non_dict_md(
    md_value: object, monkeypatch: pytest.MonkeyPatch
) -> None:
    """md_path returns None for any non-dict metrics.md (absent, the 'not available'
    placeholder, or an unexpected scalar) instead of raising AttributeError on .get.
    """
    model = Model.mace_mp_0
    monkeypatch.setattr(
        type(model), "metrics", property(lambda _self: {"md": md_value})
    )
    assert model.md_path is None


def test_model_md_path_returns_path_for_dict_md(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """md_path resolves local pred_file without downloading when no URL is available."""
    from matbench_discovery import ROOT
    from matbench_discovery import enums as enums_mod

    model = Model.mace_mp_0
    monkeypatch.setattr(
        type(model),
        "metrics",
        property(lambda _self: {"md": {"pred_file": "models/x/md.csv.gz"}}),
    )
    download_calls: list[tuple[object, ...]] = []
    monkeypatch.setattr(
        enums_mod,
        "maybe_auto_download_file",
        lambda *args, **_kwargs: download_calls.append(args),
    )
    assert model.md_path == f"{ROOT}/models/x/md.csv.gz"
    assert download_calls == []


@pytest.mark.parametrize("data_file", DataFiles)
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
        # (e.g. mp_trj_json_gz in the original MPtrj article), can't enforce md5
        pytest.skip(f"{data_file.name} file id {file_id} not in data-files article")
    declared_md5 = data_file.yaml[data_file.name].get("md5")
    assert declared_md5 == computed_md5, (
        f"data-files.yml md5 for {data_file.name} ({declared_md5}) does not match "
        f"Figshare computed_md5 ({computed_md5}) for file id {file_id}. Update the "
        "registry entry (or re-run scripts/upload_data_files_to_figshare.py)."
    )


def test_files_enum_auto_download(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture
) -> None:
    """Test auto-download behavior in Files class."""

    # Create a test Files class with our temp directory
    class TestFiles(Files, base_dir=str(tmp_path)):
        test_file = auto(), "test/file.txt"

        @property
        def url(self) -> str:
            """URL associated with the file."""
            return "https://example.com/file.txt"

        @property
        def label(self) -> str:
            """Label associated with the file."""
            return "test"

    test_file = TestFiles.test_file
    abs_path = f"{tmp_path}/test/file.txt"
    os.makedirs(os.path.dirname(abs_path), exist_ok=True)

    mock_response = make_mock_response(b"test content")

    # Test 1: Auto-download enabled (default)
    monkeypatch.setenv("MBD_AUTO_DOWNLOAD_FILES", "true")
    with patch("requests.get", return_value=mock_response):
        maybe_auto_download_file(test_file.url, abs_path, label=test_file.label)
    stdout, _ = capsys.readouterr()
    assert f"Downloading 'test' from {test_file.url!r}" in stdout
    assert os.path.isfile(abs_path)

    # Test 2: File already exists (no download attempt)
    with patch("requests.get") as mock_get:
        maybe_auto_download_file(test_file.url, abs_path, label=test_file.label)
        mock_get.assert_not_called()

    # Test 3: Auto-download disabled
    os.remove(abs_path)
    monkeypatch.setenv("MBD_AUTO_DOWNLOAD_FILES", "false")
    assert not os.path.isfile(abs_path)
    with (
        patch("requests.get", return_value=mock_response),
        patch("sys.stdin.isatty", return_value=False),
    ):
        maybe_auto_download_file(test_file.url, abs_path, label=test_file.label)
    assert os.path.isfile(abs_path)  # file should now be downloaded


def test_model_enum() -> None:
    """Test Model enum functionality."""
    # Test basic model attributes
    assert Model.alignn.name == "alignn"
    assert Model.alignn.rel_path == "alignn/alignn.yml"
    assert Model.alignn.pr_url == "https://github.com/janosh/matbench-discovery/pull/85"
    assert Model.alignn.label == "ALIGNN"

    # Test __repr__ and __str__ for Model
    assert repr(Model.alignn) == "Model.alignn"
    assert str(Model.alignn) == "alignn"

    # Test metadata property
    assert isinstance(Model.alignn.metadata, dict)

    # Test yaml_path property
    assert Model.alignn.yaml_path.endswith("alignn/alignn.yml")
    grace_kappa_path = Model.grace_2l_mptrj.kappa_103_path
    assert isinstance(grace_kappa_path, str)
    assert grace_kappa_path.endswith(
        "2024-11-20-kappa-103-FIRE-fmax=1e-4-symprec=1e-5.json.gz"
    )

    # Test Model metrics property
    metrics = Model.alignn.metrics
    assert isinstance(metrics, dict)
    assert {*metrics} >= {"discovery", "geo_opt", "phonons"}

    # Test registry-wide model properties and backing files
    for model in Model:
        assert os.path.isfile(model.yaml_path)
    for model in Model.active():
        assert "/models/" in model.discovery_path

    assert Model.mace_mp_0.label == "MACE-MP-0"
    assert Model.mace_mp_0.name == Model.mace_mp_0.value == "mace_mp_0"
    assert not Model.alphanet_mptrj.is_complete
    assert not Model.dpa_3_1_mptrj.is_complete
    model_keys = {model.key for model in Model}
    for model in Model:
        if model.metadata.get("status") == "superseded":
            assert model.metadata["superseded_by"] in model_keys


@pytest.mark.parametrize(
    "input_value, expected_model",
    [
        # Exact matches
        ("mace_mp_0", Model.mace_mp_0),
        ("eqv2_s_dens_mp", Model.eqv2_s_dens_mp),
        # Dash conversion
        ("mace-mp-0", Model.mace_mp_0),
        ("eqV2-s-dens-mp", Model.eqv2_s_dens_mp),
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
    rel_path = "models/test/md.csv.gz"
    abs_path = f"{tmp_path}/{rel_path}"
    monkeypatch.setenv("HF_TOKEN", "hf_secret")
    monkeypatch.setenv("MBD_AUTO_DOWNLOAD_FILES", "true")
    monkeypatch.setattr("matbench_discovery.enums.ROOT", str(tmp_path))
    monkeypatch.setitem(
        Model.alignn.__dict__,
        "metadata",
        {
            "model_name": "Gated MD",
            "metrics": {"md": {"pred_file": rel_path, "pred_file_url": url}},
        },
    )

    with patch("requests.get", return_value=make_mock_response(b"md")) as mock_get:
        assert Model.alignn.md_path == abs_path

    assert os.path.isfile(abs_path)
    assert mock_get.call_args.args == (url,)
    assert mock_get.call_args.kwargs["headers"] == {"Authorization": "Bearer hf_secret"}


def get_urls_from_dict(
    dct: dict[str, Any], parent_key: str = ""
) -> list[tuple[str, str]]:
    """Recursively find all keys ending in _url in a nested dictionary.
    Returns list of tuples with (dotted.path.to.key, url_value).
    """
    urls = []
    for key, val in dct.items():
        current_key = f"{parent_key}.{key}" if parent_key else key

        if key.endswith("_url") and isinstance(val, str):
            urls.append((current_key, val))
        elif isinstance(val, dict):
            urls.extend(get_urls_from_dict(val, current_key))

    return urls


TIMEOUT = 30
TRANSIENT_URL_STATUSES = (429, 500, 502, 503, 504)
VALID_URL_STATUSES = {200, 202, 403, 429}


@pytest.fixture(scope="session")
def url_session() -> requests.Session:
    """HTTP session that retries transient errors (429, 5xx) with backoff, so a
    momentary server hiccup doesn't fail URL-validation tests.
    """
    session = requests.Session()
    session.headers["User-Agent"] = "unit test"
    n_pool = 2 * (os.cpu_count() or 4)
    session.mount(
        "https://",
        requests.adapters.HTTPAdapter(
            pool_connections=n_pool,
            pool_maxsize=n_pool,
            max_retries=requests.adapters.Retry(
                total=3,
                backoff_factor=1,
                status_forcelist=TRANSIENT_URL_STATUSES,
                raise_on_status=False,
            ),
        ),
    )
    return session


def check_url(session: requests.Session, url: str) -> None:
    """Assert a model URL resolves. The session adapter retries transient errors
    (connection failures, 429, 5xx), so only persistently bad links (e.g. 404) fail.
    """
    # 200: OK, 202: figshare async, 403: restricted, 429: rate limited after retries
    status = session.head(url, allow_redirects=True, timeout=TIMEOUT).status_code
    assert status in VALID_URL_STATUSES, f"unexpected {status=} for {url}"


def test_model_prediction_urls(url_session: requests.Session) -> None:
    """Test that all model prediction file URLs are valid."""
    import concurrent.futures
    import multiprocessing as mp

    tasks: dict[str, str] = {}
    for model in Model.active():
        if model.name == Model.mace_mpa_0.name:
            continue
        tasks[model.pr_url] = model.name
        metrics = model.metrics
        if not metrics:
            continue
        for key_path, url in get_urls_from_dict(metrics):
            tasks[url] = f"{model.name}.{key_path}"

    n_workers = min(len(tasks), mp.cpu_count())
    errors: list[Exception] = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=n_workers) as executor:
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
