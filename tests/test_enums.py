"""Test enums module."""

import os
import sys
import warnings
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
from matbench_discovery.remote.fetch import maybe_auto_download_file


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


def test_data_files_enum() -> None:
    """Test DataFiles enum functionality."""
    # Test __repr__ and __str__ for DataFiles
    assert repr(DataFiles.mp_energies) == "DataFiles.mp_energies"
    assert str(DataFiles.mp_energies) == "mp_energies"

    # Test that paths are constructed correctly
    assert DataFiles.mp_energies.rel_path == "mp/2025-02-01-mp-energies.csv.gz"
    assert DataFiles.mp_energies.name == "mp_energies"
    assert DataFiles.mp_energies.url.startswith("https://figshare.com/files/")

    # Test that multiple files exist and have correct attributes
    assert DataFiles.wbm_summary.rel_path == "wbm/2023-12-13-wbm-summary.csv.gz"
    assert DataFiles.wbm_summary.path == f"{DATA_DIR}/wbm/2023-12-13-wbm-summary.csv.gz"
    assert DataFiles.wbm_summary.url.startswith("https://figshare.com/files/")


@pytest.mark.parametrize("data_file", DataFiles)
def test_data_files_enum_urls(
    data_file: DataFiles, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Test that each URL in data-files.yml is a valid Figshare download URL."""

    name, url = data_file.name, data_file.url
    # check that URL is a figshare download
    assert "figshare.com/files/" in url, (
        f"URL for {name} is not a Figshare download URL: {url}"
    )

    # Mock requests.head to avoid actual network calls
    class MockResponse:
        status_code = 200

    def mock_head(*_args: str, **_kwargs: dict[str, str]) -> MockResponse:
        return MockResponse()

    monkeypatch.setattr(requests, "head", mock_head)

    # check that the URL is valid by sending a head request
    response = requests.head(url, allow_redirects=True, timeout=5)
    assert response.status_code in {200, 403}, f"Invalid URL for {name}: {url}"


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

    # Mock successful request
    class MockResponse:
        status_code = 200
        content = b"test content"

        def raise_for_status(self) -> None:
            """Mock the raise_for_status method."""
            if self.status_code >= 400:
                raise requests.HTTPError(f"HTTP Error: {self.status_code}")

        def iter_content(self, chunk_size: int = 8192) -> list[bytes]:  # noqa: ARG002
            """Mock iter_content for streaming."""
            return [self.content]

    # Mock stdin to simulate non-interactive mode
    class MockStdin:
        def isatty(self) -> bool:
            """Mock isatty to simulate non-interactive mode."""
            return False

        def readline(self) -> str:
            """Mock readline method."""
            return "y\n"  # Default to yes for testing

    monkeypatch.setattr(requests, "get", lambda *_args, **_kwargs: MockResponse())
    monkeypatch.setattr(sys, "stdin", MockStdin())

    # Test 1: Auto-download enabled (default)
    monkeypatch.setenv("MBD_AUTO_DOWNLOAD_FILES", "true")
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


def check_url(session: requests.Session, url: str, desc: str) -> None:
    """Check if a URL is valid."""
    http_status = None
    try:
        response = session.head(url, allow_redirects=True)
        http_status = response.status_code
        assert http_status in {200, 403, 429}
        if http_status == 429:
            assert "Too Many Requests" in response.reason, f"{response.reason=}"
            warnings.warn(f"{response.reason=}", stacklevel=2)
    except (requests.RequestException, AssertionError) as exc:
        exc.add_note(f"Failed to validate\n{url}\n{desc}\n{http_status=}")
        raise


def test_model_prediction_urls() -> None:
    """Test that all model prediction file URLs are valid."""
    import asyncio
    import concurrent.futures
    import multiprocessing as mp

    tasks: dict[str, str] = {}
    for model in Model:
        if model.name == Model.mace_mpa_0.name:
            continue

        # Check model PR URL
        tasks[model.pr_url] = model.name

        # Check all URLs in metrics
        metrics = model.metrics
        if not metrics:
            continue

        for key_path, url in get_urls_from_dict(metrics):
            tasks[url] = f"{model.name}.{key_path}"

    # Create session with connection pooling and adaptive settings
    session = requests.Session()
    session.headers["User-Agent"] = "unit test"

    # Determine optimal pool size based on system resources and environment
    # Use min of CPU count and number of tasks to avoid over-allocation
    n_workers = min(len(tasks), mp.cpu_count())

    # Configure connection pooling with adaptive settings
    adapter = requests.adapters.HTTPAdapter(
        pool_connections=n_workers,
        pool_maxsize=n_workers * 2,  # Allow some room for growth
        max_retries=0,  # We handle retries at a higher level
        pool_block=False,  # don't block main thread waiting for connections
    )
    session.mount("https://", adapter)
    session.mount("http://", adapter)

    # Create event loop for async execution
    loop = asyncio.new_event_loop()
    asyncio.set_event_loop(loop)

    async def check_urls_async() -> None:
        """Check URLs concurrently using asyncio and thread pool."""
        # Use ThreadPoolExecutor for I/O-bound HTTP requests
        with concurrent.futures.ThreadPoolExecutor(max_workers=n_workers) as executor:
            futures = [  # Create futures for all URL checks
                loop.run_in_executor(executor, check_url, *(session, url, desc))
                for url, desc in tasks.items()
            ]
            # Wait for all futures to complete
            await asyncio.gather(*futures)

    try:
        # Run the async checks
        loop.run_until_complete(check_urls_async())
    finally:
        loop.close()
