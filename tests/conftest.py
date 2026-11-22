"""Shared pytest fixtures for the test suite."""

from collections.abc import Iterator

import numpy as np
import pandas as pd
import pytest
import requests
from kemio import Structure
from requests.adapters import HTTPAdapter, Retry

from matbench_discovery import cli
from matbench_discovery.metrics.diatomics import DiatomicCurves
from matbench_discovery.remote import fetch

# scripts and helpers read cli.cli_args, which strictly parses sys.argv on first use;
# give them the shared defaults instead of pytest's own flags
cli.parse_cli_args([])


@pytest.fixture(autouse=True)
def no_download_retry_delay(monkeypatch: pytest.MonkeyPatch) -> None:
    """Keep download_file's transient-error retries instant in tests."""
    monkeypatch.setattr(fetch, "HTTP_RETRY_DELAYS", (0,) * 4)


@pytest.fixture(scope="session")
def url_session() -> Iterator[requests.Session]:
    """Reuse HTTP connections and retry transient errors for live artifact checks."""
    with requests.Session() as session:
        session.headers["User-Agent"] = "unit test"
        session.mount(
            "https://",
            HTTPAdapter(
                pool_connections=16,
                pool_maxsize=16,
                max_retries=Retry(
                    total=3,
                    backoff_factor=1,
                    status_forcelist=(429, 500, 502, 503, 504),
                    raise_on_status=False,
                ),
            ),
        )
        yield session


def _monoclinic_lattice(a: float, b: float, c: float, beta_deg: float) -> np.ndarray:
    """Monoclinic lattice rows with b along y and beta between a and c."""
    beta = np.radians(beta_deg)
    return np.array([[a, 0, 0], [0, b, 0], [c * np.cos(beta), 0, c * np.sin(beta)]])


@pytest.fixture
def dummy_struct() -> Structure:
    """Cubic FeO (CsCl-like) test structure."""
    coords = ((0, 0, 0), (0.5, 0.5, 0.5))
    return Structure(4.2 * np.eye(3), ("Fe", "O"), coords, fractional=True)


@pytest.fixture
def cubic_struct() -> Structure:
    """Cubic Si2 test structure."""
    return Structure(
        4.2 * np.eye(3), ["Si", "Si"], [[0, 0, 0], [0.5, 0.5, 0.5]], fractional=True
    )


@pytest.fixture
def tetragonal_struct() -> Structure:
    """Tetragonal TiO2-like test structure."""
    coords = [[0, 0, 0], [0.5, 0.5, 0], [0, 0.5, 0.5]]
    return Structure(
        np.diag([4.0, 4.0, 6.0]), ["Ti", "O", "O"], coords, fractional=True
    )


@pytest.fixture
def monoclinic_struct() -> Structure:
    """Monoclinic PO4-like test structure."""
    x1 = [0, 0, 0]
    x2 = [0.25, 0.25, 0.25]
    x3 = [0.75, 0.75, 0.25]
    x4 = [0.25, 0.75, 0.75]
    x5 = [0.75, 0.25, 0.75]
    lattice = _monoclinic_lattice(5, 4, 6, 100)
    return Structure(
        lattice, ["P", "O", "O", "O", "O"], [x1, x2, x3, x4, x5], fractional=True
    )


@pytest.fixture
def df_float() -> pd.DataFrame:
    rng = np.random.default_rng(seed=0)
    return pd.DataFrame(rng.normal(size=(10, 5)), columns=list("ABCDE"))


@pytest.fixture
def df_mixed() -> pd.DataFrame:
    rng = np.random.default_rng(seed=0)

    floats = rng.random(size=10)
    bools = rng.choice([True, False], size=10)
    strings = rng.choice([*"abcdef"], size=10)
    return pd.DataFrame(dict(floats=floats, bools=bools, strings=strings))


@pytest.fixture
def pred_ref_diatomic_curves() -> tuple[DiatomicCurves, DiatomicCurves]:
    """Reference and predicted diatomic curves for H and He."""
    distances = np.logspace(1, -1, 40)
    # Reference: Simple Morse potential
    e_h_ref = 5 * (1 - np.exp(-2 * (distances - 1.5))) ** 2 - 5
    e_he_ref = 4 * (1 - np.exp(-2 * (distances - 1.5))) ** 2 - 4
    # Prediction: Slightly perturbed Morse potential
    e_h_pred = 5.2 * (1 - np.exp(-2.1 * (distances - 1.48))) ** 2 - 5.1
    e_he_pred = 4.2 * (1 - np.exp(-2.1 * (distances - 1.48))) ** 2 - 4.1

    np_rng = np.random.default_rng(seed=0)

    # forces have shape (2, 3), 2 atoms with 3 force components
    f_h_ref = 0.5 * np_rng.normal(size=(len(distances), 2, 3))
    f_h_pred = 0.55 * np_rng.normal(size=(len(distances), 2, 3))
    f_he_ref = 0.1 * np_rng.normal(size=(len(distances), 2, 3))
    f_he_pred = 0.15 * np_rng.normal(size=(len(distances), 2, 3))

    ref_curves = {
        "distances": distances,
        "homo-nuclear": {
            "H": {"energies": e_h_ref, "forces": f_h_ref},
            "He": {"energies": e_he_ref, "forces": f_he_ref},
        },
    }
    pred_curves = {
        "distances": distances,
        "homo-nuclear": {
            "H": {"energies": e_h_pred, "forces": f_h_pred},
            "He": {"energies": e_he_pred, "forces": f_he_pred},
        },
    }
    return DiatomicCurves.from_dict(ref_curves), DiatomicCurves.from_dict(pred_curves)
