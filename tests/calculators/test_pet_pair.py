"""Registration is dependency-free; runtime checks require the optional model."""

from pathlib import Path

import pytest

from matbench_discovery.calculators import CALCULATORS


def test_pet_pair_registration_without_optional_runtime() -> None:
    """The registry remains usable without installing every model dependency."""
    assert "pet_xs_pair_bo_v1" in CALCULATORS


def test_pet_pair_shared_registry_reaches_precision_guard(tmp_path: Path) -> None:
    """With the model installed, dropped dtype bypasses its CPU/fp32 contract."""
    pytest.importorskip("pet_pair_bo", reason="optional shared-runner model runtime")
    with pytest.raises(ValueError, match="float32"):
        CALCULATORS["pet_xs_pair_bo_v1"].make_calc(
            "cpu", checkpoint=str(tmp_path / "not-read.json"), dtype="float64"
        )
