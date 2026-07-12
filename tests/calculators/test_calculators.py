"""Tests for calculator registry and checkpoint management."""

import hashlib
import inspect
import os
import subprocess
import zipfile
from pathlib import Path
from types import SimpleNamespace

import pytest
from ase.calculators.emt import EMT

from matbench_discovery import calculators
from matbench_discovery.calculators import CALCULATORS, load_calculator
from matbench_discovery.enums import Model


def test_load_calculator(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Emt loads with no extra deps; dtype is ignored by models that don't declare it;
    unknown keys raise a helpful error.
    """
    assert isinstance(load_calculator("emt"), EMT)
    assert isinstance(load_calculator("emt", dtype="float32"), EMT)  # dtype ignored

    seen_dtype = ""
    seen_checkpoint: str | None = None

    def configurable_calc(
        device: str, dtype: str = "float64", checkpoint: str | None = None
    ) -> EMT:
        """Capture supported optional calculator arguments."""
        nonlocal seen_checkpoint, seen_dtype
        del device
        seen_dtype = dtype
        seen_checkpoint = checkpoint
        return EMT()

    monkeypatch.setitem(
        calculators.CALCULATORS, "mace_mp_0", calculators.CalcSpec(configurable_calc)
    )
    checkpoint_path = str(tmp_path / "model.pt")
    Path(checkpoint_path).write_bytes(b"checkpoint")
    load_calculator(
        "mace-mp-0", device="cpu", dtype="float32", checkpoint=checkpoint_path
    )
    assert seen_dtype == "float32"
    assert seen_checkpoint == os.path.abspath(checkpoint_path)

    def download_checkpoint_stub(model_key: str, ext: str | None = None) -> str:
        """Return the local artifact used to test registry-managed downloads."""
        assert model_key == "mace_mp_0"
        assert ext == ".pt"
        return checkpoint_path

    monkeypatch.setattr(calculators, "download_checkpoint", download_checkpoint_stub)
    monkeypatch.setitem(
        calculators.CALCULATORS,
        "mace_mp_0",
        calculators.CalcSpec(
            configurable_calc, auto_checkpoint=True, checkpoint_ext=".pt"
        ),
    )
    load_calculator("mace_mp_0", device="cpu")
    assert seen_checkpoint == checkpoint_path

    monkeypatch.setitem(
        calculators.CALCULATORS,
        "equiformer_v3_mp",
        calculators.CalcSpec(configurable_calc, requires_checkpoint=True),
    )
    with pytest.raises(ValueError, match="requires an explicit checkpoint"):
        load_calculator("equiformer_v3_mp", device="cpu")

    with pytest.raises(ValueError, match="Unknown model"):
        load_calculator("does-not-exist")


def test_equflash_uv_command_uses_isolated_override_project() -> None:
    """EquFlash relaxes fairchem's torch bound without affecting other models."""
    calc_spec = CALCULATORS["equflash_29m_oam"]
    assert calc_spec.uv_run_cmd("models/run_kappa.py") == [
        "uv",
        "run",
        "--project",
        "models/equflash/kappa-env",
        "--isolated",
        "--python",
        "3.12",
        "python",
        "models/run_kappa.py",
    ]
    assert calc_spec.project
    assert os.path.isfile(f"{calc_spec.project}/pyproject.toml")


def test_checkpoint_specs_accept_declared_overrides() -> None:
    """Managed and required checkpoints must reach their calculator factories."""
    for model_key, calc_spec in CALCULATORS.items():
        if calc_spec.auto_checkpoint or calc_spec.requires_checkpoint:
            assert "checkpoint" in inspect.signature(calc_spec.make_calc).parameters, (
                model_key
            )


@pytest.mark.parametrize(
    "model_key", ["sevennet_0", "sevennet_l3i5", "sevennet_omni_i12"]
)
def test_named_sevennet_models_keep_package_managed_checkpoints(model_key: str) -> None:
    """Name-based SevenNet loaders must not be replaced by registry downloads."""
    assert not CALCULATORS[model_key].auto_checkpoint


def test_deepmd_archive_extraction_replaces_invalid_cache(tmp_path: Path) -> None:
    """Archive extraction atomically replaces stale multi-model directories."""
    archive_path = tmp_path / "model.zip"
    with zipfile.ZipFile(archive_path, mode="w") as archive:
        archive.writestr("nested/model.pth", b"frozen model")
    extract_dir = tmp_path / "model-extracted"
    extract_dir.mkdir()
    (extract_dir / "stale.pb").write_bytes(b"stale")
    (extract_dir / "orphan.pt2").write_bytes(b"orphan")

    model_path = calculators._extract_single_deepmd_model(  # noqa: SLF001
        str(archive_path)
    )

    assert Path(model_path).read_bytes() == b"frozen model"
    assert [
        path.name
        for path in extract_dir.rglob("*")
        if path.suffix in {".pb", ".pth", ".pt2"}
    ] == ["model.pth"]
    with zipfile.ZipFile(archive_path, mode="w") as archive:
        archive.writestr("nested/model.pth", b"updated frozen model")
    updated_path = calculators._extract_single_deepmd_model(  # noqa: SLF001
        str(archive_path)
    )
    assert Path(updated_path).read_bytes() == b"updated frozen model"


def test_atomic_command_output_recovers_after_interruption(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Interrupted compilers leave no accepted partial cache and can retry."""
    destination = tmp_path / "compiled.nequip.pth"
    n_calls = 0

    def run_stub(command: list[str], *, check: bool) -> None:
        """Write a partial first result and a complete retry result."""
        nonlocal n_calls
        assert check
        n_calls += 1
        output_path = next(argument for argument in command if ".tmp." in argument)
        Path(output_path).write_bytes(b"partial" if n_calls == 1 else b"complete")
        if n_calls == 1:
            raise subprocess.CalledProcessError(1, command)

    monkeypatch.setattr(calculators.subprocess, "run", run_stub)
    with pytest.raises(subprocess.CalledProcessError):
        calculators._run_to_atomic_output(  # noqa: SLF001
            ["nequip-compile", "registry:model", "{output}"], str(destination)
        )
    assert not destination.is_file()

    calculators._run_to_atomic_output(  # noqa: SLF001
        ["nequip-compile", "registry:model", "{output}"], str(destination)
    )
    assert destination.read_bytes() == b"complete"


def test_download_checkpoint_replaces_zero_byte_cache(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """download_checkpoint redownloads zero-byte cached checkpoint files."""
    from matbench_discovery.remote import fetch

    url = "https://example.com/model.ckpt"
    monkeypatch.setattr(calculators, "CHECKPOINT_DIR", f"{tmp_path}")
    from_ref_mock = classmethod(
        lambda _model_cls, _model_key: SimpleNamespace(metadata={"checkpoint_url": url})
    )
    monkeypatch.setattr(Model, "from_ref", from_ref_mock)

    download_file_mock = (  # noqa: E731
        lambda destination_path, _url, **_kwargs: Path(destination_path).write_bytes(
            b"checkpoint"
        )
    )
    monkeypatch.setattr(fetch, "download_file", download_file_mock)
    url_hash = hashlib.sha256(url.encode()).hexdigest()[:12]
    dest = tmp_path / f"fake-{url_hash}.ckpt"
    dest.write_bytes(b"")

    actual_path = os.path.normpath(calculators.download_checkpoint("fake"))
    assert actual_path == os.path.normpath(f"{dest}")
    assert dest.read_bytes() == b"checkpoint"
