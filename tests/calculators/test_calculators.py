"""Tests for calculator registry and checkpoint management."""

import hashlib
import inspect
import os
import subprocess
import sys
import zipfile
from collections.abc import Sequence
from pathlib import Path
from types import ModuleType, SimpleNamespace

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


def test_corrected_model_factory_contracts() -> None:
    """Registry factories retain model variants, modalities, and dtype support."""
    expected_closure_values = {
        "orb_v3": ("variant", "orb-v3-conservative-inf-mpa"),
        "sevennet_mf_ompa": ("modal", "mpa"),
        "alphanet_mptrj": ("config_url", "/pretrained/MPtrj/mp.json"),
    }
    for model_key, (field_name, expected) in expected_closure_values.items():
        value = inspect.getclosurevars(CALCULATORS[model_key].make_calc).nonlocals[
            field_name
        ]
        assert value == expected or str(value).endswith(expected)
    assert "orb-v3-conservative-inf-mpa" in Model.orb_v3.metadata["checkpoint_url"]
    assert {"md", "diatomics"}.isdisjoint(Model.orb_v3.metrics)


def test_alphanet_factory_honors_requested_dtype(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """AlphaNet maps the runner dtype to its calculator precision."""

    class FakeConfig:
        """Minimal AlphaNet config loader."""

        def from_json(self, path: str) -> dict[str, str]:
            """Return the downloaded config path."""
            return {"path": path}

    config_module = ModuleType("alphanet.config")
    config_module.__dict__["All_Config"] = FakeConfig
    calculator_module = ModuleType("alphanet.infer.calc")
    calculator_module.__dict__["AlphaNetCalculator"] = SimpleNamespace
    for module_name, module in {
        "alphanet": ModuleType("alphanet"),
        "alphanet.config": config_module,
        "alphanet.infer": ModuleType("alphanet.infer"),
        "alphanet.infer.calc": calculator_module,
    }.items():
        monkeypatch.setitem(sys.modules, module_name, module)

    from matbench_discovery.remote import fetch

    def download_config(destination: str, _url: str) -> None:
        """Write a minimal AlphaNet config."""
        Path(destination).write_text("{}", encoding="utf-8")

    monkeypatch.setattr(fetch, "download_file", download_config)
    monkeypatch.setattr(calculators, "CHECKPOINT_DIR", str(tmp_path))
    factory = calculators._alphanet("fake", "https://example.com/config.json")  # noqa: SLF001
    calculator = factory("cpu", dtype="float64", checkpoint="/tmp/model.ckpt")
    assert vars(calculator)["precision"] == "64"


def test_pet_factory_casts_exported_model(monkeypatch: pytest.MonkeyPatch) -> None:
    """PET casts the loaded model and advertises the requested dtype."""
    import torch

    captured: dict[str, object] = {}
    model_capabilities = SimpleNamespace(dtype="float32")

    def capabilities_stub() -> SimpleNamespace:
        """Return mutable model capabilities."""
        return model_capabilities

    def cast_model(**kwargs: object) -> SimpleNamespace:
        """Capture the requested torch dtype and device."""
        captured.update(kwargs)
        return fake_model

    fake_model = SimpleNamespace(capabilities=capabilities_stub, to=cast_model)

    def load_model(_path: str) -> SimpleNamespace:
        """Return the fake exported model."""
        return fake_model

    def make_calculator(model: object, **kwargs: object) -> SimpleNamespace:
        """Capture the model passed to MetatomicCalculator."""
        captured.update(model=model, **kwargs)
        return SimpleNamespace()

    torch_module = ModuleType("metatomic.torch")
    torch_module.__dict__["load_atomistic_model"] = load_model
    calculator_module = ModuleType("metatomic.torch.ase_calculator")
    calculator_module.__dict__["MetatomicCalculator"] = make_calculator
    for module_name, module in {
        "metatomic": ModuleType("metatomic"),
        "metatomic.torch": torch_module,
        "metatomic.torch.ase_calculator": calculator_module,
    }.items():
        monkeypatch.setitem(sys.modules, module_name, module)

    def export_stub(
        _command: Sequence[str], _destination: str, **kwargs: Sequence[str]
    ) -> None:
        """Verify the export cache identity inputs."""
        assert kwargs["source_paths"]
        assert kwargs["tool_packages"] == ("metatrain",)

    monkeypatch.setattr(calculators, "_run_to_atomic_output", export_stub)
    calculators._pet("fake")(  # noqa: SLF001
        "cpu", dtype="float64", checkpoint="/tmp/model.ckpt"
    )
    assert fake_model.capabilities().dtype == "float64"
    assert captured == {
        "dtype": torch.float64,
        "device": "cpu",
        "model": fake_model,
        "non_conservative": False,
    }


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
    calculators._run_to_atomic_output(  # noqa: SLF001
        ["nequip-compile", "registry:model", "{output}"], str(destination)
    )
    assert n_calls == 2

    source = tmp_path / "source.ckpt"
    source.write_bytes(b"version 1")
    command = ["export", str(source), "{output}"]
    calculators._run_to_atomic_output(  # noqa: SLF001
        command, str(destination), source_paths=(str(source),)
    )
    source.write_bytes(b"version 2")
    calculators._run_to_atomic_output(  # noqa: SLF001
        command, str(destination), source_paths=(str(source),)
    )
    assert n_calls == 4

    tool_versions = {"compiler": "1.0"}

    def version_stub(package: str) -> str:
        """Return the mutable fake compiler version."""
        return tool_versions[package]

    monkeypatch.setattr(calculators, "version", version_stub)
    for expected_calls in (5, 5):
        calculators._run_to_atomic_output(  # noqa: SLF001
            command,
            str(destination),
            source_paths=(str(source),),
            tool_packages=("compiler",),
        )
        assert n_calls == expected_calls
    tool_versions["compiler"] = "2.0"
    calculators._run_to_atomic_output(  # noqa: SLF001
        command,
        str(destination),
        source_paths=(str(source),),
        tool_packages=("compiler",),
    )
    assert n_calls == 6


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
