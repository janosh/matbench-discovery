"""Tests for calculator registry and checkpoint management."""

import hashlib
import inspect
import os
import subprocess
import sys
import zipfile
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
        "alphanet_v1_mptrj": ("config_url", "/pretrained/MPtrj/mp.json"),
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
        "alphanet.config": config_module,
        "alphanet.infer.calc": calculator_module,
    }.items():
        monkeypatch.setitem(sys.modules, module_name, module)

    monkeypatch.setattr(calculators, "CHECKPOINT_DIR", str(tmp_path))
    config_url = "https://example.com/config.json"
    config_hash = hashlib.sha256(config_url.encode()).hexdigest()[:12]
    (tmp_path / f"fake-config-{config_hash}.json").write_text("{}", encoding="utf-8")
    factory = calculators._alphanet("fake", config_url)  # noqa: SLF001
    calculator = factory(
        "cpu", dtype="float64", checkpoint=str(tmp_path / "model.ckpt")
    )
    assert vars(calculator)["precision"] == "64"


def test_pet_factory_casts_exported_model(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """PET casts the loaded model and advertises the requested dtype."""
    captured: dict[str, object] = {}
    model_capabilities = SimpleNamespace(dtype="float32")
    float64_dtype = object()

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

    def make_calculator(_model: object, **_kwargs: object) -> SimpleNamespace:
        """Return a placeholder Metatomic calculator."""
        return SimpleNamespace()

    torch_module = ModuleType("torch")
    torch_module.__dict__.update(float32=object(), float64=float64_dtype)
    metatomic_torch_module = ModuleType("metatomic.torch")
    metatomic_torch_module.__dict__["load_atomistic_model"] = load_model
    calculator_module = ModuleType("metatomic.torch.ase_calculator")
    calculator_module.__dict__["MetatomicCalculator"] = make_calculator
    for module_name, module in {
        "torch": torch_module,
        "metatomic": ModuleType("metatomic"),
        "metatomic.torch": metatomic_torch_module,
        "metatomic.torch.ase_calculator": calculator_module,
    }.items():
        monkeypatch.setitem(sys.modules, module_name, module)

    def export_stub(*_args: object, **_kwargs: object) -> None:
        """Stand in for the separately tested export cache."""

    monkeypatch.setattr(calculators, "_run_to_atomic_output", export_stub)
    calculators._pet("fake")(  # noqa: SLF001
        "cpu", dtype="float64", checkpoint=str(tmp_path / "model.ckpt")
    )
    assert fake_model.capabilities().dtype == "float64"
    assert captured == {"dtype": float64_dtype, "device": "cpu"}


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


def test_atomic_command_output_recovers_after_timeout(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Timed-out compilers leave no accepted partial cache and can retry."""
    destination = tmp_path / "fresh" / "compiled.nequip.pth"
    n_calls = 0

    def run_stub(command: list[str], *, check: bool, timeout: float) -> None:
        """Write a partial first result and a complete retry result."""
        nonlocal n_calls
        assert check
        assert timeout == calculators.DERIVED_ARTIFACT_TIMEOUT_SEC
        n_calls += 1
        output_path = next(argument for argument in command if ".tmp." in argument)
        Path(output_path).write_bytes(b"partial" if n_calls == 1 else b"complete")
        if n_calls == 1:
            raise subprocess.TimeoutExpired(command, timeout)

    monkeypatch.setattr(calculators.subprocess, "run", run_stub)
    with pytest.raises(subprocess.TimeoutExpired):
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

    tool_versions: dict[str, str | None] = {"compiler": "1.0"}

    def version_stub(package: str) -> str:
        """Return the mutable fake compiler version."""
        package_version = tool_versions[package]
        if package_version is None:
            raise calculators.PackageNotFoundError(package)
        return package_version

    monkeypatch.setattr(calculators, "version", version_stub)
    seen_tool_versions: set[str | None] = set()
    for package_version in ("1.0", "1.0", "2.0", None, None):
        tool_versions["compiler"] = package_version
        seen_tool_versions.add(package_version)
        calculators._run_to_atomic_output(  # noqa: SLF001
            command,
            str(destination),
            source_paths=(str(source),),
            tool_packages=("compiler",),
        )
        assert n_calls == 4 + len(seen_tool_versions)


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
