"""Tests for calculator registry and checkpoint management."""

import hashlib
import inspect
import os
import subprocess
import sys
import tomllib
import zipfile
from collections.abc import Callable
from pathlib import Path
from types import ModuleType, SimpleNamespace
from unittest.mock import Mock

import pytest
from ase.calculators.calculator import Calculator
from ase.calculators.emt import EMT
from packaging.requirements import Requirement
from packaging.version import Version

from matbench_discovery import ROOT, calculators, hpc
from matbench_discovery.calculators import CALCULATORS, load_calculator
from matbench_discovery.enums import Model
from tests.utils import patch_model_metadata


def test_load_calculator(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """EMT loads with no deps, undeclared non-default dtype and unknown keys raise."""
    # emt has no model YAML: must stay an eager CalcSpec, not _runtime_calc_spec.
    assert isinstance(CALCULATORS._data["emt"], calculators.CalcSpec)  # noqa: SLF001
    assert isinstance(load_calculator("emt"), EMT)
    assert isinstance(load_calculator("emt", dtype="float64"), EMT)
    with pytest.raises(ValueError, match="emt calculator has no dtype option"):
        load_calculator("emt", dtype="float32")

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


def test_load_calculator_with_runtime(monkeypatch: pytest.MonkeyPatch) -> None:
    """Runtime records weights' device/dtype and load time, rejecting ignored args."""
    calc, runtime = calculators.load_calculator_with_runtime("emt")
    assert isinstance(calc, EMT)
    cpu = hpc.detect_hardware("cpu")
    labels = ("backend", "device", "hardware", "acceleration", "n_processes")
    assert [runtime[key] for key in labels] == ["none", "cpu", cpu, [], 1]
    assert "dtype" not in runtime  # unknown for non-torch calculators w/o dtype arg
    assert runtime["load_time_sec"] >= 0
    assert "ase" in runtime["versions"]

    torch = pytest.importorskip("torch")

    class TorchCalc(EMT):
        """EMT carrying a float32 torch module, like an MLIP calculator."""

        def __init__(self) -> None:
            super().__init__()
            self.model = torch.nn.Linear(2, 2).to(torch.float32)

    def float64_factory(device: str, dtype: str = "float64") -> EMT:
        """Factory that ignores dtype (weights stay float32)."""
        del device, dtype
        return TorchCalc()

    def no_dtype_factory(device: str) -> Calculator:
        """Factory without a dtype option: the detected dtype is recorded."""
        del device
        return calculators._with_acceleration(TorchCalc(), "cuequivariance")  # noqa: SLF001

    monkeypatch.setitem(
        calculators.CALCULATORS, "mace_mp_0", calculators.CalcSpec(no_dtype_factory)
    )
    _, runtime = calculators.load_calculator_with_runtime("mace_mp_0", device="cpu")
    assert (runtime["backend"], runtime["dtype"]) == ("torch", "float32")
    assert runtime["acceleration"] == ["cuequivariance"]
    with pytest.raises(ValueError, match="factory ignores the device"):
        calculators.load_calculator_with_runtime("mace_mp_0", device="cuda")
    monkeypatch.setitem(
        calculators.CALCULATORS, "mace_mp_0", calculators.CalcSpec(float64_factory)
    )
    with pytest.raises(ValueError, match="weights are float32 but dtype='float64'"):
        calculators.load_calculator_with_runtime("mace_mp_0", device="cpu")


def test_equflash_uv_command_uses_isolated_override_project() -> None:
    """EquFlash installs via project env; YAML deps are provenance-only."""
    calc_spec = CALCULATORS["equflash_29m_oam"]
    # exact equality also pins that project mode injects no --with, which would fight
    # the project's override-dependencies
    assert calc_spec.uv_run_cmd("models/run_kappa.py") == [
        "uv",
        "run",
        "--project",
        "models/equflash/kappa-env",
        "--isolated",
        "python",
        "models/run_kappa.py",
    ]
    assert calc_spec.deps  # still exposed for site provenance
    assert calc_spec.project


def test_uv_command_layers_model_packages_on_runner_env(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """--with mode drops YAML Python and runner-provided pins, keeps model packages."""
    environment = {
        "python_version": "3.11",  # provenance only
        "dependencies": ["ase==3.22.1", "Pymatgen>=2024", "torch==2.0.1", "fake-mlip"],
        "find_links": ["https://data.pyg.org/whl/torch-2.0.1.html"],
        "extra_index_urls": ["https://example.org/simple"],
    }
    patch_model_metadata(monkeypatch, Model["mace_mp_0"], {"environment": environment})
    calc_spec = calculators._runtime_calc_spec("mace_mp_0", EMT)()  # noqa: SLF001
    assert calc_spec.deps == tuple(environment["dependencies"])
    assert calc_spec.uv_run_cmd("models/run_md.py", "--model", "mace_mp_0") == [
        *("uv", "run", "--no-project"),
        *("--with", "torch==2.0.1", "--with", "fake-mlip"),
        *("--find-links", "https://data.pyg.org/whl/torch-2.0.1.html"),
        *("--extra-index-url", "https://example.org/simple"),
        *("models/run_md.py", "--model", "mace_mp_0"),
    ]
    provided = calculators.runner_provided_packages()
    assert {"ase", "numpy", "pymatgen", "phonopy", "matbench-discovery"} <= provided
    assert "torch" not in provided


def test_project_backed_yaml_deps_match_install_project() -> None:
    """Provenance pins must match the authoritative project (no silent drift)."""
    for model_key, calc_spec in CALCULATORS.items():
        if not calc_spec.project:
            continue
        pyproject_path = f"{ROOT}/{calc_spec.project}/pyproject.toml"
        assert os.path.isfile(pyproject_path), model_key
        with open(pyproject_path, mode="rb") as file:
            pyproject = tomllib.load(file)
        project_deps = {
            *pyproject.get("project", {}).get("dependencies", []),
            *pyproject.get("tool", {}).get("uv", {}).get("override-dependencies", []),
        }
        missing = [dep for dep in calc_spec.deps if dep not in project_deps]
        assert not missing, (
            f"{model_key}: YAML environment.dependencies {missing} are not in "
            f"{calc_spec.project}; keep display pins identical to the install project"
        )


def test_checkpoint_specs_accept_declared_overrides() -> None:
    """Managed/required checkpoints are a required factory arg (no inner download)."""
    for model_key, calc_spec in CALCULATORS.items():
        if calc_spec.auto_checkpoint or calc_spec.requires_checkpoint:
            params = inspect.signature(calc_spec.make_calc).parameters
            assert "checkpoint" in params, model_key
            assert params["checkpoint"].default is inspect.Parameter.empty, model_key


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
    assert "orb-v3-conservative-inf-mpa" in Model["orb_v3"].metadata["checkpoint_url"]
    # TACE's OEQ/AOTI export needs torch>=2.11, so TACE factories may only accelerate
    # when their YAML pins a new enough torch (older pins must stay eager)
    tace_keys = [
        key
        for key in CALCULATORS
        if getattr(CALCULATORS[key].make_calc, "__qualname__", "").startswith("_tace.")
    ]
    assert len(tace_keys) >= 2, tace_keys
    for model_key in tace_keys:
        calc_spec = CALCULATORS[model_key]
        accelerate = inspect.getclosurevars(calc_spec.make_calc).nonlocals["accelerate"]
        torch_pin = next(
            Requirement(dep) for dep in calc_spec.deps if dep.startswith("torch==")
        )
        torch_version = Version(next(iter(torch_pin.specifier)).version)
        assert accelerate is (torch_version >= Version("2.11")), model_key


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


@pytest.mark.parametrize("device", ["cpu", "cuda", "cuda:1"])
def test_prophet_factory_honors_requested_device(
    device: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Prophet receives the requested device before allocating its model."""
    constructor = Mock()
    monkeypatch.setitem(
        sys.modules, "prophet", SimpleNamespace(KairosCalculator=constructor)
    )
    checkpoint = "model.pt"
    calc = CALCULATORS["prophet_oame_mbd"].make_calc(device, checkpoint=checkpoint)
    constructor.assert_called_once_with(
        model_path=checkpoint,
        use_kernel=device.startswith("cuda"),
        use_compile=False,
        device=device,
    )
    assert calc is constructor.return_value


def test_pet_factory_casts_exported_model(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """PET casts the loaded model and advertises the requested dtype."""
    captured: dict[str, object] = {}
    model_capabilities = SimpleNamespace(dtype="float32")
    float64_dtype = object()

    def cast_model(**kwargs: object) -> SimpleNamespace:
        """Capture the requested torch dtype and device."""
        captured.update(kwargs)
        return fake_model

    fake_model = SimpleNamespace(capabilities=lambda: model_capabilities, to=cast_model)

    torch_module = ModuleType("torch")
    torch_module.__dict__.update(float32=object(), float64=float64_dtype)
    metatomic_torch_module = ModuleType("metatomic.torch")
    metatomic_torch_module.__dict__["load_atomistic_model"] = lambda _path: fake_model
    calculator_module = ModuleType("metatomic.torch.ase_calculator")
    calculator_module.__dict__["MetatomicCalculator"] = lambda _model, **_kwargs: (
        SimpleNamespace()
    )
    for module_name, module in {
        "torch": torch_module,
        "metatomic": ModuleType("metatomic"),
        "metatomic.torch": metatomic_torch_module,
        "metatomic.torch.ase_calculator": calculator_module,
    }.items():
        monkeypatch.setitem(sys.modules, module_name, module)

    # Stand in for the separately tested export cache.
    monkeypatch.setattr(
        calculators, "_run_to_atomic_output", lambda *_args, **_kwargs: None
    )
    calculators._pet(  # noqa: SLF001
        "cpu", dtype="float64", checkpoint=str(tmp_path / "model.ckpt")
    )
    assert fake_model.capabilities().dtype == "float64"
    assert captured == {"dtype": float64_dtype, "device": "cpu"}


def test_tace_factory_acceleration_paths(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """TACE exports a cached AOTI package on CUDA; v1 and CPU stay eager."""
    calc_calls: list[tuple[str, dict[str, object]]] = []
    export_calls: list[tuple[list[str], str, dict[str, str]]] = []
    ase_module = ModuleType("tace.interface.ase")
    ase_module.__dict__["TACEAseCalc"] = lambda model, **kwargs: (
        calc_calls.append((model, kwargs)) or SimpleNamespace()
    )
    torch_module = ModuleType("torch")
    torch_module.__dict__["cuda"] = SimpleNamespace(
        get_device_capability=lambda _device: (8, 9)
    )
    for module_name, module in {
        "tace": ModuleType("tace"),
        "tace.interface": ModuleType("tace.interface"),
        "tace.interface.ase": ase_module,
        "torch": torch_module,
    }.items():
        monkeypatch.setitem(sys.modules, module_name, module)

    accel_flags = {"TACE_USE_OEQ": "1", "TACE_USE_COMPILE": "1"}
    for flag in accel_flags:
        monkeypatch.delenv(flag, raising=False)
    monkeypatch.setattr(calculators, "CHECKPOINT_DIR", str(tmp_path))
    monkeypatch.setattr(
        calculators,
        "_run_to_atomic_output",
        lambda command, dest, **kwargs: export_calls.append(
            (list(command), dest, dict(kwargs["env_overrides"]))
        ),
    )
    checkpoint = str(tmp_path / "model.pt")
    Path(checkpoint).write_bytes(b"checkpoint")

    calc = calculators._tace("tace_oam_l")("cuda", checkpoint=checkpoint)  # noqa: SLF001
    # the enabled acceleration paths are recorded as cost labels
    assert getattr(calc, calculators.ACCELERATION_ATTR) == (
        "openequivariance",
        "torch_aot_inductor",
    )

    compiled = f"{tmp_path}/tace_oam_l-cuda-sm89-oeq-aoti.pt2"
    command, dest, env_overrides = export_calls[0]
    assert command[:3] == [sys.executable, "-m", "tace.scripts.export_eval"]
    assert command[-4:] == ["--backend", "aoti", "--device", "cuda"]
    assert (dest, env_overrides) == (compiled, accel_flags)
    # TACE never unsets these, so leaking them would accelerate later eager calcs
    assert accel_flags.keys().isdisjoint(os.environ)
    accelerated = {"device": "cuda", "enable_oeq": True, "enable_compile": True}
    assert calc_calls == [(compiled, accelerated)]

    # a pre-compiled .pt2 and both eager paths (v1, CPU) skip the export
    packaged = str(tmp_path / "ready.pt2")
    Path(packaged).write_bytes(b"package")
    calculators._tace("tace_oam_l")("cuda", checkpoint=packaged)  # noqa: SLF001
    calculators._tace("tace_v1_oam_m", accelerate=False)(  # noqa: SLF001
        "cuda", checkpoint=checkpoint
    )
    eager = calculators._tace("tace_oam_l")("cpu", checkpoint=checkpoint)  # noqa: SLF001
    assert not hasattr(eager, calculators.ACCELERATION_ATTR)
    assert len(export_calls) == 1
    assert calc_calls[1:] == [
        (packaged, accelerated),
        (checkpoint, {"device": "cuda"}),
        (checkpoint, {"device": "cpu"}),
    ]


@pytest.mark.parametrize(
    "model_key", ["sevennet_0", "sevennet_l3i5", "sevennet_omni_i12"]
)
def test_named_sevennet_models_keep_package_managed_checkpoints(model_key: str) -> None:
    """Name-based SevenNet loaders must not be replaced by registry downloads."""
    assert not CALCULATORS[model_key].auto_checkpoint


@pytest.mark.parametrize(
    ("find_models", "model_file"),
    [
        (calculators._deepmd_model_paths, "nested/model.pth"),  # noqa: SLF001
        (calculators._alignn_model_dirs, "nested/best_model.pt"),  # noqa: SLF001
    ],
)
def test_archive_extraction_replaces_invalid_cache(
    tmp_path: Path, find_models: Callable[[str], list[str]], model_file: str
) -> None:
    """Extraction refuses folders it did not create, replaces its own stale caches,
    reuses a matching one, re-extracts a changed archive and rejects members escaping
    the extraction directory.
    """
    archive_path = tmp_path / "model.zip"

    def write_archive(weights: bytes) -> None:
        with zipfile.ZipFile(archive_path, mode="w") as archive:
            archive.writestr(model_file, weights)
            archive.writestr("nested/config.json", b"{}")

    write_archive(b"weights")
    extract_dir = tmp_path / "model-extracted"
    (extract_dir / "stale").mkdir(parents=True)
    for name in ("stale.pb", "stale/config.json", "stale/best_model.pt"):
        (extract_dir / name).write_bytes(b"stale")

    extract = calculators._extract_single_model  # noqa: SLF001
    for unmarked_dir in (extract_dir, tmp_path / "model-extracted.tmp"):
        (unmarked_dir / "user-file").parent.mkdir(exist_ok=True)
        (unmarked_dir / "user-file").write_bytes(b"keep")
        with pytest.raises(ValueError, match="Refusing to replace"):
            extract(str(archive_path), find_models)
        assert (unmarked_dir / "user-file").read_bytes() == b"keep"
        (unmarked_dir / "user-file").unlink()
        if unmarked_dir != extract_dir:
            unmarked_dir.rmdir()
    # an empty leftover temporary folder (run stopped before marking it) is cleared
    (tmp_path / "model-extracted.tmp").mkdir()
    (extract_dir / ".archive-stat").write_text("old archive")
    model_path = extract(str(archive_path), find_models)

    # finders join with "/" but the extractor returns OS-normalized paths
    assert [os.path.normpath(path) for path in find_models(str(extract_dir))] == [
        model_path
    ]
    assert (extract_dir / model_file).read_bytes() == b"weights"
    sentinel = extract_dir / "sentinel"
    sentinel.touch()
    assert extract(str(archive_path), find_models) == model_path
    assert sentinel.is_file()  # unchanged archive: cache hit, no re-extraction
    write_archive(b"updated weights")
    assert extract(str(archive_path), find_models) == model_path
    assert (extract_dir / model_file).read_bytes() == b"updated weights"
    assert not sentinel.exists()

    with zipfile.ZipFile(archive_path, mode="w") as archive:
        archive.writestr("../evil", b"escape")
    with pytest.raises(ValueError, match=r"Unsafe path '\.\./evil'"):
        extract(str(archive_path), find_models)
    assert not (tmp_path.parent / "evil").exists()


@pytest.mark.parametrize(
    ("model_key", "eager"), [("nequip_oam_l_0_1", False), ("nequip_oam_xl_0_1", True)]
)
def test_nequip_factory_load_path(
    model_key: str, eager: bool, monkeypatch: pytest.MonkeyPatch
) -> None:
    """NequIP/Allegro load AOTI .nequip.pt2s, except eager NequIP-OAM-XL (SIGFPEs)."""
    compile_stub = Mock()
    compiled_loader, eager_loader = Mock(return_value=EMT()), Mock(return_value=EMT())
    nequip_ase = ModuleType("nequip.ase")
    nequip_ase.NequIPCalculator = SimpleNamespace(  # ty: ignore[unresolved-attribute]
        from_compiled_model=compiled_loader, _from_saved_model=eager_loader
    )
    monkeypatch.setitem(sys.modules, "nequip", ModuleType("nequip"))
    monkeypatch.setitem(sys.modules, "nequip.ase", nequip_ase)
    torch = SimpleNamespace(
        device=lambda name: SimpleNamespace(type=name.split(":")[0]),
        cuda=SimpleNamespace(get_device_capability=lambda _device: (9, 0)),
    )
    monkeypatch.setitem(sys.modules, "torch", torch)
    monkeypatch.setattr(calculators, "_run_to_atomic_output", compile_stub)

    calc = CALCULATORS[model_key].make_calc("cuda")
    if eager:
        compile_stub.assert_not_called()
        (registry,), kwargs = eager_loader.call_args
        assert registry.startswith("nequip.net:mir-group/")
        assert kwargs == {"device": "cuda", "compile_mode": "eager"}
        assert not hasattr(calc, calculators.ACCELERATION_ATTR)
        return
    (command, output_path), compile_kwargs = compile_stub.call_args
    assert command[command.index("--mode") + 1] == "aotinductor"
    # AOTI kernels are arch- and torch-specific, so both key the cached package
    assert output_path.endswith(f"{model_key}-cuda-sm90.nequip.pt2")
    assert compile_kwargs["tool_packages"] == ("nequip", "torch")
    compiled_loader.assert_called_once_with(compile_path=output_path, device="cuda")
    assert getattr(calc, calculators.ACCELERATION_ATTR) == ("torch_aot_inductor",)


def test_dpa_4_0_1_loads_eager_checkpoint(monkeypatch: pytest.MonkeyPatch) -> None:
    """DPA-4.0.1 loads its .pt eagerly (the frozen .pt2 breaks on isolated dimers)."""
    deepmd_calculator = ModuleType("deepmd.calculator")
    deepmd_calculator.DP = Mock(return_value=EMT())  # ty: ignore[unresolved-attribute]
    monkeypatch.setitem(sys.modules, "deepmd", ModuleType("deepmd"))
    monkeypatch.setitem(sys.modules, "deepmd.calculator", deepmd_calculator)

    spec = CALCULATORS["dpa_4_0_1_pro_mptrj"]
    assert spec.make_calc is calculators._deepmd  # noqa: SLF001
    assert isinstance(spec.make_calc("cuda", checkpoint="/ckpt/dpa.pt"), EMT)
    deepmd_calculator.DP.assert_called_once_with("/ckpt/dpa.pt")  # no freeze, no .pt2


def test_atomic_command_output_recovers_after_timeout(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Timed-out compilers leave no accepted partial cache and can retry."""
    destination = tmp_path / "fresh" / "compiled.nequip.pth"
    n_calls = 0
    seen_envs: list[dict[str, str] | None] = []

    def run_stub(
        command: list[str],
        *,
        check: bool,
        timeout: float,
        env: dict[str, str] | None,
    ) -> None:
        """Write a partial first result and a complete retry result."""
        nonlocal n_calls
        assert check
        assert timeout == calculators.DERIVED_ARTIFACT_TIMEOUT_SEC
        seen_envs.append(env)
        n_calls += 1
        output_path = next(argument for argument in command if ".tmp." in argument)
        assert output_path.endswith(".nequip.pth")
        assert Path(output_path).parent == destination.parent
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

    # uv gives each run a fresh interpreter path, which must not bust the cache
    n_before = n_calls
    for python in ("/uv/run-a/bin/python", "/uv/run-b/bin/python"):
        monkeypatch.setattr(calculators.sys, "executable", python)
        calculators._run_to_atomic_output(  # noqa: SLF001
            [python, "-m", "exporter", "{output}"], str(destination)
        )
    assert n_calls == n_before + 1

    # env_overrides reach only the child process, merged over the inherited env
    assert all(env is None for env in seen_envs)
    source.write_bytes(b"version 3")
    calculators._run_to_atomic_output(  # noqa: SLF001
        command, str(destination), source_paths=(str(source),), env_overrides={"X": "1"}
    )
    assert seen_envs[-1] == {**os.environ, "X": "1"}
    assert "X" not in os.environ


@pytest.mark.parametrize(
    ("url", "resolved_url"),
    [
        ("https://example.com/model.ckpt", "https://example.com/model.ckpt"),
        (
            "https://huggingface.co/org/model/blob/main/model.ckpt",
            "https://huggingface.co/org/model/resolve/main/model.ckpt",
        ),
        (
            "https://huggingface.co.example.org/blob/model.ckpt",
            "https://huggingface.co.example.org/blob/model.ckpt",
        ),
        (
            "https://github.com/org/model/blob/main/model.ckpt",
            "https://github.com/org/model/raw/main/model.ckpt",
        ),
        (
            "https://uni.sciebo.de/s/token?key=abc",
            "https://uni.sciebo.de/s/token/download?key=abc",
        ),
    ],
)
def test_download_checkpoint_replaces_zero_byte_cache(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    url: str,
    resolved_url: str,
) -> None:
    """download_checkpoint redownloads zero-byte cached checkpoint files."""
    from matbench_discovery.remote import fetch

    monkeypatch.setenv("HF_TOKEN", "hf_test")
    monkeypatch.setattr(calculators, "CHECKPOINT_DIR", f"{tmp_path}")
    from_ref_mock = classmethod(
        lambda _model_cls, _model_key: SimpleNamespace(metadata={"checkpoint_url": url})
    )
    monkeypatch.setattr(Model, "from_ref", from_ref_mock)

    def download_file_mock(destination_path: str, download_url: str) -> None:
        """Verify normalized URLs; the shared downloader owns authentication."""
        assert download_url == resolved_url
        Path(destination_path).write_bytes(b"checkpoint")

    monkeypatch.setattr(fetch, "download_file", download_file_mock)
    url_hash = hashlib.sha256(resolved_url.encode()).hexdigest()[:12]
    dest = tmp_path / f"fake-{url_hash}.ckpt"
    dest.write_bytes(b"")

    actual_path = os.path.normpath(calculators.download_checkpoint("fake", ext=".ckpt"))
    assert actual_path == os.path.normpath(f"{dest}")
    assert dest.read_bytes() == b"checkpoint"


@pytest.mark.parametrize(
    ("model_key", "head", "checkpoint"),
    [
        ("mace_soyo_mptrj", "mptraj", "model.pth"),
        ("mace_soyo_oam", "salexmp", "model.pth"),
        ("mace_soyo_oam", "salexmp", "exported.pt2"),
    ],
)
def test_mace_soyo_exports_per_gpu_arch(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    model_key: str,
    head: str,
    checkpoint: str,
) -> None:
    """A .pth checkpoint exports to a per-arch AOTI package; a .pt2 loads as is."""
    calc_cls, exporter = Mock(return_value=SimpleNamespace()), Mock()
    cuda = SimpleNamespace(get_device_capability=lambda _device: (9, 0))
    torch = SimpleNamespace(
        device=lambda name: SimpleNamespace(type=name.split(":")[0]), cuda=cuda
    )
    monkeypatch.setitem(sys.modules, "torch", torch)
    monkeypatch.setitem(
        sys.modules, "mace_soyo", SimpleNamespace(MACESoyoCalculator=calc_cls)
    )
    monkeypatch.setattr(calculators, "CHECKPOINT_DIR", f"{tmp_path}/fresh-cache")
    monkeypatch.setattr(calculators, "_run_to_atomic_output", exporter)
    path = f"{tmp_path}/{checkpoint}"
    with open(path, mode="wb") as file:
        file.write(b"weights")
    make_calc = CALCULATORS[model_key].make_calc

    make_calc("cuda", path)
    kwargs = calc_cls.call_args.kwargs
    assert (kwargs["head"], kwargs["device"], kwargs["use_d3"]) == (head, "cuda", False)
    if checkpoint.endswith(".pt2"):
        exporter.assert_not_called()
        assert kwargs["package_path"] == path
    else:
        command, package = exporter.call_args.args
        assert package == kwargs["package_path"]
        digest = hashlib.sha256(b"weights").hexdigest()[:16]
        assert package.endswith(f"-{digest}-sm90-fp32.pt2")
        assert command[command.index("--ckpt-path") + 1] == path
        assert os.path.isfile(command[command.index("--structure-file") + 1])
    with pytest.raises(ValueError, match="requires CUDA"):
        make_calc("cpu", path)
