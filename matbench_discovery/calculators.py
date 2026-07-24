"""Registry of MLIP ASE calculators shared across benchmark tasks.

Runners resolve each model's YAML environment via CalcSpec.uv_run_cmd
(uv run --with or --project for override-heavy stacks like EquFlash).
Factories import packages lazily so listing models needs only core deps.

Registry keys are Model enum names so metrics write to the right YAML.
"""

import hashlib
import inspect
import os
import shutil
import subprocess
import zipfile
from collections.abc import Callable, Iterator, Sequence
from dataclasses import dataclass
from importlib.metadata import PackageNotFoundError, version
from typing import TYPE_CHECKING

from filelock import FileLock

from matbench_discovery import DEFAULT_CACHE_DIR

if TYPE_CHECKING:
    import argparse
    from collections.abc import Mapping

    from ase.calculators.calculator import Calculator

CHECKPOINT_DIR = f"{DEFAULT_CACHE_DIR}/md-checkpoints"
DERIVED_ARTIFACT_TIMEOUT_SEC = 30 * 60


def _is_non_empty_file(path: str) -> bool:
    """Return whether path exists and has non-zero size."""
    return os.path.isfile(path) and os.path.getsize(path) > 0


def _file_sha256(path: str) -> str:
    """Return a file's SHA-256 digest."""
    digest = hashlib.sha256()
    with open(path, mode="rb") as file:
        for chunk in iter(lambda: file.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _file_state(path: str) -> tuple[int, int, int]:
    """Return cheap fields that change whenever normal file writes occur."""
    file_stat = os.stat(path)
    return file_stat.st_size, file_stat.st_mtime_ns, file_stat.st_ctime_ns


def _stable_file_sha256(path: str) -> tuple[str, tuple[int, int, int]]:
    """Hash a file outside locks and reject concurrent mutation."""
    state = _file_state(path)
    digest = _file_sha256(path)
    if _file_state(path) != state:
        raise RuntimeError(f"File changed while hashing: {path}")
    return digest, state


def _cached_file_sha256(path: str) -> tuple[str, tuple[int, int, int]] | None:
    """Hash an existing cache entry, tolerating a concurrent atomic replacement."""
    try:
        return _stable_file_sha256(path) if _is_non_empty_file(path) else None
    except (FileNotFoundError, RuntimeError):
        return None


def download_checkpoint(model_key: str, ext: str | None = None) -> str:
    """Download a model's weights from its YAML checkpoint_url to a local cache
    and return the path. Normalizes HuggingFace /blob/ to /resolve/ and
    sciebo share links to direct downloads; figshare URLs are handled by
    download_file. Cached after the first call.

    Args:
        model_key: Model enum name whose YAML carries checkpoint_url.
        ext: Force this file extension (e.g. '.pth' for deepmd, whose loader picks
            its backend by suffix). Defaults to the extension parsed from the URL.
    """
    from matbench_discovery.enums import Model
    from matbench_discovery.remote.fetch import download_file

    url = Model.from_ref(model_key).metadata["checkpoint_url"]
    if not url:
        raise ValueError(f"{model_key} has no checkpoint_url in its YAML")
    headers = None
    if "huggingface.co" in url:
        url = url.replace("/blob/", "/resolve/")
        # gated repos (e.g. fairchem OMAT24) need a bearer token + license acceptance
        token = os.getenv("HF_TOKEN") or os.getenv("HUGGING_FACE_HUB_TOKEN")
        if token:
            headers = {"Authorization": f"Bearer {token}"}
    if "sciebo" in url and not url.endswith("/download"):
        url = f"{url}/download"
    if "github.com" in url and "/blob/" in url:  # serve the raw file, not the HTML page
        url = url.replace("/blob/", "/raw/")

    ext = ext or os.path.splitext(url.split("?")[0])[1] or ".ckpt"
    url_hash = hashlib.sha256(url.encode()).hexdigest()[:12]
    dest = f"{CHECKPOINT_DIR}/{model_key}-{url_hash}{ext}"
    os.makedirs(CHECKPOINT_DIR, exist_ok=True)
    # serialize concurrent same-model downloads: parallel array tasks otherwise race on
    # download_file's temp-file rename, leaving losers with a vanished .part file
    with FileLock(f"{dest}.lock"):
        if os.path.isfile(dest) and not _is_non_empty_file(dest):
            os.remove(dest)
        if not _is_non_empty_file(dest):
            download_file(dest, url, headers=headers)
    if not _is_non_empty_file(dest):
        raise RuntimeError(
            f"Failed to download {model_key} checkpoint from {url}. If the repo is "
            "gated (e.g. fairchem OMAT24), accept its license on HuggingFace and set "
            "HF_TOKEN in the environment."
        )
    return dest


def _stage_checkpoint(
    model_key: str,
    dest: str,
    *,
    ext: str | None = None,
    source_path: str | None = None,
) -> None:
    """Download a checkpoint and atomically stage it in a framework-specific cache.

    This bypasses framework downloaders that a WAF may serve as empty files.
    """
    source = source_path or download_checkpoint(model_key, ext=ext)
    source_hash, source_state = _stable_file_sha256(source)
    cached_dest = _cached_file_sha256(dest)
    os.makedirs(os.path.dirname(dest) or ".", exist_ok=True)
    with FileLock(f"{dest}.lock"):
        if _file_state(source) != source_state:
            raise RuntimeError(f"Source changed while staging {dest}")
        if (
            cached_dest is not None
            and _is_non_empty_file(dest)
            and _file_state(dest) == cached_dest[1]
            and cached_dest[0] == source_hash
        ):
            return
        tmp_dest = f"{dest}.tmp"
        shutil.copy(source, tmp_dest)
        if _file_state(source) != source_state:
            os.remove(tmp_dest)
            raise RuntimeError(f"Source changed while staging {dest}")
        os.replace(tmp_dest, dest)


def _run_to_atomic_output(
    command: Sequence[str],
    dest: str,
    *,
    source_paths: Sequence[str] = (),
    tool_packages: Sequence[str] = (),
) -> None:
    """Atomically cache output by command, source contents, and tool versions."""
    if command.count("{output}") != 1:
        raise ValueError("Atomic output command must contain one {output} placeholder")
    os.makedirs(os.path.dirname(dest) or ".", exist_ok=True)
    dest_base, dest_ext = os.path.splitext(dest)
    tmp_dest = f"{dest_base}.tmp{dest_ext}"
    identity_path = f"{dest}.sha256"
    source_identities = tuple(
        (path, *_stable_file_sha256(path)) for path in source_paths
    )
    tool_versions = []
    for package in tool_packages:
        try:
            tool_versions.append((package, version(package)))
        except PackageNotFoundError:
            tool_versions.append((package, "unknown"))
    recipe = (
        tuple(command),
        tuple(digest for _path, digest, _state in source_identities),
        tuple(tool_versions),
    )
    recipe_hash = hashlib.sha256(repr(recipe).encode()).hexdigest()
    cached_dest = _cached_file_sha256(dest)
    with FileLock(f"{dest}.lock"):
        if any(
            _file_state(path) != state for path, _digest, state in source_identities
        ):
            raise RuntimeError(f"Source changed while creating {dest}")
        if (
            cached_dest is not None
            and _is_non_empty_file(dest)
            and _file_state(dest) == cached_dest[1]
            and os.path.isfile(identity_path)
        ):
            with open(identity_path, encoding="utf-8") as file:
                if file.read() == f"{recipe_hash}\n{cached_dest[0]}":
                    return
        if os.path.isfile(tmp_dest):
            os.remove(tmp_dest)
        try:
            rendered_command = [
                tmp_dest if argument == "{output}" else argument for argument in command
            ]
            subprocess.run(
                rendered_command, check=True, timeout=DERIVED_ARTIFACT_TIMEOUT_SEC
            )
            if not _is_non_empty_file(tmp_dest):
                command_text = " ".join(rendered_command)
                raise RuntimeError(f"{command_text=} wrote no output to {tmp_dest}")
            if any(
                _file_state(path) != state for path, _digest, state in source_identities
            ):
                raise RuntimeError(f"Source changed while creating {dest}")
            output_hash = _file_sha256(tmp_dest)
            os.replace(tmp_dest, dest)
            with open(identity_path, mode="w", encoding="utf-8") as file:
                file.write(f"{recipe_hash}\n{output_hash}")
        finally:
            if os.path.isfile(tmp_dest):
                os.remove(tmp_dest)


@dataclass(frozen=True)
class CalcSpec:
    """A registered MLIP calculator and its YAML-owned runtime environment."""

    make_calc: Callable[..., "Calculator"]
    deps: tuple[str, ...] = ()  # extra uv requirements beyond CORE_DEPS
    find_links: tuple[str, ...] = ()  # uv --find-links (e.g. PyG/dgl wheel pages)
    extra_index_url: tuple[str, ...] = ()  # uv --extra-index-url entries
    # pin uv's Python (e.g. HIENet's torch 2.1.2 has no cp312 wheels, so needs 3.11);
    # None lets uv pick the default
    python_version: str | None = None
    project: str | None = None
    requires_checkpoint: bool = False
    auto_checkpoint: bool = False
    checkpoint_ext: str | None = None

    def uv_run_cmd(self, script: str, *args: str) -> list[str]:
        """Uv run command that resolves this model's env and runs the script.

        With project, that pyproject owns resolution (incl. override-dependencies);
        YAML dependencies stay for display only. Otherwise deps / indexes apply.
        """
        py_args = ["--python", self.python_version] if self.python_version else []
        if self.project:
            mode = ["--project", self.project, "--isolated", *py_args, "python"]
        else:
            pairs = (
                *(("--with", dep) for dep in self.deps),
                *(("--find-links", url) for url in self.find_links),
                *(("--extra-index-url", url) for url in self.extra_index_url),
            )
            mode = ["--no-project", *py_args, *(tok for pair in pairs for tok in pair)]
        return ["uv", "run", *mode, script, *args]


def _runtime_calc_spec(
    model_key: str,
    make_calc: Callable[..., "Calculator"],
    *,
    requires_checkpoint: bool = False,
    auto_checkpoint: bool = False,
    checkpoint_ext: str | None = None,
) -> Callable[[], CalcSpec]:
    """Lazy CalcSpec factory that reads environment from the model YAML.

    model_key must be a Model enum name with an environment block.
    Debug-only keys like emt (no YAML) must register an eager CalcSpec.
    """

    def build() -> CalcSpec:
        from matbench_discovery.enums import Model

        environment = Model[model_key].metadata["environment"]
        return CalcSpec(
            make_calc,
            deps=tuple(environment["dependencies"]),
            find_links=tuple(environment.get("find_links", ())),
            extra_index_url=tuple(environment.get("extra_index_urls", ())),
            python_version=environment["python_version"],
            project=environment.get("project"),
            requires_checkpoint=requires_checkpoint,
            auto_checkpoint=auto_checkpoint,
            checkpoint_ext=checkpoint_ext,
        )

    return build


def _checkpoint_spec(
    model_key: str,
    make_calc: Callable[..., "Calculator"],
    *,
    ext: str | None = None,
    requires_checkpoint: bool = False,
) -> Callable[[], CalcSpec]:
    """Lazy YAML-backed spec with auto-download or required-checkpoint flags."""
    return _runtime_calc_spec(
        model_key,
        make_calc,
        requires_checkpoint=requires_checkpoint,
        auto_checkpoint=not requires_checkpoint,
        checkpoint_ext=ext,
    )


def _named_spec(
    factory: Callable[[str], Callable[..., "Calculator"]],
    model_key: str,
    *,
    checkpoint: bool = False,
    ext: str | None = None,
    requires_checkpoint: bool = False,
) -> Callable[[], CalcSpec]:
    """Register factory(model_key) under the same YAML-backed key."""
    if checkpoint:
        return _checkpoint_spec(
            model_key,
            factory(model_key),
            ext=ext,
            requires_checkpoint=requires_checkpoint,
        )
    return _runtime_calc_spec(model_key, factory(model_key))


class _CalcRegistry:
    """Resolve lazy CalcSpec builders on first [] / .get / iteration."""

    def __init__(self, mapping: dict[str, CalcSpec | Callable[[], CalcSpec]]) -> None:
        self._data = dict(mapping)

    def __getitem__(self, key: str) -> CalcSpec:
        value = self._data[key]
        if not isinstance(value, CalcSpec):
            self._data[key] = value = value()
        return value

    def __setitem__(self, key: str, value: CalcSpec | Callable[[], CalcSpec]) -> None:
        self._data[key] = value

    def __delitem__(self, key: str) -> None:
        del self._data[key]

    def __contains__(self, key: object) -> bool:
        return key in self._data

    def __iter__(self) -> Iterator[str]:
        return iter(self._data)

    def get(self, key: str, default: CalcSpec | None = None) -> CalcSpec | None:
        """Return a resolved CalcSpec, or default when the key is missing."""
        if key not in self._data:
            return default
        return self[key]

    def items(self) -> Iterator[tuple[str, CalcSpec]]:
        return ((key, self[key]) for key in self._data)


def resolve_checkpoint(model_key: str, checkpoint: str | None = None) -> str | None:
    """Resolve and validate an explicit or registry-managed checkpoint file."""
    spec = CALCULATORS[resolve_calculator_key(model_key)]
    if checkpoint is not None:
        checkpoint = os.path.abspath(checkpoint)
        if not _is_non_empty_file(checkpoint):
            raise FileNotFoundError(f"Checkpoint file not found or empty: {checkpoint}")
        return checkpoint
    if spec.requires_checkpoint:
        raise ValueError(f"{model_key} requires an explicit checkpoint")
    if spec.auto_checkpoint:
        return download_checkpoint(model_key, ext=spec.checkpoint_ext)
    return None


def _detect_device() -> str:
    """Return 'cuda' if a GPU is visible to torch, else 'cpu'. Non-torch backends
    (TensorFlow GRACE, JAX nequix) have no torch installed and ignore the device, so
    fall back to 'cpu' rather than failing on the import.
    """
    try:
        import torch
    except ImportError:
        return "cpu"
    return "cuda" if torch.cuda.is_available() else "cpu"


def resolve_device(model_key: str, device: str | None = None) -> str:
    """Resolve an explicit or auto-detected calculator device."""
    if device is not None:
        return device
    return "cpu" if model_key == "emt" else _detect_device()


def _mace(
    checkpoint: str, head: str | None = None
) -> Callable[[str, str], "Calculator"]:
    def make_calc(device: str, dtype: str = "float64") -> "Calculator":
        from mace.calculators import mace_mp

        kwargs = {"head": head} if head else {}
        return mace_mp(
            model=checkpoint,
            device=device,
            default_dtype=dtype,
            enable_cueq=device == "cuda",
            **kwargs,
        )

    return make_calc


def _orb(variant: str) -> Callable[[str], "Calculator"]:
    def make_calc(device: str) -> "Calculator":
        from orb_models.forcefield.calculator import ORBCalculator
        from orb_models.forcefield.pretrained import ORB_PRETRAINED_MODELS

        model = ORB_PRETRAINED_MODELS[variant]()
        model.to(device)
        return ORBCalculator(model, device=device)

    return make_calc


def _mattersim(checkpoint: str) -> Callable[[str], "Calculator"]:
    def make_calc(device: str) -> "Calculator":
        from mattersim.forcefield import MatterSimCalculator

        return MatterSimCalculator(load_path=checkpoint, device=device)

    return make_calc


def _sevennet(
    model: str, modal: str | None = None, *, from_url: bool = False
) -> Callable[..., "Calculator"]:
    def make_calc(device: str, checkpoint: str | None = None) -> "Calculator":
        """Construct SevenNet, locking package downloads or resolving a YAML URL."""
        from sevenn.calculator import SevenNetCalculator

        kwargs = {"modal": modal} if modal else {}
        if from_url:
            checkpoint = checkpoint or download_checkpoint(model)
            return SevenNetCalculator(model=checkpoint, device=device, **kwargs)
        os.makedirs(CHECKPOINT_DIR, exist_ok=True)
        with FileLock(f"{CHECKPOINT_DIR}/sevennet-{model}.lock"):
            return SevenNetCalculator(
                model=checkpoint or model, device=device, **kwargs
            )

    return make_calc


def _grace(model_name: str) -> Callable[[str], "Calculator"]:
    def make_calc(device: str) -> "Calculator":  # noqa: ARG001 - TF picks the device
        from tensorpotential.calculator import grace_fm

        return grace_fm(model_name)

    return make_calc


def _fairchem(model_key: str) -> Callable[..., "Calculator"]:
    def make_calc(device: str, checkpoint: str | None = None) -> "Calculator":
        from fairchem.core import OCPCalculator

        if checkpoint is None and model_key.startswith("equiformer_v3"):
            raise ValueError(
                f"{model_key} requires --checkpoint because its YAML links to the "
                "checkpoint repository, not a model artifact"
            )
        checkpoint = checkpoint or download_checkpoint(model_key)
        return OCPCalculator(checkpoint_path=checkpoint, cpu=device == "cpu", seed=0)

    return make_calc


def _deepmd(model_key: str) -> Callable[..., "Calculator"]:
    def make_calc(
        device: str,  # noqa: ARG001 - DP picks the device
        checkpoint: str | None = None,
    ) -> "Calculator":
        from deepmd.calculator import DP

        # DP selects its backend by file suffix, so force '.pth' (figshare URLs are
        # extensionless and would otherwise save as '.ckpt')
        return DP(checkpoint or download_checkpoint(model_key, ext=".pth"))

    return make_calc


def _deepmd_model_paths(extract_dir: str) -> list[str]:
    """Return frozen DeePMD artifacts recursively contained in a directory."""
    return [
        f"{root_dir}/{file_name}"
        for root_dir, _dir_names, file_names in os.walk(extract_dir)
        for file_name in file_names
        if os.path.splitext(file_name)[1] in {".pb", ".pth", ".pt2"}
    ]


def _extract_single_deepmd_model(archive_path: str) -> str:
    """Safely and atomically extract one frozen model from a DeePMD ZIP."""
    extract_dir = f"{os.path.splitext(archive_path)[0]}-extracted"
    tmp_extract_dir = f"{extract_dir}.tmp"
    archive_identity = ":".join(map(str, _file_state(archive_path)))
    identity_path = f"{extract_dir}/.archive-stat"
    with FileLock(f"{extract_dir}.lock"):
        existing_paths = _deepmd_model_paths(extract_dir)
        if len(existing_paths) == 1 and os.path.isfile(identity_path):
            with open(identity_path, encoding="utf-8") as file:
                if file.read() == archive_identity:
                    return existing_paths[0]

        shutil.rmtree(tmp_extract_dir, ignore_errors=True)
        os.makedirs(tmp_extract_dir)
        try:
            with zipfile.ZipFile(archive_path) as archive:
                real_extract_dir = os.path.realpath(tmp_extract_dir)
                for member in archive.infolist():
                    target_path = os.path.realpath(
                        f"{tmp_extract_dir}/{member.filename}"
                    )
                    common_path = os.path.commonpath((real_extract_dir, target_path))
                    if common_path != real_extract_dir:
                        raise ValueError(
                            f"Unsafe path {member.filename!r} in {archive_path}"
                        )
                archive.extractall(tmp_extract_dir)

            model_paths = _deepmd_model_paths(tmp_extract_dir)
            if len(model_paths) != 1:
                raise ValueError(
                    f"Expected one frozen DeePMD model in {archive_path}, "
                    f"got {model_paths}"
                )
            relative_model_path = os.path.relpath(model_paths[0], tmp_extract_dir)
            with open(
                f"{tmp_extract_dir}/.archive-stat", mode="w", encoding="utf-8"
            ) as file:
                file.write(archive_identity)
            shutil.rmtree(extract_dir, ignore_errors=True)
            os.replace(tmp_extract_dir, extract_dir)
            return f"{extract_dir}/{relative_model_path}"
        finally:
            shutil.rmtree(tmp_extract_dir, ignore_errors=True)


def _deepmd_archive(model_key: str) -> Callable[..., "Calculator"]:
    """Build a DeePMD calculator from a downloadable model archive."""

    def make_calc(
        device: str,  # noqa: ARG001 - DP picks the device
        checkpoint: str | None = None,
    ) -> "Calculator":
        """Extract the model archive when needed and load its frozen artifact."""
        from deepmd.calculator import DP

        checkpoint = checkpoint or download_checkpoint(model_key, ext=".zip")
        if zipfile.is_zipfile(checkpoint):
            checkpoint = _extract_single_deepmd_model(checkpoint)
        return DP(checkpoint)

    return make_calc


def _deepmd_freeze(model_key: str) -> Callable[..., "Calculator"]:
    def make_calc(
        device: str,  # noqa: ARG001 - DP picks the device
        checkpoint: str | None = None,
    ) -> "Calculator":
        from deepmd.calculator import DP

        # figshare ships a training checkpoint (state dict), so freeze it to the
        # torch-export .pt2 artifact that DP(frozen) loads.
        ckpt = checkpoint or download_checkpoint(model_key, ext=".pt")
        if ckpt.endswith(".pt2"):
            return DP(ckpt)
        # deepmd 3.2's pt backend freezes to a torch-export .pt2 (not TorchScript .pth)
        frozen = f"{os.path.splitext(ckpt)[0]}-frozen.pt2"
        _run_to_atomic_output(
            ["dp", "--pt", "freeze", "-c", ckpt, "-o", "{output}"],
            frozen,
            source_paths=(ckpt,),
            tool_packages=("deepmd-kit",),
        )
        return DP(frozen)

    return make_calc


def _tace(model_key: str) -> Callable[..., "Calculator"]:
    def make_calc(device: str, checkpoint: str | None = None) -> "Calculator":
        import sys

        from tace.interface.ase import TACEAseCalc

        checkpoint = checkpoint or download_checkpoint(model_key)
        if model_key == "tace_v1_oam_m":
            return TACEAseCalc(checkpoint, device=device)

        # Configure acceleration before loading the source checkpoint. The export
        # subprocess inherits these variables, so its graph contains OEQ operators
        # and uses TACE's compile-compatible model path.
        os.environ["TACE_USE_OEQ"] = "1"
        os.environ["TACE_USE_COMPILE"] = "1"

        # A user-supplied .pt2 is already deployable. Otherwise export from scratch.
        # Cache the package outside benchmark timing and rebuild it when the source
        # or toolchain changes.
        if checkpoint.endswith(".pt2"):
            compiled = checkpoint
        else:
            import torch

            torch_device = torch.device(device)
            device_tag = torch_device.type
            if torch_device.type == "cuda":
                major, minor = torch.cuda.get_device_capability(torch_device)
                device_tag = f"cuda-sm{major}{minor}"
            compiled = f"{CHECKPOINT_DIR}/{model_key}-{device_tag}-oeq-aoti.pt2"
            _run_to_atomic_output(
                [
                    sys.executable,
                    "-m",
                    "tace.scripts.export_eval",
                    "-m",
                    checkpoint,
                    "-o",
                    "{output}",
                    "--backend",
                    "aoti",
                    "--device",
                    device,
                ],
                compiled,
                source_paths=(checkpoint,),
                tool_packages=("tace", "torch", "openequivariance"),
            )
        return TACEAseCalc(
            compiled,
            device=device,
            enable_oeq=True,
            enable_compile=True,
        )

    return make_calc


def _hienet(model_key: str) -> Callable[..., "Calculator"]:
    def make_calc(device: str, checkpoint: str | None = None) -> "Calculator":
        from hienet.hienet_calculator import HIENetCalculator

        checkpoint = checkpoint or download_checkpoint(model_key)
        return HIENetCalculator(model=checkpoint, device=device)

    return make_calc


def _nequip(model_key: str) -> Callable[[str], "Calculator"]:
    def make_calc(device: str) -> "Calculator":
        from nequip.ase import NequIPCalculator

        from matbench_discovery.enums import Model

        # checkpoint_url is a nequip.net registry model; nequip-compile fetches it and
        # compiles it (needs a GPU) to a .nequip.pth that from_compiled_model loads.
        # torchscript not aotinductor: aotinductor's missing c-shim for
        # aten._linalg_det incorrectly compiles atomic_energy on the OAM models.
        # Cache + lock so it builds once.
        url = Model.from_ref(model_key).metadata["checkpoint_url"]
        registry = f"nequip.net:{url.split('nequip.net/models/')[-1]}"
        compiled = f"{CHECKPOINT_DIR}/{model_key}-{device}.nequip.pth"
        _run_to_atomic_output(
            [
                "nequip-compile",
                registry,
                "{output}",
                "--mode",
                "torchscript",
                "--device",
                device,
                "--target",
                "ase",
            ],
            compiled,
            tool_packages=("nequip",),
        )
        return NequIPCalculator.from_compiled_model(
            compile_path=compiled, device=device
        )

    return make_calc


def _eqnorm(model_key: str, model_variant: str) -> Callable[..., "Calculator"]:
    def make_calc(device: str, checkpoint: str | None = None) -> "Calculator":
        from eqnorm.calculator import EqnormCalculator

        # EqnormCalculator downloads via the wget package, which figshare's WAF serves
        # as 0 bytes (-> torch.load EOFError), so stage our checkpoint where it looks
        dest = os.path.expanduser(f"~/.cache/eqnorm/{model_variant}.pt")
        _stage_checkpoint(model_key, dest, ext=".pt", source_path=checkpoint)
        return EqnormCalculator(
            model_name="eqnorm", model_variant=model_variant, device=device
        )

    return make_calc


def _nequix(model_key: str) -> Callable[..., "Calculator"]:
    def make_calc(
        device: str,  # noqa: ARG001 - jax picks the device
        checkpoint: str | None = None,
    ) -> "Calculator":
        from nequix.calculator import NequixCalculator

        # NequixCalculator auto-downloads by name via figshare (WAF can serve 0 bytes);
        # stage our checkpoint and pass model_path so it loads ours directly.
        # use_kernel=False avoids the openequivariance extension (needs a separate pip)
        dest = os.path.expanduser(f"~/.cache/nequix/{model_key}.nqx")
        _stage_checkpoint(model_key, dest, ext=".nqx", source_path=checkpoint)
        return NequixCalculator(model_path=dest, backend="jax", use_kernel=False)

    return make_calc


def _matris(model: str, cache_name: str) -> Callable[..., "Calculator"]:
    def make_calc(device: str, checkpoint: str | None = None) -> "Calculator":
        from matris.applications import MatRISCalculator

        # MatRIS.load() only takes a registered model name and resolves it to
        # ~/.cache/matris/<cache_name> (its built-in figshare downloader serves 0 bytes
        # behind a WAF), so stage our YAML checkpoint there before constructing the calc
        _stage_checkpoint(
            model,
            os.path.expanduser(f"~/.cache/matris/{cache_name}"),
            source_path=checkpoint,
        )
        return MatRISCalculator(model=model, device=device)

    return make_calc


def _alphanet(model_key: str, config_url: str) -> Callable[..., "Calculator"]:
    def make_calc(
        device: str, dtype: str = "float64", checkpoint: str | None = None
    ) -> "Calculator":
        from alphanet.config import All_Config
        from alphanet.infer.calc import AlphaNetCalculator

        from matbench_discovery.remote.fetch import download_file

        # AlphaNet needs an architecture config json (not bundled in the weights) that
        # matches the checkpoint; use the model-specific commit-pinned upstream config
        config_hash = hashlib.sha256(config_url.encode()).hexdigest()[:12]
        config_path = f"{CHECKPOINT_DIR}/{model_key}-config-{config_hash}.json"
        os.makedirs(CHECKPOINT_DIR, exist_ok=True)
        with FileLock(f"{config_path}.lock"):
            if not _is_non_empty_file(config_path):
                download_file(config_path, config_url)
        config = All_Config().from_json(config_path)
        return AlphaNetCalculator(
            ckpt_path=checkpoint or download_checkpoint(model_key),
            device=device,
            precision={"float32": "32", "float64": "64"}[dtype],
            config=config,
        )

    return make_calc


def _pet(model_key: str) -> Callable[..., "Calculator"]:
    def make_calc(
        device: str, dtype: str = "float64", checkpoint: str | None = None
    ) -> "Calculator":
        import torch
        from metatomic.torch import load_atomistic_model
        from metatomic.torch.ase_calculator import MetatomicCalculator

        # PET ships a metatrain .ckpt that must be exported to a TorchScript .pt before
        # metatomic can load it; cache the export and lock it against parallel tasks
        ckpt = checkpoint or download_checkpoint(model_key)
        pt_file = f"{os.path.splitext(ckpt)[0]}.pt"
        # The temp name must end in .pt: mtt appends .pt otherwise, writing elsewhere.
        _run_to_atomic_output(
            ["mtt", "export", ckpt, "-o", "{output}"],
            pt_file,
            source_paths=(ckpt,),
            tool_packages=("metatrain", "upet"),
        )
        model = load_atomistic_model(pt_file)
        model.capabilities().dtype = dtype
        model = model.to(
            dtype={"float32": torch.float32, "float64": torch.float64}[dtype],
            device=device,
        )
        return MetatomicCalculator(model, device=device, non_conservative=False)

    return make_calc


def _chgnet(device: str) -> "Calculator":
    from chgnet.model.dynamics import CHGNetCalculator

    return CHGNetCalculator(use_device=device)


def _m3gnet(device: str) -> "Calculator":  # noqa: ARG001 - matgl manages device
    import matgl
    from matgl.ext.ase import PESCalculator

    # M3GNet-MP-2021.2.8-PES ships in DGL format; matgl>=2 defaults to a non-DGL
    # backend and removed this checkpoint entirely in v4, so pin matgl<4 + DGL.
    # matgl 3.x set_backend expects the uppercase "DGL"/"PYG" literal.
    matgl.set_backend("DGL")
    return PESCalculator(matgl.load_model("M3GNet-MP-2021.2.8-PES"))


def _equflash(model_key: str) -> Callable[..., "Calculator"]:
    def make_calc(device: str, checkpoint: str | None = None) -> "Calculator":
        from GGNN.common.calculator import UCalculator

        # UCalculator is an ASE Calculator (the kappa task drives phonopy with it);
        # figshare ships a torch .pt checkpoint that it loads by path
        ckpt = checkpoint or download_checkpoint(model_key, ext=".pt")
        return UCalculator(checkpoint_path=ckpt, cpu=device == "cpu")

    return make_calc


def _emt(device: str) -> "Calculator":  # noqa: ARG001 - CPU only, debug model
    from ase.calculators.emt import EMT

    return EMT()


# Executable pins live in each model YAML's environment. Factories above keep
# checkpoint behavior; _runtime_calc_spec supplies the uv env. Compatibility notes:
# - fairchem-core v1: torch 2.4, matching PyG wheels, numpy<2, scipy<1.15
# - HIENet: torch 2.1.2 needs Python 3.11 + matching PyG wheels
# - NequIP/Allegro: torch<2.10 (TorchScript GPU path)
# - MACE on CUDA: cuEquivariance CUDA 12 extras
# - EqNorm: torch 2.2.2 + vesin 0.3.2 + matching PyG wheels
# - M3GNet: matgl<4 + DGL backend + dgl 2.4 wheel
# - TACE: install from upstream git (PyPI tace is unrelated). Exact torch and
#   OpenEquivariance versions are test pins, not TACE compatibility bounds.

_ALPHANET_PRETRAINED = (
    "https://raw.githubusercontent.com/zmyybc/AlphaNet/"
    "65f8ea9330459e0106867d1c694aec4139c6cb19/pretrained"
)


def _deepmd_spec(
    model_key: str,
    factory: Callable[[str], Callable[..., "Calculator"]] = _deepmd,
    *,
    ext: str = ".pth",
) -> Callable[[], CalcSpec]:
    """YAML-backed DeePMD spec with managed checkpoint download."""
    return _named_spec(factory, model_key, checkpoint=True, ext=ext)


CALCULATORS: _CalcRegistry = _CalcRegistry(
    {
        "mace_mp_0": _runtime_calc_spec("mace_mp_0", _mace("medium")),
        "mace_mpa_0": _runtime_calc_spec("mace_mpa_0", _mace("medium-mpa-0")),
        "orb_v2": _runtime_calc_spec("orb_v2", _orb("orb-v2")),
        "orb_v3": _runtime_calc_spec("orb_v3", _orb("orb-v3-conservative-inf-mpa")),
        "orb_v2_mptrj": _runtime_calc_spec("orb_v2_mptrj", _orb("orb-mptraj-only-v2")),
        "mattersim_v1_5m": _runtime_calc_spec(
            "mattersim_v1_5m", _mattersim("mattersim-v1.0.0-5m.pth")
        ),
        "sevennet_0": _runtime_calc_spec("sevennet_0", _sevennet("7net-0")),
        "sevennet_l3i5": _runtime_calc_spec("sevennet_l3i5", _sevennet("7net-l3i5")),
        "sevennet_mf_ompa": _checkpoint_spec(
            "sevennet_mf_ompa",
            _sevennet("sevennet_mf_ompa", modal="mpa", from_url=True),
        ),
        "sevennet_omni_i12": _runtime_calc_spec(
            "sevennet_omni_i12", _sevennet("7net-omni-i12", modal="mpa")
        ),
        "grace_2l_oam": _runtime_calc_spec("grace_2l_oam", _grace("GRACE-2L-OAM")),
        "grace_1l_oam": _runtime_calc_spec("grace_1l_oam", _grace("GRACE-1L-OAM")),
        "grace_2l_oam_l": _runtime_calc_spec(
            "grace_2l_oam_l", _grace("GRACE-2L-OMAT-large-ft-AM")
        ),
        "grace_3l_oam_l": _runtime_calc_spec(
            "grace_3l_oam_l", _grace("GRACE-3L-OMAT-large-ft-AM")
        ),
        "grace_2l_mptrj": _runtime_calc_spec(
            "grace_2l_mptrj", _grace("GRACE-2L-MP-r6")
        ),
        "chgnet_0_3_0": _runtime_calc_spec("chgnet_0_3_0", _chgnet),
        "hienet": _named_spec(_hienet, "hienet", checkpoint=True),
        "nequip_mp_l_0_1": _named_spec(_nequip, "nequip_mp_l_0_1"),
        "nequip_oam_l_0_1": _named_spec(_nequip, "nequip_oam_l_0_1"),
        "nequip_oam_xl_0_1": _named_spec(_nequip, "nequip_oam_xl_0_1"),
        "allegro_mp_l_0_1": _named_spec(_nequip, "allegro_mp_l_0_1"),
        "allegro_oam_l_0_1": _named_spec(_nequip, "allegro_oam_l_0_1"),
        "matris_10m_oam": _checkpoint_spec(
            "matris_10m_oam", _matris("matris_10m_oam", "MatRIS_10M_OAM.pth.tar")
        ),
        "matris_10m_mp": _checkpoint_spec(
            "matris_10m_mp", _matris("matris_10m_mp", "MatRIS_10M_MP.pth.tar")
        ),
        "eqnorm_mptrj": _checkpoint_spec(
            "eqnorm_mptrj", _eqnorm("eqnorm_mptrj", "eqnorm-mptrj"), ext=".pt"
        ),
        "nequix_mp_1": _named_spec(_nequix, "nequix_mp_1", checkpoint=True, ext=".nqx"),
        "nequix_mp_1_pft": _named_spec(
            _nequix, "nequix_mp_1_pft", checkpoint=True, ext=".nqx"
        ),
        "pet_oam_xl_1_0_0": _named_spec(_pet, "pet_oam_xl_1_0_0", checkpoint=True),
        "alphanet_v1_mptrj": _checkpoint_spec(
            "alphanet_v1_mptrj",
            _alphanet("alphanet_v1_mptrj", f"{_ALPHANET_PRETRAINED}/MPtrj/mp.json"),
        ),
        "alphanet_v1_oam": _checkpoint_spec(
            "alphanet_v1_oam",
            _alphanet("alphanet_v1_oam", f"{_ALPHANET_PRETRAINED}/OMA/oma.json"),
        ),
        "m3gnet": _runtime_calc_spec("m3gnet", _m3gnet),
        "eqv2_s_dens_mp": _named_spec(_fairchem, "eqv2_s_dens_mp", checkpoint=True),
        "eqv2_m_omat_salex_mp": _named_spec(
            _fairchem, "eqv2_m_omat_salex_mp", checkpoint=True
        ),
        "esen_30m_oam": _named_spec(_fairchem, "esen_30m_oam", checkpoint=True),
        "esen_30m_mp": _named_spec(_fairchem, "esen_30m_mp", checkpoint=True),
        "equiformer_v3_mp": _named_spec(
            _fairchem, "equiformer_v3_mp", checkpoint=True, requires_checkpoint=True
        ),
        "equiformer_v3_oam": _named_spec(
            _fairchem, "equiformer_v3_oam", checkpoint=True, requires_checkpoint=True
        ),
        "tace_v1_oam_m": _named_spec(_tace, "tace_v1_oam_m", checkpoint=True),
        "tace_oam_l": _named_spec(_tace, "tace_oam_l", checkpoint=True),
        "tace_oam_rra_preview": _named_spec(
            _tace, "tace_oam_rra_preview", checkpoint=True
        ),
        "tece_oam_rra_1_0": _named_spec(_tace, "tece_oam_rra_1_0", checkpoint=True),
        "dpa_3_1_mptrj": _deepmd_spec("dpa_3_1_mptrj"),
        "dpa_3_1_3m_ft": _deepmd_spec("dpa_3_1_3m_ft"),
        "dpa3_v2_mptrj": _deepmd_spec("dpa3_v2_mptrj"),
        "dpa3_v2_openlam": _deepmd_spec("dpa3_v2_openlam"),
        "dpa3_v1_mptrj": _deepmd_spec("dpa3_v1_mptrj", _deepmd_archive, ext=".zip"),
        "dpa3_v1_openlam": _deepmd_spec("dpa3_v1_openlam", _deepmd_archive, ext=".zip"),
        "dpa_4_0_pro_mptrj": _deepmd_spec(
            "dpa_4_0_pro_mptrj", _deepmd_freeze, ext=".pt"
        ),
        "dpa_4_0_1_pro_mptrj": _deepmd_spec(
            "dpa_4_0_1_pro_mptrj", _deepmd_freeze, ext=".pt"
        ),
        "equflash_29m_oam": _named_spec(
            _equflash, "equflash_29m_oam", checkpoint=True, ext=".pt"
        ),
        "equflashv2_45m_oam": _named_spec(
            _equflash, "equflashv2_45m_oam", checkpoint=True, ext=".pt"
        ),
        # Eager: no Model YAML / environment; do not wrap with _runtime_calc_spec.
        "emt": CalcSpec(_emt),
    }
)


def _model_ref_key(model_ref: str) -> str:
    """Resolve a Model ref to its registry key, preserving debug-only keys."""
    from matbench_discovery.enums import Model

    try:
        return Model.from_ref(model_ref).name
    except ValueError:
        return model_ref


def resolve_calculator_key(model_ref: str) -> str:
    """Resolve a Model ref or debug key to a registered calculator key."""
    model_key = _model_ref_key(model_ref)
    if model_key not in CALCULATORS:
        raise ValueError(
            f"Unknown model {model_ref!r}, pick from {sorted(CALCULATORS)} or "
            "register it in matbench_discovery/calculators.py"
        )
    return model_key


def resolve_cli_calculator(
    parser: "argparse.ArgumentParser",
    model_ref: str | None,
    *,
    list_models: bool = False,
    archived_reasons: "Mapping[str, str] | None" = None,
    task: str = "task",
) -> str | None:
    """Shared-runner CLI preamble: print the registry or resolve --model.

    Returns None after printing one key: deps line per registered calculator
    for --list-models. Otherwise resolves the model reference to a calculator
    key, exiting via parser.error when the model is missing, archived for this
    task (per archived_reasons), or not registered.
    """
    if list_models:
        for model_key, spec in CALCULATORS.items():
            if archived_reasons and model_key in archived_reasons:
                continue
            print(f"{model_key}: {', '.join(spec.deps) or '(core deps only)'}")
        return None
    if not model_ref:
        parser.error("--model is required (or pass --list-models)")

    model_key = _model_ref_key(model_ref)
    if archived_reasons and (reason := archived_reasons.get(model_key)):
        parser.error(f"{model_key} {task} is archived: {reason}")
    try:
        return resolve_calculator_key(model_key)
    except ValueError as exc:
        parser.error(f"{exc}, see --list-models")


def load_calculator(
    model_key: str,
    device: str | None = None,
    dtype: str = "float64",
    checkpoint: str | None = None,
) -> "Calculator":
    """Instantiate the ASE calculator for a registered model.

    Args:
        model_key: Key into CALCULATORS (a Model enum name).
        device: 'cuda' or 'cpu'. Defaults to auto-detection (cuda if torch sees a GPU,
            except the CPU-only 'emt' debug model).
        dtype: Floating-point precision ('float64' or 'float32'). Passed to MACE,
            AlphaNet, and PET; other calculators keep their package defaults.
        checkpoint: Explicit local model artifact for factories that support one.

    Returns:
        Calculator: The model's ASE calculator.
    """
    model_key = resolve_calculator_key(model_key)
    device = resolve_device(model_key, device)
    calc_spec = CALCULATORS[model_key]
    checkpoint = resolve_checkpoint(model_key, checkpoint)
    factory_params = inspect.signature(calc_spec.make_calc).parameters
    kwargs: dict[str, str] = {"dtype": dtype} if "dtype" in factory_params else {}
    # Pass optional runtime values only to factories that explicitly declare them.
    if checkpoint is not None:
        if "checkpoint" not in factory_params:
            raise ValueError(
                f"{model_key} does not support an explicit checkpoint override"
            )
        kwargs["checkpoint"] = checkpoint
    return calc_spec.make_calc(device, **kwargs)
