"""Registry of MLIP ASE calculators for the molecular dynamics benchmark.

Every model the leaderboard runs is launched through the single ``models/run_md.py``
script. Because MLIP dependency trees conflict (torch vs jax vs tensorflow, mutually
exclusive CUDA builds), they cannot share one environment. Instead each model declares
its own ``uv`` requirements here; the launcher resolves a per-model environment
on the fly with ``uv run --with`` (see ``MdModel.uv_run_cmd``), on top of the core
dependencies declared in ``models/run_md.py``'s inline script metadata. Calculator
construction is lazy (imports happen inside the factory) so listing models and
printing dependencies work with only the core dependencies installed.

Registry keys are ``Model`` enum names so metrics can be written to the right YAML.
"""

import inspect
import os
from collections.abc import Callable
from dataclasses import dataclass
from typing import TYPE_CHECKING

from matbench_discovery import DEFAULT_CACHE_DIR

if TYPE_CHECKING:
    from ase.calculators.calculator import Calculator

MD_CHECKPOINT_DIR = f"{DEFAULT_CACHE_DIR}/md-checkpoints"


def download_checkpoint(model_key: str, ext: str | None = None) -> str:
    """Download a model's weights from its YAML ``checkpoint_url`` to a local cache
    and return the path. Normalizes HuggingFace ``/blob/`` to ``/resolve/`` and
    sciebo share links to direct downloads; figshare URLs are handled by
    ``download_file``. Cached after the first call.

    Args:
        model_key: Model enum name whose YAML carries ``checkpoint_url``.
        ext: Force this file extension (e.g. '.pth' for deepmd, whose loader picks
            its backend by suffix). Defaults to the extension parsed from the URL.
    """
    from matbench_discovery.enums import Model
    from matbench_discovery.remote.fetch import download_file

    url = Model.from_ref(model_key).metadata.get("checkpoint_url")
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
    dest = f"{MD_CHECKPOINT_DIR}/{model_key}{ext}"
    os.makedirs(MD_CHECKPOINT_DIR, exist_ok=True)
    # serialize concurrent same-model downloads: parallel array tasks otherwise race on
    # download_file's temp-file rename, leaving losers with a vanished .part file
    from filelock import FileLock

    with FileLock(f"{dest}.lock"):
        if not os.path.isfile(dest):
            download_file(dest, url, headers=headers)
    if not os.path.isfile(dest):  # download_file prints instead of raising
        raise RuntimeError(
            f"Failed to download {model_key} checkpoint from {url}. If the repo is "
            "gated (e.g. fairchem OMAT24), accept its license on HuggingFace and set "
            "HF_TOKEN in the environment."
        )
    return dest


@dataclass(frozen=True)
class MdModel:
    """A registered MLIP: how to build its calculator and its uv requirements."""

    make_calc: Callable[..., "Calculator"]
    deps: tuple[str, ...] = ()  # extra uv requirements beyond CORE_DEPS
    find_links: tuple[str, ...] = ()  # uv --find-links (e.g. PyG/dgl wheel pages)
    # pin uv's Python (e.g. HIENet's torch 2.1.2 has no cp312 wheels, so needs 3.11);
    # None lets uv pick the default
    python_version: str | None = None

    def uv_run_cmd(self, script: str, *args: str) -> list[str]:
        """``uv run`` command that resolves this model's env and runs the script."""
        py_args = ["--python", self.python_version] if self.python_version else []
        with_args = [tok for dep in self.deps for tok in ("--with", dep)]
        link_args = [tok for url in self.find_links for tok in ("--find-links", url)]
        return ["uv", "run", *py_args, *with_args, *link_args, script, *args]


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


def _mace_mp(checkpoint: str) -> Callable[[str, str], "Calculator"]:
    def make_calc(device: str, dtype: str = "float64") -> "Calculator":
        from mace.calculators import mace_mp

        # enable_cueq left off: the cuequivariance fast path is version-brittle across
        # mace/cueq releases; correctness is identical, only throughput differs
        return mace_mp(model=checkpoint, device=device, default_dtype=dtype)

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
    model_name: str, modal: str | None = None
) -> Callable[[str], "Calculator"]:
    def make_calc(device: str) -> "Calculator":
        from sevenn.calculator import SevenNetCalculator

        kwargs = {"modal": modal} if modal else {}
        return SevenNetCalculator(model=model_name, device=device, **kwargs)

    return make_calc


def _grace(model_name: str) -> Callable[[str], "Calculator"]:
    def make_calc(device: str) -> "Calculator":  # noqa: ARG001 - TF picks the device
        from tensorpotential.calculator import grace_fm

        return grace_fm(model_name)

    return make_calc


def _fairchem(model_key: str) -> Callable[[str], "Calculator"]:
    def make_calc(device: str) -> "Calculator":
        from fairchem.core import OCPCalculator

        checkpoint = download_checkpoint(model_key)
        return OCPCalculator(checkpoint_path=checkpoint, cpu=device == "cpu", seed=0)

    return make_calc


def _deepmd(model_key: str) -> Callable[[str], "Calculator"]:
    def make_calc(device: str) -> "Calculator":  # noqa: ARG001 - DP picks the device
        from deepmd.calculator import DP

        # DP selects its backend by file suffix, so force '.pth' (figshare URLs are
        # extensionless and would otherwise save as '.ckpt')
        return DP(download_checkpoint(model_key, ext=".pth"))

    return make_calc


def _tace(model_key: str) -> Callable[[str], "Calculator"]:
    def make_calc(device: str) -> "Calculator":
        from tace.interface.ase import TACEAseCalc

        return TACEAseCalc(download_checkpoint(model_key), use_ema=True, device=device)

    return make_calc


def _hienet(model_key: str) -> Callable[[str], "Calculator"]:
    def make_calc(device: str) -> "Calculator":
        from hienet.hienet_calculator import HIENetCalculator

        return HIENetCalculator(model=download_checkpoint(model_key), device=device)

    return make_calc


def _matris(model: str, cache_name: str) -> Callable[[str], "Calculator"]:
    def make_calc(device: str) -> "Calculator":
        import shutil

        from filelock import FileLock
        from matris.applications import MatRISCalculator

        # MatRIS.load() only takes a registered model name and resolves it to
        # ~/.cache/matris/<cache_name> (its built-in figshare downloader serves 0 bytes
        # behind a WAF), so stage our YAML checkpoint there before constructing the calc
        cache_dir = os.path.expanduser("~/.cache/matris")
        os.makedirs(cache_dir, exist_ok=True)
        dest = f"{cache_dir}/{cache_name}"
        # lock + atomic replace so parallel tasks don't race and an interrupted copy
        # can't leave a truncated checkpoint that later runs treat as valid cache
        with FileLock(f"{dest}.lock"):
            if not os.path.isfile(dest):
                tmp_dest = f"{dest}.tmp"
                shutil.copy(download_checkpoint(model), tmp_dest)
                os.replace(tmp_dest, dest)
        return MatRISCalculator(model=model, device=device)

    return make_calc


def _alphanet(model_key: str, config_url: str) -> Callable[[str], "Calculator"]:
    def make_calc(device: str) -> "Calculator":
        from alphanet.config import All_Config
        from alphanet.infer.calc import AlphaNetCalculator

        from matbench_discovery.remote.fetch import download_file

        # AlphaNet needs an architecture config json (not bundled in the weights) that
        # matches the checkpoint; fetch the OMA config that ships in the AlphaNet repo
        cfg_path = f"{MD_CHECKPOINT_DIR}/{model_key}-config.json"
        if not os.path.isfile(cfg_path):
            download_file(cfg_path, config_url)
        return AlphaNetCalculator(
            ckpt_path=download_checkpoint(model_key),
            device=device,
            precision="32",
            config=All_Config().from_json(cfg_path),
        )

    return make_calc


def _pet(model_key: str) -> Callable[[str], "Calculator"]:
    def make_calc(device: str) -> "Calculator":
        import subprocess

        from filelock import FileLock
        from metatomic.torch.ase_calculator import MetatomicCalculator

        # PET ships a metatrain .ckpt that must be exported to a TorchScript .pt before
        # metatomic can load it; cache the export and lock it against parallel tasks
        ckpt = download_checkpoint(model_key)
        pt_file = f"{os.path.splitext(ckpt)[0]}.pt"
        with FileLock(f"{pt_file}.lock"):
            if not os.path.isfile(pt_file):
                # export to a temp file and atomically rename so a failed `mtt export`
                # (raises via check=True) can't leave a partial .pt that the next run
                # would treat as a valid cached export
                tmp_pt = f"{pt_file}.tmp"
                subprocess.run(["mtt", "export", ckpt, "-o", tmp_pt], check=True)  # noqa: S607
                os.replace(tmp_pt, pt_file)
        return MetatomicCalculator(pt_file, device=device)

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


def _emt(device: str) -> "Calculator":  # noqa: ARG001 - CPU only, debug model
    from ase.calculators.emt import EMT

    return EMT()


# key = Model enum name; deps are uv requirement strings (https://peps.python.org/pep-0508)
# Models get their weights either from their calculator's own auto-download (mace, orb,
# sevennet, grace, chgnet, mattersim) or via download_checkpoint() from the YAML
# checkpoint_url (fairchem, deepmd, tace). Gated repos (fairchem OMAT24) need HF_TOKEN.

# fairchem-core v1 pins torch~=2.4 and imports torch_geometric (needs torch-scatter).
# PyG ships prebuilt wheels per torch+CUDA build; PYG_LINKS must match the torch pin.
# scipy<1.15 because fairchem v1 imports scipy.special.sph_harm (removed in 1.15)
FAIRCHEM_DEPS = (
    "fairchem-core[torch-extras]==1.10.0",
    "torch==2.4.0",
    "numpy<2",
    "scipy<1.15",
)
PYG_LINKS = ("https://data.pyg.org/whl/torch-2.4.0+cu121.html",)

# HIENet pins torch 2.1.2 (no cp312 wheels -> python_version="3.11") and pulls compiled
# torch-scatter/torch-geometric from the torch-2.1.2 PyG wheel page.
HIENET_DEPS = (
    "hienet @ git+https://github.com/divelab/AIRS.git#subdirectory=OpenMat/HIENet",
    "torch==2.1.2",
    "torch-geometric==2.6.1",
    "torch-scatter==2.1.2",
    "e3nn==0.5.6",
    "braceexpand",
    "numpy<2",
    # this py311 resolution pulls a pymatviz (imported by matbench_discovery/__init__)
    # that still imports plotly.validators.scatter, removed in plotly 6
    "plotly<6",
)

MD_MODELS: dict[str, MdModel] = {
    "mace_mp_0": MdModel(_mace_mp("medium"), deps=("mace-torch>=0.3.6",)),
    "mace_mpa_0": MdModel(_mace_mp("medium-mpa-0"), deps=("mace-torch>=0.3.6",)),
    "orb_v2": MdModel(_orb("orb-v2"), deps=("orb-models==0.4.3",)),
    "orb_v3": MdModel(
        _orb("orb-v3-conservative-inf-omat"), deps=("orb-models==0.5.4",)
    ),
    "orb_v2_mptrj": MdModel(_orb("orb-mptraj-only-v2"), deps=("orb-models==0.4.3",)),
    "mattersim_v1_5m": MdModel(
        _mattersim("mattersim-v1.0.0-5m.pth"), deps=("mattersim",)
    ),
    "sevennet_l3i5": MdModel(_sevennet("7net-l3i5"), deps=("sevenn",)),
    "sevennet_omni_i12": MdModel(
        _sevennet("7net-omni-i12", modal="mpa"), deps=("sevenn",)
    ),
    "grace_2l_oam": MdModel(_grace("GRACE-2L-OAM"), deps=("tensorpotential",)),
    "grace_1l_oam": MdModel(_grace("GRACE-1L-OAM"), deps=("tensorpotential",)),
    "grace_2l_oam_l": MdModel(
        _grace("GRACE-2L-OMAT-large-ft-AM"), deps=("tensorpotential",)
    ),
    # grace_2l_mptrj == registry name "GRACE-2L-MP-r6" (sciebo 42Ivgi3eaLCynwC), but
    # tensorpotential>=0.5 (all that's on PyPI) dropped the MP-r6 models, so pin the
    # 0.4.4-era commit from git that still registers it
    "grace_2l_mptrj": MdModel(
        _grace("GRACE-2L-MP-r6"),
        deps=(
            "tensorpotential @ git+https://github.com/ICAMS/grace-tensorpotential@3115a9314",
        ),
    ),
    "chgnet_030": MdModel(_chgnet, deps=("chgnet",)),
    # HIENet (e3nn-based, github checkpoint). torch 2.1.2 needs python 3.11.
    "hienet": MdModel(
        _hienet("hienet"),
        deps=HIENET_DEPS,
        find_links=("https://data.pyg.org/whl/torch-2.1.2+cu121.html",),
        python_version="3.11",
    ),
    # MatRIS (figshare checkpoint staged into MatRIS's own ~/.cache/matris dir)
    "matris_10m_oam": MdModel(
        _matris("matris_10m_oam", "MatRIS_10M_OAM.pth.tar"),
        deps=(
            "matris @ git+https://github.com/HPC-AI-Team/MatRIS",
            "torch==2.6.0",
            "numpy<3",
        ),
    ),
    # PET (metatrain/metatomic): download .ckpt, `mtt export` to .pt, then load
    "pet_oam_xl_1_0_0": MdModel(_pet("pet_oam_xl_1_0_0"), deps=("upet==0.1.0",)),
    # AlphaNet (config json from repo + figshare weights, torch-scatter from PyG links)
    "alphanet_v1_oam": MdModel(
        _alphanet(
            "alphanet_v1_oam",
            "https://raw.githubusercontent.com/zmyybc/AlphaNet/main/pretrained/OMA/oma.json",
        ),
        deps=(
            "alphanet @ git+https://github.com/zmyybc/AlphaNet",
            "torch==2.5.1",
            "torch-geometric==2.6.1",
            "torch-scatter",
            "numpy<2",
        ),
        find_links=("https://data.pyg.org/whl/torch-2.5.1+cu121.html",),
    ),
    # matgl<4 + DGL backend; DGL wheels come from its torch-matched find-links page.
    # Pin dgl==2.4.0: data.dgl.ai serves the index but 403s the dgl-2.5.0 cu121 wheel,
    # so an unpinned dgl resolves to 2.5.0 and fails to download.
    "m3gnet_ms": MdModel(
        _m3gnet,
        deps=("matgl==3.0.5", "dgl==2.4.0", "torch==2.4.0"),
        find_links=("https://data.dgl.ai/wheels/torch-2.4/cu121/repo.html",),
    ),
    # fairchem family (OCPCalculator + HF .pt checkpoints, fairchem-core v1 API).
    # torch~=2.4 + torch-geometric pulls compiled PyG extensions (torch-scatter etc.);
    # their prebuilt wheels resolve from PYG_LINKS matching torch 2.4.0+cu121.
    "eqv2_s_dens_mp": MdModel(
        _fairchem("eqv2_s_dens_mp"), deps=FAIRCHEM_DEPS, find_links=PYG_LINKS
    ),
    "eqv2_m_omat_salex_mp": MdModel(
        _fairchem("eqv2_m_omat_salex_mp"), deps=FAIRCHEM_DEPS, find_links=PYG_LINKS
    ),
    # NOTE: esen_30m_mp omitted - its YAML checkpoint_url wrongly points to
    # esen_30m_oam.pt (the MP-trained file is esen_30m_mptrj.pt); fix the YAML to add it
    "esen_30m_oam": MdModel(
        _fairchem("esen_30m_oam"), deps=FAIRCHEM_DEPS, find_links=PYG_LINKS
    ),
    # TACE: install from the model's github repo (PyPI 'tace' is a different/stale
    # package that can't instantiate the checkpoint config); pins per its YAML
    "tace_oam_l": MdModel(
        _tace("tace_oam_l"),
        deps=(
            "tace @ git+https://github.com/xvzemin/tace",
            "torch==2.9.1",
            "torch-geometric==2.7.0",
            "pytorch-lightning==2.5.5",
        ),
    ),
    # deepmd DPA (figshare frozen .pth checkpoints, DP() loads by suffix).
    # NOTE: dpa_4_0_pro_mptrj omitted - its figshare file is a training checkpoint
    # (state dict), not a frozen TorchScript model; needs `dp --pt freeze` first
    "dpa_3_1_3m_ft": MdModel(_deepmd("dpa_3_1_3m_ft"), deps=("deepmd-kit[torch]",)),
    # CPU-only debug model for smoke-testing the pipeline without heavy installs
    "emt": MdModel(_emt),
}


def load_calculator(
    model_key: str, device: str | None = None, dtype: str = "float64"
) -> "Calculator":
    """Instantiate the ASE calculator for a registered model.

    Args:
        model_key: Key into MD_MODELS (a Model enum name).
        device: 'cuda' or 'cpu'. Defaults to auto-detection (cuda if torch sees a GPU,
            except the CPU-only 'emt' debug model).
        dtype: Floating-point precision ('float64' or 'float32'). Only MACE models
            honor it; other calculators keep their package defaults.

    Returns:
        Calculator: The model's ASE calculator.
    """
    if model_key not in MD_MODELS:
        raise ValueError(
            f"Unknown {model_key=}, pick from {sorted(MD_MODELS)} or register it "
            "in matbench_discovery/md_models.py"
        )
    if device is None:
        device = "cpu" if model_key == "emt" else _detect_device()
    make_calc = MD_MODELS[model_key].make_calc
    # pass dtype only to factories that declare it (currently MACE), so a new
    # dtype-aware model is honored automatically without editing a hardcoded key set
    if "dtype" in inspect.signature(make_calc).parameters:
        return make_calc(device, dtype=dtype)
    return make_calc(device)
