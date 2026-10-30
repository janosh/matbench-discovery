# MACE-SOYO

A multi-head interatomic potential with experimental compiled training, AOTInductor (AOTI) inference for
ASE/TorchSim, and a LAMMPS ML-IAP interface(in test). The current branch predicts energy,
forces, and stress; molecular charge/spin conditioning is optional.
This project aims to unlock the full potential of MACE-mh-1. (https://github.com/ACEsuit/mace)

The main artitecture and idea is the same as MACE-mh-1, while simple attention for edge_weight, Equivariant RMSnorm, code modification for speed,
Spin&Charge embedding, heavier readout, Muon(optimizer) and WSD(schedule from DPA4) are added.
Broadly speaking, MACE-SOYO offers an accuracy–speed trade-off between DPA4-Neo and DPA4-Mini.
Training, AOTI export and inference require an NVIDIA CUDA GPU. CPU and other
device types are rejected explicitly due to the use of cuequivariance and nvalchemi-toolkit.

This is an independently developed project with assistance from GPT-6 Astra to formalize the code.
Although all code except for LAMMPS ML-IAP is double checked by both Astra and me, and the matbench
results are encouraging, this is the first time for me to complete such formal project in my spare time.
Some minor bugs may still exist, it would be very appreciated if you report them. Thank You!

Also, this code is inspired or inherited from several other superior potentials
like Nequip, Nequix, DPA4 and EquflashV2. I give my special thank to those MLIAPs' creators.

## Installation

Python >=3.12,<3.15 is required by the pinned dependencies.
From the repository root:
I recommend installing PyTorch(torch 2.12) separately first, choosing a CUDA build compatible with your NVIDIA driver.

nvcc is NEEDED for AOTI export. Make sure you have it in environment.
Type nvcc to check whether there is an output.

```bash
python -m pip install -r requirements.txt
```

## Checkpoints

The evaluated MPtraj checkpoint has only one head, `mptraj`, and does not use
charge/spin conditioning. Download
[`mace_soyo_mptraj.pth`](https://github.com/PrincessHakuryu/MACE-SOYO/releases/download/v0.1/mace_soyo_mptraj.pth)
from the [v0.1 release](https://github.com/PrincessHakuryu/MACE-SOYO/releases/tag/v0.1)

The provided training checkpoint is a `.pth` file. ASE and TorchSim load an AOTI `.pt2`
package; follow the export example below on the intended inference GPU.

A nine-head model trained on the full datasets may be released in the future:
OMat24, sAlex + MPTrj, OC20, OC22, OC25, OMol25 + OPoly26, ODAC25, OMC25, and
MatPES-r2SCAN.

## Matbench-discovery results

Discovery metrics use the `unique_prototypes` subset; normalized RMSD uses
`symprec=1e-2`. Combined scores use the current Matbench default weights.

| Metric | MACE-MH-1 reference | MACE-SOYO MPtraj |
| --- | ---: | ---: |
| CPS (higher is better) | 0.735 | 0.763 |
| F1 | 0.7549 | 0.7766 |
| Accuracy | 0.9172 | 0.9264 |
| DAF | 4.5574 | 4.7949 |
| Precision | 0.6967 | 0.7330 |
| Discovery MAE (eV/atom; lower is better) | 0.04449 | 0.04136 |
| Discovery R2 | 0.7849 | 0.7987 |
| kappa_SRME (lower is better) | 0.4197 | 0.3498 |
| RMSD (normalized; lower is better) | 0.0880 | 0.0832 |
| CMDS (higher is better) | Not measured | 0.750 |
| CDS (higher is better) | Not measured | 0.790 |

Both models used the same fixed MPtraj train/validation split, 100 epochs, and E:F:S loss weights of 1:1:1.
The MACE-MH-1 inference reference used cuEquivariance and `torch.compile`.
Its training configuration is
[`example/mace-mh1-mptraj.yml`](example/mace-mh1-mptraj.yml).
These benchmark results DO NOT establish that MACE-SOYO is superior to MACE-MH-1 in all respects; broader evaluation is needed.

Evaluation settings:

- Discovery: ASE FIRE with FrechetCellFilter, `fmax=0.02` eV/Angstrom and
  `max_steps=500`.
- Phonons: `phonondb-v1`, displacement 0.03 Angstrom, relaxation
  `force_max=1e-4` eV/Angstrom and `max_steps=300`.
- MACE-SOYO hardware: RTX 4090 for discovery/phonons; H200 NVL for MD/diatomics.


## Training

### Recommended PyTorch version

I strongly recommend **PyTorch 2.12.0** for compiled training. Constraints on
Pytorch's version may be looser if you only want to export and inference.

Other PyTorch versions may have incompatibilities with `make_fx`, compiled
force/stress training, or dynamic shapes. Disable compile_training in config.yml
to avoid any incompatibilities if the automatic gradcheck fails.

Enable compiled training would speed up training for about 30% in MACE-SOYO.

### Time and VRAM
Training Time and VRAM cost:
Taking MPtraj as an example, the total training time is 5.0 GPU-days on 4 A100
with compiled_training on.(100 epochs)The batchsize is 48 per GPU. The peak VRAM is up to 52GB.
(Refernece: DPA4-Mini:5.8(150 epochs), DPA4-Neo:8.6(100 epochs)
 https://arxiv.org/pdf/2606.02419 Table1, MACE-MH-1-cueq:8.3(100 epochs, test by us))

Inference Time:
Using the supplied 576-atom [example/md.arc](example/md.arc), the average time per NVT step is listed below.
(DPA4-Neo:20.5ms DPA4-Mini:10.2ms MACE-SOYO: 11.2ms MACE-mh-1:22.7ms)
The inference timing benchmark was using the ASE calculator of these models,
DP_TRITON_INFER=2, DP_TF32_INFER=0, DP_AMP_INFER=0 and vesin on for exporting and using DPA4 models,
CUEQ and torch.compile on for MACE-mh-1.
It needs to be mentioned that building neighbour matrix is also included in the time.

### Data and configuration

Training reads ASE LMDB (`.aselmdb`) files. Each entry in `datasets` defines one
named head; a directory may contain multiple shards. Copy the configuration
template before your first run:

```bash
cp config/config.example.yml config/config.yml
```

For the evaluated single-head MPtraj settings, copy
`config/config-mptraj.example.yml` instead. A nine-head fine-tuning template is
provided as `config/config-full.example.yml`. These public templates use
placeholder paths. Local `config.yml`, `config-mptraj.yml` and `config-full.yml`
are excluded from Git so local training paths can be kept without republishing them.

Edit the dataset paths and training settings in `config/config.yml`. This local
configuration is excluded from Git:

```yaml
num_heads: 2
readout_hidden: 64
max_correlations: 3
use_spin_charge: false
datasets:
  - name: omat_pbe
    train_path: /data/omat24/train
    valid_path: /data/omat24/val
    pbc: true
    e0_yaml_path: ./config/e0_omat.yaml
  - name: salexmp
    train_path: /data/salex-mptraj/train
    valid_path: /data/salex-mptraj/val
    pbc: true
    e0_yaml_path: ./config/e0_salex.yaml
```

For every run, `num_heads` must equal the number of dataset
entries. `datasets` declares the target model's head names and their order;
checkpoint heads omitted from this list are not saved in the fine-tuned model.
Head names are stored in
checkpoints and AOTI packages. `pbc: true` means TTT; `false` means FFF and allows
molecules without a cell. Periodic structures require a valid three-dimensional
cell. An empty `valid_path` splits that dataset using `train_percent`,
`valid_percent`, and `split_seed`.


Energy and forces are required, in eV and eV/Å. Stress uses ASE's eV/Å³ convention;
missing stress labels and nonperiodic structures are excluded from stress loss.
Floating-point training data, E0 references, total-energy accumulation, and
E/F/stress losses use FP64. Set `model_dtype: float32` (default) or `float64`
to select network precision. Predicted forces/stress, and parameter
gradients follow the network dtype.
Prepare an E0 YAML for each head before training:

```yaml
atom_energies:
  1: -13.6
  8: -400.0
```

These values illustrate the format only; use references appropriate to the actual
dataset and DFT settings. Every element in the dataset must be present. Training
does not fit E0 automatically. To fit E0 externally using only training data:

```bash
python -m mace_soyo.utils.sample_aselmdb --input /data/train --fit-e0 --e0-output /data/e0.yaml
```

Training and the E0/sampling tool accept one file or discover `.aselmdb` shards
recursively under a directory; the tool's JSON report records the input file list
and shard count.

To convert extxyz data into ASE LMDB:

```bash
python -m mace_soyo.utils.load_data \
  --input_extxyz /data/input.extxyz \
  --output_lmdb /data/train \
  --pbc True \
  --ignorewarn
```

Use `--input_extxyz /data/` to process a directory of extxyz files.
Use `--pbc False` for nonperiodic molecules. For periodic inputs, `--ignorewarn`
keeps frames without stress; without it, those frames are skipped.

### Launch

Run from the repository root. GPU and node counts are controlled by `torchrun`,
not by configuration fields:

```bash
CUDA_VISIBLE_DEVICES=0,1 torchrun --standalone --nproc_per_node=2 run_ddp.py
```

Default checkpoints are saved in `model_pth/` beside `run_ddp.py`, regardless of
the launch directory. The directory is created automatically on the first save.

Use `--output-dir` (alias `--output_dir`) to save checkpoints for different runs
in separate directories. Relative output paths are resolved from the launch
directory, and `~` is expanded. Per-rank logs are also saved under that directory:

```bash
torchrun --nproc_per_node=2 run_ddp.py --config config/my_run.yml --output-dir ./runs/my_run
```

When resuming with `--output-dir`, `resume_path` in the YAML selects the input
checkpoint; new best/last checkpoints are saved in the selected output directory.
Use distinct output directories for runs that should not overwrite each other.
With `resume: true`, `resume_path` must explicitly name the input checkpoint;
a missing, empty or blank value causes an error and exits.Without `--output-dir`,
all last/best checkpoints are saved in `model_pth/` beside `run_ddp.py`; logs use `./logs/debug`.
To continue saving in another run directory, select it with `--output-dir`.

For a new run, set `resume: false` and `finetune: false`. For continuation, set
`resume: true` and `resume_path`, keeping the same architecture and heads.
Resume requires exactly the checkpoint's head count, names and dataset order;
adding, removing, renaming or reordering heads requires fine-tuning.
For fine-tuning, set `finetune: true` and `finetune_path` with `resume: false`:

- Keep the shared backbone and initialize every target readout from a checkpoint
  head. Matching names automatically inherit their own readout. With a single-head
  checkpoint, any new target head automatically inherits that sole readout.
- With a multihead checkpoint, each new name must specify `finetune_head` inside
  its dataset entry to select the source head. This field can also override the
  source for an existing name. Several target heads can copy the same source;
  their readout weights can then train separately.
- Keep `readout_hidden` unchanged. Transfer accounts for the change in Linear
  input normalization when the head count changes.
- `finetune_reset_e0: false` copies each target's E0 from its selected source head.
  SET IT to `true` to replace E0 using each target dataset's prepared YAML in config.
  Fine-tuning starts a new optimizer/schedule; inherited weights remain trainable.

For example, to fine-tune B and D from a checkpoint with heads `[A, B, C, D, E]`,
set `num_heads: 2` and put B and D in `datasets`. They inherit the original B and D
readouts; the saved model contains only `[B, D]`, at head IDs 0 and 1. To expand
a single-head checkpoint `[A]` into `[A, B, C, D, E]`, declare those five dataset
entries and set `num_heads: 5`. All five readouts start from A's weights.

For a new name with a multihead source, choose its initializer explicitly:

```yaml
- name: new_task
  finetune_head: omat_pbe  # Copy this checkpoint head's readout and, unless reset, E0.
  train_path: /data/new_task/train
  valid_path: /data/new_task/val
  pbc: true
  e0_yaml_path: ./config/e0_new_task.yaml
```
For molecular conditioning, set `use_spin_charge: true`. The loader reads graph
labels `charge` and `spin` from ASE row data. Spin means multiplicity (singlet = 1).
Conditioning is enabled only when both labels are present; it does not predict
atomic charges or use atomic magnetic moments as molecular spin.
Charge must be an integer in [-100, 100] and spin in [1, 100].

## AOTI export

Export a trained checkpoint using structures which can read by ASE.
export.arc is provided in example folder.

```bash
python -m mace_soyo.export.AOTI_export \
  --ckpt-path checkpoints/mace_soyo_mptraj.pth \
  --structure-file example/export.arc \
  --output-path model.pt2 \
  --device cuda:0 \
  --dtype float32 \
  --max-graphs 20000
```

Use a structure with at least two atoms and a nonempty neighbor list for tracing.
The tracing example always uses B=2. Every new package supports variable batch
sizes, including B=1 for ASE. `--max-graphs` defaults to 20000 and must be at
least 2. Optional `--max-nodes` and `--max-edges` set upper bounds on atoms and
directed edges.

One `.pt2` contains all heads and their names. Charge/spin inputs are included
automatically when the checkpoint enables them. D3 is not embedded in the package.
Export `--dtype` is independent of training `model_dtype`; its default remains
`float32`, while E0 and output energies remain FP64 for large energy of OMol dataset.
AOTI packages target a particular hardware/runtime environment. Export on the
intended inference GPU and use compatible PyTorch/CUDA/CUEQ versions. DO NOT
assume that a package exported for RTX 4090 is portable to H200; re-export the
`.pth` checkpoint for the target GPU and check the E/F/S consistency results.

## ASE

```python
from ase.io import read, write
from ase.optimize import FIRE
from mace_soyo import MACESoyoCalculator

atoms = read("input.extxyz")
atoms.calc = MACESoyoCalculator(
    package_path="model.pt2",
    head="mptraj",
    device="cuda:0",
    use_d3=False,
)

print(atoms.get_potential_energy())
print(atoms.get_forces())
if atoms.pbc.all():
    print(atoms.get_stress())

FIRE(atoms).run(fmax=0.05, steps=1000)
write("relaxed.extxyz", atoms)
```

This optimizes atomic positions at fixed cell. For atomic and cell relaxation,
pass `ase.filters.FrechetCellFilter(atoms)` to FIRE instead (periodic structures
with a stress-capable head only).

A single-head package may omit `head`. `calc.head_names` lists the names. Each
calculator selects one head at construction. Boundary conditions come from
`atoms.pbc`, not from the head name. Non-TTT stress is a zero placeholder.
Outputs use eV, eV/Å, and eV/Å³.

For a charge/spin-conditioned molecular package, set labels on each structure
before calculation:

```python
atoms.info.update(charge=1, spin=1)
```

Both labels are needed to enable conditioning. Missing either disables it for
that structure; it does not imply a neutral singlet. Do not pass valid electronic
labels to an unconditioned package. And do not use this embedding on PBC structures,
which is unphysical and unreliable.

## TorchSim

```python
import torch_sim as ts
from ase.io import read, write
from mace_soyo import MACESoyoTorchSimAOTIModel
from mace_soyo.torchsim import atoms_to_state

model = MACESoyoTorchSimAOTIModel(
    package_path="model.pt2",
    head="mptraj",
    device="cuda:0",
    use_d3=False,
)
atoms_list = read("input.extxyz", index=":")
state = atoms_to_state(atoms_list, device=model.device, dtype=model.dtype)
results = model(state)
print(results["energy"], results["forces"], results["stress"])

final_state = ts.optimize(
    state,
    model,
    optimizer=ts.Optimizer.fire,
    convergence_fn=ts.generate_force_convergence_fn(force_tol=0.05),
    max_steps=1000,
)
write("relaxed.extxyz", final_state.to_atoms())
```

Use `atoms_to_state` to retain `atoms.info['charge']` and `atoms.info['spin']` when
present. The batch can mix labeled and unlabeled systems; one model instance uses
one head throughout. The example relaxes positions at fixed cell. Match
`model.device` and `model.dtype` when constructing a state. TorchSim 0.6.1 requires
the same PBC flags for all systems passed together to atoms_to_state.

ASE and TorchSim support optional ordinary D3(BJ) with `use_d3=True`. Only enable
it when the selected head's labels DO NOT already include that dispersion
correction. Both interfaces also support `use_laspd3=True` after building
`other-modules/LASP-D3-torchsim` and installing its shared library in `mace_soyo/utils/`;
LASP-D3 requires TTT(CELL IS NEEDED). Do not enable both D3 backends together.

`d3_cutoff_radius` is in Å for both backends, independently of the neural-network
cutoff. The default is 24.59394 Å. LASP-D3 is up to 3x faster than D3 in nvalchemi-toolkit.
But it needs to be manually installed.

For LASP-D3(BJ) in ASE:

```python
calc = MACESoyoCalculator(
    package_path="model.pt2",
    head="mptraj",
    use_laspd3=True,
    use_laspd3_BJ=True,
    functional_type=0,  # LASP-D3's PBE parameter set
)
atoms.calc = calc
FIRE(atoms).run(fmax=0.05, steps=1000)
calc.close()
```

LASP-D3 uses zero damping unless `use_laspd3_BJ=True`. Ordinary D3 parameters
`a1`, `a2`, `s8`, `s6`, and `d3_params_path` apply only to `use_d3=True`;
LASP-D3 selects its parameters with `functional_type` and `use_laspd3_BJ`.

With CUDA (`nvcc`), a C++ compiler, a Fortran compiler and CMake >= 3.25.2:

```bash
cmake -S other-modules/LASP-D3-torchsim -B build/lasp-d3 -DBUILD_D3_TESTS=OFF
cmake --build build/lasp-d3 --target d3_shared -j
cp build/lasp-d3/libd3.so.1.0.0 mace_soyo/utils/
pip install --no-deps .
```

Build before installation, or reinstall after adding the shared library.

LASP-D3 source: [LipidL/LASP-D3](https://github.com/LipidL/LASP-D3).

## LAMMPS ML-IAP

This interface is experimental and completely written by GPT-6 Astra;
its correctness has not yet been fully validated.

This exports a fixed-head Python/PyTorch object (`lammps.pt`), not AOTI or
TorchScript. Use a LAMMPS build with PYTHON, ML-IAP unified Python support, and
CUDA Kokkos. The embedded Python environment must contain MACE-SOYO, PyTorch,
CUEQ, and CuPy:

```bash
python -m pip install cupy-cuda12x
python -m mace_soyo.export.export_lammps_mliap \
  --ckpt-path /path/best.pth \
  --head mptraj \
  --output-path lammps.pt
```

The selected head is fixed at export. Compilation occurs on first use; add
`--no-compile` when exporting to use eager execution. This interface does not
add D3 or expose molecular charge/spin conditioning.

Use `units metal`, `newton on`, and element names matching the LAMMPS atom-type
order. For a C/H system with a 6 Å model cutoff and 1 Å neighbor skin:

```text
units metal
atom_style atomic
newton on
read_data structure.data
neighbor 1.0 bin
comm_modify cutoff 7.0
pair_style mliap unified lammps.pt 0
pair_coeff * * C H
```

Keep the final `0` in `pair_style`; the adapter exchanges local/ghost node
features between layers. Communication cutoff must cover the model cutoff plus
neighbor skin. See `example/in.lammps_mliap` for a complete input.

Single GPU:

```bash
CUDA_VISIBLE_DEVICES=0 OMP_NUM_THREADS=1 lmp -k on g 1 -sf kk \
  -pk kokkos newton on neigh half comm device gpu/aware on -in in.lammps
```

Two GPUs, one MPI rank per GPU (CUDA-aware MPI required):

```bash
CUDA_VISIBLE_DEVICES=0,1 OMP_NUM_THREADS=1 mpirun -np 2 lmp -k on g 2 -sf kk \
  -pk kokkos newton on neigh half comm device gpu/aware on -in in.lammps
```
## Future plan

Keep formalizing the code, and add long-range interactions with multipole_particle_mesh_ewald.

## License

MACE-SOYO's original code is licensed under the [MIT License](LICENSE).
The v0.1 MPtraj model checkpoint is also released under the MIT License.

Identified third-party code adaptations and their original MIT notices are listed
in [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md). These notices are included in
source distributions and wheels.

The vendored sources in `other-modules/LASP-D3-torchsim/` are adapted from
[LipidL/LASP-D3](https://github.com/LipidL/LASP-D3) for GPU-accelerated
dispersion calculations in ASE and TorchSim.
