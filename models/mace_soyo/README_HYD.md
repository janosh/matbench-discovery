# MACE-SOYO-HYD radial v8

This submission targets the **Diatomics CDS subranking**, not the composite CPS
leaderboard. The underlying MACE-SOYO-OAM weights are unchanged. HYD fits 16
energy/force-consistent radial correction parameters on independently generated
heteronuclear PBE/def2-SVP curves. The correction applies to isolated dimers only;
ordinary periodic structures retain the original `salexmp` potential, without D3.

The original MACE-SOYO model and implementation are by Qi-Ming Liang, Jia-Jie Du,
and Zhi-Pan Liu. Its source revision is
`387b5114c974fb04b8aceffac82d1d8f10ea4ac0`. The foundation checkpoint has SHA256
`180784bad2fdee6aec58dae2cdb0d581ad1c4f3d2272d2857bce6315ea43b82b`.
HYD does not claim to have pretrained or fine-tuned those foundation weights.

`requirements_HYD.txt` records the evaluated Linux/Python 3.14/CUDA 12.6
environment, including unused development packages. It is a reproducibility
record rather than a minimal dependency list. Its local source paths are replaced
by the exact upstream source revision and an instruction to install this checkout.

## Training and model selection

The two bundled fit files contain 233 attempted geometries from 28
heteronuclear molecule groups, including unsuccessful reference calculations.
The separate molecule-disjoint validation file adds 64 attempted geometries
from eight groups, bringing the released dataset to 297 geometries in 36 groups.
The fit uses the 24 groups with at least five converged finite reference points.
The reference calculations use neutral RKS/UKS PBE/def2-SVP with level-3 grids.
They are separate from the official homonuclear benchmark. Absolute molecular
energy offsets are centered within each species; the force residual is fitted
jointly with the centered energy residual.

Molecule-grouped development cross-validation selected energy weight 16 and
ridge penalty 0.85. Separate molecule-disjoint development validation and a
broader independent confirmation preceded official evaluation. The frozen
candidate's SHA256 is
`cf601a2cda2510069f9efc118c2f7ce47f6201734e5b44480d3ed547af33f3d9`.
Official labels, public predictions, and official accuracy results were not
used to fit or select these parameters.

For transparency, the submission conservatively flags `diatomics` as
test-set-informed: canonical batching was implemented and verified after the
serial official run. This flag discloses the execution-path work; the scientific
weights and hyperparameters remain the pre-evaluation frozen model.

Refit the selected parameters from the bundled independent training data:

```sh
uv run models/mace_soyo/train_mace_soyo_hyd.py --output refit-candidate.json
```

This recomputes the chosen fit using the checksum-verified `frozen-source.zip`
in the model artifact folder. It also contains the independent fit/validation
data and development evidence. The archive preserves the original grid-search,
reference-generation, and inference sources byte for byte. The production
calculator packages only the required frozen method bodies, with type annotations
and documentation added. Original method bodies and the coefficient file are
verified separately; the archive is not imported during ordinary inference.

## Inference through the shared runners

Register `mace_soyo_hyd_radial_v8` with the shared calculator registry. The
registry factory downloads the original foundation checkpoint and reuses the
upstream MACE-SOYO exporter. `--checkpoint` also accepts an existing AOTInductor
package exported for the current GPU architecture. The evaluated A100/SM80
package SHA256 is
`5ccecb3d856181c57ac4c080618b98e51cbc860a71978c7bccc3f3aaba872d31`.

```sh
uv run models/run_diatomics.py --print-cmd --model mace_soyo_hyd_radial_v8 --dry-run
eval "$(uv run models/run_diatomics.py --print-cmd --model mace_soyo_hyd_radial_v8)"
```

Isolated-dimer energy/force calls prefetch the fixed 119-point 0.1–6 Å geometric
distance grid in batches of 32. Each system retains the original scalar ASE
fractional-coordinate and neighbor-list conversion. This is input batching,
not a lookup of reference or previously published predictions. An initial off-grid
input is included as an additional system; later off-grid inputs for a cached
species use scalar ASE, retaining the grid cache instead of rebuilding it.
Crystals and stress requests use the
unchanged scalar ASE path with the same frozen correction rule.

The evaluated numerical repair was first validated on 144 independent off-grid
geometries. All 10,948 official points subsequently satisfied the already-fixed
combined absolute/relative tolerances against the frozen serial predictions:
energy `5e-4 eV + 1e-5 * abs(serial_energy)` and forces
`1e-3 eV/Å + 1e-5 * abs(serial_force)`. This is tolerance equivalence, not
bit-identical output. The recorded canonical sweep took 69.26 seconds on one
NVIDIA A100-SXM4-80GB; cold export, environment installation, and platform
archival are not part of that sweep time. Failed acceleration attempts remain
in the private audit record and were not substituted for successful outputs.
The native vibrational-frequency diagnostic has 80 valid fits among 81 eligible
reference species, with the K fit unsuccessful and no species missing. This
failure is retained in the metadata; all 92 predicted curves remain present.

The portable adapter's compiled and normal raw-checkpoint/export paths passed
all 10,948 diatomic points against the frozen outputs under the original
tolerances, plus 32 WBM structures and seven independent synthetic structures.
This is numerical reproduction, not model selection or replacement of the
recorded canonical sweep timing.

## WBM protocol

WBM uses **zero-step single-point initial-structure energy predictions** for all
256,963 materials: `max_steps=0`, `test_task=IS2RE`. No structure relaxation was
performed. Initial structures saved for provenance are not presented as a
completed geometry-optimization task. Phonons and molecular dynamics have not
been evaluated. The zero-step protocol was frozen before any WBM prediction.

Predictions and explicit failure records cover all 256,963 inputs and are
verified and published at the metadata URLs. Mark the submission ready for review only after the actual
native ingestion command has completed. A local score or a pull request
does not establish an official public ranking.

## Partial-task ingestion

The accompanying ingestion change skips task-specific reference analysis and
figure generation when every explicitly selected model declares that optional
task pending or unavailable. This avoids requesting phonon and relaxed-geometry
reference files for this single-point submission. Discovery and diatomic
validation still run. A mixed-model refresh retains any task submitted by one
selected model; a full-roster refresh keeps its original behavior. The native
ingestion tests cover these cases, independently of benchmark labels or scores.

## Attribution and license

HYD additions are MIT licensed. The upstream MACE-SOYO code and checkpoint are
MIT licensed; see its original license and third-party notices:
https://github.com/PrincessHakuryu/MACE-SOYO/blob/387b5114c974fb04b8aceffac82d1d8f10ea4ac0/LICENSE
and `THIRD_PARTY_NOTICES.md` in that repository. The upstream training sets are
MPtrj, OMat24, and sAlex. Independent HYD references and their generating code
are released separately from the foundation data and checkpoint.
The independent HYD reference data are released under CC-BY-4.0. The DOI in the
dataset entry cites the PySCF method, not a DOI minted for the HYD dataset;
no dataset-specific DOI is claimed.
