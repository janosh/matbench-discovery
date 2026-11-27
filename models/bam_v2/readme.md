# BAM-V2-Athena-OAM

BAM-V2-Athena-OAM (internal name RACE-V2.1-9M-OAM-Replay) is an equivariant interatomic potential of the RACE v2.1 architecture in
[BAM_V2](https://github.com/myung-group/BAM_V2) (branch `main`). It predicts energy, forces and stress;
forces and stress are derivatives of the energy.

## Architecture

- Layer-wise equivariant cluster expansion with per-layer energy strata: 5 layers, `max_ell=3`, 128 channels,
  cutoff 7.0 Å, 89 species, 9.0M parameters (a single model, not an ensemble)
- Changes with respect to RACE v2.0: bilinear layer update, 64-channel bottleneck of the edge tensor product, low-rank
  linear maps, C4 cutoff envelope `(1 − (r/r_c)^4)^5`
- Physical priors: UFF repulsion and element shell-radius channels

## Training data and schedule

- OMat24 pretraining from scratch (2 data epochs)
- Fine-tuning on MPtrj ×8 + sAlex ×1 (6.2 data epochs)
- Replay stage: about 1 data epoch on the fine-tuning stream plus 20 % of the OMat24 training rows (forces and stress
  only)
- The submitted checkpoint is the last epoch of this schedule, which was fixed before the evaluation; no
  checkpoint was selected on a benchmark score. Configs: `configs/race_v2.1_oam/` in the repository.

## Inference

- The model was trained with the JAX code (`race_v2`), and every benchmark run uses that JAX model: the registered
  calculator wraps `race_v2.inference.relax_batched.RaceEvaluator`, the evaluator of the authors' relaxations, MD and
  diatomic curves (float64 inputs). A PyTorch port (`race_v2_torch`) is in the same repository but was not used here.
- float64 for every task (the thermal conductivity needs it).
- Core exclusion (no parameters, fixed before the evaluation): below the shortest distance of an element pair in the
  training data (capped at half the sum of the covalent radii) the learned edges are switched off with a C4 switch of
  width 0.1 Å and replaced by ZBL. No atom pair of the WBM initial structures, the relaxed WBM structures or the
  PhononDB structures is inside this range, so discovery and κ are unaffected; it acts on the diatomic curves.
- On GPUs with compute capability ≥ 8.0 the edge tensor products and the symmetric contraction use the fused
  cuEquivariance CUDA kernel; older GPUs (V100) run the same polynomials in pure JAX (`RACE_V2_TP_METHOD=naive`). The
  two agree to rounding, but the pure-JAX path needs more than 32 GB of GPU memory for the κ supercells of dense
  cells, so κ was computed on A100.

## Benchmark runs

| Task | How it was run |
| --- | --- |
| Molecular dynamics, diatomic curves | the benchmark's own code (`matbench_discovery.md.run_md_benchmark`, `predict_curves` / `write_results` of `models/run_diatomics.py`, commit `6d587c1`) called on 2026-11-22 from the authors' script `scripts/eval/mbd/md_diatomics.py` with an ASE calculator around the JAX model (float64, core exclusion on, cost label `tight_padding`), NVIDIA GH200; not through the registered calculator |
| Discovery, geometry optimization | relaxed on 2026-11-22 with the JAX model and a batched FIRE + FrechetCellFilter of the authors (fmax 0.05 eV/Å, at most 500 steps, float64: the settings of the shared runner), converted to the canonical files with `matbench_discovery.discovery.write_discovery_artifacts` |
| Thermal conductivity | shared runner `models/run_kappa.py` (matbench-discovery `6d587c1`, Python 3.14.10, jax 0.12.0.1, cuequivariance-jax 0.13.0) with the registered JAX calculator (float64, fused cuEquivariance kernel), NVIDIA A100-SXM4-80GB, 8 shards, 2026-11-26/27 |
