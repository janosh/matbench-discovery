# PET

PET (Potential Energy Transformer) metatrain models.

## PET-XS Pair-BO v1

[PET-XS Pair-BO](https://github.com/nkwork9999/pet-xs-pair-bo) combines the frozen
PET-OMat-XS v1.0.0 backbone with a 36-coefficient conservative radial correction.
A Gaussian-process expected-improvement search selects the correction's fitting
hyperparameters on independently generated heteronuclear validation groups.

This entry is an isolated one/two-atom specialist for the **diatomics** task.
One atom receives no correction; a pair receives the correction only when both
atomic numbers lie in 1–17. Other elements retain PET's output, and systems with
more than two atoms are rejected. Periodic inputs must have sufficient vacuum;
the shared runner's 50 Å cubic cell is supported. Stress, bulk relaxation,
phonons and general materials dynamics are outside its scope.

The PET weights run on CPU in float32; the analytic residual uses float64.
The adapter keeps all three Cartesian force components from the base model and
adds the residual energy gradient. It does not remove PET's rotational or
transverse-force errors, and it applies no reference-dependent energy shift.

Use the registered shared runner with the declared precision and device:

```sh
uv run models/run_diatomics.py --print-cmd --model pet-xs-pair-bo-v1 --device cpu --dtype float32
```

The checkpoint URL resolves to a small frozen JSON manifest. The adapter fetches
and SHA256-verifies the separately hosted upstream PET checkpoint and caches it
under `~/.cache/pet-xs-pair-bo/checkpoints/`. Both the new code/residual and the
upstream PET weights use BSD-3-Clause licensing, with upstream attribution
preserved in the author repository.

Training and reference-generation code, configuration hashes, data splits and
quality checks are in the author repository. `train_pet_xs_pair_bo.py` delegates
to that package's single training implementation; pass `--project-root` pointing
to a checkout of the author repository so provenance is resolved there. Training
requires the versioned independent reference archive and its quality manifest.

The failed v1 pilot and revision 2.1 reference generation are retained as failed
artifacts. Reference protocol revision 2.2 repairs the generic solver on training
chemistry while keeping the same 16 chemical groups, 176 geometries and quality
thresholds. LiH/BeS/MgS remain the auxiliary partition, without replacement; their
QA observations have informed solver repair, so their later metrics are
descriptive, not an unseen confirmatory audit. Candidate training/search does
not use official homonuclear labels or feedback, and the first official run
follows a source/state/runtime freeze. This is not a guarantee of complete
benchmark non-exposure in the upstream backbone's pretraining history.
