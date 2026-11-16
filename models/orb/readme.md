# Orb Models from [Orbital Materials](https://orbitalmaterials.com)

All models described here are available for use via our Python package, [`orb-models`].

[`orb-models`]: https://github.com/orbital-materials/orb-models

## Models

The benchmark submission contains two models. These models are architecturally identical and use the same training procedure, but differ in the training data used to train them.

_Note: The models in the benchmark submissions are not trained with the D3 dispersion correction. If you are interested in models with integrated D3 corrections, please see [`orb-models`]._

- **ORB**: `orb-v2` - This is our full model, pretrained as a de-noising model on a wide variety of collected crystal structures. It is then finetuned using [`MPtrj`](https://figshare.com/articles/dataset/23713842) and [`Alexandria`](https://alexandria.icams.rub.de/) jointly.

- **ORB MPtrj**: `orb-v2-mptrj-only` - This is our model which fits the data requirements of the Matbench Discovery benchmark - namely that it is trained only on the MPtrj dataset. This model uses the full MPtrj dataset for both the pretraining steps and the finetuning steps. It uses no other data sources.

For more information on the models, please refer to our [Github Repository][`orb-models`]. We are working on a more complete technical report which will be available soon.

## **ORB**: `orb-v2` - Full dataset pretraining, MPtrj + Alexandria finetuning

Public model weights:
https://orbitalmaterials-public-models.s3.us-west-1.amazonaws.com/forcefields/orb-v2-20241011.ckpt

## **ORB MPtrj**: `orb-v2-mptrj-only` - MPtrj pretraining, MPtrj finetuning

Public model weights:
https://orbitalmaterials-public-models.s3.us-west-1.amazonaws.com/forcefields/orb-mptraj-only-v2-20241014.ckpt

## V1 Models

- **ORB**: `orb-v1` - This is our full model, pretrained as a de-noising model on a wide variety of collected crystal structures. It is then finetuned using [`MPtrj`](https://figshare.com/articles/dataset/23713842) and [`Alexandria`](https://alexandria.icams.rub.de/) jointly.

- **ORB MPtrj**: `orb-v1-mptrj-only` - This is our model which fits the data requirements of the Matbench Discovery benchmark - namely that it is trained only on the MPtrj dataset. This model uses the full MPtrj dataset for both the pretraining steps and the finetuning steps. It uses no other data sources.

## **ORB**: `orb-v1` - Full dataset pretraining, MPtrj + Alexandria finetuning

Public model weights:
https://orbitalmaterials-public-models.s3.us-west-1.amazonaws.com/forcefields/orbff-v1-20240827.ckpt

## **ORB MPtrj**: `orb-v1-mptrj-only` - MPtrj pretraining, MPtrj finetuning

Public model weights:
https://orbitalmaterials-public-models.s3.us-west-1.amazonaws.com/forcefields/orbff-mptraj-only-v1-20240827.ckpt

## ORB-relaxed structures

ORB-relaxed structures can be downloaded from:

```bash
gsutil ls gs://orbitalmaterials-public-models/matbench-evals/
```

E.g. to download structures for a specific model, run:

```bash
gsutil cp gs://orbitalmaterials-public-models/matbench-evals/orb-mptraj-only-v2-20241014.json.gz models/orb/
```

## ORB-G2 GroupCV v1 (independent calibration submission)

[Source, fixed data, plans and reproduction](https://github.com/nkwork9999/orb-g2-groupcv) / [release v0.1.3](https://github.com/nkwork9999/orb-g2-groupcv/releases/tag/v0.1.3). The frozen ORB-v3 conservative MPA backbone is credited separately above; this calibration is by nk (Enuke).

A 36-coefficient conservative pair-energy head is fitted on 83 multi-atom ASE G2 molecules (68 formula groups). Training-only five-fold formula-group CV selects width/ridge/shrinkage with a fixed 16-trial GP-EI budget. A metadata-only split reserves 19 molecules in 18 disjoint formula groups; all one/two-atom molecules are excluded. Targets are experimental atomization energies extrapolated with ZPE/thermal corrections, not PBE force labels. The fixed molecular holdout macro-MAE changes from 0.4933 to 0.2175 eV (15 groups improve, 3 worsen). Upstream pretraining overlap has not been proven absent.

The head supplies general multi-atom/periodic energies, forces and stress. Only Diatomics is benchmarked here; no bulk/force/MD accuracy improvement is claimed. Historical Diatomics feedback informed redesign, so `test_set_informed: [diatomics]` is mandatory. Official labels and outer scores never select this trained state. Runtime-only v0.1.1 preserves all learned coefficients, selecting exact fast neighbor enumeration and four CPU threads using training structures only. Compilation failed the prospective speed threshold and is disabled.

The model-specific uv project documents the Python 3.14 dm-tree/filelock compatibility overrides. A portable training wrapper is `models/orb/orb-g2-groupcv-v1/train.py`; run it against the unpacked, SHA256-verified release data archive. It reproduces the fixed model using training rows only. Upstream weights remain separately downloaded and hash-verified.

Deployment v0.1.3 uses MLX float32 geometry, ZBL and dense inference on Apple GPU, with CPU graph construction and conservative edge-vector derivatives. All learned weights stay unchanged.183 training/synthetic cases meet atol1e-4+rtol1e-5. Compilation/warmup is included in load time. The common CLI adds `mps` to match the existing schema, and run-cost instrumentation measures real MLX GPU memory. The official sampling/scoring code is unchanged. Run `uv run --project models/orb/orb-g2-groupcv-v1 --isolated python models/run_diatomics.py --model orb_g2_groupcv_v1 --device mps --dtype float32`; use `--device cpu` for the numerically equivalent CPU reference. Apple GPU evaluation requires an Apple Silicon Mac; CPU reproduction is portable.
