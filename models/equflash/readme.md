# EquFlash

This folder contains the EquFlash model submission support code and metadata for Matbench Discovery.

## WBM discovery

The published EquFlash predictions use the retained custom batched relaxation
pipeline in `test_equflash_discovery.py`. It packs structures by atom count and
calls Fair-Chem's `ml_relax`, so it is not equivalent to the shared per-structure
ASE relaxation loop.

```bash
./models/equflash/install.sh
python models/equflash/test_equflash_discovery.py \
  --checkpoint CHECKPOINT --out OUT_DIR \
  --init-structs-dir WBM_INITIAL_STRUCTURES \
  --wbm-metadata-file WBM_METADATA_CSV
```

`install.sh` intentionally installs `fairchem-core==1.10.0` without its pinned
Torch dependencies so EquFlash can use Torch 2.9.1. This environment therefore
cannot be expressed by the shared runner's ordinary `uv run --with` resolution.

The original submission metadata bundle remains archived at
https://figshare.com/articles/dataset/Equflash_Matbench_Discovery_submission_Metadata/32625693?file=65423763.
