# EquFlash

## WBM discovery

The published EquFlash predictions were produced with a batched Fair-Chem
`ml_relax` pipeline. The original submission metadata bundle remains archived at
https://figshare.com/articles/dataset/Equflash_Matbench_Discovery_submission_Metadata/32625693?file=65423763.

The legacy one-off discovery executable was retired when
[PR #375](https://github.com/janosh/matbench-discovery/pull/375) centralized
discovery execution. EquFlash is explicitly marked as archived because the
exercised Fair-Chem path requires installing `fairchem-core==1.10.0` with `--no-deps`:
Fair-Chem declares Torch `<2.5`, while EquFlash requires Torch 2.9.1. Adding that
special-case environment and batched relaxation to the shared runner would make
the common path less reliable.
