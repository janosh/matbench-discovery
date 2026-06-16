# EquFlash

This folder contains the EquFlash model submission scripts and metadata for Matbench Discovery.

## WBM discovery input files

`test_equflash_discovery.py` requires two WBM input paths:

```bash
--init-structs-dir
--wbm-metadata-file
```

These files can be downloaded from the Figshare metadata bundle:

https://figshare.com/articles/dataset/Equflash_Matbench_Discovery_submission_Metadata/32625693?file=65423763

After downloading and extracting the files, run the discovery script by passing the extracted WBM initial structures directory to `--init-structs-dir` and the WBM metadata file path to `--wbm-metadata-file`.

Example:

```bash
python models/equflash/test_equflash_discovery.py \
  --checkpoint /path/to/checkpoint.pt \
  --out /path/to/output_dir \
  --init-structs-dir /path/to/wbm_initial_structures \
  --wbm-metadata-file /path/to/wbm_metadata.csv
```

Additional optional arguments include:

```bash
--rank 0
--worldsize 1
--fmax 0.02
--max-atoms 1000
```
