 #!/bin/bash

export PYTHONPATH=/data/fywang/code/fairchem-main/src:$PYTHONPATH

#  python models/eqV2/join_fairchem_preds.py \
#     --input-path /data/fywang/code/matbench-discovery/pred_results/eqV2-small-dens/test/ \
#     --out-path /data/fywang/code/matbench-discovery/pred_results/eqV2-small-dens/test/ \
#     --model-name eqV2-s-dens-mp \
#     --apply-mp-corrections

 python join_fairchem_preds.py \
    --apply-mp-corrections \
    --model-name eqV2-small-dens-tracegrad \
    --input-path /data/fywang/code/matbench-discovery/pred_results/eqV2-small-dens-tracegrad/2025-06-19@03-43-53/ \
    --out-path /data/fywang/code/matbench-discovery/pred_results/eqV2-small-dens-tracegrad/2025-06-19@03-43-53/ \
