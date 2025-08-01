 #!/bin/bash
export WANDB_MODE=offline
export PYTHONPATH=/data/fywang/code/fairchem-tracegrad-hanhai/src:$PYTHONPATH

 python test_fairchem_tracegrad_discovery.py \
    --out-path /data/fywang/code/matbench-discovery/pred_results/ \
    --model-name eqV2-small-dens-tracegrad \
    --optimizer FIRE \
    --cell-filter unit \
    --force-max 0.02 \
    --max-steps 500 \
    --num-jobs 10 \
    --device cuda \
    --checkpoint-path /data/fywang/code/fairchem-tracegrad-hanhai/checkpoints/2025-06-16-14-11-12/best_checkpoint.pt \

#  python models/eqV2_tracegrad/test_fairchem_tracegrad_discovery_debug.py \
#     --checkpoint-path /data/fywang/code/fairchem-tracegrad-hanhai/checkpoints/2025-06-16-14-11-12/best_checkpoint.pt \
#     --out-path /data/fywang/code/matbench-discovery/pred_results/ \
#     --model-name eqV2-small-dens-tracegrad \
#     --optimizer FIRE \
#     --cell-filter unit \
#     --force-max 0.02 \
#     --max-steps 500 \
#     --num-jobs 1\
#     --debug \
