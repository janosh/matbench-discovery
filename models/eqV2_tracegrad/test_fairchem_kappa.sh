export PYTHONPATH=/data/fywang/code/fairchem-tracegrad-hanhai/src:$PYTHONPATH

python test_fairchem_kappa.py \
  --checkpoint-path /data/fywang/code/fairchem-tracegrad-hanhai/checkpoints/2025-06-16-14-11-12/best_checkpoint.pt \
  --out-path /data/fywang/code/matbench-discovery/pred_results/ \
  --model-name eqV2-small-dens-tracegrad \
  --identifier kappa-103 \
  --atom-disp 0.01 \
  --seed 42 \
  --num-jobs 10 \
  # --checkpoint-path /data/fywang/code/fairchem-main/checkpoints/2025-01-20-10-54-56/best_checkpoint.pt \

# python models/eqV2_tracegrad/test_fairchem_kappa_debug.py \
#   --checkpoint-path /data/fywang/code/fairchem-tracegrad-hanhai/checkpoints/2025-04-01-10-12-16/best_checkpoint.pt \
#   --out-path /data/fywang/code/matbench-discovery/pred_results/ \
#   --model-name eqV2-small-dens-tracegrad \
#   --identifier kappa-103 \
#   --atom-disp 0.01 \
#   --seed 42 \
#   --num-jobs 1\
#   --debug \
