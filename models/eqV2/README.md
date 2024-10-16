# FAIRChem EquiformerV2 models

Model checkpoints for the two models submitted and additional models trained are available [here](https://huggingface.co/fairchem/OMAT24).

The OMat24 training and validation dataset splits used to pre-train non-compliant models can be downloaded [here](https://huggingface.co/datasets/fairchem/OMAT24)
The OMat24 train and val splits are fully compatible with the Matbench-Discovery test set.

1. The splits do not contain any structure that has a protostructure label present in the initial or relaxed structures of the WBM dataset.
2. The splits do not include any structure that was generated starting from an Alexandria relaxed structure with protostructure label in the initial or relaxed structures of the
   WBM dataset.

## Models

We submitted two equiformerV2 models:

1. **eqV2-S DeNS**
2. **eqV2-L**
