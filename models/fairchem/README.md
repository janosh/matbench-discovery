FAIRChem models
===============

Model checkpoints for the two models submitted and additional models trained are avaiable [here](https://huggingface.co/fairchem/OMAT24)

The OMat24 training and validation dataset splits used to pretrain non-compliant models can be downloaded [here](https://fair-chem.github.io/core/datasets/omat24.html)
The OMat24 train and val splits are fully compatible with the Matbench-Discovery test set.
   1. The splits do not contain any structure that has a protostructure label present in the initial or relaxed structures of the WBM dataset.
   2. The splits do not include any structure that was generated starting from an Alexandria relaxed structure with protostructure lable in the intitial or relaxed structures of the
      WBM datset.

The code and training config files necessary to train and evaluate the models is available in the official [FAIRChem repo](https://github.com/FAIR-Chem/fairchem)
- Currently the functionality to reproduce training results and run relaxations with stress enabled is on this [branch](https://github.com/FAIR-Chem/fairchem/tree/stress-relaxations), which will be merged into the main codebase in the coming weeks.

If you use any of our models or the OMat24 dataset in your work please cite our pre-print,


## Models

We submitted two equiformerV2 models:

1. **eqV2-31M-dens-mp** - A compliant model trained on MPTrj only
2. **eqV2-86M-omat-mp** - A non-compliant model pretrained with the OMat24 dataset and finetuned with MPTrj

### eqV2-31M-mp metrics
```txt

```

### eqV2-86M-dens-mp metrics
```txt
```