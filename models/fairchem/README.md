FAIRChem models
===============

Model checkpoints for the two models submitted and additional models trained are avaiable at the [FAIRChem HF org](https://huggingface.co/fairchem)

The OMat24 training and validation dataset splits used to pretrain non-compliant models can be downloaded [here](https://fair-chem.github.io/core/datasets/omat24.html)

The code and training config files necessary to train and evaluate the models is available in the official [FAIRChem repo](https://github.com/FAIR-Chem/fairchem)
- Currently the functionality to reproduce training results and run relaxations with stress enabled is on this [branch](https://github.com/FAIR-Chem/fairchem/tree/stress-relaxations), which we plan to merge into the main codebase in the coming weeks.

If you use any of our models or the OMat24 dataset in your work please cite our pre-print,


## Models

We submitted two equiformerV2 models:

1. **eqV2-31M-MP** - A compliant model trained on MPTrj only
2. **eqV2-31M-OMat-MP** - A non-compliant model pretrained with the OMat24 dataset and finetuned with MPTrj

### eqV2-31M-MP metrics
```txt

```

### eqV2-31M-MP metrics
```txt
```