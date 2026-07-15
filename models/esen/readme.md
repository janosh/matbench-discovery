# FAIRChem eSEN models

## Models

Two equivariant Smooth Energy Network (eSEN) models are submitted:

1. **eSEN-30M-MP** - A compliant model trained on MPtrj only
2. **eSEN-30M-OAM** - A non-compliant model pretrained with the OMat24 dataset and finetuned with MPtrj + sAlex

The code and training config is currently WIP and will be available in the coming weeks in the official [FAIRChem repo](https://github.com/FAIR-Chem/fairchem). The checkpoints will be made available in the official FAIRChem [HuggingFace site](https://huggingface.co/fairchem).

Additional information about eSEN can be found in our paper:

```bib
@article{fu2025learning,
  title={Learning Smooth and Expressive Interatomic Potentials for Physical Property Prediction},
  author={Fu, Xiang and Wood, Brandon M and Barroso-Luque, Luis and Levine, Daniel S and Gao, Meng and Dzamba, Misko and Zitnick, C Lawrence},
  journal={arXiv preprint arXiv:2502.12147},
  year={2025}
}
```
