# FAIRChem EquiformerV2 models

## Models

Two equiformerV2 models submitted:
1. **eqV2 DeNS** - A compliant model trained on MPTrj only
2. **eqV2** - A non-compliant model pretrained with the OMat24 dataset and finetuned with MPTrj + sAlex

Model checkpoints for the two models submitted and additional models of different size are available [here](https://huggingface.co/fairchem/OMAT24).

## Datasets

The OMat24 training and validation dataset splits used to pre-train non-compliant models can be downloaded [here](https://huggingface.co/datasets/fairchem/OMAT24)
1. The splits do not contain any structure that has a protostructure label present in the initial or relaxed structures of the WBM dataset.
2. The splits do not include any structure that was generated starting from an Alexandria relaxed structure with protostructure label in the initial or relaxed structures of the
   WBM dataset.

We processed the following datasets for finetuning as follows
1. MPTrj - We use uncorrected energies and removed all structures with atoms at least 12 Angstroms apart.
2. sAlex - We removed all trajectories from the Alexandria 3D PBE dataset if any structure had a structure prototype that is in the set of all structure prototypes for all WBM initial and final structures. We then filtered the remaining trajectories by removing all structures with energy > 0 eV, forces norm > 50 eV/A, absolute stress > 80GPa. Finally we subsampled each trajectory by taking structures with energies at least 10 meV/atom apart.

The sAlex training and validation dataset splits used for fine-tuning are also available [here](https://huggingface.co/datasets/fairchem/OMAT24).

The code and training config files necessary to train and evaluate the models is available in the official [FAIRChem repo](https://github.com/FAIR-Chem/fairchem)

- Currently the functionality to reproduce training results and run relaxations with stress enabled is on this [branch](https://github.com/FAIR-Chem/fairchem/tree/stress-relaxations), which will be merged into the main codebase in the coming weeks.

Additional information about the datasets and models can be found in our pre-print:
```bib
@article{barroso_omat24,
  title={Open Materials 2024 (OMat24) Inorganic Materials Dataset and Models},
  author={Barroso-Luque, Luis and Muhammed, Shuaibi and Fu, Xiang and Wood, Brandon, Dzamba, Misko, and Gao, Meng and Rizvi, Ammar and  Zitnick, C. Lawrence and Ulissi, Zachary W.},
  journal={arXiv preprint arXiv:2410.12771},
  year={2024}
}
```