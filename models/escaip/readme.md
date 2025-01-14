# EScAIP: Efficiently Scaled Attention Interatomic Potential

EScAIP is an NNIP architecture designed to be scalable and efficient. It leverages a novel multi-head self-attention mechanism within graph neural networks, to enhance expressivity while maintaining computational efficiency. We submitted the model trained on MPTrj, without using denoising objectives.

The code and checkpoints can be found in [our repo](https://github.com/ASK-Berkeley/EScAIP). Our model is heavily based on [FAIRChem](https://github.com/FAIR-Chem/fairchem). All testing and evaluation scripts follows the FAIRChem settings ([ref](/models/eqV2)).

```bib
@inproceedings{
qu2024the,
title={The Importance of Being Scalable: Improving the Speed and Accuracy of Neural Network Interatomic Potentials Across Chemical Domains},
author={Eric Qu and Aditi S. Krishnapriyan},
booktitle={The Thirty-eighth Annual Conference on Neural Information Processing Systems},
year={2024},
url={https://openreview.net/forum?id=Y4mBaZu4vy}
}
```
