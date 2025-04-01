# EScAIP: Efficiently Scaled Attention Interatomic Potential

EScAIP is an NNIP architecture designed to be scalable and efficient. It leverages a novel multi-head self-attention mechanism within graph neural networks, to enhance expressivity while maintaining computational efficiency. We submitted the model trained only on MPTrj, without using any denoising objectives during training.

The code and checkpoints can be found in [our repo](https://github.com/ASK-Berkeley/EScAIP). Our model is based on the infrastructure in [FAIRChem](https://github.com/FAIR-Chem/fairchem). All testing and evaluation scripts follow the FAIRChem settings ([ref](/models/eqV2)).

```bib
@inproceedings{
  qu_2024_importance,
  title={The Importance of Being Scalable: Improving the Speed and Accuracy of Neural Network Interatomic Potentials Across Chemical Domains},
  author={Eric Qu and Aditi S. Krishnapriyan},
  booktitle={The Thirty-eighth Annual Conference on Neural Information Processing Systems},
  year={2024},
  url={https://openreview.net/forum?id=Y4mBaZu4vy}
}
```
