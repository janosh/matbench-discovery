
## Orb Models from [Orbital Materials](https://orbitalmaterials.com)


All models described here are available for use via our Python package, [`orbitalmaterials`](https://github.com/orbital-materials/orbitalmaterials). 


### Models

The benchmark submission contains two models. These models are architecturally identical and use the same training proceedure, but differ in the training data used to train them.


*Note: The models in the benchmark submissions are not trained with the D3 dispersion correction. If you are interested in models with integrated D3 corrections, please see our github link above.*

- ORB - This is our full model, pretrained as a de-noising model on a wide variety of collected crystal structures. It is then finetuned using [`MPtraj`](https://figshare.com/articles/dataset/23713842) and [`Alexandria`](https://alexandria.icams.rub.de/) jointly.

- ORB-mptraj - This is our model which fits the data requirements of the Matbench Discovery benchmark - namely that it is trained only on the MPtraj dataset. This model uses the full MPTraj dataset for both the pretraining steps and the finetuning steps. It uses no other data sources.


For more information on the models, please refer to our [Github Repository](https://github.com/orbital-materials/orbitalmaterials) and technical blog post [here](TODO: Add link). We are working on a more complete technical report which will be available soon.
