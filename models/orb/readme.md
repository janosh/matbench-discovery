# Orb Models from [Orbital Materials](https://orbitalmaterials.com)

All models described here are available for use via our Python package, [`orb-models`].

[`orb-models`]: https://github.com/orbital-materials/orb-models

## Models

The benchmark submission contains two models. These models are architecturally identical and use the same training procedure, but differ in the training data used to train them.

_Note: The models in the benchmark submissions are not trained with the D3 dispersion correction. If you are interested in models with integrated D3 corrections, please see [`orb-models`]._

- **ORB**: `orb-v2` - This is our full model, pretrained as a de-noising model on a wide variety of collected crystal structures. It is then finetuned using [`MPtrj`](https://figshare.com/articles/dataset/23713842) and [`Alexandria`](https://alexandria.icams.rub.de/) jointly.

- **ORB MPtrj**: `orb-v2-mptrj-only` - This is our model which fits the data requirements of the Matbench Discovery benchmark - namely that it is trained only on the MPtrj dataset. This model uses the full MPtrj dataset for both the pretraining steps and the finetuning steps. It uses no other data sources.

For more information on the models, please refer to our [Github Repository][`orb-models`]. We are working on a more complete technical report which will be available soon.

## **ORB**: `orb-v2` - Full dataset pretraining, MPtrj + Alexandria finetuning

Public model weights:
https://orbitalmaterials-public-models.s3.us-west-1.amazonaws.com/forcefields/orb-v2-20241011.ckpt

## **ORB MPtrj**: `orb-v2-mptrj-only` - MPtrj pretraining, MPtrj finetuning

Public model weights:
https://orbitalmaterials-public-models.s3.us-west-1.amazonaws.com/forcefields/orb-mptraj-only-v2-20241014.ckpt

## V1 Models

- **ORB**: `orb-v1` - This is our full model, pretrained as a de-noising model on a wide variety of collected crystal structures. It is then finetuned using [`MPtrj`](https://figshare.com/articles/dataset/23713842) and [`Alexandria`](https://alexandria.icams.rub.de/) jointly.

- **ORB MPtrj**: `orb-v1-mptrj-only` - This is our model which fits the data requirements of the Matbench Discovery benchmark - namely that it is trained only on the MPtrj dataset. This model uses the full MPtrj dataset for both the pretraining steps and the finetuning steps. It uses no other data sources.

## **ORB**: `orb-v1` - Full dataset pretraining, MPtrj + Alexandria finetuning

Public model weights:
https://orbitalmaterials-public-models.s3.us-west-1.amazonaws.com/forcefields/orbff-v1-20240827.ckpt

## **ORB MPtrj**: `orb-v1-mptrj-only` - MPtrj pretraining, MPtrj finetuning

Public model weights:
https://orbitalmaterials-public-models.s3.us-west-1.amazonaws.com/forcefields/orbff-mptraj-only-v1-20240827.ckpt

## ORB-relaxed structures

ORB-relaxed structures can be downloaded from:

```bash
gsutil ls gs://orbitalmaterials-public-models/matbench-evals/
```

E.g. to download structures for a specific model, run:

```bash
gsutil cp gs://orbitalmaterials-public-models/matbench-evals/orb-mptraj-only-v2-20241014.json.gz models/orb/
```
