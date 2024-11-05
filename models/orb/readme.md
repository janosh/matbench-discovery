# Orb Models from [Orbital Materials](https://orbitalmaterials.com)

All models described here are available for use via our Python package, [`orb-models`].

[`orb-models`]: https://github.com/orbital-materials/orb-models

## Models

The benchmark submission contains two models. These models are architecturally identical and use the same training procedure, but differ in the training data used to train them.

_Note: The models in the benchmark submissions are not trained with the D3 dispersion correction. If you are interested in models with integrated D3 corrections, please see [`orb-models`]._

- **ORB**: `orb-v2` - This is our full model, pretrained as a de-noising model on a wide variety of collected crystal structures. It is then finetuned using [`MPtrj`](https://figshare.com/articles/dataset/23713842) and [`Alexandria`](https://alexandria.icams.rub.de/) jointly.

- **ORB MPtrj**: `orb-v2-mptrj-only` - This is our model which fits the data requirements of the Matbench Discovery benchmark - namely that it is trained only on the MPtrj dataset. This model uses the full MPTrj dataset for both the pretraining steps and the finetuning steps. It uses no other data sources.

For more information on the models, please refer to our [Github Repository][`orb-models`]. We are working on a more complete technical report which will be available soon.

## **ORB**: `orb-v2` - Full dataset pretraining, MPtrj + Alexandria finetuning

| Item                 | Link/Url                                                                                       |
| -------------------- | ---------------------------------------------------------------------------------------------- |
| Results CSV          | [orbff-v2-20241011.csv.gz](./orbff-v2-20241011.csv.gz)                                         |
| Public Model Weights | https://storage.googleapis.com/orbitalmaterials-public-models/forcefields/orb-v2-20241011.ckpt |

```txt
                     orb          10k         unique
F1              0.859644     0.984978       0.880185
DAF             5.413980     6.347810       6.035050
Precision       0.902284     0.970400       0.922588
Recall          0.820852     1.000000       0.841509
Accuracy        0.955328     0.970400       0.964977
TPR             0.820852     1.000000       0.841509
FPR             0.017778     1.000000       0.012742
TNR             0.982222     0.000000       0.987258
FNR             0.179148     0.000000       0.158491
TP          35153.000000  9704.000000   27721.000000
FP           3807.000000   296.000000    2326.000000
TN         210331.000000     0.000000  180220.000000
FN           7672.000000     0.000000    5221.000000
MAE             0.028499     0.019454       0.028244
RMSE            0.079011     0.068514       0.077203
R2              0.808497     0.896680       0.824233
```

## **ORB MPtrj**: `orb-v2-mptrj-only` - MPTrj pretraining, MPTrj finetuning

| Item                 | Link/URL                                                                                                   |
| -------------------- | ---------------------------------------------------------------------------------------------------------- |
| Results CSV          | [orbff-mptrj-only-v2-20241014.csv.gz](./orbff-mptrj-only-v2-20241014.csv.gz)                               |
| Public Model Weights | https://storage.googleapis.com/orbitalmaterials-public-models/forcefields/orb-mptraj-only-v2-20241014.ckpt |

```txt
                     orb          10k         unique
F1              0.753764     0.970929       0.764406
DAF             4.279928     6.171845       4.694615
Precision       0.713285     0.943500       0.717673
Recall          0.799113     1.000000       0.817649
Accuracy        0.912987     0.943500       0.922952
TPR             0.799113     1.000000       0.817649
FPR             0.064239     1.000000       0.058046
TNR             0.935761     0.000000       0.941954
FNR             0.200887     0.000000       0.182351
TP          34222.000000  9435.000000   26935.000000
FP          13756.000000   565.000000   10596.000000
TN         200382.000000     0.000000  171950.000000
FN           8603.000000     0.000000    6007.000000
MAE             0.043234     0.037052       0.044555
RMSE            0.090746     0.097382       0.090975
R2              0.747388     0.801666       0.755931
```

## V1 Models

- **ORB**: `orb-v1` - This is our full model, pretrained as a de-noising model on a wide variety of collected crystal structures. It is then finetuned using [`MPtrj`](https://figshare.com/articles/dataset/23713842) and [`Alexandria`](https://alexandria.icams.rub.de/) jointly.

- **ORB MPtrj**: `orb-v1-mptrj-only` - This is our model which fits the data requirements of the Matbench Discovery benchmark - namely that it is trained only on the MPtrj dataset. This model uses the full MPTrj dataset for both the pretraining steps and the finetuning steps. It uses no other data sources.

For more information on the models, please refer to our [Github Repository][`orb-models`]. We are working on a more complete technical report which will be available soon.

## **ORB**: `orb-v1` - Full dataset pretraining, MPtrj + Alexandria finetuning

| Item                 | Link/Url                                                                                         |
| -------------------- | ------------------------------------------------------------------------------------------------ |
| Results CSV          | [orbff-v1-20240827.csv.gz](./orbff-v1-20240827.csv.gz)                                           |
| Public Model Weights | https://storage.googleapis.com/orbitalmaterials-public-models/forcefields/orbff-v1-20230827.ckpt |

```txt
                      orb          10k         unique
F1              0.846577     0.988213       0.867282
DAF             5.394101     6.389021       6.015771
Precision       0.898971     0.976700       0.919641
Recall          0.799953     1.000000       0.820563
Accuracy        0.951678     0.976700       0.961608
TPR             0.799953     1.000000       0.820563
FPR             0.017979     1.000000       0.012939
TNR             0.982021     0.000000       0.987061
FNR             0.200047     0.000000       0.179437
TP          34258.000000  9767.000000   27031.000000
FP           3850.000000   233.000000    2362.000000
TN         210288.000000     0.000000  180184.000000
FN           8567.000000     0.000000    5911.000000
MAE             0.030884     0.019012       0.030589
RMSE            0.080986     0.064470       0.079003
R2              0.798803     0.907903       0.815941
```

## **ORB MPtrj**: `orb-v1-mptrj-only` - MPTrj pretraining, MPTrj finetuning

| Item                 | Link/URL                                                                                                     |
| -------------------- | ------------------------------------------------------------------------------------------------------------ |
| Results CSV          | [orbff-mptrj-only-v1-20240827.csv.gz](./orbff-mptrj-only-v1-20240827.csv.gz)                                 |
| Public Model Weights | https://storage.googleapis.com/orbitalmaterials-public-models/forcefields/orbff-mptraj-only-v1-20230827.ckpt |

```txt
                     orb          10k         unique
F1              0.752143     0.963193       0.761336
DAF             4.267540     6.076994       4.667345
Precision       0.711221     0.929000       0.713505
Recall          0.798062     1.000000       0.816040
Accuracy        0.912341     0.929000       0.921787
TPR             0.798062     1.000000       0.816040
FPR             0.064804     1.000000       0.059130
TNR             0.935196     0.000000       0.940870
FNR             0.201938     0.000000       0.183960
TP          34177.000000  9290.000000   26882.000000
FP          13877.000000   710.000000   10794.000000
TN         200261.000000     0.000000  171752.000000
FN           8648.000000     0.000000    6060.000000
MAE             0.044745     0.040998       0.046230
RMSE            0.093426     0.102950       0.093919
R2              0.732243     0.780546       0.739879
```

## ORB-relaxed structures

ORB-relaxed structures can be downloaded from:

```bash
~ $ gsutil ls gs://orbitalmaterials-public-models/matbench-evals/
gs://orbitalmaterials-public-models/matbench-evals/
gs://orbitalmaterials-public-models/matbench-evals/orb-mptraj-only-v2-20241014.json.gz
gs://orbitalmaterials-public-models/matbench-evals/orb-v2-20241011.json.gz
gs://orbitalmaterials-public-models/matbench-evals/orbff-mptraj-only-v1-20240827.json.gz
gs://orbitalmaterials-public-models/matbench-evals/orbff-v1-20240827.json.gz
```

E.g. to download structures for a specific model, run:

```bash
gsutil cp gs://orbitalmaterials-public-models/matbench-evals/orb-mptraj-only-v2-20241014.json.gz matbench-discovery/models/orb/
```
