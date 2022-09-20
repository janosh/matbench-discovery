# File Index

The `block_{0,1}.p` files in this directory were obtained on 2022-08-05 from [figshare](https://figshare.com/articles/dataset/MPF_2021_2_8/19470599) as referenced in the [M3GNet repo](https://github.com/materialsvirtuallab/m3gnet/#datasets).

## Figshare description

Dataset posted on 30.03.2022, 21:14 authored by Chi ChenChi Chen, Shyue Ping Ong
This dataset contains the MPF.2021.2.8 data used to train the M3GNet model reported in `https://arxiv.org/abs/2202.02450`

I have split the dataset into two pickle files. To load the data, you can use example code as below.

```py
import pickle

with open('block_0.p', 'rb') as f:
    data = pickle.load(f)

with open('block_1.p', 'rb') as f:
    data.update(pickle.load(f))
```

where `data` will be a dictionary with `material_id` as the key and an inner dictionary as the value.

The inner dictionary contains the snapshots of this `material_id`, with the following keys.

- structure
- energy
- force
- stress
- id

Each id in the `id` list is of format `material_id-calc_id-ionic_step_id`, where `calc_id` is 0 (second) or 1 (first) in the double relaxation process.

The `stress` here is the raw output from VASP, meaning that it is really the negative stress using the convention in our paper. Hence to train the model, please multiply stress with -0.1 (kBa to GPa and change sign)

The units for energy, force and stress in the data are eV, eV/A, and kBa. Remember to convert the stress to GPa and take the negative sign to work with `m3gnet` training.
