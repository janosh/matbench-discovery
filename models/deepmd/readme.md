# DeePMD-DPA3

## Model

- dpa3-openlam：model pretrained on the `OpenLAM` datasest, and finetuned with `MPtraj` & `Sub-Alex` datasets.

  ```bash
  wget https://figshare.com/files/52989059
  ```

- dpa3-mptrj：model trained only on the `MPtraj` dataset.

  ```bash
  wget https://figshare.com/files/52989056

  ```

### How to install

```bash
pip install git+https://github.com/deepmodeling/deepmd-kit.git@devel
```

### How to use

```py
from ase import Atoms
from deepmd.calculator import DP

water = Atoms(
   "H2O",
   positions=[(0.7601, 1.9270, 1), (1.9575, 1, 1), (1.0, 1.0, 1.0)],
   cell=[100, 100, 100],
   calculator=DP(model="2025-03-14-dpa3-mptrj.pth"),
)

print(water.get_potential_energy())
print(water.get_forces())
```

## OPENLAM Data

- To be released soon.

| Dataset Name               | Amount                                                                                      | Downsample                                                                                        | Data Cleaning                                                                                    | Weights | References                                                                                                                                                                                                                                                                                  |
| -------------------------- | ------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------ | ------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Domains_Alloy              | trn: 71482; val: 1240                                                                       |                                                                                                   |                                                                                                  | 1.0     | -                                                                                                                                                                                                                                                                                           |
| Domains_Anode              | trn: 88692; val: 9695                                                                       |                                                                                                   |                                                                                                  | 2.46    | -                                                                                                                                                                                                                                                                                           |
| Domains_Cluster            | trn: 143418; val: 15331                                                                     |                                                                                                   |                                                                                                  | 1.6     | -                                                                                                                                                                                                                                                                                           |
| Domains_Drug               | trn: 1379956; val: 24257                                                                    |                                                                                                   |                                                                                                  | 7.01    | -                                                                                                                                                                                                                                                                                           |
| Domains_FerroEle           | trn: 14487; val: 1357                                                                       |                                                                                                   |                                                                                                  | 1.15    | Jing Wu, Jiyuan Yang, Yuan-Jinsheng Liu, Duo Zhang, Yudi Yang, Yuzhi Zhang, Linfeng Zhang, Shi Liu, et al. Universal interatomic potential for perovskite oxides. Physical Review B, 108(18):L180104, 2023.                                                                                 |
| Domains_SSE_PBE            | trn: 17582; val: 886                                                                        |                                                                                                   |                                                                                                  | 2.58    | Jianxing Huang, Linfeng Zhang, Han Wang, Jinbao Zhao, Jun Cheng, and Weinan E. Deep potential generation scheme and simulation protocol for the li10gep2s12-type superionic conductors. The Journal of Chemical Physics, 154(9):094703, 2021.                                               |
| Domains_SemiCond           | trn: 215481; val: 23343                                                                     |                                                                                                   |                                                                                                  | 5.82    | Jianchuan Liu, Xingchen Zhang, Yuzhi Zhang, Duo Zhang, Linfeng Zhang, and Mohan Chen. Machine-learning-based interatomic potentials for group iib to via semiconductors: A comparative study of universal and independent models. arXiv preprint arXiv:2311.11305, 2023.                    |
| H2O_H2O_PD                 | trn: 46077; val: 2342                                                                       |                                                                                                   |                                                                                                  | 6.59    | Linfeng Zhang, Han Wang, Roberto Car, and Weinan E. Phase diagram of a deep potential water model. Physical review letters, 126(23):236001, 2021.                                                                                                                                           |
| Metals_AlMgCu              | trn: 138194; val: 3965                                                                      |                                                                                                   |                                                                                                  | 0.77    | Wanrun Jiang, Yuzhi Zhang, Linfeng Zhang, and Han Wang. Accurate deep potential model for the al–cu–mg alloy in the full concentration space. Chinese Physics B, 30(5):050706, 2021.                                                                                                        |
| Metals_Sn                  | trn: 6449; val: 276                                                                         |                                                                                                   |                                                                                                  | 0.05    | Tao Chen, Fengbo Yuan, Jianchuan Liu, Huayun Geng, Linfeng Zhang, Han Wang, and Mohan Chen. Modeling the high-pressure solid and liquid phases of tin from deep potentials with ab initio accuracy. Physical Review Materials, 7(5):053603, 2023.                                           |
| Metals_Ti                  | trn: 10054; val: 474                                                                        |                                                                                                   |                                                                                                  | 0.07    | Tongqi Wen, Rui Wang, Lingyu Zhu, Linfeng Zhang, Han Wang, David J Srolovitz, and Zhaoxuan Wu. Specialising neural network potentials for accurate properties and application to the mechanical response of titanium. npj Computational Materials, 7(1):206, 2021.                          |
| Metals_V                   | trn: 14935; val: 738                                                                        |                                                                                                   |                                                                                                  | 0.11    | Rui Wang, Xiaoxiao Ma, Linfeng Zhang, Han Wang, David J Srolovitz, Tongqi Wen, and Zhaoxuan Wu. Classical and machine learning interatomic potentials for bcc vanadium. Physical Review Materials, 6(11):113603, 2022.                                                                      |
| Metals_W                   | trn: 42297; val: 2100                                                                       |                                                                                                   |                                                                                                  | 0.18    | Xiaoyang Wang, Yinan Wang, Linfeng Zhang, Fuzhi Dai, and Han Wang. A tungsten deep neural-network potential for simulating mechanical property degradation under fusion service environment. Nuclear Fusion, 2022.                                                                          |
| Others_HfO2                | trn: 27660; val: 917                                                                        |                                                                                                   |                                                                                                  | 1.42    | Jing Wu, Yuzhi Zhang, Linfeng Zhang, and Shi Liu. Deep learning of accurate force field of ferroelectric hfo 2. Physical Review B, 103(2):024108, 2021.                                                                                                                                     |
| Domains_SSE_PBESol         | trn: 7502; val: 384                                                                         |                                                                                                   |                                                                                                  | 1.24    | Jianxing Huang, Linfeng Zhang, Han Wang, Jinbao Zhao, Jun Cheng, and Weinan E. Deep potential generation scheme and simulation protocol for the li10gep2s12-type superionic conductors. The Journal of Chemical Physics, 154(9):094703, 2021.                                               |
| Domains_Transition1x       | trn: 7632328; val: 967454                                                                   |                                                                                                   |                                                                                                  | 3.42    | Mathias Schreiner, Arghya Bhowmik, Tejs Vegge, Jonas Busk, and Ole Winther. Transition1x a dataset for building generalizable reactive machine learning potentials. Scientific Data, 9(1):779, 2022.                                                                                        |
| Metals_AgAu_PBED3          | trn: 64239; val: 2256                                                                       |                                                                                                   |                                                                                                  | 0.39    | YiNan Wang, LinFeng Zhang, Ben Xu, XiaoYang Wang, and Han Wang. A generalizable machine learning potential of ag–au nanoalloys and its application to surface reconstruction, segregation and diffusion. Modelling and Simulation in Materials Science and Engineering, 30(2):025003, 2021. |
| Others_In2Se3              | trn: 11621; val: 568                                                                        |                                                                                                   |                                                                                                  | 0.21    | Jing Wu, Liyi Bai, Jiawei Huang, Liyang Ma, Jian Liu, and Shi Liu. Accurate force field of two-dimensional ferroelectrics from deep learning. Physical Review B, 104(17):174107, 2021.                                                                                                      |
| MP_traj_v024_alldata_mixu  | trn: 1401956; val: 110918                                                                   |                                                                                                   |                                                                                                  | 6.99    | https://www.nature.com/articles/s42256-023-00716-3                                                                                                                                                                                                                                          |
| Alloy_tongqi               | trn: 24097; val: 100                                                                        |                                                                                                   |                                                                                                  | 0.63    | -                                                                                                                                                                                                                                                                                           |
| SSE_ABACUS                 | trn: 125083; val: 6587                                                                      |                                                                                                   |                                                                                                  | 3.05    | -                                                                                                                                                                                                                                                                                           |
| Hybrid_Perovskite          | trn: 48078; val: 2530                                                                       |                                                                                                   |                                                                                                  | 0.88    | Tuo P, Li L, Wang X, et al. Hybrid nano-domain structures of organic-inorganic perovskites from molecule-cage coupling effects\[J]. arXiv preprint arXiv:2209.12445, 2022                                                                                                                   |
| solvated_protein_fragments | trn: 2594609; val: 136571                                                                   |                                                                                                   |                                                                                                  | 6.14    | https://pubs.acs.org/doi/10.1021/acs.jctc.9b00181                                                                                                                                                                                                                                           |
| Electrolyte                | trn: 65393; val: 3438                                                                       |                                                                                                   |                                                                                                  | 8.78    | https://www.aissquare.com/datasets/detail?name=Electrolyte\&id=216\&pageType=datasets                                                                                                                                                                                                       |
| ODAC23_new                 | trn: 2682332Vaid: 63623                                                                     | The first three frames, the last three frames of the track, and every 20 frames of the trajectory | Removed frames with energy/atom > 0.5eV/atom, energy/atom < -0.2eV/atom, max_abs_force > 25eV/A  | 1.78    | https://pubs.acs.org/doi/10.1021/acscentsci.3c01629                                                                                                                                                                                                                                         |
| Alex2D_new                 | trn: 1223831 val: 135810                                                                    | The first three frames, the last three frames of the track, and every 10 frames of the trajectory | Removed frames with energy/atom > 0eV/atom, max_abs_force > 5eV/A,max_abs_virial /atom> 8eV/atom | 2.75    | https://iopscience.iop.org/article/10.1088/2053-1583/accc43                                                                                                                                                                                                                                 |
| OMAT24                     | trn: 100,568,820 (removed sub-alex, trn: 10M, val 0.5M） val: 1,074,647 (removed sub-alex） | -                                                                                                 | Removed frames with energy/atom > 0eV/atom, energy/atom < -25eV/atom, max_abs_force > 50eV/A     | 33.24   | https://arxiv.org/html/2410.12771v1                                                                                                                                                                                                                                                         |
| SPICE2_new                 | trn: 1,621,168 val: 180089                                                                  | -                                                                                                 | Removed frames with energy/atom <-10000eV/atom, max_abs_force > 15eV/A                           | 8.08    | https://arxiv.org/pdf/2406.13112                                                                                                                                                                                                                                                            |
| OC20M                      | trn: 20,000,000 val: 999,866                                                                | -                                                                                                 | -                                                                                                | 17.48   | Lowik Chanussot, Abhishek Das, Siddharth Goyal, Thibaut Lavril, Muhammed Shuaibi, Morgane Riviere, Kevin Tran, Javier Heras-Domingo, Caleb Ho, Weihua Hu, et al. Open catalyst 2020 (oc20) dataset and community challenges. ACS Catalysis, 11(10):6059–6072, 2021.                         |
| OC22                       | trn: 8,194,770 val: 394,727                                                                 | -                                                                                                 | Removed frames with energy/atom > 0eV/atom, max_abs_force > 10eV/A,                              | 12.21   | https://pubs.acs.org/doi/10.1021/acscatal.2c05426                                                                                                                                                                                                                                           |
| Organic_Reactions          | trn: 14,024,587 val: 1,558,260                                                              | -                                                                                                 | Removed frames with max_abs_force > 20eV/A                                                       | 11.93   | Schreiner, Mathias, et al. "Transition1x - a dataset for building generalizable reactive machine learning potentials." _Scientific Data_ 9.1 (2022): 779.（Active learning based on Transition-1x）                                                                                         |

&#x20;

## Results

### dpa3-openlam

`2025-01-10-dpa3-openlam-preds.csv`

```txt
    Full-set    Unique  10K
F1              0.869284       0.890341     0.986469
DAF             4.976639       5.747061     6.366780
Precision       0.853936       0.878562     0.973300
Recall          0.885195       0.902439     1.000000
Accuracy        0.954320       0.965571     0.973300
TPR             0.885195       0.902439     1.000000
FPR             0.031362       0.022859     1.000000
TNR             0.968638       0.977141     0.000000
FNR             0.114805       0.097561     0.000000
TP          39030.000000   30118.000000  9733.000000
FP           6676.000000    4163.000000   267.000000
TN         206195.000000  177951.000000     0.000000
FN           5062.000000    3256.000000     0.000000
MAE             0.021902       0.022137     0.018143
RMSE            0.068287       0.066571     0.062670
R2              0.856952       0.869308     0.912547
```

### dpa3-mptrj

`2025-01-10-dpa3-mptrj-pred.csv`

```txt
    Full-set    Unique  10K
F1              0.774286       0.785753     0.979592
DAF             4.249362       4.822293     6.279779
Precision       0.729143       0.737192     0.960000
Recall          0.825388       0.841164     1.000000
Accuracy        0.917428       0.928957     0.960000
TPR             0.825388       0.841164     1.000000
FPR             0.063508       0.054955     1.000000
TNR             0.936492       0.945045     0.000000
FNR             0.174612       0.158836     0.000000
TP          36393.000000   28073.000000  9600.000000
FP          13519.000000   10008.000000   400.000000
TN         199352.000000  172106.000000     0.000000
FN           7699.000000    5301.000000     0.000000
MAE             0.037755       0.039197     0.032419
RMSE            0.081581       0.081490     0.078476
R2              0.795834       0.804168     0.865694
```

### Relaxed Structure

```sh
# dpa3-openlam
wget https://figshare.com/files/52988759

# dpa3-mptrj
wget https://figshare.com/files/52988747
```
