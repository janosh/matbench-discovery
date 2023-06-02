#%% Imports
### ----------------------------------------------------------------------------

import torch
import os
import csv

from monty.serialization import dumpfn
from tqdm import tqdm

from alignn.models.alignn import ALIGNN
from alignn.config import TrainingConfig

from jarvis.core.atoms import Atoms
from jarvis.core.graphs import Graph
from jarvis.db.jsonutils import loadjson

from matbench_discovery.data import df_wbm

#%% Cuda setup
### ----------------------------------------------------------------------------

device = "cpu"
if torch.cuda.is_available():
    device = torch.device("cuda")

#%% ALIGNN config
### ----------------------------------------------------------------------------

config = loadjson("alignn_config.json")
config = TrainingConfig(**config)

#%% Import an ALIGNN model
### ----------------------------------------------------------------------------

def load_model(basename):

    filename = os.path.join(basename, 'model_best.pth')

    model = ALIGNN(config.model)
    model.load_state_dict(torch.load(filename, map_location=device))
    model = model.to(device)
    model = model.eval()

    return model

#%% Load test data
### ----------------------------------------------------------------------------

def load_data_directory(basename):

    id_prop_dat = os.path.join(basename, "id_prop.csv")

    with open(id_prop_dat, "r") as f:
        reader = csv.reader(f)
        data = [row for row in reader]

    dataset = []
    for i in data:
        filename = os.path.join(basename, i[0])

        info           = {}
        info["atoms"]  = Atoms.from_poscar(filename)
        info["jid"]    = i[0]
        info['target'] = float(i[1])

        dataset.append(info)
    
    return dataset

#%% Get predictions
### ----------------------------------------------------------------------------

def get_predictions(model, dataset):
    model.eval()
    targets     = []
    predictions = []
    with torch.no_grad():
        for datum in tqdm(dataset):
            g, lg = Graph.atom_dgl_multigraph(
                    datum['atoms'],
                    cutoff        = config.cutoff,
                    atom_features = config.atom_features,
                    max_neighbors = config.max_neighbors)
            y_hat = model([g.to(device), lg.to(device)]).item()

            predictions.append(y_hat)
            targets    .append(datum['target'])

    return targets, predictions

#%% Compute test result
### ----------------------------------------------------------------------------

if __name__ == "__main__":

    model   = load_model('data_train_result')
    dataset = load_data_directory('data_test_wbm')

    targets, predictions = get_predictions(model, dataset)

    mae = torch.nn.L1Loss()(torch.tensor(predictions), torch.tensor(targets)).item()

    print(f'Test MAE: {mae}')

    # %%
    dumpfn({'material_id': list(df_wbm['material_id']),
            'predictions': predictions,
            'targets'    : targets,
            'mae'        : mae },
            'test_alignn_result.json')
