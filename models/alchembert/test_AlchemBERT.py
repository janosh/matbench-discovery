import lightning as l
import pandas as pd
import torch
from fire import Fire
from pymatgen.core import Structure
from torch.utils.data import DataLoader
from train_AlchemBERT import MatBert, MyDataset, task

torch.manual_seed(42)


def get_test_data(only_y=False):
    id_col = "material_id"
    input_col = "initial_structure"
    target_col = "e_form_per_atom_mp2020_corrected"
    data_path = "2022-10-19-wbm-init-structs.json"
    df_wbm = pd.read_csv("2022-10-19-wbm-summary.csv")
    if only_y is False:
        df_in = pd.read_json(data_path).set_index(id_col)
        X = pd.Series(
            [Structure.from_dict(x) for x in df_in[input_col]], index=df_in.index
        )
        y = pd.Series(df_wbm[target_col])
        return X[y.index], y
    y = pd.Series(df_wbm[target_col])
    return y


def get_train_data(only_y=False):
    target_col = "formation_energy_per_atom"
    input_col = "structure"
    id_col = "material_id"
    if only_y is False:
        df_cse = pd.read_json(
            "2023-02-07-mp-computed-structure-entries.json"
        ).set_index(id_col)
        df_eng = pd.read_csv("2023-01-10-mp-energies.csv").set_index(id_col)
        X = pd.Series(
            [Structure.from_dict(cse[input_col]) for cse in df_cse.entry],
            index=df_cse.index,
        )
        y = pd.Series(df_eng[target_col], index=df_eng.index)
        return X[y.index], y
    df_eng = pd.read_csv("2023-01-10-mp-energies.csv").set_index(id_col)
    y = pd.Series(df_eng[target_col], index=df_eng.index)
    return y


bert_path = "bert-base-cased"
predictions_path = "2024-12-25-alchembert-wbm-IS2RE.csv.gz"
test_pad_cased_path = f"test_{task}_pad_cased_inputs.json"


# %%
def main(best_epoch, val_mae):
    best_model = f"epoch={best_epoch}_val_MAE={val_mae}_best_model.ckpt"
    best_model_path = f"checkpoints/model_epoch5000_{task}/{best_model}"
    test_inputs = pd.read_json(test_pad_cased_path)
    test_outputs = get_test_data(only_y=True)

    best_model = torch.load(best_model_path, weights_only=True)
    model = MatBert(bert_path)
    model.load_state_dict(best_model["state_dict"])
    model.eval()

    # %% test
    trainer = l.Trainer(accelerator="gpu", devices=[0])

    test_input_ids = torch.tensor(test_inputs["input_ids"])
    test_attention_mask = torch.tensor(test_inputs["attention_mask"])
    test_outputs = torch.tensor(test_outputs.values)
    test_dataset = MyDataset(test_input_ids, test_attention_mask, test_outputs)
    test_loader = DataLoader(test_dataset, batch_size=1, shuffle=False)
    predictions = trainer.predict(model, test_loader)

    predictions = [tensor.cpu().item() for tensor in predictions]
    results = {"e_form_per_atom_alchembert": predictions}
    results = pd.DataFrame(results)
    print(results)
    results.to_csv(predictions_path, index=False, compression="gzip")
    print(predictions_path)


# %%
if __name__ == "__main__":
    Fire(main)
