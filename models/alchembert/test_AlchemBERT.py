import lightning as l
import pandas as pd
import torch
from fire import Fire
from torch.utils.data import DataLoader
from train_AlchemBERT import MatBert, MyDataset, task


torch.manual_seed(42)


def get_test_data() -> pd.Series:
    target_col = "e_form_per_atom_mp2020_corrected"
    df_wbm = pd.read_csv("2022-10-19-wbm-summary.csv")
    return pd.Series(df_wbm[target_col])


def get_train_data() -> pd.Series:
    target_col = "formation_energy_per_atom"
    id_col = "material_id"
    df_eng = pd.read_csv("2023-01-10-mp-energies.csv").set_index(id_col)
    return pd.Series(df_eng[target_col], index=df_eng.index)


bert_path = "bert-base-cased"
predictions_path = "2024-12-25-alchembert-wbm-IS2RE.csv.gz"
test_pad_cased_path = f"test_{task}_pad_cased_inputs.json"


# %%

def main(
    best_epoch: int,
    val_mae: float
) -> None:

    best_model = f"epoch={best_epoch}_val_MAE={val_mae}_best_model.ckpt"
    best_model_path = f"checkpoints/model_epoch5000_{task}/{best_model}"
    test_inputs = pd.read_json(test_pad_cased_path)
    test_outputs = get_test_data()

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
