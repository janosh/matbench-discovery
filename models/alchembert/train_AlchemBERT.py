import os
import sys
import warnings

import lightning
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as f
from lightning.pytorch.callbacks import EarlyStopping, ModelCheckpoint
from torch.utils.data import DataLoader, Dataset, random_split
from transformers import BertConfig, BertModel

seed = 42
torch.manual_seed(seed)
if torch.cuda.is_available():
    torch.cuda.manual_seed_all(seed)

max_length = 512
train_batch_size = 32
val_batch_size = 32
epoch = 5000
patience = 200
log_every_n_steps = 50
save_top_k = 1
l_r = 1e-5

task = "nl"
bert_path = "bert-base-cased"

train_pad_cased_path = "train_nl_pad_cased_inputs.json"
test_pad_cased_path = "test_nl_pad_cased_inputs.json"


def get_test_data() -> pd.Series:
    target_col = "e_form_per_atom_mp2020_corrected"
    df_wbm = pd.read_csv("2022-10-19-wbm-summary.csv")
    return pd.Series(df_wbm[target_col])


def get_train_data() -> pd.Series:
    target_col = "formation_energy_per_atom"
    id_col = "material_id"
    df_eng = pd.read_csv("2023-01-10-mp-energies.csv").set_index(id_col)
    return pd.Series(df_eng[target_col], index=df_eng.index)


# %%
class MyDataset(Dataset):
    def __init__(
        self,
        input_ids: torch.Tensor,
        attention_mask: torch.Tensor,
        labels: torch.Tensor,
    ) -> None:
        self.input_ids = input_ids
        self.attention_mask = attention_mask
        self.labels = labels

    def __len__(self) -> int:
        return len(self.labels)

    def __getitem__(
        self, index: int
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        return self.input_ids[index], self.attention_mask[index], self.labels[index]


# %%
class MatBert(lightning.LightningModule):
    def __init__(self, b_path: str) -> None:
        super().__init__()
        self.bert = BertModel.from_pretrained(b_path, output_hidden_states=True)
        self.config = BertConfig.from_pretrained(bert_path)
        self.linear = nn.Linear(self.config.hidden_size, 1)

    def forward(
        self, input_ids: torch.Tensor, attention_mask: torch.Tensor
    ) -> torch.Tensor:
        outputs = self.bert(input_ids=input_ids, attention_mask=attention_mask)
        cls_representation = outputs.last_hidden_state[:, 0, :]
        return self.linear(cls_representation).squeeze(-1)

    def training_step(
        self, batch: tuple[torch.Tensor, torch.Tensor, torch.Tensor]
    ) -> torch.Tensor:
        input_ids, attention_mask, y = batch
        input_ids.cuda()
        attention_mask.cuda()
        y.cuda()
        y_hat = self(input_ids, attention_mask)
        loss = f.mse_loss(y_hat.float(), y.float())
        self.log("train_mse_loss", loss, on_epoch=True, sync_dist=True)
        return loss

    def validation_step(
        self, batch: tuple[torch.Tensor, torch.Tensor, torch.Tensor]
    ) -> dict[str, torch.Tensor]:
        input_ids, attention_mask, y = batch
        input_ids.cuda()
        attention_mask.cuda()
        y.cuda()
        y_hat = self(input_ids, attention_mask)
        loss = nn.functional.mse_loss(y_hat.float(), y.float())
        mae = torch.mean(torch.absolute(y_hat - y))
        self.log("val_MAE", mae, on_epoch=True, sync_dist=True)
        return {"val_loss": loss, "val_MAE": mae}

    def predict_step(
        self, batch: tuple[torch.Tensor, torch.Tensor, torch.Tensor]
    ) -> torch.Tensor:
        input_ids, attention_mask, y = batch
        return self(input_ids, attention_mask)

    def configure_optimizers(self) -> torch.optim.Optimizer:
        return torch.optim.Adam(self.parameters(), lr=l_r)


# %% data
def main() -> None:
    if os.path.exists(train_pad_cased_path):
        print(f"file {train_pad_cased_path} exists")
        train_inputs = pd.read_json(train_pad_cased_path)
        train_outputs = get_train_data()

        input_ids = torch.tensor(train_inputs["input_ids"])
        attention_mask = torch.tensor(train_inputs["attention_mask"])
        train_outputs = torch.tensor(train_outputs.values)
    else:
        warnings.warn("file doesn't exist", UserWarning, stacklevel=2)
        sys.exit()

    dataset = MyDataset(input_ids, attention_mask, train_outputs)
    train_set, val_set = random_split(dataset, [0.9, 0.1])

    train_loader = DataLoader(
        train_set, batch_size=train_batch_size, shuffle=True, num_workers=2
    )
    val_loader = DataLoader(
        val_set, batch_size=val_batch_size, shuffle=False, num_workers=2
    )

    # %% train

    model = MatBert(bert_path)
    model.cuda()
    early_stopping = EarlyStopping(
        monitor="val_MAE", patience=patience, verbose=True, mode="min"
    )
    check_point = ModelCheckpoint(
        monitor="val_MAE",
        save_top_k=save_top_k,
        dirpath=f"checkpoints/model_epoch{epoch}_{task}",
        filename="{epoch}_{val_MAE:.4f}_best_model",
        mode="min",
    )
    trainer = lightning.Trainer(
        max_epochs=epoch,
        accelerator="gpu",
        callbacks=[check_point, early_stopping],
        log_every_n_steps=log_every_n_steps,
        devices=-1,
        strategy="ddp_find_unused_parameters_true",
    )
    model.train()
    trainer.fit(model=model, train_dataloaders=train_loader, val_dataloaders=val_loader)


# %%
if __name__ == "__main__":
    main()
