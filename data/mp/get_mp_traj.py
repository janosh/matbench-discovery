"""Download all MP ionic steps using direct read-access to the mp_core DB.

Gzipped JSON is ~15GB.
On a good connection, takes about 15 min per batch * 140 batches = 35 h to download
all 1.6M task docs.
"""

# %%
import os
import subprocess
from glob import glob

import pandas as pd
from emmet.core.tasks import TaskDoc
from pymatviz.enums import Key
from pymongo import MongoClient
from pymongo.database import Database
from tqdm import tqdm, trange

from matbench_discovery import ROOT, today

__author__ = "Janosh Riebesell"
__date__ = "2023-03-15"

module_dir = os.path.dirname(__file__)


# %% access mp_core database directly through pymongo instead of API for speed
host = "knowhere.lbl.gov"
db_name = "mp_core"

with open(f"{ROOT}/site/.env") as file:
    text = file.read()
    user = text.split("user=")[1].split("\n")[0]
    password = text.split("password=")[1].split("\n")[0]

uri = f"mongodb://{user}:{password}@{host}/?authSource={db_name}"
db: Database[TaskDoc] = MongoClient(uri)[db_name]


# %%
ids_path = f"{module_dir}/2023-03-15-mp-task-ids.csv.bz2"
fields = (
    Key.task_id,
    "formula_pretty",
    "run_type",
    "nsites",
    Key.task_type,
    "tags",
    "completed_at",
)

if os.path.isfile(ids_path):
    print(f"Found existing list of task IDs to query at\n{ids_path=}")
    df_tasks = pd.read_csv(ids_path, low_memory=False).set_index(Key.task_id)
else:
    print(f"Querying all task docs from {db_name}\n{fields=}.\nThis takes a while...")
    task_docs = sorted(
        db["tasks"].find({}, fields),
        key=lambda doc: int(doc[Key.task_id].split("-")[1]),
    )

    print(f"{today}: {len(task_docs)=:,}")

    df_tasks = pd.DataFrame(task_docs).drop(columns=["_id"]).set_index(Key.task_id)
    df_tasks.task_type.value_counts(dropna=False).plot.pie()

    df_tasks.to_csv(f"{module_dir}/{today}-mp-task-ids.csv.bz2")


# %% inspect schema of a single task doc
doc = db.tasks.find_one({Key.task_id: "mp-288"})
# the most relevant task data is found in the 1st calc's ionic steps which are
# the relaxation trajectory frames with the highest rate of change
# see docs[0]["calcs_reversed"][-1]["output"]["ionic_steps"]


# %%
batch_size = 10_000
task_ids = df_tasks.index.tolist()

os.makedirs(f"{module_dir}/mp-tasks", exist_ok=True)
# Iterate over task_ids in batches
desc = "Fetching MP task docs..."
pbar = trange(0, len(task_ids), batch_size, desc=desc, unit_scale=batch_size)
for start_idx in pbar:
    # Define start and end indices for batch
    end_idx = min(start_idx + batch_size, len(task_ids))
    start_id = task_ids[start_idx]
    end_id = task_ids[end_idx - 1]
    batch_ids = task_ids[start_idx:end_idx]
    pbar.set_postfix_str(f"{start_id} to {end_id}")

    out_path = f"{module_dir}/mp-tasks/{start_id}__{end_id}.json.gz"

    # Check if output file for batch already exists
    if os.path.isfile(out_path):
        continue

    # query batch of task docs
    batch_docs = list(
        db["tasks"].find(
            {Key.task_id: {"$in": batch_ids}},
            [*fields, "calcs_reversed.output.ionic_steps"],
        )
    )

    # Convert documents to DataFrame and save to file
    df_batch = pd.DataFrame(batch_docs).set_index(Key.task_id).drop(columns=["_id"])
    # handler=str needed since MongoDB ObjectId is not JSON serializable
    df_batch.reset_index().to_json(out_path, default_handler=str)
    # don't store df_batch to save memory


# %% inspect saved task docs for expected data
df_batch = pd.read_json(
    f"{module_dir}/mp-tasks/mp-531529__mp-568116.json.gz"
).set_index(Key.task_id)

print(f"{len(df_batch)=}")
df_batch.head()


# %% use gzip CLI to check all files for archive corruption
for path in tqdm(glob(f"{module_dir}/mp-tasks/*.json.gz")):
    try:
        subprocess.run(["gzip", "--test", path], check=True)
    except subprocess.CalledProcessError as exc:
        print(f"{path} raised {exc.stderr}")
        # ask user to delete corrupted file
        if input("Delete corrupted file? [y/N] ").lower() == "y":
            os.remove(path)
