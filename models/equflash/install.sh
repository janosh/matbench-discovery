#!/usr/bin/bash
uv venv --python 3.12
source .venv/bin/activate
uv pip install torch==2.9.1+cu126 --extra-index-url https://download.pytorch.org/whl/cu126
uv pip install -r requirements.txt

# Intentional workaround: only the exercised fairchem code path is used.
uv pip install fairchem-core==1.10.0 scipy==1.16.1 --no-deps
