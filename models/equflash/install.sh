#!/usr/bin/bash
uv venv --python 3.12
source .venv/bin/activate
uv pip install -r requirements.txt --index-strategy unsafe-best-match
uv pip install fairchem-core==1.10.0 scipy==1.16.1 --no-deps