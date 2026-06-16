# Justfile for Matbench Discovery
# https://github.com/casey/just
#
# Recipes are thin entry points; the ingestion logic lives in
# scripts/ingest_model.py (tested in tests/test_ingest_model.py).

set dotenv-load := false

# List available commands
default:
    @just --list

# Prepare a model submission: PR checklist + evals + per-model figures + refresh of
# all multi-model site figure payloads
prepare-model-submission model_name overwrite="false":
    uv run python scripts/ingest_model.py {{ model_name }} {{ if overwrite == "true" { "--overwrite" } else { "" } }}

# Fully ingest a model submission: prepare-model-submission plus archiving the
# author's prediction files to the project's figshare articles (rewrites the YAML's
# *_url keys for longevity) and publishing parity assets to the GitHub release the
# site build downloads from. Requires FIGSHARE_TOKEN and gh auth.
ingest-model model_name overwrite="false":
    uv run python scripts/ingest_model.py {{ model_name }} --archive {{ if overwrite == "true" { "--overwrite" } else { "" } }}

# Refresh the multi-model site figure payloads (.jsonl) so pages like /models/tmi and
# /tasks/geo-opt include every model, then run the payload shape tests. Run after
# adding/updating a model, then commit the changed payloads.
# `just update-site-figs <model>` merges only that model's entries into the committed
# payloads (no other model's prediction files needed locally); without a model, all
# payloads are regenerated from the full active roster.
update-site-figs model_name="":
    uv run python scripts/ingest_model.py {{ model_name }} --payloads-only
