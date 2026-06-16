"""Summarize per-payload model-roster changes for the update-site-figs PR body.

Diffs each regenerated JSONL figure payload (site/src/figs/*.jsonl + the route-local
per-element-errors.jsonl) against its committed (HEAD) version and prints a markdown
bullet list of which payloads grew and which models were added. The update-site-figs
workflow captures stdout into the PR body. Run from the repo root.
"""

import json
import os
import subprocess
from glob import glob

# per-element columns that are metadata, not models (excluded from roster counts)
METADATA_KEYS = {"MP Occurrences", "Test set standard deviation"}


def roster(text: str) -> set[str]:
    """Set of model ids in one .jsonl payload (one JSON object per line; the lone
    ``{"_base": ...}`` line and metadata-only columns are not models).
    """
    ids: set[str] = set()
    for line in text.splitlines():
        if not (line := line.strip()):
            continue
        entry = json.loads(line)
        if set(entry) != {"_base"}:  # the lone _base line isn't a model
            ids.add(str(entry.get("key") or entry.get("label")))
    return ids - METADATA_KEYS


def summarize() -> str:
    """Markdown summary of roster growth across all committed JSONL payloads."""
    paths = [
        *glob("site/src/figs/*.jsonl"),
        "site/src/routes/models/per-element-each-errors.jsonl",
    ]
    rows: list[str] = []
    for path in sorted(paths):
        with open(path) as file:
            now = roster(file.read())
        # command via a variable so ruff can't flag the partial git path (S607)
        git_show = ["git", "show", f"HEAD:{path}"]
        head = subprocess.run(git_show, capture_output=True, text=True, check=False)
        was = roster(head.stdout) if head.returncode == 0 else set()
        if was == now:
            continue
        added = sorted(now - was)
        name = os.path.basename(path).removesuffix(".jsonl")
        tail = f" (+{len(added)}: {', '.join(added)})" if added else ""
        rows.append(f"- `{name}`: {len(was)} → {len(now)}{tail}")

    body = "\n".join(rows) or "_No roster changes (data-only refresh)._"
    return f"### Model-count changes\n\n{body}"


if __name__ == "__main__":
    print(summarize())
