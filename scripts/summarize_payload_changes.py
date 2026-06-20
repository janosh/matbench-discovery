"""Summarize per-payload model-roster changes for the update-site-figs PR body.

Diffs each regenerated JSONL figure payload (site/src/figs/*.jsonl + the route-local
per-element-errors.jsonl) against its committed (HEAD) version and prints a markdown
bullet list of which payloads changed and which models were added/removed. The
update-site-figs workflow captures stdout into the PR body. Run from the repo root.
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
    entries = (json.loads(line) for line in text.splitlines() if line.strip())
    # the lone {"_base": ...} line isn't a model; drop metadata-only columns too
    return {
        str(entry.get("key") or entry.get("label"))
        for entry in entries
        if set(entry) != {"_base"}
    } - METADATA_KEYS


def format_change(name: str, was: set[str], now: set[str]) -> str | None:
    """Markdown bullet for one payload's roster change, or None if unchanged. Lists
    both additions (+N) and removals (-N) so the count delta is always explained.
    """
    if was == now:
        return None
    # was != now guarantees at least one non-empty delta, so changes is never empty
    deltas = {"+": sorted(now - was), "-": sorted(was - now)}
    changes = [
        f"{sign}{len(ids)}: {', '.join(ids)}" for sign, ids in deltas.items() if ids
    ]
    return f"- `{name}`: {len(was)} → {len(now)} ({'; '.join(changes)})"


def summarize() -> str:
    """Markdown summary of roster changes across all committed JSONL payloads."""
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
        name = os.path.basename(path).removesuffix(".jsonl")
        if (row := format_change(name, was, now)) is not None:
            rows.append(row)

    body = "\n".join(rows) or "_No roster changes (data-only refresh)._"
    return f"### Model-count changes\n\n{body}"


if __name__ == "__main__":
    print(summarize())
