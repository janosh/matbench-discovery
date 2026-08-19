"""Summarize shared and keyed JSONL payload changes for CI commit messages."""

from __future__ import annotations

import json
import os
import subprocess
from typing import TypeGuard


def _is_object(value: object) -> TypeGuard[dict[str, object]]:
    """Return whether a value is a JSON object."""
    return isinstance(value, dict)


def _is_number(value: object) -> TypeGuard[int | float]:
    """Return whether a value is a JSON number rather than a Boolean."""
    return isinstance(value, (int, float)) and not isinstance(value, bool)


def canonical_json(value: object) -> str:
    """Return compact deterministic JSON for exact value comparisons."""
    return json.dumps(value, allow_nan=False, separators=(",", ":"), sort_keys=True)


def _git_output(*args: str) -> str:
    """Return text output from a checked Git command."""
    command = ["git", *args]
    return subprocess.check_output(command, text=True)


def _parse_payload(
    text: str,
) -> tuple[dict[str, object], dict[str, dict[str, object]]]:
    """Parse one schema-v2 JSONL payload indexed by immutable model key."""
    base: dict[str, object] = {}
    models: dict[str, dict[str, object]] = {}
    for record in map(json.loads, filter(str.strip, text.splitlines())):
        if set(record) == {"_base"}:
            base = record["_base"]
        elif not isinstance(model_key := record.get("model_key"), str):
            raise TypeError(f"Payload record has no model_key: {record!r}")
        elif model_key in models:
            raise ValueError(f"Duplicate payload model_key: {model_key!r}")
        else:
            models[model_key] = record
    return base, models


def _numeric_delta(old: object, new: object) -> str:
    """Format maximum absolute and symmetric-relative numeric differences."""
    stack = [(old, new)]
    max_absolute = max_relative = 0.0
    while stack:
        old_value, new_value = stack.pop()
        if _is_number(old_value) and _is_number(new_value):
            absolute = abs(float(old_value) - float(new_value))
            denominator = max(abs(float(old_value)), abs(float(new_value)))
            max_absolute = max(max_absolute, absolute)
            max_relative = max(
                max_relative, absolute / denominator if denominator else 0
            )
        elif _is_object(old_value) and _is_object(new_value):
            stack.extend(
                (old_value[key], new_value[key])
                for key in old_value.keys() & new_value.keys()
                if key != "input_artifacts"
            )
        elif isinstance(old_value, list) and isinstance(new_value, list):
            stack.extend(zip(old_value, new_value, strict=False))
    return (
        f" (max |Δ|={max_absolute:.6g}, max rel={max_relative:.6g})"
        if max_absolute
        else ""
    )


def _change_summary(
    name: str,
    old_base: dict[str, object],
    new_base: dict[str, object],
    old_models: dict[str, dict[str, object]],
    new_models: dict[str, dict[str, object]],
) -> str | None:
    """Return one compact Markdown summary, or None when nothing changed."""
    old_keys, new_keys = set(old_models), set(new_models)
    updated = [
        model_key + _numeric_delta(old_models[model_key], new_models[model_key])
        for model_key in sorted(old_keys & new_keys)
        if canonical_json(old_models[model_key])
        != canonical_json(new_models[model_key])
    ]
    changes: list[str] = []
    if (
        old_base
        and new_base
        and old_base.get("schema_version") != new_base.get("schema_version")
    ):
        changes.append("schema")
    old_identity, new_identity = old_base.get("identity"), new_base.get("identity")
    if canonical_json(old_identity) != canonical_json(new_identity):
        identity_change = "identity"
        if _is_object(old_identity) and _is_object(new_identity):
            components = ", ".join(
                key
                for key in sorted(old_identity.keys() | new_identity.keys())
                if canonical_json(old_identity.get(key))
                != canonical_json(new_identity.get(key))
            )
            identity_change += f" ({components})"
        changes.append(identity_change)
    old_derived, new_derived = old_base.get("derived"), new_base.get("derived")
    if canonical_json(old_derived) != canonical_json(new_derived):
        changes.append("shared data" + _numeric_delta(old_derived, new_derived))
    if canonical_json(old_base.get("audit")) != canonical_json(new_base.get("audit")):
        changes.append("audit")
    if added := new_keys - old_keys:
        changes.append(f"added: {', '.join(sorted(added))}")
    if removed := old_keys - new_keys:
        changes.append(f"removed: {', '.join(sorted(removed))}")
    if updated:
        changes.append(f"updated: {', '.join(updated)}")
    return f"- `{name}`: {'; '.join(changes)}" if changes else None


def summarize() -> str:
    """Return a Markdown summary of all changed multi-model payload records."""
    changes = {
        entry[3:]: entry[:2]
        for entry in _git_output(
            "status",
            "--short",
            "--no-renames",
            "-z",
            "--untracked-files=all",
            "--",
            "site/src/figs",
            "site/src/routes/models/per-element-each-errors.jsonl",
        ).split("\0")
        if entry and entry[3:].endswith(".jsonl")
    }
    rows: list[str] = []
    for path, status in sorted(changes.items()):
        new_text = ""
        if os.path.isfile(path):
            with open(path, encoding="utf-8") as file:
                new_text = file.read()
        old_text = (
            _git_output("show", f"HEAD:{path.replace(os.sep, '/')}")
            if status != "??" and "A" not in status
            else ""
        )
        old_base, old_models = _parse_payload(old_text)
        new_base, new_models = _parse_payload(new_text)
        name = os.path.basename(path).removesuffix(".jsonl")
        if row := _change_summary(name, old_base, new_base, old_models, new_models):
            rows.append(row)

    body = "\n".join(rows) or "_No payload records changed._"
    return f"### Payload changes\n\n{body}"


if __name__ == "__main__":
    print(summarize())
