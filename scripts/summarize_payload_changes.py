"""Report structural, roster, identity, and numerical JSONL payload changes.

This script is stdlib-only because trusted CI runs it from the pull-request worktree.
"""

from __future__ import annotations

import hashlib
import json
import os
import subprocess
from dataclasses import dataclass, field
from glob import glob
from typing import TypeGuard

METADATA_FIELDS = frozenset(
    {"schema_version", "identity", "audit", "model_key", "label", "input_artifacts"}
)


@dataclass
class Delta:
    """Recursive comparison statistics for one shared or model record."""

    structural_paths: list[str] = field(default_factory=list)
    numeric_leaf_count: int = 0
    max_absolute_error: float = 0
    max_relative_error: float = 0


def canonical_json(value: object) -> str:
    """Return deterministic compact JSON for equality checks."""
    return json.dumps(value, allow_nan=False, separators=(",", ":"), sort_keys=True)


def _git_output(*args: str) -> str:
    """Return text output from a checked Git command."""
    git_command = ["git", *args]
    return subprocess.check_output(git_command, text=True)


def _is_number(value: object) -> TypeGuard[int | float]:
    """Return whether value is a JSON number, excluding booleans."""
    return isinstance(value, (int, float)) and not isinstance(value, bool)


def _is_object(value: object) -> TypeGuard[dict[object, object]]:
    """Return whether value is a JSON object."""
    return isinstance(value, dict)


def compare_values(
    old: object,
    new: object,
    *,
    path: str = "$",
    delta: Delta | None = None,
    exclude_numeric: bool = False,
) -> Delta:
    """Recursively compare JSON values and collect structural and numeric deltas."""
    delta = delta or Delta()
    if not exclude_numeric and _is_number(old) and _is_number(new):
        absolute_error = abs(float(old) - float(new))
        denominator = max(abs(float(old)), abs(float(new)))
        relative_error = absolute_error / denominator if denominator else 0
        delta.numeric_leaf_count += 1
        delta.max_absolute_error = max(delta.max_absolute_error, absolute_error)
        delta.max_relative_error = max(delta.max_relative_error, relative_error)
        return delta
    if type(old) is not type(new) and not (_is_number(old) and _is_number(new)):
        delta.structural_paths.append(
            f"{path} (type {type(old).__name__} -> {type(new).__name__})"
        )
        return delta
    if _is_object(old) and _is_object(new):
        for key in sorted(old.keys() | new.keys(), key=str):
            child_path = f"{path}.{key}"
            if key not in old:
                delta.structural_paths.append(f"{child_path} (added)")
            elif key not in new:
                delta.structural_paths.append(f"{child_path} (removed)")
            else:
                compare_values(
                    old[key],
                    new[key],
                    path=child_path,
                    delta=delta,
                    exclude_numeric=exclude_numeric or key in METADATA_FIELDS,
                )
        return delta
    if isinstance(old, list) and isinstance(new, list):
        if len(old) != len(new):
            delta.structural_paths.append(f"{path} (length {len(old)} -> {len(new)})")
        for index, (old_item, new_item) in enumerate(zip(old, new, strict=False)):
            compare_values(
                old_item,
                new_item,
                path=f"{path}[{index}]",
                delta=delta,
                exclude_numeric=exclude_numeric,
            )
        return delta
    if old != new:
        delta.structural_paths.append(path)
    return delta


def _parse_payload(
    text: str, label_keys: dict[str, str] | None = None
) -> tuple[dict[str, object], dict[str, dict[str, object]]]:
    """Parse one schema-v2 JSONL payload indexed by immutable model key."""
    base: dict[str, object] = {}
    models: dict[str, dict[str, object]] = {}
    for record in map(json.loads, filter(str.strip, text.splitlines())):
        if set(record) == {"_base"}:
            base = record["_base"]
            continue
        # The reporter reads historical Git payloads while reviewing the schema-v2
        # migration. Runtime payload readers remain schema-v2-only.
        model_key = record.get("model_key", record.get("key"))
        if not isinstance(model_key, str) and isinstance(record.get("label"), str):
            model_key = (label_keys or {}).get(record["label"], record["label"])
        if not isinstance(model_key, str):
            raise TypeError(f"Payload record has no model_key: {record!r}")
        models[model_key] = record
    return base, models


def _format_paths(paths: list[str]) -> str:
    """Format a bounded structural-path list for a Markdown table cell."""
    visible = paths[:6]
    suffix = f"; +{len(paths) - len(visible)} more" if len(paths) > len(visible) else ""
    return "; ".join(f"`{path}`" for path in visible) + suffix if paths else "—"


def _identity_hash(identity: object) -> str:
    """Hash complete computation identity for compact reporting."""
    if identity is None:
        return "none"
    return hashlib.sha256(canonical_json(identity).encode()).hexdigest()


def _changed_identity_components(old: object, new: object) -> str:
    """List changed top-level computation-identity components."""
    if not isinstance(old, dict) or not isinstance(new, dict):
        return "record"
    changed = [
        str(key)
        for key in sorted(old.keys() | new.keys(), key=str)
        if old.get(key) != new.get(key)
    ]
    return ", ".join(changed) or "none"


def _changed_row(
    payload_name: str,
    record_key: str,
    old: object | None,
    new: object | None,
) -> str | None:
    """Return one Markdown delta row, or None for byte-equivalent JSON values."""
    if (
        old is not None
        and new is not None
        and canonical_json(old) == canonical_json(new)
    ):
        return None
    if old is None:
        delta = Delta(structural_paths=["$ (record added)"])
    elif new is None:
        delta = Delta(structural_paths=["$ (record removed)"])
    else:
        delta = compare_values(old, new)
    return (
        f"| `{payload_name}` | `{record_key}` | "
        f"{_format_paths(delta.structural_paths)} | {delta.numeric_leaf_count} | "
        f"{delta.max_absolute_error:.6g} | {delta.max_relative_error:.6g} |"
    )


def summarize() -> str:
    """Return a complete Markdown review report for all changed JSONL payloads."""
    working_paths = {
        *glob("site/src/figs/*.jsonl"),
        *glob("site/src/routes/models/per-element-each-errors.jsonl"),
    }
    head_listing = _git_output(
        "ls-tree",
        "-r",
        "--name-only",
        "HEAD",
        "--",
        "site/src/figs",
        "site/src/routes/models/per-element-each-errors.jsonl",
    )
    head_paths = {path for path in head_listing.splitlines() if path.endswith(".jsonl")}
    delta_rows: list[str] = []
    summaries: list[str] = []
    for path in sorted(working_paths | head_paths):
        if os.path.isfile(path):
            with open(path, encoding="utf-8") as file:
                new_text = file.read()
        else:
            new_text = ""
        old_text = (
            _git_output("show", "HEAD:" + path.replace("\\", "/"))
            if path in head_paths
            else ""
        )
        new_base, new_models = _parse_payload(new_text)
        label_keys = {
            label: model_key
            for model_key, model in new_models.items()
            if isinstance(label := model.get("label"), str)
        }
        old_base, old_models = _parse_payload(old_text, label_keys)
        payload_name = os.path.basename(path).removesuffix(".jsonl")
        if row := _changed_row(payload_name, "<shared>", old_base, new_base):
            delta_rows.append(row)
        delta_rows.extend(
            row
            for model_key in sorted(set(old_models) | set(new_models))
            if (
                row := _changed_row(
                    payload_name,
                    model_key,
                    old_models.get(model_key),
                    new_models.get(model_key),
                )
            )
        )

        old_keys, new_keys = set(old_models), set(new_models)
        if old_keys != new_keys:
            added = ", ".join(sorted(new_keys - old_keys)) or "—"
            removed = ", ".join(sorted(old_keys - new_keys)) or "—"
            summaries.append(
                f"- `{payload_name}` roster {len(old_keys)} -> {len(new_keys)}; "
                f"added: {added}; removed: {removed}"
            )
        old_identity = old_base.get("identity")
        new_identity = new_base.get("identity")
        if canonical_json(old_identity) != canonical_json(new_identity):
            summaries.append(
                f"- `{payload_name}` identity `{_identity_hash(old_identity)}` "
                f"-> `{_identity_hash(new_identity)}`; changed: "
                f"{_changed_identity_components(old_identity, new_identity)}"
            )

    summary = "\n".join(summaries) or "_No roster or identity changes._"
    rows = "\n".join(delta_rows) or "_No payload records changed._"
    return (
        f"### Roster and identity changes\n\n{summary}\n\n"
        "### Payload record deltas\n\n"
        "| Payload | Record | Structural paths | Numeric leaves | Max abs error | "
        "Max symmetric relative error |\n"
        "| --- | --- | --- | ---: | ---: | ---: |\n"
        f"{rows}"
    )


if __name__ == "__main__":
    print(summarize())
