"""Overlay changed model YAMLs from an untrusted PR without executing its code.

Used by .github/workflows/update-site-figs.yml ingest mode: the workflow checks out
main (trusted code) plus the PR head into a side directory, then runs this script to
apply canonical ``models/<arch>/<model>.yml`` changes reported by the GitHub API.

This makes secrets exposure during ingestion safe by construction: all *executed*
code comes from main; the PR contributes data only.
"""

import argparse
import os
import re
import shutil
from collections.abc import Sequence

MODEL_YAML_PATTERN = re.compile(
    r"^models/[A-Za-z0-9][A-Za-z0-9_.+@-]*/"
    r"[A-Za-z0-9][A-Za-z0-9_.+@-]*\.yml$"
)


def _canonical_model_yaml_path(path: str) -> str | None:
    """Return a POSIX models/<family>/<key>.yml path, or None if unsafe."""
    normalized = path.replace("\\", "/")
    if not MODEL_YAML_PATTERN.fullmatch(normalized):
        return None
    return normalized


def _canonicalize_paths(paths: Sequence[str]) -> list[str] | str:
    """Validate and POSIX-normalize model YAML paths, or return an error.

    Canonicalizes before deduplicating so separator-equivalent inputs collapse
    to one entry while preserving first-seen order.
    """
    canonical_paths: list[str] = []
    seen: set[str] = set()
    for relative_path in paths:
        canonical = _canonical_model_yaml_path(relative_path)
        if canonical is None:
            return f"::error::Unsafe model YAML path: {relative_path!r}"
        if canonical in seen:
            continue
        seen.add(canonical)
        canonical_paths.append(canonical)
    return canonical_paths


def _validate_overlay(
    pr_tree: str,
    yaml_paths: Sequence[str],
    removed_paths: Sequence[str],
) -> tuple[list[str], list[str]] | str:
    """Return canonical path lists, or an error message string."""
    canonical_yaml_paths = _canonicalize_paths(yaml_paths)
    if isinstance(canonical_yaml_paths, str):
        return canonical_yaml_paths
    canonical_removed_paths = _canonicalize_paths(removed_paths)
    if isinstance(canonical_removed_paths, str):
        return canonical_removed_paths

    checks = (
        (
            not canonical_yaml_paths and not canonical_removed_paths,
            "::error::No changed model YAML paths supplied",
        ),
        (
            os.path.islink(f"{pr_tree}/models"),
            f"::error::Unsafe submitted models directory: {pr_tree}/models",
        ),
    )
    for failed, message in checks:
        if failed:
            return message

    for relative_path in canonical_removed_paths:
        if not (os.path.isfile(relative_path) or os.path.islink(relative_path)):
            return f"::error::Missing trusted YAML to remove: {relative_path}"

    for relative_path in canonical_yaml_paths:
        source_path = f"{pr_tree}/{relative_path}"
        if (
            os.path.islink(os.path.dirname(source_path))
            or not os.path.isfile(source_path)
            or os.path.islink(source_path)
        ):
            return f"::error::Missing or unsafe submitted YAML: {source_path}"

    return canonical_yaml_paths, canonical_removed_paths


def main(
    pr_tree: str,
    yaml_paths: Sequence[str],
    removed_paths: Sequence[str] = (),
) -> int:
    """Apply explicit, canonical model YAML changes to the trusted checkout."""
    validated = _validate_overlay(pr_tree, yaml_paths, removed_paths)
    if isinstance(validated, str):
        print(validated)
        return 1

    canonical_yaml_paths, canonical_removed_paths = validated
    for relative_path in canonical_removed_paths:
        os.remove(relative_path)
        print(f"Removed {relative_path}")

    for relative_path in canonical_yaml_paths:
        source_path = f"{pr_tree}/{relative_path}"
        os.makedirs(os.path.dirname(relative_path), exist_ok=True)
        shutil.copy2(source_path, relative_path)
        print(f"Overlaid {relative_path}")
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pr_tree")
    parser.add_argument("yaml_paths", nargs="*")
    parser.add_argument("--remove", nargs="*", default=())
    args = parser.parse_args()
    raise SystemExit(main(args.pr_tree, args.yaml_paths, args.remove))
