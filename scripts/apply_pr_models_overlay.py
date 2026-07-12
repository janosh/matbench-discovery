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


def main(
    pr_tree: str,
    yaml_paths: Sequence[str],
    removed_paths: Sequence[str] = (),
) -> int:
    """Apply explicit, canonical model YAML changes to the trusted checkout."""
    yaml_paths = tuple(dict.fromkeys(yaml_paths))
    removed_paths = tuple(dict.fromkeys(removed_paths))
    all_paths = (*yaml_paths, *removed_paths)
    if not all_paths:
        print("::error::No changed model YAML paths supplied")
        return 1
    if os.path.islink(f"{pr_tree}/models"):
        print(f"::error::Unsafe submitted models directory: {pr_tree}/models")
        return 1

    for relative_path in all_paths:
        if not MODEL_YAML_PATTERN.fullmatch(relative_path):
            print(f"::error::Unsafe model YAML path: {relative_path!r}")
            return 1

    for relative_path in removed_paths:
        if not (os.path.isfile(relative_path) or os.path.islink(relative_path)):
            print(f"::error::Missing trusted YAML to remove: {relative_path}")
            return 1

    for relative_path in yaml_paths:
        source_path = f"{pr_tree}/{relative_path}"
        if (
            os.path.islink(os.path.dirname(source_path))
            or not os.path.isfile(source_path)
            or os.path.islink(source_path)
        ):
            print(f"::error::Missing or unsafe submitted YAML: {source_path}")
            return 1

    for relative_path in removed_paths:
        os.remove(relative_path)
        print(f"Removed {relative_path}")

    for relative_path in yaml_paths:
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
