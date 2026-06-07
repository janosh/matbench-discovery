"""Overlay a model-submission PR's *data* onto a trusted checkout without ever
executing PR code.

Used by .github/workflows/update-site-figs.yml ingest mode: the workflow checks out
main (trusted code) plus the PR head into a side directory, then runs this script to
copy across only what the ingestion pipeline needs from the PR:

1. ``models/`` - model YAMLs (parsed with safe loaders only) and never-executed
   author files (test scripts are checked for existence by the PR checklist, not run)
2. ``matbench_discovery/enums.py`` - ONLY if the diff vs the trusted copy consists
   purely of added lines matching the strict ``Model`` member pattern (plus blank or
   comment lines). Anything else fails closed so a human runs ingestion locally after
   real review.

This makes secrets exposure during ingestion safe by construction: all *executed*
code comes from main; the PR contributes data only.
"""

import difflib
import re
import shutil
import sys
from pathlib import Path

# `model_name = auto(), "arch/model-file.yml"  # optional comment` (4-space indent)
MEMBER_LINE_RE = re.compile(
    r'^ {4}[a-z0-9_]+ = auto\(\), "(?P<path>[\w./@+-]+\.yml)"(\s*#.*)?$'
)
HARMLESS_LINE_RE = re.compile(r"^\s*(#.*)?$")  # blank or comment-only lines


def is_safe_member_line(line: str) -> bool:
    """Match a Model member line whose YAML path can't escape the models/ dir."""
    match = MEMBER_LINE_RE.match(line)
    return bool(match) and ".." not in match["path"]


ENUMS_REL_PATH = "matbench_discovery/enums.py"


def validate_enums_diff(trusted: str, submitted: str) -> list[str]:
    """Return a list of violations (empty = diff is safe to apply).

    Safe means: every change is an *insertion* of a Model-member, blank or comment
    line. Deletions or modifications of existing lines are never safe.
    """
    violations: list[str] = []
    trusted_lines, submitted_lines = trusted.splitlines(), submitted.splitlines()
    matcher = difflib.SequenceMatcher(
        None, trusted_lines, submitted_lines, autojunk=False
    )
    for op, idx1, idx2, jdx1, jdx2 in matcher.get_opcodes():
        if op == "equal":
            continue
        if op in ("delete", "replace"):
            removed = trusted_lines[idx1:idx2]
            violations += [f"removed/modified line: {line!r}" for line in removed]
        if op in ("insert", "replace"):
            for line in submitted_lines[jdx1:jdx2]:
                if not (is_safe_member_line(line) or HARMLESS_LINE_RE.match(line)):
                    violations.append(f"added non-member line: {line!r}")
    return violations


def main(pr_tree: str) -> int:
    """Overlay models/ + validated enums.py from ``pr_tree`` onto the cwd checkout."""
    pr_root = Path(pr_tree)
    if not (pr_root / "models").is_dir():
        print(f"::error::{pr_tree}/models not found - is this a repo checkout?")
        return 1

    # 1. models/ wholesale (YAML + author files; nothing in there is executed)
    shutil.copytree(pr_root / "models", "models", dirs_exist_ok=True)
    print(f"Overlaid {pr_tree}/models/ onto trusted checkout")

    # 2. enums.py only if the diff is provably additive Model members
    trusted_src = Path(ENUMS_REL_PATH).read_text()
    submitted_src = (pr_root / ENUMS_REL_PATH).read_text()
    if submitted_src == trusted_src:
        print(f"{ENUMS_REL_PATH} unchanged")
        return 0

    violations = validate_enums_diff(trusted_src, submitted_src)
    if violations:
        print(
            f"::error::{ENUMS_REL_PATH} diff goes beyond adding Model members - "
            "refusing to run PR code in CI. Review the PR and run "
            "`just ingest-model <model>` locally instead. Violations:"
        )
        for violation in violations:
            print(f"  {violation}")
        return 1

    Path(ENUMS_REL_PATH).write_text(submitted_src)
    print(f"Applied additive-only {ENUMS_REL_PATH} diff from PR")
    return 0


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(f"usage: {sys.argv[0]} <path-to-pr-checkout>")
    raise SystemExit(main(sys.argv[1]))
