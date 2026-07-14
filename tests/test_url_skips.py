"""Require network URL skips to declare metadata in tests/url-skip-registry.yml."""

from __future__ import annotations

import ast
import os
import re
from glob import glob
from typing import Any, cast

import pytest
import yaml

from matbench_discovery import ROOT

SKIP_CALL_PATTERN = re.compile(r"pytest\.skip\(|@pytest\.mark\.skipif\(")
URL_SKIP_PATTERNS = (
    "figshare",
    "GitHub API",
    "not in data-files article",
    "rate-limited",
    "unreachable",
)


def _iter_test_python_files() -> list[str]:
    """Return pytest modules under tests/."""
    return sorted(glob(f"{ROOT}/tests/**/test_*.py", recursive=True))


def _authorized_skip_rules() -> dict[str, dict[str, list[str] | str]]:
    """Load and sanity-check authorized skip metadata."""
    with open(f"{ROOT}/tests/url-skip-registry.yml", encoding="utf-8") as file:
        registry = cast("dict[str, Any]", yaml.safe_load(file))
    rules = cast("dict[str, dict[str, list[str] | str]]", registry["skips"])
    for rule_key, rule in rules.items():
        tests = rule["tests"]
        match_substrings = rule["match_substrings"]
        reason = rule["reason"]
        assert isinstance(tests, list), f"{rule_key=} missing tests"
        assert tests, f"{rule_key=} missing tests"
        assert isinstance(match_substrings, list), (
            f"{rule_key=} missing match_substrings"
        )
        assert match_substrings, f"{rule_key=} missing match_substrings"
        assert isinstance(reason, str), f"{rule_key=} missing reason"
        assert reason, f"{rule_key=} missing reason"
    return rules


def _skip_literals(node: ast.AST) -> list[str]:
    """Extract string literals passed to pytest.skip/skipif in an AST subtree."""
    literals: list[str] = []
    for child in ast.walk(node):
        if (
            isinstance(child, ast.Call)
            and isinstance(child.func, ast.Attribute)
            and child.func.attr in {"skip", "skipif"}
        ):
            literals.extend(
                arg.value
                for arg in child.args
                if isinstance(arg, ast.Constant) and isinstance(arg.value, str)
            )
            literals.extend(
                keyword.value.value
                for keyword in child.keywords
                if keyword.arg == "reason"
                and isinstance(keyword.value, ast.Constant)
                and isinstance(keyword.value.value, str)
            )
    return literals


def test_url_skip_registry_entries_are_well_formed() -> None:
    """Authorized skip registry entries include tests, matchers, and reasons."""
    _authorized_skip_rules()


@pytest.mark.parametrize("test_file", _iter_test_python_files())
def test_url_related_skips_declare_registry_metadata(test_file: str) -> None:
    """Network URL skips must cite an entry in tests/url-skip-registry.yml."""
    with open(test_file, encoding="utf-8") as file:
        source = file.read()
    if not SKIP_CALL_PATTERN.search(source):
        return
    url_literals = [
        literal
        for literal in _skip_literals(ast.parse(source))
        if any(
            pattern.casefold() in literal.casefold() for pattern in URL_SKIP_PATTERNS
        )
    ]
    if not url_literals:
        return

    rules = _authorized_skip_rules()
    relative_test_file = os.path.relpath(test_file, ROOT)
    for literal in url_literals:
        matching_rules = [
            rule
            for rule in rules.values()
            if isinstance(rule["tests"], list)
            and relative_test_file in rule["tests"]
            and isinstance(rule["match_substrings"], list)
            and any(
                substring.casefold() in literal.casefold()
                for substring in rule["match_substrings"]
            )
        ]
        assert matching_rules, (
            f"{relative_test_file} has undocumented URL skip {literal!r}. "
            "Add metadata to tests/url-skip-registry.yml."
        )
