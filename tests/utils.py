"""Shared helpers for repository script tests."""

import importlib.util
import sys
from types import ModuleType

from matbench_discovery import ROOT


def import_repo_script(module_name: str, rel_path: str) -> ModuleType:
    """Import a repository-local script without package-name collisions."""
    spec = importlib.util.spec_from_file_location(module_name, f"{ROOT}/{rel_path}")
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot import {module_name} from {rel_path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    try:
        spec.loader.exec_module(module)
    except BaseException:
        sys.modules.pop(module_name, None)
        raise
    return module
