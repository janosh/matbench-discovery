"""Generate the site's API reference."""

# /// script
# dependencies = ["lazydocs @ git+https://github.com/ml-tooling/lazydocs"]
# ///

import annotationlib
import dataclasses
import functools
import importlib
import inspect
import json
import os
import pkgutil
import re
import shutil
import subprocess
import sys
from types import ModuleType, SimpleNamespace

import matbench_discovery

ROOT = os.path.normpath(f"{os.path.dirname(__file__)}/..")
OUT_DIR = f"{ROOT}/site/src/routes/api/modules"
PREFIX = f"{matbench_discovery.__name__}."


def import_api_module(module_name: str) -> ModuleType | None:
    """Omit optional inference modules when CPU documentation lacks their extras.

    Missing core dependencies and unexpected transitive import failures remain
    errors; only this model package's declared inference extras are optional.
    """
    try:
        return importlib.import_module(module_name)
    except ModuleNotFoundError as exc:
        optional_model = module_name.startswith(f"{PREFIX}mace_soyo_hyd.")
        if not optional_model or exc.name not in {"torch", "mace_soyo"}:
            raise
        print(
            f"Skipping optional inference API module {module_name}: missing {exc.name}",
            file=sys.stderr,
        )
        return None


def main() -> None:
    """Generate core API docs without requiring optional inference environments."""
    from lazydocs.generation import MarkdownGenerator, to_md_file

    inspect.signature = functools.partial(
        inspect.signature, annotation_format=annotationlib.Format.STRING
    )

    with open(f"{ROOT}/site/package.json") as file:
        repo_url = json.load(file)["repository"]
    commit = subprocess.check_output(
        ["git", "rev-parse", "HEAD"],  # noqa: S607
        cwd=ROOT,
        text=True,
    ).strip()

    shutil.rmtree(OUT_DIR, ignore_errors=True)
    os.makedirs(OUT_DIR)
    generator = MarkdownGenerator(ROOT, f"{repo_url}/blob/{commit}")
    walked_modules = pkgutil.walk_packages(matbench_discovery.__path__, PREFIX)
    modules = [matbench_discovery.__name__, *(mod.name for mod in walked_modules)]

    written = 0
    for mod_name in modules:
        module = import_api_module(mod_name)
        if module is None:
            continue
        for obj in vars(module).values():
            if dataclasses.is_dataclass(obj):
                for field in dataclasses.fields(obj):
                    field.type = SimpleNamespace(
                        __name__=inspect.formatannotation(
                            field.type, quote_annotation_strings=False
                        )
                    )

        md_text = generator.module2md(module).replace("<b>", "").replace("</b>", "")
        md_text = md_text.replace(".py#L0", ".py")
        md_text = re.sub(r" at 0x[0-9a-f]+(?=>)", "", md_text)
        md_text = re.sub(r'<a href="[^"]*%3Cstring%3E"><img[^>]*></a>\n*', "", md_text)
        md_text = md_text.replace(
            'src="https://img.shields.io/badge/-source-cccccc?style=flat-square"',
            'src="https://img.shields.io/badge/source-blue?style=flat" '
            'alt="source link"',
        )
        file_name = mod_name.removeprefix(PREFIX).replace(".", "-")
        to_md_file(md_text, file_name, OUT_DIR, watermark=False)
        written += 1

    print(f"Wrote docs for {written} modules to {OUT_DIR}")


if __name__ == "__main__":
    main()
