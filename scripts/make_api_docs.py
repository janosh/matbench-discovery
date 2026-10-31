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
from types import SimpleNamespace

from lazydocs.generation import MarkdownGenerator, to_md_file

import matbench_discovery

ROOT = os.path.normpath(f"{os.path.dirname(__file__)}/..")
OUT_DIR = f"{ROOT}/site/src/routes/api/modules"
PREFIX = f"{matbench_discovery.__name__}."


def escape_svelte_braces(md_text: str) -> str:
    """Entity-escape braces outside code, which the site's markdown preprocessor
    would otherwise parse as Svelte expressions (e.g. "({} if absent)" in a docstring).
    """
    lines, in_fence = [], False
    for line in md_text.split("\n"):
        if line.lstrip().startswith("```"):
            in_fence = not in_fence
        elif not in_fence:
            # odd indices are inline code spans (``double`` or `single` backticks)
            segments = re.split(r"(``.+?``|`[^`]+`)", line)
            segments[::2] = [
                seg.replace("{", "&#123;").replace("}", "&#125;")
                for seg in segments[::2]
            ]
            lines.append("".join(segments))
            continue
        lines.append(line)
    return "\n".join(lines)


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

for mod_name in modules:
    module = importlib.import_module(mod_name)
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
        'src="https://img.shields.io/badge/source-blue?style=flat" alt="source link"',
    )
    md_text = escape_svelte_braces(md_text)
    file_name = mod_name.removeprefix(PREFIX).replace(".", "-")
    to_md_file(md_text, file_name, OUT_DIR, watermark=False)

print(f"Wrote docs for {len(modules)} modules to {OUT_DIR}")
