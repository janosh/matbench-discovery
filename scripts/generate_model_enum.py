"""Generate the ``Model`` enum member block from model YAML metadata."""

import keyword
import os
import re
from glob import glob

import yaml

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ENUM_PATH = f"{ROOT}/matbench_discovery/enums.py"
BEGIN_MARKER = "    # BEGIN GENERATED MODEL MEMBERS"
END_MARKER = "    # END GENERATED MODEL MEMBERS"


def generate_source(source: str) -> str:
    """Replace the generated member block with entries from model YAMLs."""
    if (source.count(BEGIN_MARKER), source.count(END_MARKER)) != (1, 1):
        raise ValueError("Expected one generated Model marker pair")
    prefix, remainder = source.split(BEGIN_MARKER)
    _, suffix = remainder.split(END_MARKER)
    members: list[tuple[str, str]] = []
    seen_names: set[str] = set()
    for yaml_path in glob(f"{ROOT}/models/[!_]*/[!_]*.yml"):
        with open(yaml_path, encoding="utf-8") as file:
            metadata = yaml.safe_load(file)
        if not isinstance(metadata, dict):
            raise TypeError(f"{yaml_path} must contain a YAML mapping")
        if metadata.get("status") == "aborted":
            continue
        model_key = metadata.get("model_key")
        if not isinstance(model_key, str):
            raise TypeError(f"{yaml_path} has invalid {model_key=}")
        name = re.sub(r"[^a-z0-9]+", "_", model_key.casefold()).strip("_")
        if not name.isidentifier() or keyword.iskeyword(name) or name in seen_names:
            raise ValueError(f"{yaml_path} has invalid or duplicate enum name {name!r}")
        seen_names.add(name)
        rel_path = os.path.relpath(yaml_path, f"{ROOT}/models").replace("\\", "/")
        members.append((name, rel_path))
    member_source = "".join(
        f'    {name} = auto(), "{path}"\n' for name, path in sorted(members)
    )
    return f"{prefix}{BEGIN_MARKER}\n{member_source}{END_MARKER}{suffix}"


def main() -> int:
    """Update the generated enum member block."""
    with open(ENUM_PATH, encoding="utf-8") as file:
        source = file.read()
    generated_source = generate_source(source)
    if generated_source == source:
        return 0
    with open(ENUM_PATH, mode="w", encoding="utf-8") as file:
        file.write(generated_source)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
