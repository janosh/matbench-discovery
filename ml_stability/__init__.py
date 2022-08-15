import os
from os.path import dirname
from typing import Any, Generator, Sequence


PKG_DIR = dirname(__file__)
ROOT = dirname(PKG_DIR)

os.makedirs(f"{PKG_DIR}/plots", exist_ok=True)


def chunks(xs: Sequence[Any], n: int) -> Generator[Sequence[Any], None, None]:
    return (xs[i : i + n] for i in range(0, len(xs), n))


def as_dict_handler(obj: Any) -> dict[str, Any] | None:
    """Use as default_handler kwarg to json.dump() or pandas.to_json()."""
    try:
        return obj.as_dict()  # all MSONable objects implement as_dict()
    except AttributeError:
        return None  # replace unhandled objects with None in serialized data
