from __future__ import annotations

import os
from collections.abc import Generator, Sequence
from typing import Any

PKG_DIR = os.path.dirname(__file__)
ROOT = os.path.dirname(PKG_DIR)


def chunks(xs: Sequence[Any], n: int) -> Generator[Sequence[Any], None, None]:
    return (xs[i : i + n] for i in range(0, len(xs), n))


def as_dict_handler(obj: Any) -> dict[str, Any] | None:
    """Use as default_handler kwarg to json.dump() or pandas.to_json()."""
    try:
        return obj.as_dict()  # all MSONable objects implement as_dict()
    except AttributeError:
        return None  # replace unhandled objects with None in serialized data
        # removes e.g. non-serializable AseAtoms from M3GNet relaxation trajectories
