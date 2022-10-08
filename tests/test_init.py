import os
from typing import Any

from mb_discovery import as_dict_handler, chunks


def test_has_root_pkg_dir() -> None:

    from mb_discovery import PKG_DIR, ROOT

    assert os.path.isdir(ROOT)
    assert os.path.isdir(PKG_DIR)


def test_chunks() -> None:

    assert list(chunks([], 1)) == []
    assert list(chunks([1], 1)) == [[1]]
    assert list(chunks([1, 2], 1)) == [[1], [2]]
    assert list(chunks([1, 2, 3], 1)) == [[1], [2], [3]]
    assert list(chunks([1, 2, 3], 2)) == [[1, 2], [3]]
    assert list(chunks(range(1, 4), 2)) == [range(1, 3), range(3, 4)]
    assert list(chunks(range(1, 5), 2)) == [range(1, 3), range(3, 5)]
    assert list(chunks(range(1, 5), 3)) == [range(1, 4), range(4, 5)]


def test_as_dict_handler() -> None:
    class C:
        def as_dict(self) -> dict[str, Any]:
            return {"foo": "bar"}

    assert as_dict_handler(C()) == {"foo": "bar"}
    assert as_dict_handler(1) is None
    assert as_dict_handler("foo") is None
    assert as_dict_handler([1, 2, 3]) is None
    assert as_dict_handler({"foo": "bar"}) is None
