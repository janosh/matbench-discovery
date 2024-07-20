import inspect
from enum import StrEnum

import matbench_discovery.enums
from matbench_discovery.enums import LabelEnum


def test_label_enum() -> None:
    # assert all classes in enums are LabelEnums
    for enum in dir(matbench_discovery.enums):
        if inspect.isclass(enum) and issubclass(enum, StrEnum):
            assert issubclass(enum, LabelEnum)
            val_dict = enum.val_dict()
            assert isinstance(val_dict, dict)
            label_dict = enum.label_dict()
            assert isinstance(label_dict, dict)
            assert val_dict != label_dict
            assert len(val_dict) == len(label_dict)
