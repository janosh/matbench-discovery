"""Regress optional inference imports without suppressing real documentation errors."""

from types import ModuleType

import pytest

import scripts.make_api_docs as docs


@pytest.mark.parametrize("dependency", ["torch", "mace_soyo"])
@pytest.mark.parametrize("module", ["adapter", "calculator"])
def test_optional_model_dependency_missing(
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
    dependency: str,
    module: str,
) -> None:
    """CPU docs omit the optional GPU adapter while reporting its missing extra."""

    def missing(_name: str) -> ModuleType:
        raise ModuleNotFoundError("optional inference extra missing", name=dependency)

    monkeypatch.setattr(docs.importlib, "import_module", missing)
    name = f"matbench_discovery.mace_soyo_hyd.{module}"
    assert docs.import_api_module(name) is None
    assert name in capsys.readouterr().err


@pytest.mark.parametrize(
    ("module", "dependency"),
    [
        ("matbench_discovery.data", "torch"),
        ("matbench_discovery.mace_soyo_hyd.adapter", "ase"),
        ("matbench_discovery.mace_soyo_hyd.calculator", "mace_soyo.broken"),
    ],
)
def test_unexpected_import_failure_is_fatal(
    monkeypatch: pytest.MonkeyPatch, module: str, dependency: str
) -> None:
    """Core errors and broken installed inference packages must still fail CI."""
    error = ModuleNotFoundError("unexpected dependency failure", name=dependency)

    def missing(_name: str) -> ModuleType:
        raise error

    monkeypatch.setattr(docs.importlib, "import_module", missing)
    with pytest.raises(ModuleNotFoundError) as caught:
        docs.import_api_module(module)
    assert caught.value is error


def test_available_optional_module_is_documented(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """An installed optional module is returned normally rather than excluded."""
    module = ModuleType("matbench_discovery.mace_soyo_hyd.calculator")
    monkeypatch.setattr(docs.importlib, "import_module", lambda _name: module)
    assert docs.import_api_module(module.__name__) is module
