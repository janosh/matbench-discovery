"""Test enums module."""

from matbench_discovery.enums import MbdKey, ModelType, Open, Task, TestSubset


def test_mbd_key() -> None:
    """Test MbdKey enum."""
    # Test basic enum functionality
    assert MbdKey.e_form_dft == "e_form_per_atom_mp2020_corrected"
    # HTML tags are part of the label, so we need to test the exact string
    assert MbdKey.e_form_dft.label == (
        "DFT E<sub>form</sub> "
        "<span style='font-size: 0.8em; font-weight: lighter;'>(eV/atom)</span>"
    )
    assert MbdKey.each_true == "e_above_hull_mp2020_corrected_ppd_mp"
    assert MbdKey.each_true.label == "E<sub>MP hull dist</sub>"

    # Test that all keys have values and labels
    for key in MbdKey:
        assert isinstance(key.value, str), f"{key=}"
        assert isinstance(key.label, str), f"{key=}"
        assert key.value != ""
        assert key.label != ""

    # Test uniqueness of values and labels
    values = [key.value for key in MbdKey]
    labels = [key.label for key in MbdKey]
    assert len(values) == len(set(values)), "Values must be unique"
    assert len(labels) == len(set(labels)), "Labels must be unique"


def test_task() -> None:
    """Test Task enum."""
    # Test basic enum functionality
    assert Task.S2E == "S2E"
    assert Task.S2E.label == "structure to energy"
    assert Task.S2EFS == "S2EFS"
    assert Task.S2EFS.label == "structure to energy, force, stress"

    # Test that all tasks have values and labels
    for task in Task:
        assert isinstance(task.value, str)
        assert isinstance(task.label, str)
        assert task.value != ""
        assert task.label != ""

    # Test task descriptions make sense
    assert "energy" in (Task.S2E.label or "")
    assert "force" in (Task.S2EF.label or "")
    assert "stress" in (Task.S2EFS.label or "")
    assert "magmoms" in (Task.S2EFSM.label or "")


def test_model_type() -> None:
    """Test ModelType enum."""
    # Test basic enum functionality
    assert ModelType.GNN == "GNN"
    assert ModelType.GNN.label == "Graph Neural Network"
    assert ModelType.RF == "RF"
    assert ModelType.RF.label == "Random Forest"

    # Test that all model types have values and labels
    for model_type in ModelType:
        assert isinstance(model_type.value, str)
        assert isinstance(model_type.label, str)
        assert model_type.value != ""
        assert model_type.label != ""

    # Test model type descriptions make sense
    assert "Neural" in (ModelType.GNN.label or "")
    assert "Forest" in (ModelType.RF.label or "")
    assert "Transformer" in (ModelType.Transformer.label or "")
    assert "Fingerprint" in (ModelType.Fingerprint.label or "")


def test_open() -> None:
    """Test Open enum."""
    # Test basic enum functionality
    assert Open.OSOD == "OSOD"
    assert Open.OSOD.label == "open source, open data"
    assert Open.CSCD == "CSCD"
    assert Open.CSCD.label == "closed source, closed data"

    # Test that all openness types have values and labels
    for open_type in Open:
        assert isinstance(open_type.value, str)
        assert isinstance(open_type.label, str)
        assert open_type.value != ""
        assert open_type.label != ""

    # Test openness descriptions make sense
    assert "open source" in Open.OSOD.label
    assert "open data" in Open.OSOD.label
    assert "closed source" in Open.CSCD.label
    assert "closed data" in Open.CSCD.label


def test_test_subset() -> None:
    """Test TestSubset enum."""
    # Test basic enum functionality
    assert TestSubset.uniq_protos == "unique_prototypes"
    assert TestSubset.uniq_protos.label == "Unique Structure Prototypes"
    assert TestSubset.most_stable_10k == "most_stable_10k"
    assert TestSubset.most_stable_10k.label == "10k Most Stable Materials"

    # Test that all test subsets have values and labels
    for subset in TestSubset:
        assert isinstance(subset.value, str)
        assert isinstance(subset.label, str)
        assert subset.value != ""
        assert subset.label != ""

    # Test subset descriptions make sense
    assert "Unique" in (TestSubset.uniq_protos.label or "")
    assert "Stable" in (TestSubset.most_stable_10k.label or "")
    assert "Full" in (TestSubset.full_test_set.label or "")
