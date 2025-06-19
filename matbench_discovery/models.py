"""Model utilities for matbench-discovery."""

from typing import Any

from matbench_discovery.enums import ModelType, Open


def model_is_compliant(metadata: dict[str, str | list[str]]) -> bool:
    """Check if model complies with benchmark restrictions on allowed training sets and
    open model code and weights.

    Args:
        metadata: model metadata dictionary

    Returns:
        bool: True if model is compliant
    """
    openness = metadata.get("openness", Open.OSOD)
    if openness != Open.OSOD:
        return False
    training_sets = metadata.get("training_set")

    if not isinstance(training_sets, list):
        model_name = metadata.get("model_name")
        raise TypeError(
            f"{model_name}: expected list of training sets, got {training_sets=}"
        )

    return set(training_sets) <= {"MP 2022", "MPtrj", "MPF", "MP Graphs"}


def validate_model_metadata(metadata: dict[str, Any], metadata_file: str) -> None:
    """Validate model metadata for required fields and types.

    Args:
        metadata: model metadata dictionary
        metadata_file: path to metadata file for error reporting

    Raises:
        ValueError: if validation fails
    """
    if metadata.get("status", "complete") != "complete":  # skip incomplete models
        raise ValueError(f"Model {metadata_file} has status != 'complete'")

    try:  # check if model_type is valid
        ModelType(metadata.get("model_type"))  # type: ignore[arg-type]
    except ValueError as exc:
        exc.add_note(f"{metadata_file=}\nPick from {', '.join(ModelType)}")
        raise
