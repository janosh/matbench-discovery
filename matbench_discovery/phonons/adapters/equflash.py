"""EquFlash two-stage relaxation and batched force evaluation."""

from __future__ import annotations

from collections.abc import Sequence
from typing import TYPE_CHECKING, Any

import numpy as np

from matbench_discovery.phonons.adapters.standard import StandardKappaAdapter

if TYPE_CHECKING:
    from kemio import Structure
    from kemio.model import Model as ForceModel

    from matbench_discovery.phonons.pipeline import KappaSettings


class EquFlashBatchedModel:
    """EquFlash forces for a batch of structures through its trainer's graph API;
    other property requests go through the wrapped ASE calculator.
    """

    def __init__(self, calculator: object) -> None:
        """Wrap an EquFlash ASE calculator exposing ``trainer``."""
        trainer = getattr(calculator, "trainer", None)
        if trainer is None:
            raise TypeError("EquFlash calculator has no trainer for batched forces")
        self.calculator = calculator
        self.trainer = trainer

    def compute(
        self, structures: Sequence[Structure], properties: Sequence[str]
    ) -> list[dict[str, Any]]:
        """Evaluate a batch; forces-only requests run as one GPU graph batch."""
        from kemio.interop import AseCalculatorModel, structure_to_ase

        if list(properties) != ["forces"]:
            return AseCalculatorModel(self.calculator).compute(structures, properties)
        from fairchem.core.preprocessing import AtomsToGraphs
        from torch_geometric.data import Batch

        device = getattr(self.trainer, "device", None)
        if device is None:
            import torch

            device = "cuda" if torch.cuda.is_available() else "cpu"
        converter = AtomsToGraphs(r_edges=False)
        graphs = [converter.convert(structure_to_ase(struct)) for struct in structures]
        graph_batch = Batch.from_data_list(graphs).to(device)
        outputs = self.trainer.predict(graph_batch, per_image=False, disable_tqdm=True)
        forces = outputs["forces"].detach().cpu().numpy()
        offsets = np.cumsum([0, *(len(struct) for struct in structures)])
        return [
            {"forces": forces[offsets[idx] : offsets[idx + 1]]}
            for idx in range(len(structures))
        ]


class EquFlashKappaAdapter(StandardKappaAdapter):
    """Use the shared two-stage relaxation and EquFlash's batched forces."""

    name = "equflash"

    def prepare_model(
        self,
        model: ForceModel,
        settings: KappaSettings,  # noqa: ARG002
    ) -> ForceModel:
        """Route force batches through the EquFlash trainer."""
        return EquFlashBatchedModel(getattr(model, "calculator", model))
