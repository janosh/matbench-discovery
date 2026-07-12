"""EquFlash two-stage relaxation and batched FC3 evaluation."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

from matbench_discovery.phonons.adapters.standard import StandardKappaAdapter
from matbench_discovery.phonons.thermal_conductivity import batched_displacement_forces

if TYPE_CHECKING:
    import numpy as np
    from ase import Atoms
    from ase.calculators.calculator import Calculator
    from phono3py.api_phono3py import Phono3py

    from matbench_discovery.phonons.pipeline import KappaSettings


class EquFlashKappaAdapter(StandardKappaAdapter):
    """Use the shared two-stage relaxation and EquFlash's batched FC3."""

    name = "equflash"

    def calculate_fc3(
        self,
        phono3py: Phono3py,
        calculator: Calculator,
        settings: KappaSettings,
        *,
        progress: dict[str, Any] | None = None,  # noqa: ARG002
        max_evaluations: int | None = None,
    ) -> np.ndarray:
        """Evaluate EquFlash FC3 displacement graphs in GPU batches."""
        from fairchem.core.preprocessing import AtomsToGraphs
        from torch_geometric.data import Batch

        n_atoms = len(phono3py.supercell)
        graph_converter = AtomsToGraphs(r_edges=False)
        batch_size = settings.batch_size
        if settings.max_atoms_per_batch is not None:
            batch_size = max(settings.max_atoms_per_batch // n_atoms, 1)

        trainer = getattr(calculator, "trainer", None)
        if trainer is None:
            raise TypeError("EquFlash calculator has no trainer for batched FC3")
        device = getattr(trainer, "device", None)
        if device is None:
            import torch

            device = "cuda" if torch.cuda.is_available() else "cpu"

        def evaluate_batch(batch_atoms: list[Atoms]) -> np.ndarray:
            """Evaluate one EquFlash displacement batch through its trainer."""
            graphs = [
                graph_converter.convert(displaced_atoms)
                for displaced_atoms in batch_atoms
            ]
            graph_batch = Batch.from_data_list(graphs).to(device)  # ty: ignore[unresolved-attribute]
            outputs = trainer.predict(graph_batch, per_image=False, disable_tqdm=True)
            return outputs["forces"].detach().cpu().numpy().reshape(-1, n_atoms, 3)

        force_set = batched_displacement_forces(
            phono3py.supercells_with_displacements,
            evaluate_batch,
            batch_size=batch_size,
            n_atoms=n_atoms,
            max_evaluations=max_evaluations,
        )
        phono3py.forces = force_set
        return force_set
