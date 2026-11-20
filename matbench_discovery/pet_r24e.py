"""PET-OMat-XS with tetrahedral averaging and fitted elemental references."""

from __future__ import annotations

import json
from collections.abc import Sequence
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.calculators.calculator import (
    Calculator,
    PropertyNotImplementedError,
    all_changes,
)
from ase.stress import full_3x3_to_voigt_6_stress, voigt_6_to_full_3x3_stress
from scipy.spatial.transform import Rotation

from matbench_discovery import ROOT, file_digest

BASE_SHA256 = "64235c15f9351577db36068a833c36ba2bc0682f77121345926668c58c582299"
PARAMETERS = Path(ROOT) / "models/pet/pet-omat-xs-r24e/frozen-model.json"
PARAMETERS_SHA256 = "feaf9c6e41c932a56f8ac0e8baa9c993f6a399a10f48188ae9993591a8ad5301"


def rotations() -> np.ndarray:
    """The 12 proper tetrahedral rotations and their spatial inverses, in fit order."""
    proper = Rotation.create_group("T").as_matrix()
    return np.concatenate((proper, -proper))


class R24E(Calculator):
    """Average conservative PET energies before adding composition-only references."""

    implemented_properties = ("energy", "forces", "stress")

    def __init__(
        self,
        base: Calculator,
        *,
        intercept: float,
        offsets: Sequence[float] | np.ndarray,
        native_batch: bool = False,
    ) -> None:
        """Attach the parent calculator and fitted per-atom energy references."""
        super().__init__()
        self.base = base
        self.intercept = intercept
        self.offsets = np.asarray(offsets, dtype=float)
        if self.offsets.shape != (119,) or not np.isfinite(self.offsets).all():
            raise ValueError("Expected 119 finite elemental offsets")
        self.rotations = rotations()
        self.native_batch = native_batch

    def calculate(
        self,
        atoms: Atoms | None = None,
        properties: Sequence[str] = ("energy",),
        system_changes: Sequence[str] = all_changes,
    ) -> None:
        """Compute the requested ASE properties in the input frame."""
        super().calculate(atoms, properties, system_changes)
        atoms = self.atoms
        periodic = atoms.pbc.all() and atoms.cell.rank == 3
        if "stress" in properties and not periodic:
            raise PropertyNotImplementedError(
                "Bulk stress requires a fully periodic cell"
            )
        batch_size = min(12, max(1, 384 // len(atoms)))
        if self.native_batch and periodic and batch_size > 1:
            self.results = self._batched_bulk(atoms, batch_size)
        else:
            energy = 0.0
            forces = np.zeros((len(atoms), 3))
            stress = np.zeros((3, 3))
            props = ["energy", "forces", "stress"] if periodic else ["energy", "forces"]
            for rotation in self.rotations:
                moved = atoms.copy()
                moved.set_cell(atoms.cell.array @ rotation.T, apply_constraint=False)
                moved.set_positions(
                    atoms.positions @ rotation.T, apply_constraint=False
                )
                self.base.calculate(moved, props, all_changes)
                value = self.base.results
                energy += float(value["energy"])
                forces += value["forces"] @ rotation
                if periodic:
                    tensor = voigt_6_to_full_3x3_stress(value["stress"])
                    stress += rotation.T @ tensor @ rotation
            self.results = {"energy": energy / 24, "forces": forces / 24}
            if periodic:
                self.results["stress"] = full_3x3_to_voigt_6_stress(stress / 24)
        self.results["energy"] += (
            self.intercept * len(atoms) + self.offsets[atoms.numbers].sum()
        )

    def _batched_bulk(
        self, atoms: Atoms, batch_size: int
    ) -> dict[str, float | np.ndarray]:
        """Differentiate the averaged energy in the common, unrotated ASE frame."""
        import torch
        from metatomic.torch import ModelEvaluationOptions, ModelOutput, System
        from metatomic_ase._calculator import _ase_to_torch_data, _get_ase_input

        base = self.base.calculator
        model = base.model()
        key = base._energy_key  # noqa: SLF001
        dtype, device = base._dtype, base._device  # noqa: SLF001
        options = ModelEvaluationOptions(
            length_unit="angstrom",
            outputs={key: ModelOutput(unit="eV", sample_kind="system")},
            selected_atoms=None,
        )
        types, pos, cell, pbc = _ase_to_torch_data(atoms, dtype=dtype, device=device)
        matrices = torch.as_tensor(self.rotations, dtype=dtype, device=device)
        energy, force, stress = 0.0, np.zeros((len(atoms), 3)), np.zeros((3, 3))
        for chunk in matrices.split(batch_size):
            coords = pos.detach().requires_grad_(requires_grad=True)
            strain = torch.eye(3, dtype=dtype, device=device, requires_grad=True)
            systems = [
                System(types, (coords @ strain) @ r.T, (cell @ strain) @ r.T, pbc)
                for r in chunk
            ]
            systems = base._nl_calculators.compute(systems=systems)  # noqa: SLF001
            for system in systems:
                for name, spec in model.requested_inputs(use_new_names=True).items():
                    system.add_data(
                        name,
                        _get_ase_input(atoms, name, spec, dtype=dtype, device=device),
                    )
            values = (
                model(systems, options, check_consistency=False)[key].block().values
            )
            if len(values) != len(chunk):
                raise ValueError("Incomplete PET orientation batch")
            value = values.sum() / len(self.rotations)
            dp, ds = torch.autograd.grad(value, (coords, strain))
            energy += value.item()
            force -= dp.detach().double().cpu().numpy()
            stress += ds.detach().double().cpu().numpy() / atoms.get_volume()
        return {
            "energy": energy,
            "forces": force,
            "stress": full_3x3_to_voigt_6_stress((stress + stress.T) / 2),
        }


def load(device: str, checkpoint: str, dtype: str = "float32") -> Calculator:
    """Load the exact retained checkpoint and elemental coefficients."""
    from upet.ase import UPETCalculator

    from matbench_discovery.calculators import CHECKPOINT_DIR, _stage_checkpoint

    if (
        file_digest(checkpoint) != BASE_SHA256
        or file_digest(str(PARAMETERS)) != PARAMETERS_SHA256
    ):
        raise ValueError(
            "Frozen PET-OMat-XS checkpoint or R24E coefficients do not match"
        )
    # UPET 0.3.1 parses model, size and version from the basename, including for
    # local overrides. Shared runners store downloads under a hash-based filename.
    canonical = "pet-omat-xs-v1.0.0.ckpt"
    if Path(checkpoint).name != canonical:
        target = str(Path(CHECKPOINT_DIR) / "pet-omat-xs-r24e" / canonical)
        _stage_checkpoint(checkpoint, target)
        checkpoint = target
    fit = json.loads(PARAMETERS.read_text(encoding="utf-8"))
    base = UPETCalculator(
        checkpoint_path=checkpoint, device=device, dtype=dtype, non_conservative=False
    )
    return R24E(
        base, intercept=fit["intercept"], offsets=fit["offsets"], native_batch=True
    )
