"""
Copyright (c) Meta, Inc. and its affiliates.
Modifications Copyright (c) 2025 Samsung Electronics Co., Ltd.

This source code is licensed under the MIT license.
You may obtain a copy of the license at:
https://opensource.org/licenses/MIT

Code based on ase.optimize
"""

from __future__ import annotations

from types import MappingProxyType, SimpleNamespace
from typing import TYPE_CHECKING, ClassVar

import numpy as np
import torch
from ase import Atoms
from ase.calculators.calculator import PropertyNotImplementedError
from ase.calculators.singlepoint import SinglePointCalculator
from ase.geometry import wrap_positions
from ase.stress import voigt_6_to_full_3x3_stress
from torch_scatter import scatter

if TYPE_CHECKING:
    from collections.abc import Generator, Sequence

    from GGNN.trainer.utrainer import OCPTrainer
    from torch_geometric.data import Batch


ASE_PROP_RESHAPE = MappingProxyType(
    {"stress": (-1, 3, 3), "dielectric_tensor": (-1, 3, 3)}
)
ALL_CHANGES: set[str] = {"pos"}


def batch_to_atoms(
    batch: Batch,
    results: dict[str, torch.Tensor] | None = None,
    *,
    wrap_pos: bool = True,
    eps: float = 1e-7,
) -> list[Atoms]:
    """Convert a data batch to ase Atoms

    Args:
        batch: data batch
        results: dictionary with predicted result tensors that
        will be added to a SinglePointCalculator. If no results
            are given no calculator will be added to the atoms objects.
        wrap_pos: wrap positions back into the cell.
        eps: Small number to prevent slightly negative coordinates from being wrapped.

    Returns:
        list of Atoms
    """
    n_systems = batch.natoms.shape[0]
    n_atoms = batch.natoms.tolist()
    atomic_nums = torch.split(batch.atomic_numbers, n_atoms)
    bs = int((batch.batch.max() + 1).detach().cpu())
    if results is not None:
        results = {
            key: (
                val.view(ASE_PROP_RESHAPE.get(key, -1)).tolist()
                if len(val) == bs
                else [v.cpu().detach().numpy() for v in torch.split(val, n_atoms)]
            )
            for key, val in results.items()
        }

    positions = torch.split(batch.pos, n_atoms)

    atoms_objects = []
    for idx in range(n_systems):
        pos = positions[idx].cpu().detach().numpy()
        cell = batch.cell[idx].cpu().detach().numpy()

        # TODO take pbc from data
        if wrap_pos:
            pos = wrap_positions(pos, cell, pbc=[True, True, True], eps=eps)

        pbc = [True, True, True]
        atoms = Atoms(
            numbers=atomic_nums[idx].tolist(), cell=cell, positions=pos, pbc=pbc
        )

        if results is not None:
            calc = SinglePointCalculator(
                atoms=atoms, **{key: val[idx] for key, val in results.items()}
            )
            atoms.set_calculator(calc)

        atoms_objects.append(atoms)

    return atoms_objects


def compare_batches(
    batch1: Batch | SimpleNamespace | None,
    batch2: Batch,
    tol: float = 1e-6,
    excluded_properties: set[str] | None = None,
) -> bool:
    """Compare properties between two batches

    Args:
        batch1: atoms batch
        batch2: atoms batch
        tol: tolerance used to compare equility of floating point properties
        excluded_properties: list of properties to exclude from comparison

    Returns:
        list of system changes,
          property names that are difference between batch1 and batch2
    """
    system_changes = []

    if batch1 is None:
        system_changes = ALL_CHANGES
    else:
        properties_to_check = set(ALL_CHANGES)
        if excluded_properties:
            properties_to_check -= set(excluded_properties)

        # Check properties that aren't
        for prop in ALL_CHANGES:
            if prop in properties_to_check:
                properties_to_check.remove(prop)
                if not torch.allclose(
                    getattr(batch1, prop), getattr(batch2, prop), atol=tol
                ):
                    system_changes.append(prop)

    return len(system_changes) > 0


class OptimizableBatch:
    """A Batch version of ase Optimizable Atoms

    This class can be used with ML relaxations
    in fairchem.core.relaxations.ml_relaxation
    or in ase relaxations classes, i.e. ase.optimize.lbfgs
    """

    ignored_changes: ClassVar[set[str]] = set()

    def __init__(
        self,
        batch: Batch,  # list of ase atoms
        trainer: OCPTrainer,
        transform: torch.nn.Module | None = None,
        *,
        mask_converged: bool = True,
        numpy: bool = False,
        masked_eps: float = 1e-8,
        device: torch.device | str = "cuda",
    ) -> None:
        """Initialize Optimizable Batch

        Args:
            batch: A batch of atoms graph data
            model: An instance of a BaseTrainer derived class
            transform: graph transform
            mask_converged: if true will mask systems in batch
                that are already converged
            numpy: whether to cast results to numpy arrays
            masked_eps: masking systems that are converged
                when using ASE optimizers results in divisions by zero
                from zero differences in masked positions at future steps,
                we add a small number to prevent this.
        """

        self.device = device
        torch.set_default_dtype(torch.float32)
        self.trainer = trainer
        torch.set_default_dtype(torch.float64)

        self.batch = batch.to(self.device)
        self.transform = transform
        self.numpy = numpy
        self.mask_converged = mask_converged
        self._cached_batch: Batch | SimpleNamespace | None = None
        self._update_mask = None
        self.torch_results = {}
        self.results = {}
        self._eps = masked_eps

        self.otf_graph = True  # trainer._unwrapped_model.otf_graph
        if not self.otf_graph and "edge_index" not in self.batch:
            self.update_graph()

        self.cell_factor = self.batch.natoms

    @property
    def batch_indices(self) -> torch.Tensor:
        """Get the batch indices specifying
        which position/force corresponds to which batch."""
        return self.batch.batch

    @property
    def converged_mask(self) -> torch.Tensor | None:
        if self._update_mask is not None:
            return torch.logical_not(self._update_mask)
        return None

    @property
    def update_mask(self) -> torch.Tensor:
        if self._update_mask is None:
            return torch.ones(len(self.batch)).bool()
        return self._update_mask

    def check_state(self, batch: Batch, tol: float = 1e-12) -> bool:
        """Check for any system changes since last calculation."""
        return compare_batches(
            self._cached_batch,
            batch,
            tol=tol,
            excluded_properties=set(self.ignored_changes),
        )

    def _predict(self) -> None:
        """Run prediction if batch has any changes."""
        system_changes = self.check_state(self.batch)
        if system_changes:
            # convert batch to fp32
            batch_fp32 = self.batch.clone()
            for k, v in batch_fp32:
                if isinstance(v, torch.Tensor) and v.dtype == torch.float64:
                    batch_fp32[k] = v.to(torch.float32)
            res = self.trainer.predict(batch_fp32, per_image=False, disable_tqdm=True)
            stress = res["stress"].reshape(-1, 3, 3)

            self.torch_results = {
                "forces": res["forces"].to(torch.float64),
                "energy": res["energy"].to(torch.float64),
                "stress": stress.to(torch.float64),
            }  # reorder stress to voigt notation}

            # set batch key w.r.t ocp convention
            # save only subset of props in simple namespace instead of
            # cloning the whole batch to save memory
            changes = ALL_CHANGES - set(self.ignored_changes)
            self._cached_batch = SimpleNamespace(
                **{prop: self.batch[prop].clone() for prop in changes}
            )

    def get_property(self, name: str, *, no_numpy: bool = False) -> torch.Tensor:
        """Get a predicted property by name."""
        self._predict()
        if self.numpy:
            self.results = {
                key: pred.item() if pred.numel() == 1 else pred.cpu().numpy()
                for key, pred in self.torch_results.items()
            }
        else:
            self.results = self.torch_results

        if name not in self.results:
            raise PropertyNotImplementedError(f"{name} not present in this calculation")

        return (
            self.results[name]
            if no_numpy is False
            else self.torch_results[name].detach()
        )

    def get_positions(self) -> torch.Tensor:
        """Get the batch positions"""
        pos = self.batch.pos.clone()
        if self.numpy:
            if self.mask_converged:
                pos[~self.update_mask[self.batch.batch]] = self._eps
            pos = pos.cpu().numpy()

        return pos

    def set_positions(self, positions: torch.Tensor) -> None:
        """Set the atom positions in the batch."""
        if isinstance(positions, np.ndarray):
            positions = torch.tensor(positions)

        positions = positions.to(dtype=torch.float64, device=self.device)

        if self.mask_converged and self._update_mask is not None:
            mask = self.update_mask[self.batch.batch]
            self.batch.pos[mask] = positions[mask]
        else:
            self.batch.pos = positions

        self.update_graph()

    def get_forces(self, *, no_numpy: bool = False) -> torch.Tensor:
        """Get predicted batch forces."""
        return self.get_property("forces", no_numpy=no_numpy)

    def get_potential_energy(self) -> torch.Tensor:
        """Get predicted energy as the sum of all batch energies."""
        if (
            len(self.batch) == 1
        ):  # unfortunately batch size 1 returns a float, not a tensor
            return self.get_property("energy")
        return self.get_property("energy").sum()

    def get_potential_energies(self) -> torch.Tensor:
        """Get the predicted energy for each system in batch."""
        return self.get_property("energy")

    def get_cells(self) -> torch.Tensor:
        """Get batch crystallographic cells."""
        return self.batch.cell

    def set_cells(self, cells: torch.Tensor) -> None:
        """Set batch cells."""
        if self.batch.cell.shape != cells.shape:
            raise ValueError("Cell shape mismatch")
        if isinstance(cells, np.ndarray):
            cells = torch.tensor(cells, dtype=torch.float64, device=self.device)
        cells = cells.to(dtype=torch.float64, device=self.device)
        self.batch.cell[self.update_mask] = cells[self.update_mask]

    def get_volumes(self) -> torch.Tensor:
        """Get a tensor of volumes for each cell in batch"""
        cells = self.get_cells()
        return torch.linalg.det(cells)

    def iterimages(self) -> Generator[Batch, None, None]:
        yield self.batch

    def get_max_forces(self, forces: torch.Tensor | None = None) -> torch.Tensor:
        """Get the maximum forces per structure in batch"""
        if forces is None:
            forces = self.get_forces(no_numpy=True)
        return scatter((forces**2).sum(dim=1).sqrt(), self.batch_indices, reduce="max")

    def converged(
        self,
        forces: torch.Tensor | None,
        fmax: float,
        max_forces: torch.Tensor | None = None,
    ) -> bool:
        """Check if norm of all predicted forces are below fmax"""
        if forces is not None:
            if isinstance(forces, np.ndarray):
                forces = torch.tensor(forces, device=self.device, dtype=torch.float64)
            max_forces = self.get_max_forces(forces)
        elif max_forces is None:
            max_forces = self.get_max_forces()

        update_mask = max_forces.ge(fmax)
        # update cached mask
        if self.mask_converged:
            if self._update_mask is None:
                self._update_mask = update_mask
            else:
                # some models can have random noise in their predictions,
                # so the mask is updated by keeping all previously
                # converged structures masked even if new force predictions
                # push it slightly above threshold
                self._update_mask = torch.logical_and(self._update_mask, update_mask)
            update_mask = self._update_mask

        return not torch.any(update_mask).item()

    def get_atoms_list(self) -> list[Atoms]:
        """Get ase Atoms objects corresponding to the batch"""
        # self.update_graph()
        # self._predict()  # in case no predictions have been run
        return batch_to_atoms(self.batch, results=self.torch_results)

    def update_graph(self) -> None:
        """Update the graph if model does not use otf_graph."""
        return

    def __len__(self) -> int:
        # TODO: this might be changed in ASE to be 3 * len(self.atoms)
        return len(self.batch.pos)


class OptimizableUnitCellBatch(OptimizableBatch):
    """Modify the supercell and the atom positions in relaxations.

    Based on ase UnitCellFilter to work on data batches
    """

    def __init__(
        self,
        batch: Batch,  # list of ase atoms
        trainer: OCPTrainer,
        transform: torch.nn.Module | None = None,
        *,
        numpy: bool = False,
        mask_converged: bool = True,
        mask: torch.Tensor | Sequence[bool] | None = None,
        cell_factor: float | torch.Tensor | None = None,
        hydrostatic_strain: bool = False,
        constant_volume: bool = False,
        scalar_pressure: float = 0.0,
        masked_eps: float = 1e-8,
        device: torch.device | str = "cuda",
    ) -> None:
        super().__init__(
            batch=batch,
            trainer=trainer,
            transform=transform,
            numpy=numpy,
            mask_converged=mask_converged,
            masked_eps=masked_eps,
            device=device,
        )

        self.orig_cells = self.get_cells().clone()
        self.stress = None

        if mask is None:
            mask_tensor = torch.ones((3, 3), device=self.device)
        else:
            mask_tensor = torch.tensor(mask, device=self.device)

        # TODO make sure mask is on GPU
        if mask_tensor.shape == (6,):
            self.mask = torch.tensor(
                voigt_6_to_full_3x3_stress(mask_tensor.detach().cpu()),
                device=self.device,
            )
        elif mask_tensor.shape == (3, 3):
            self.mask = mask_tensor
        else:
            raise ValueError("shape of mask should be (3,3) or (6,)")

        if cell_factor is None:
            cell_factor = self.batch.natoms.repeat_interleave(3).unsqueeze(dim=1)

        self.hydrostatic_strain = hydrostatic_strain
        self.constant_volume = constant_volume
        self.pressure = scalar_pressure * torch.eye(3, device=self.device)
        self.cell_factor = cell_factor
        self.stress = None
        self._batch_trace = torch.vmap(torch.trace)
        self._batch_diag = torch.vmap(lambda x: x * torch.eye(3, device=x.device))
        self.deform_grad_ = (
            torch.stack([torch.eye(3) for i in range(self.orig_cells.shape[0])])
            .to(self.device)
            .to(torch.float64)
        )

    @property
    def batch_indices(self) -> torch.Tensor:
        """Get the batch indices specifying
        which position/force corresponds to which batch.

        We augment this to specify the batch indices for augmented positions and forces.
        """
        augmented_batch = torch.repeat_interleave(
            torch.arange(
                len(self.batch), dtype=self.batch.batch.dtype, device=self.device
            ),
            3,
        )
        return torch.cat([self.batch.batch, augmented_batch])

    @property
    def deform_grad(self) -> torch.Tensor:
        """Get the cell deformation matrix"""
        return self.deform_grad_
        # return torch.transpose(
        #     torch.linalg.solve(self.orig_cells, self.get_cells()), 1, 2
        # )

    def get_positions(self) -> torch.Tensor:
        """Get positions and cell deformation gradient."""
        cur_deform_grad = self.deform_grad
        n_atoms = self.batch.num_nodes
        pos = torch.zeros(
            (n_atoms + 3 * len(self.get_cells()), 3),
            dtype=self.batch.pos.dtype,
            device=self.device,
        )

        # Augmented positions are the self.atoms.positions
        #  but without the applied deformation gradient

        pos[:n_atoms] = torch.linalg.solve(
            cur_deform_grad[self.batch.batch, :, :],
            self.batch.pos.view(-1, 3, 1).detach(),
        ).view(-1, 3)
        # cell DOFs are the deformation gradient times a scaling factor

        pos[n_atoms:] = self.cell_factor * cur_deform_grad.view(-1, 3)
        return pos.cpu().numpy() if self.numpy else pos

    def set_positions(self, positions: torch.Tensor) -> None:
        """Set positions and cell.

        positions has shape (n_atoms + n_cells * 3, 3).
        the first n_atoms rows are the positions of the atoms,
        the last n_cells * three rows are the deformation tensor
        for each cell.
        """
        if isinstance(positions, np.ndarray):
            positions = torch.tensor(positions)

        # fp64
        # positions = positions.to(dtype=torch.float64, device=self.device)
        n_atoms = self.batch.num_nodes
        new_atom_positions = positions[:n_atoms]
        new_deform_grad = (positions[n_atoms:] / self.cell_factor).view(-1, 3, 3)
        self.deform_grad_ = new_deform_grad

        # TODO check that in fact symmetry is preserved setting cells and positions
        # Set the new cell from the original cell and the new deformation gradient.
        # Both current and final structures should preserve symmetry.

        new_cells = torch.bmm(
            self.orig_cells.to(torch.float64), torch.transpose(new_deform_grad, 1, 2)
        )
        self.set_cells(new_cells.to(torch.float64))

        # Set the positions from the ones passed in
        # (which are without the deformation gradient applied) and the new
        # deformation gradient. This should also preserve symmetry
        new_atom_positions = torch.bmm(
            new_atom_positions.view(-1, 1, 3),
            torch.transpose(
                new_deform_grad[self.batch.batch, :, :].view(-1, 3, 3), 1, 2
            ),
        )

        super().set_positions(new_atom_positions.view(-1, 3))

    def get_potential_energy(
        self,
    ) -> torch.Tensor:
        """
        returns potential energy including enthalpy PV term.
        """
        atoms_energy = super().get_potential_energy()
        return atoms_energy + self.pressure[0, 0] * self.get_volumes().sum()

    def get_forces(self, *, no_numpy: bool = False) -> torch.Tensor:
        """Get forces and unit cell stress."""
        stress = self.get_property("stress", no_numpy=True).view(-1, 3, 3).detach()
        atom_forces = self.get_property("forces", no_numpy=True).detach()

        volumes = self.get_volumes().view(-1, 1, 1)
        virial = -volumes * stress + self.pressure.view(-1, 3, 3)
        cur_deform_grad = self.deform_grad
        atom_forces = torch.bmm(
            atom_forces.view(-1, 1, 3),
            cur_deform_grad[self.batch.batch, :, :].view(-1, 3, 3),
        )

        virial = torch.linalg.solve(
            cur_deform_grad, torch.transpose(virial, dim0=1, dim1=2)
        )
        # enforce symmetry
        virial = 0.5 * (torch.transpose(virial, dim0=1, dim1=2) + virial)
        # TODO this does not work yet! maybe _batch_trace gives an issue
        if self.hydrostatic_strain:
            virial = self._batch_diag(self._batch_trace(virial) / 3.0)

        # Zero out components corresponding to fixed lattice elements
        if (self.mask != 1.0).any():
            virial *= self.mask.view(-1, 3, 3)

        if self.constant_volume:
            virial[:, range(3), range(3)] -= self._batch_trace(virial).view(3, -1) / 3.0

        n_atoms = self.batch.num_nodes
        augmented_forces = torch.zeros(
            (n_atoms + 3 * len(self.get_cells()), 3),
            device=self.device,
            dtype=atom_forces.dtype,
        )
        augmented_forces[:n_atoms] = atom_forces.view(-1, 3)
        augmented_forces[n_atoms:] = virial.view(-1, 3) / self.cell_factor

        self.stress = -virial.view(-1, 9) / volumes.view(-1, 1)

        if self.numpy and not no_numpy:
            augmented_forces = augmented_forces.cpu().numpy()

        return augmented_forces

    def __len__(self) -> int:
        return len(self.batch.pos) + 3 * len(self.batch)


class OptimizableFrechetBatch(OptimizableUnitCellBatch):
    """Modify the supercell and the atom positions in relaxations.

    Based on ase FretchetCellFilter to work on data batches
    """

    def __init__(
        self,
        batch: Batch,  # list of ase atoms
        trainer: OCPTrainer,
        transform: torch.nn.Module | None = None,
        *,
        numpy: bool = False,
        mask_converged: bool = True,
        mask: torch.Tensor | Sequence[bool] | None = None,
        cell_factor: float | torch.Tensor | None = None,
        hydrostatic_strain: bool = False,
        constant_volume: bool = False,
        scalar_pressure: float = 0.0,
        masked_eps: float = 1e-8,
        device: torch.device | str = "cuda",
    ) -> None:
        super().__init__(
            batch=batch,
            trainer=trainer,
            transform=transform,
            numpy=numpy,
            mask_converged=mask_converged,
            mask=mask,
            cell_factor=cell_factor,
            hydrostatic_strain=hydrostatic_strain,
            constant_volume=constant_volume,
            scalar_pressure=scalar_pressure,
            masked_eps=masked_eps,
            device=device,
        )
        self.exp_cell_factor = self.cell_factor

    def get_positions(self) -> torch.Tensor:
        pos = OptimizableUnitCellBatch.get_positions(self)
        n_atoms = self.batch.num_nodes
        cells = pos[n_atoms:].reshape(-1, 3, 3)
        cells = self.logm(cells)

        pos[n_atoms:] = self.exp_cell_factor * cells.reshape(-1, 3)
        return pos

    def set_positions(self, new: torch.Tensor) -> None:  # type: ignore[override]
        n_atoms = self.batch.num_nodes
        new2 = new.clone()
        batched_cell = (new[n_atoms:] / self.exp_cell_factor).reshape(-1, 3, 3)
        new2[n_atoms:] = self.expm(batched_cell).reshape(-1, 3)
        OptimizableUnitCellBatch.set_positions(self, new2)

    def get_forces(self, *, no_numpy: bool = False) -> torch.Tensor:
        # forces on atoms are same as UnitCellFilter, we just
        # need to modify the stress contribution
        stress = self.get_property("stress", no_numpy=True).view(-1, 3, 3)
        atom_forces = self.get_property("forces", no_numpy=True)
        volumes = self.get_volumes().view(-1, 1, 1)
        virial = -volumes * stress + self.pressure.view(-1, 3, 3)
        cur_deform_grad = self.deform_grad

        cur_deform_grad_log = self.logm(cur_deform_grad)
        atom_forces = torch.bmm(
            atom_forces.view(-1, 1, 3),
            cur_deform_grad[self.batch.batch, :, :].view(-1, 3, 3),
        )

        if self.hydrostatic_strain:
            virial = self._batch_diag(self._batch_trace(virial) / 3.0)

        # Zero out components corresponding to fixed lattice elements
        if (self.mask != 1.0).any():
            virial *= self.mask.view(-1, 3, 3)

        # Cell gradient for UnitCellFilter
        ucf_cell_grad = torch.bmm(
            virial, torch.linalg.inv(cur_deform_grad.transpose(1, 2))
        )
        ucf_cell_grad = 0.5 * (ucf_cell_grad + ucf_cell_grad.transpose(1, 2))
        # Cell gradient for FrechetCellFilter
        deform_grad_log_force = self.expm_frechet(cur_deform_grad_log, ucf_cell_grad)
        # Cauchy stress used for convergence testing
        stress = -(virial.view(-1, 9) / volumes.view(-1, 1))

        if self.constant_volume:
            raise NotImplementedError

        n_atoms = self.batch.num_nodes
        augmented_forces = torch.zeros(
            (n_atoms + 3 * len(self.get_cells()), 3),
            device=self.device,
            dtype=atom_forces.dtype,
        )
        augmented_forces[:n_atoms] = atom_forces.view(-1, 3)
        augmented_forces[n_atoms:] = (
            deform_grad_log_force.reshape(-1, 3) / self.exp_cell_factor
        )

        if self.numpy and not no_numpy:
            augmented_forces = augmented_forces.cpu().numpy()
        return augmented_forces

    def logm(self, a: torch.Tensor) -> torch.Tensor:
        # ensure A is symmetric

        #
        s, V = torch.linalg.eigh(a)
        return torch.bmm(
            torch.bmm(V, torch.diag_embed(torch.log(torch.abs(s)))), V.transpose(-1, -2)
        )

    def expm(self, a: torch.Tensor) -> torch.Tensor:
        return torch.linalg.matrix_exp(a)

    def expm_frechet(self, a: torch.Tensor, h: torch.Tensor) -> torch.Tensor:
        Z = torch.zeros(a.shape[0], 6, 6, dtype=torch.float64, device=self.device)
        Z[:, 0:3, 0:3] = a.detach()
        Z[:, 3:6, 3:6] = a.detach()
        Z[:, 0:3, 3:6] = h.detach()

        return torch.linalg.matrix_exp(Z)[:, 0:3, 3:6]
