"""
Copyright (c) Meta, Inc. and its affiliates.
Modifications Copyright (c) 2025 Samsung Electronics Co., Ltd.

This source code is licensed under the MIT license.
You may obtain a copy of the license at:
https://opensource.org/licenses/MIT
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import ase
import torch
from torch_scatter import scatter

if TYPE_CHECKING:
    from collections.abc import Callable

    from ..optimizable import OptimizableBatch  # noqa: TID252


class FIRE:
    """FIRE optimizer for batch ML relaxations."""

    def __init__(
        self,
        optimizable_batch: OptimizableBatch,
        dt: float = 0.1,
        maxstep: float | None = None,
        dtmax: float = 1.0,
        nmin: int = 5,
        finc: float = 1.1,
        fdec: float = 0.5,
        astart: float = 0.1,
        fa: float = 0.99,
        a: float = 0.1,
        *,
        downhill_check: bool = False,
        position_reset_callback: Callable | None = None,
        save_full_traj: bool = True,
        traj_dir: Path | None = None,
        traj_names: list[str] | None = None,
    ) -> None:
        """
        Args:
            optimizable_batch: an optimizable batch which
            includes a model and a batch of data
            maxstep: largest step that any atom is allowed to move
            memory: Number of steps to be stored in memory
            damping: The calculated step is
            multiplied with this number before added to the positions.
            alpha: Initial guess for the Hessian (curvature of energy surface)
            save_full_traj: whether to save full trajectory
            traj_dir: path to save trajectories in
            traj_names: list of trajectory files names
        """
        self.optimizable = optimizable_batch
        self.maxstep = maxstep
        self.save_full = save_full_traj
        self.traj_dir = traj_dir
        self.traj_names = traj_names
        self.trajectories = None

        self.fmax = None
        self.steps = None

        self.v = None
        batch_size = len(self.optimizable.batch.natoms)
        device = self.optimizable.batch.batch.device
        self.dt = dt * torch.ones(
            batch_size,
            dtype=torch.float64,
            device=device,
        )
        self.Nsteps = torch.zeros(
            batch_size,
            dtype=torch.int,
            device=device,
        )

        if maxstep is not None:
            self.maxstep = maxstep
        else:
            self.maxstep = 0.2

        self.dtmax = torch.tensor(dtmax)
        self.Nmin = nmin
        self.finc = finc
        self.fdec = fdec
        self.astart = astart

        self.fa = fa

        self.a = a * torch.ones(
            batch_size,
            dtype=torch.float64,
            device=device,
        )

        self.downhill_check = downhill_check
        self.position_reset_callback = position_reset_callback
        if self.traj_dir and (not traj_dir or not len(traj_names)):
            raise ValueError(
                "Trajectory names should be specified to save trajectories"
            )

    def run(self, fmax: float, steps: int) -> tuple[torch.Tensor, torch.Tensor]:
        self.fmax = fmax
        self.steps = steps

        self.trajectories = None
        if self.traj_dir:
            self.traj_dir.mkdir(exist_ok=True, parents=True)
            self.trajectories = [
                ase.io.Trajectory(self.traj_dir / f"{name}.traj_tmp", mode="w")
                for name in self.traj_names
            ]

        iteration = 0
        max_forces = self.optimizable.get_max_forces()
        n_traj = torch.zeros(self.optimizable.batch.batch.max() + 1, dtype=torch.int)
        while iteration < steps and not self.optimizable.converged(
            forces=None, fmax=self.fmax, max_forces=max_forces
        ):
            self.iteration = iteration
            if self.trajectories is not None and (
                self.save_full is True or iteration == 0
            ):
                self.write()
            self.step()
            max_forces = self.optimizable.get_max_forces()
            iteration += 1
            n_traj = (
                n_traj
                + torch.ones_like(n_traj)
                * self.optimizable.update_mask.detach().cpu()
            )

        # save after converged or all iterations ran
        if iteration > 0 and self.trajectories is not None:
            self.write(force_write=True)

        # GPU memory usage as per nvidia-smi seems to gradually build up as
        # batches are processed. This releases unoccupied cached memory.
        torch.cuda.empty_cache()

        if self.trajectories is not None:
            for traj in self.trajectories:
                traj.close()
            for name in self.traj_names:
                traj_fl = Path(self.traj_dir / f"{name}.traj_tmp", mode="w")
                traj_fl.rename(traj_fl.with_suffix(".traj"))

        # set predicted values to batch
        for name, value in self.optimizable.results.items():
            setattr(self.optimizable.batch, name, value)

        return n_traj, self.optimizable.converged(
            forces=None, fmax=self.fmax, max_forces=max_forces
        )

    def determine_step(self, dr: torch.Tensor) -> torch.Tensor:
        steplengths = torch.norm(dr, dim=1)
        longest_steps = scatter(
            steplengths, self.optimizable.batch_indices, reduce="max"
        )
        longest_steps = longest_steps[self.optimizable.batch_indices]
        maxstep = longest_steps.new_tensor(self.maxstep)
        scale = (longest_steps + 1e-12).reciprocal() * torch.min(longest_steps, maxstep)
        dr *= scale.unsqueeze(1)
        return dr * self.damping

    def _batched_dot(self, x: torch.Tensor, y: torch.Tensor) -> torch.Tensor:
        return scatter(
            (x * y).sum(dim=-1), self.optimizable.batch_indices, reduce="sum"
        )

    def step(self) -> None:
        f = self.optimizable.get_forces()
        if self.v is None:
            self.v = torch.zeros(
                (len(self.optimizable), 3), device=f.device, dtype=f.dtype
            )
            if self.downhill_check:
                self.e_last = self.optimizable.get_potential_energy()
                self.r_last = self.optimizable.get_positions().copy()
                self.v_last = self.v.copy()
        else:
            is_uphill = False
            if self.downhill_check:
                raise NotImplementedError

                e = self.optimizable.get_potential_energy()
                # Check if the energy actually decreased
                if e > self.e_last:
                    # If not, reset to old positions...
                    if self.position_reset_callback is not None:
                        self.position_reset_callback(
                            self.optimizable, self.r_last, e, self.e_last
                        )
                    self.optimizable.set_positions(self.r_last)
                    is_uphill = True
                self.e_last = self.optimizable.get_potential_energy()
                self.r_last = self.optimizable.get_positions().copy()
                self.v_last = self.v.clone()

            vf = self._batched_dot(f, self.v)

            # batch index for vf>0.0 and not is_uphill
            batch_idx_ = (not is_uphill) * torch.ones(
                self.dt.shape[0], dtype=torch.bool, device=self.dt.device
            )
            batch_idx_ = batch_idx_ * (vf > 0.0)

            # update v
            denorm = (
                torch.sqrt(self._batched_dot(f, f))
                / torch.sqrt(self._batched_dot(self.v, self.v))
            )[self.optimizable.batch_indices].unsqueeze(1)

            coef = self.a[self.optimizable.batch_indices].unsqueeze(1)
            self.v = (1 - coef) * self.v + coef * f / denorm
            self.v = self.v * batch_idx_[self.optimizable.batch_indices].unsqueeze(1)

            self.a[~batch_idx_] = self.astart
            self.dt[~batch_idx_] = self.dt[~batch_idx_] * self.fdec
            self.Nsteps[~batch_idx_] = 0

            update_dt_idx = batch_idx_ * self.Nsteps > self.Nmin

            self.a[update_dt_idx] = self.a[update_dt_idx] * self.fa
            self.dt[update_dt_idx] = torch.minimum(
                self.dt[update_dt_idx] * self.finc, self.dtmax
            )
            self.Nsteps[batch_idx_] = self.Nsteps[batch_idx_] + 1

            # if vf > 0.0 and not is_uphill:
            #     self.v = (1.0 - self.a) * self.v + self.a * f / torch.sqrt(
            #         torch.vdot(f, f)
            #     ) * torch.sqrt(torch.vdot(self.v, self.v))
            #     if self.Nsteps > self.Nmin:
            #         self.dt = min(self.dt * self.finc, self.dtmax)
            #         self.a *= self.fa
            #     self.Nsteps += 1
            # else:
            #     self.v[:] *= 0.0
            #     self.a = self.astart
            #     self.dt *= self.fdec
            #     self.Nsteps = 0
        batch_to_node_idx = self.optimizable.batch_indices

        self.v += self.dt[batch_to_node_idx].unsqueeze(1) * f

        dr = self.dt[batch_to_node_idx].unsqueeze(1) * self.v
        normdr = torch.sqrt(self._batched_dot(dr, dr))

        batch_idx_ = normdr > self.maxstep
        dr_update = torch.ones_like(self.dt)
        dr_update[batch_idx_] = self.maxstep / normdr[batch_idx_]
        dr = dr * dr_update[batch_to_node_idx].unsqueeze(1)

        r = self.optimizable.get_positions()
        self.optimizable.set_positions(r + dr)

    def write(self, *, force_write: bool = False) -> None:
        atoms_objects = self.optimizable.get_atoms_list()
        # import pdb;pdb.set_trace()
        for atm, traj, mask in zip(
            atoms_objects, self.trajectories, self.optimizable.update_mask, strict=False
        ):
            if mask or force_write:
                traj.write(atm)
