# models/equflash/relaxation/optimizers/fire.py
from __future__ import annotations

from typing import TYPE_CHECKING

import torch

from .base_optimizer import BaseBatchOptimizer

if TYPE_CHECKING:
    from collections.abc import Callable
    from pathlib import Path

    from ..optimizable import OptimizableBatch  # noqa: TID252


class FIRE(BaseBatchOptimizer):
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
        super().__init__(
            optimizable_batch,
            save_full_traj=save_full_traj,
            traj_dir=traj_dir,
            traj_names=traj_names,
        )
        self.maxstep = 0.2 if maxstep is None else maxstep
        self.dtmax = torch.tensor(dtmax)
        self.Nmin = nmin
        self.finc = finc
        self.fdec = fdec
        self.astart = astart
        self.fa = fa

        self.v = None
        batch_size = len(self.optimizable.batch.natoms)
        device = self.optimizable.batch.batch.device
        self.dt = dt * torch.ones(batch_size, dtype=torch.float64, device=device)
        self.Nsteps = torch.zeros(batch_size, dtype=torch.int, device=device)
        self.a = a * torch.ones(batch_size, dtype=torch.float64, device=device)

        self.downhill_check = downhill_check
        self.position_reset_callback = position_reset_callback

    def run(self, fmax: float, steps: int) -> tuple[torch.Tensor, bool]:
        self.fmax = fmax
        self.steps = steps

        self._open_trajs()

        iteration = 0
        max_forces = self.optimizable.get_max_forces()
        n_traj = torch.zeros(self.optimizable.batch.batch.max() + 1, dtype=torch.int)
        while iteration < steps and not self.optimizable.converged(
            forces=None, fmax=self.fmax, max_forces=max_forces
        ):
            self.iteration = iteration
            if self._should_write(iteration):
                self.write()
            self.step()
            max_forces = self.optimizable.get_max_forces()
            iteration += 1

            n_traj = (
                n_traj
                + torch.ones_like(n_traj) * self.optimizable.update_mask.detach().cpu()
            )

        if iteration > 0 and self.trajectories is not None:
            self.write(force_write=True)

        self._teardown_cuda_cache()
        self._finalize_trajs()
        self._apply_results_to_batch()

        return n_traj, self.optimizable.converged(
            forces=None, fmax=self.fmax, max_forces=max_forces
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

            vf = self._batched_dot(f, self.v)

            batch_idx_ = (not is_uphill) * torch.ones(
                self.dt.shape[0], dtype=torch.bool, device=self.dt.device
            )
            batch_idx_ = batch_idx_ * (vf > 0.0)

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
