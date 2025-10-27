# models/equflash/relaxation/optimizers/lbfgs_torch.py
from __future__ import annotations

from collections import deque
from typing import TYPE_CHECKING

import torch
from torch_scatter import scatter

from .base_optimizer import BaseBatchOptimizer

if TYPE_CHECKING:
    from pathlib import Path

    from ..optimizable import OptimizableBatch  # noqa: TID252


class LBFGS(BaseBatchOptimizer):
    """Limited memory BFGS optimizer for batch ML relaxations."""

    def __init__(
        self,
        optimizable_batch: OptimizableBatch,
        maxstep: float = 0.2,
        memory: int = 100,
        damping: float = 1.2,
        alpha: float = 100.0,
        *,
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
        self.maxstep = maxstep
        self.memory = memory
        self.damping = damping
        self.alpha = alpha
        self.H0 = 1.0 / self.alpha

        self.s = deque(maxlen=self.memory)
        self.y = deque(maxlen=self.memory)
        self.rho = deque(maxlen=self.memory)
        self.r0 = None
        self.f0 = None

    def run(self, fmax: float, steps: int) -> tuple[torch.Tensor, torch.Tensor]:
        self.fmax = fmax
        self.steps = steps

        # reset memory
        self.s.clear()
        self.y.clear()
        self.rho.clear()
        self.r0 = self.f0 = None

        # open trajectories
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
            self.step(iteration)
            max_forces = self.optimizable.get_max_forces()
            iteration += 1
            n_traj = (
                n_traj
                + torch.ones_like(n_traj) * self.optimizable.update_mask.detach().cpu()
            )

        # save after converged or all iterations ran
        if iteration > 0 and self.trajectories is not None:
            self.write(force_write=True)

        # teardown & finalize
        self._teardown_cuda_cache()
        self._finalize_trajs()

        # set predicted values to batch
        self._apply_results_to_batch()

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

    def step(self, iteration: int) -> None:
        forces = self.optimizable.get_forces().to(dtype=torch.float64)
        pos = self.optimizable.get_positions().to(dtype=torch.float64)
        if iteration > 0:
            s0 = pos - self.r0
            self.s.append(s0)

            y0 = -(forces - self.f0)
            self.y.append(y0)

            self.rho.append(1.0 / self._batched_dot(y0, s0))

        loopmax = min(self.memory, iteration)
        alpha = forces.new_empty(loopmax, self.optimizable.batch.natoms.shape[0])
        q = -forces
        for i in range(loopmax - 1, -1, -1):
            alpha[i] = self.rho[i] * self._batched_dot(self.s[i], q)  # b
            q -= alpha[i][self.optimizable.batch_indices, ..., None] * self.y[i]

        z = (1.0 / self.alpha) * q  # self.H0 * q
        for i in range(loopmax):
            beta = self.rho[i] * self._batched_dot(self.y[i], z)
            z += self.s[i] * (
                alpha[i][self.optimizable.batch_indices, ..., None]
                - beta[self.optimizable.batch_indices, ..., None]
            )

        p = -z
        dr = self.determine_step(p)
        if torch.abs(dr).max() < 1e-7:
            return
        self.optimizable.set_positions(pos + dr)
        self.r0 = pos
        self.f0 = forces
