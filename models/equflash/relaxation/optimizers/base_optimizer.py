# models/equflash/relaxation/optimizers/base_optimizer.py
from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import ase
import torch
from torch_scatter import scatter

if TYPE_CHECKING:
    from ..optimizable import OptimizableBatch  # noqa: TID252


class BaseBatchOptimizer:
    """base class for optimizer."""

    def __init__(
        self,
        optimizable_batch: OptimizableBatch,
        *,
        save_full_traj: bool = True,
        traj_dir: Path | None = None,
        traj_names: list[str] | None = None,
    ) -> None:
        self.optimizable = optimizable_batch
        self.save_full = save_full_traj
        self.traj_dir = traj_dir
        self.traj_names = traj_names or []
        self.trajectories: list[ase.io.trajectory.Trajectory] | None = None

        if self.traj_dir and (not traj_dir or not self.traj_names):
            raise ValueError(
                "Trajectory names should be specified to save trajectories"
            )

        self.fmax: float | None = None
        self.steps: int | None = None
        self.iteration: int | None = None

    # ---------- Trajectory helpers ----------
    def _open_trajs(self) -> None:
        if not self.traj_dir:
            self.trajectories = []
            return
        self.traj_dir.mkdir(exist_ok=True, parents=True)
        self.trajectories = [
            ase.io.Trajectory(self.traj_dir / f"{name}.traj_tmp", mode="w")
            for name in self.traj_names
        ]

    def _close_trajs(self) -> None:
        if self.trajectories is not None:
            for traj in self.trajectories:
                traj.close()

    def _finalize_trajs(self) -> None:
        self._close_trajs()
        if self.traj_dir is not None:
            for name in self.traj_names:
                tmp_path = self.traj_dir / f"{name}.traj_tmp"
                Path(tmp_path).rename(Path(tmp_path).with_suffix(".traj"))

    def _should_write(self, iteration: int) -> bool:
        return self.trajectories is not None and (self.save_full or iteration == 0)

    def write(self, *, force_write: bool = False) -> None:
        if self.trajectories is None:
            return
        atoms_objects = self.optimizable.get_atoms_list()

        mask_iter = self.optimizable.update_mask

        for atm, traj, mask in zip(
            atoms_objects, self.trajectories, mask_iter, strict=False
        ):
            if mask or force_write:
                traj.write(atm)

    # ---------- Utilities ----------
    def _batched_dot(self, x: torch.Tensor, y: torch.Tensor) -> torch.Tensor:
        return scatter(
            (x * y).sum(dim=-1), self.optimizable.batch_indices, reduce="sum"
        )

    def _apply_results_to_batch(self) -> None:
        for name, value in self.optimizable.results.items():
            setattr(self.optimizable.batch, name, value)

    def _teardown_cuda_cache(self) -> None:
        torch.cuda.empty_cache()
