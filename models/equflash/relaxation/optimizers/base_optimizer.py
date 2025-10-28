"""Base optimizer for batched ML relaxations."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import ase
import ase.io.trajectory
from torch_scatter import scatter

if TYPE_CHECKING:
    from ase.io.trajectory import TrajectoryWriter
    from torch import Tensor

    from models.equflash.relaxation.optimizable import OptimizableBatch


class BaseBatchOptimizer:
    """Base class for batch optimizers."""

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
        self.trajectories: list[TrajectoryWriter] | None = None

        if self.traj_dir and not self.traj_names:
            raise ValueError(
                "Trajectory names should be specified to save trajectories"
            )

        self.fmax: float | None = None
        self.steps: int | None = None
        self.iteration: int | None = None

    def _open_trajs(self) -> None:
        """Open trajectory files for writing."""
        if not self.traj_dir:
            return
        self.traj_dir.mkdir(exist_ok=True, parents=True)
        self.trajectories = [
            ase.io.Trajectory(self.traj_dir / f"{name}.traj_tmp", mode="w")
            for name in self.traj_names
        ]

    def _finalize_trajs(self) -> None:
        """Finalize trajectory files by renaming from tmp."""
        if self.trajectories is not None:
            for traj in self.trajectories:
                traj.close()
        if self.traj_dir is not None:
            for name in self.traj_names:
                tmp_path = self.traj_dir / f"{name}.traj_tmp"
                Path(tmp_path).rename(Path(tmp_path).with_suffix(".traj"))

    def write(self, *, force_write: bool = False) -> None:
        """Write current state to trajectory file."""
        if self.trajectories is None:
            return
        atoms_objects = self.optimizable.get_atoms_list()
        mask_iter = self.optimizable.update_mask

        for atm, traj, mask in zip(
            atoms_objects, self.trajectories, mask_iter, strict=True
        ):
            if mask or force_write:
                traj.write(atm)

    def _batched_dot(self, x: Tensor, y: Tensor) -> Tensor:
        """Compute batched dot product."""
        return scatter(
            (x * y).sum(dim=-1), self.optimizable.batch_indices, reduce="sum"
        )
