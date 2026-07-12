"""Standard ASE implementation of the shared phonon/kappa adapter hooks."""

from __future__ import annotations

import os
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, cast

import numpy as np
from ase.constraints import FixSymmetry

from matbench_discovery.ase_relax import resolve_cell_filter, resolve_optimizer

if TYPE_CHECKING:
    from ase import Atoms
    from ase.calculators.calculator import Calculator
    from ase.filters import Filter
    from phono3py.api_phono3py import Phono3py

    from matbench_discovery.phonons.pipeline import KappaSettings

NO_TILT_MASK = [True, True, True, False, False, False]


@dataclass(frozen=True)
class RelaxationResult:
    """A relaxed structure and its symmetry/convergence metadata."""

    atoms: Atoms
    initial_spg_num: int
    final_spg_num: int
    max_stress: np.ndarray | None
    reached_max_steps: bool
    n_steps: int
    redirected_to_symmetry: bool = False


def space_group_number(atoms: Atoms, symprec: float) -> int:
    """Detect an ASE structure's space-group number with Moyo."""
    from moyopy import MoyoDataset
    from moyopy.interface import MoyoAdapter

    return int(MoyoDataset(MoyoAdapter.from_atoms(atoms), symprec=symprec).number)


class StandardKappaAdapter:
    """Run the common ASE relaxation and serial phono3py force evaluation."""

    name = "standard"

    def prepare_calculator(
        self,
        calculator: Calculator,
        settings: KappaSettings,  # noqa: ARG002
    ) -> Calculator:
        """Return a standard calculator unchanged."""
        return calculator

    def relax(
        self,
        atoms: Atoms,
        calculator: Calculator,
        settings: KappaSettings,
        *,
        log_file: str | None,
    ) -> RelaxationResult:
        """Relax atoms once with the configured ASE optimizer and filter."""
        initial_spg_num = space_group_number(atoms, settings.symprec)
        if settings.relaxation_mode == "none" or settings.max_steps == 0:
            return RelaxationResult(
                atoms=atoms,
                initial_spg_num=initial_spg_num,
                final_spg_num=initial_spg_num,
                max_stress=None,
                reached_max_steps=False,
                n_steps=0,
            )
        initial_info = dict(atoms.info)
        if settings.relaxation_mode == "two-stage":
            return self._relax_two_stage(
                atoms,
                calculator,
                settings,
                initial_spg_num=initial_spg_num,
                initial_info=initial_info,
                log_file=log_file,
            )

        atoms.calc = calculator
        filter_class = resolve_cell_filter(settings.ase_filter)
        if settings.enforce_relax_symm:
            atoms.set_constraint(FixSymmetry(atoms, symprec=settings.relax_symprec))
        if filter_class is None:
            relax_target: Atoms | Filter = atoms
        elif settings.enforce_relax_symm:
            relax_target = filter_class(atoms, mask=NO_TILT_MASK)
        else:
            relax_target = filter_class(atoms)

        optimizer_class = resolve_optimizer(settings.ase_optimizer)
        if log_file and (log_parent := os.path.dirname(log_file)):
            os.makedirs(log_parent, exist_ok=True)
        optimizer = optimizer_class(relax_target, logfile=log_file)  # ty: ignore[invalid-argument-type]
        optimizer.run(fmax=settings.force_max, steps=settings.max_steps)
        n_steps = int(optimizer.get_number_of_steps())
        max_stress = atoms.get_stress().reshape((2, 3), order="C").max(axis=1)

        atoms.calc = None
        atoms.constraints = None
        atoms.info = initial_info | atoms.info
        final_spg_num = space_group_number(atoms, settings.symprec)
        return RelaxationResult(
            atoms=atoms,
            initial_spg_num=initial_spg_num,
            final_spg_num=final_spg_num,
            max_stress=np.asarray(max_stress),
            reached_max_steps=n_steps >= settings.max_steps,
            n_steps=n_steps,
        )

    def _relax_two_stage(
        self,
        atoms: Atoms,
        calculator: Calculator,
        settings: KappaSettings,
        *,
        initial_spg_num: int,
        initial_info: dict[str, Any],
        log_file: str | None,
    ) -> RelaxationResult:
        """Run constrained relaxation, then test unconstrained symmetry stability."""
        filter_class = resolve_cell_filter(settings.ase_filter)
        optimizer_class = resolve_optimizer(settings.ase_optimizer)
        if log_file and (log_parent := os.path.dirname(log_file)):
            os.makedirs(log_parent, exist_ok=True)
        log_stem, log_suffix = os.path.splitext(log_file or "")

        atoms.calc = calculator
        if settings.enforce_relax_symm:
            atoms.set_constraint(FixSymmetry(atoms, symprec=settings.relax_symprec))
        stage1_target: Atoms | Filter = (
            atoms if filter_class is None else filter_class(atoms, mask=NO_TILT_MASK)
        )
        stage1_optimizer = optimizer_class(
            cast("Atoms", stage1_target),
            logfile=f"{log_stem}-stage1{log_suffix}" if log_file else None,
        )
        stage1_optimizer.run(fmax=settings.force_max, steps=settings.max_steps)
        stage1_steps = int(stage1_optimizer.get_number_of_steps())
        stage1_spg_num = space_group_number(atoms, settings.symprec)
        stage1_stress = atoms.get_stress().reshape((2, 3), order="C").max(axis=1)
        stage1_atoms = atoms.copy()

        atoms.constraints = None
        stage2_target: Atoms | Filter = (
            atoms if filter_class is None else filter_class(atoms, mask=NO_TILT_MASK)
        )
        stage2_optimizer = optimizer_class(
            cast("Atoms", stage2_target),
            logfile=f"{log_stem}-stage2{log_suffix}" if log_file else None,
        )
        stage2_optimizer.run(fmax=settings.force_max, steps=settings.max_steps)
        stage2_steps = int(stage2_optimizer.get_number_of_steps())
        stage2_spg_num = space_group_number(atoms, settings.symprec)
        stage2_stress = atoms.get_stress().reshape((2, 3), order="C").max(axis=1)

        redirected = settings.enforce_relax_symm and stage1_spg_num != stage2_spg_num
        relaxed_atoms = stage1_atoms if redirected else atoms
        relaxed_atoms.calc = None
        relaxed_atoms.constraints = None
        relaxed_atoms.info = initial_info | relaxed_atoms.info
        return RelaxationResult(
            atoms=relaxed_atoms,
            initial_spg_num=initial_spg_num,
            final_spg_num=stage1_spg_num if redirected else stage2_spg_num,
            max_stress=np.asarray(stage1_stress if redirected else stage2_stress),
            reached_max_steps=(
                stage1_steps >= settings.max_steps or stage2_steps >= settings.max_steps
            ),
            n_steps=stage1_steps + stage2_steps,
            redirected_to_symmetry=redirected,
        )

    def init_phono3py(self, atoms: Atoms, settings: KappaSettings) -> Phono3py:
        """Initialize phono3py from canonical PhononDB metadata keys."""
        from matbench_discovery.phonons import thermal_conductivity as ltc

        q_point_mesh = atoms.info.get("q_point_mesh", atoms.info.get("q_mesh"))
        if q_point_mesh is None:
            raise ValueError("Structure has no q_point_mesh or q_mesh metadata")
        return ltc.init_phono3py(
            atoms,
            fc2_supercell=np.asarray(atoms.info["fc2_supercell"]),
            fc3_supercell=np.asarray(atoms.info["fc3_supercell"]),
            q_point_mesh=tuple(q_point_mesh),
            displacement_distance=settings.displacement_distance,
            symprec=settings.symprec,
            is_plusminus=settings.is_plusminus,
        )

    def calculate_fc2(
        self,
        phono3py: Phono3py,
        calculator: Calculator,
        settings: KappaSettings,  # noqa: ARG002
        *,
        progress: dict[str, Any] | None = None,
    ) -> tuple[Phono3py, np.ndarray, np.ndarray]:
        """Calculate FC2 and frequencies with serial ASE force calls."""
        from matbench_discovery.phonons import thermal_conductivity as ltc

        return ltc.get_fc2_and_freqs(
            phono3py, calculator, pbar_kwargs=progress or {"disable": True}
        )

    def calculate_fc3(
        self,
        phono3py: Phono3py,
        calculator: Calculator,
        settings: KappaSettings,  # noqa: ARG002
        *,
        progress: dict[str, Any] | None = None,
        max_evaluations: int | None = None,
    ) -> np.ndarray:
        """Calculate FC3 with serial ASE force calls."""
        from matbench_discovery.phonons import thermal_conductivity as ltc

        return ltc.calculate_fc3_set(
            phono3py,
            calculator,
            pbar_kwargs=progress or {"disable": True},
            max_evaluations=max_evaluations,
        )


class FairchemKappaAdapter(StandardKappaAdapter):
    """Disable training-only force scaling before standard ASE evaluation."""

    name = "fairchem"

    def prepare_calculator(
        self,
        calculator: Calculator,
        settings: KappaSettings,  # noqa: ARG002
    ) -> Calculator:
        """Remove a checkpoint trainer's scaler when the API exposes it."""
        if hasattr(trainer := getattr(calculator, "trainer", None), "scaler"):
            trainer.scaler = None
        return calculator
