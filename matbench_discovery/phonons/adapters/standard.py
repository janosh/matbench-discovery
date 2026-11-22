"""Standard implementation of the shared phonon/kappa adapter hooks."""

from __future__ import annotations

import os
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

import numpy as np
from kemio import RelaxResult, Structure, symmetry_dataset
from kemio.model import matrix_to_voigt

from matbench_discovery.relaxation import relax_structures

if TYPE_CHECKING:
    from kemio.model import Model as ForceModel

    from matbench_discovery.phonons.pipeline import KappaSettings
    from matbench_discovery.phonons.thermal_conductivity import KappaCalculation

NO_TILT_MASK = (True, True, True, False, False, False)


@dataclass(frozen=True)
class RelaxationResult:
    """A relaxed structure and its symmetry/convergence metadata."""

    atoms: Structure
    initial_spg_num: int
    final_spg_num: int
    max_stress: np.ndarray | None
    reached_max_steps: bool
    n_steps: int
    redirected_to_symmetry: bool = False


def space_group_number(structure: Structure, symprec: float) -> int:
    """Detect a structure's space-group number with moyo."""
    return symmetry_dataset(structure, symprec=symprec).number


def supercell_batch_size(settings: KappaSettings, n_atoms: int) -> int:
    """Supercells per model call: ``batch_size``, or as many ``n_atoms``-atom cells as
    fit in ``max_atoms_per_batch`` (at least one: every kappa-103 FC2 supercell has 300
    or 512 atoms, beyond EquFlash's 128).
    """
    if settings.max_atoms_per_batch is None:
        return settings.batch_size
    return max(settings.max_atoms_per_batch // n_atoms, 1)


def _max_stress(stress: np.ndarray | None) -> np.ndarray:
    """Largest absolute normal and shear Voigt stress components (eV/Å³) of a 3x3
    stress tensor (xx, yy, zz | yz, xz, xy).

    Absolute values so a large compressive (negative) stress isn't reported as ~0.
    """
    if stress is None:
        raise RuntimeError("relaxation reported no stress although stress=True")
    return np.abs(matrix_to_voigt(stress)).reshape((2, 3), order="C").max(axis=1)


def _write_relax_log(log_file: str | None, stage: str, result: RelaxResult) -> None:
    """Append one line summarizing a relaxation stage."""
    if not log_file:
        return
    if log_parent := os.path.dirname(log_file):
        os.makedirs(log_parent, exist_ok=True)
    with open(log_file, mode="a", encoding="utf-8") as file:
        file.write(
            f"{stage}: n_steps={result.n_steps} converged={result.converged} "
            f"energy={result.energy:.10f} eV\n"
        )


class StandardKappaAdapter:
    """Run the common relaxation and batched displacement force evaluation."""

    name = "standard"

    def prepare_model(
        self,
        model: ForceModel,
        settings: KappaSettings,  # noqa: ARG002
    ) -> ForceModel:
        """Return a standard model unchanged."""
        return model

    def _relax_stage(
        self,
        structure: Structure,
        model: ForceModel,
        settings: KappaSettings,
        *,
        fix_symmetry: bool,
        no_tilt: bool,
    ) -> RelaxResult:
        """One relaxation with the configured optimizer, filter and symmetry (always
        reporting the final, symmetrized stress for ``max_stress``).
        """
        (result,) = relax_structures(
            [structure],
            model,
            optimizer=settings.ase_optimizer,
            cell_filter=settings.ase_filter,
            max_force=settings.force_max,
            max_steps=settings.max_steps,
            fix_symmetry=settings.relax_symprec if fix_symmetry else None,
            cell_mask=NO_TILT_MASK if no_tilt else None,
            stress=True,
        )
        if result.error is not None:
            raise RuntimeError(f"relaxation failed: {result.error}")
        return result

    def relax(
        self,
        atoms: Structure,
        model: ForceModel,
        settings: KappaSettings,
        *,
        log_file: str | None,
    ) -> RelaxationResult:
        """Relax a structure once with the configured optimizer and filter."""
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
        if settings.relaxation_mode == "two-stage":
            return self._relax_two_stage(
                atoms,
                model,
                settings,
                initial_spg_num=initial_spg_num,
                log_file=log_file,
            )

        result = self._relax_stage(
            atoms,
            model,
            settings,
            fix_symmetry=settings.enforce_relax_symm,
            no_tilt=settings.enforce_relax_symm,
        )
        _write_relax_log(log_file, "relax", result)
        relaxed = result.structure
        n_steps = result.n_steps
        return RelaxationResult(
            atoms=relaxed,
            initial_spg_num=initial_spg_num,
            final_spg_num=space_group_number(relaxed, settings.symprec),
            max_stress=np.asarray(_max_stress(result.stress)),
            reached_max_steps=n_steps >= settings.max_steps,
            n_steps=n_steps,
        )

    def _relax_two_stage(
        self,
        atoms: Structure,
        model: ForceModel,
        settings: KappaSettings,
        *,
        initial_spg_num: int,
        log_file: str | None,
    ) -> RelaxationResult:
        """Run constrained relaxation, then test unconstrained symmetry stability."""
        stage1 = self._relax_stage(
            atoms,
            model,
            settings,
            fix_symmetry=settings.enforce_relax_symm,
            no_tilt=True,
        )
        _write_relax_log(log_file, "stage1", stage1)
        stage1_structure = stage1.structure
        stage1_spg_num = space_group_number(stage1_structure, settings.symprec)
        stage2 = self._relax_stage(
            stage1_structure, model, settings, fix_symmetry=False, no_tilt=True
        )
        _write_relax_log(log_file, "stage2", stage2)
        stage2_spg_num = space_group_number(stage2.structure, settings.symprec)

        redirected = settings.enforce_relax_symm and stage1_spg_num != stage2_spg_num
        chosen = stage1 if redirected else stage2
        n_steps = (stage1.n_steps, stage2.n_steps)
        return RelaxationResult(
            atoms=chosen.structure,
            initial_spg_num=initial_spg_num,
            final_spg_num=stage1_spg_num if redirected else stage2_spg_num,
            max_stress=np.asarray(_max_stress(chosen.stress)),
            reached_max_steps=max(n_steps) >= settings.max_steps,
            n_steps=sum(n_steps),
            redirected_to_symmetry=redirected,
        )

    def init_phonons(
        self, atoms: Structure, settings: KappaSettings
    ) -> KappaCalculation:
        """Set up fc2/fc3 displacements from canonical PhononDB metadata keys."""
        from matbench_discovery.phonons import thermal_conductivity as ltc

        return ltc.init_phonons(
            atoms,
            fc2_supercell=np.asarray(atoms.info["fc2_supercell"]),
            fc3_supercell=np.asarray(atoms.info["fc3_supercell"]),
            q_point_mesh=tuple(np.asarray(atoms.info["q_point_mesh"]).ravel()),
            displacement_distance=settings.displacement_distance,
            symprec=settings.symprec,
            is_plusminus=settings.is_plusminus,
        )

    def calculate_fc2(
        self,
        calc: KappaCalculation,
        model: ForceModel,
        settings: KappaSettings,
    ) -> np.ndarray:
        """Fit symmetrized FC2 to the model's displacement forces (kept as
        ``fc2_forces``) in batches of ``supercell_batch_size`` supercells and return
        the frequencies (n_grid, n_bands) in THz on each regular grid point once
        (Gamma first), the uniform-weight counterpart of the weighted irreducible mesh.
        """
        batch_size = supercell_batch_size(settings, len(calc.phonons.phonon_supercell))
        calc.phonons.produce_fc2(
            model, symmetrize=True, batch={"structures": batch_size}
        )
        mesh = calc.phonons.harmonic.mesh(list(calc.mesh))
        return mesh.frequencies[mesh.ir_of_grid]

    def harmonic_phonons(
        self,
        calc: KappaCalculation,
        settings: KappaSettings,  # noqa: ARG002
    ) -> dict[str, Any]:
        """Derive band path, eigenvectors, DOS and thermal properties from FC2.

        Runs after calculate_fc2 for every material, so it must not depend on FC3 or
        the conductivity solver.
        """
        from matbench_discovery.phonons.harmonic import harmonic_phonon_data

        return harmonic_phonon_data(calc.phonons.harmonic, mesh=calc.mesh)

    def calculate_fc3(
        self,
        calc: KappaCalculation,
        model: ForceModel,
        settings: KappaSettings,
        *,
        max_evaluations: int | None = None,
    ) -> None:
        """Fit FC3 (with sum rules and permutation symmetry) to the model's
        displacement forces (kept as ``fc3_forces``), batched as ``calculate_fc2``.
        ``max_evaluations`` caps the evaluated supercells, the rest taking zero forces
        (dry runs).
        """
        batch_size = supercell_batch_size(settings, len(calc.phonons.supercell))
        calc.phonons.produce_fc3(
            model,
            symmetrize=True,
            batch={"structures": batch_size},
            max_evaluations=max_evaluations,
        )


class FairchemKappaAdapter(StandardKappaAdapter):
    """Disable training-only force scaling before standard evaluation."""

    name = "fairchem"

    def prepare_model(
        self,
        model: ForceModel,
        settings: KappaSettings,  # noqa: ARG002
    ) -> ForceModel:
        """Remove a checkpoint trainer's scaler when the wrapped calculator has one."""
        calculator = getattr(model, "calculator", None)
        if hasattr(trainer := getattr(calculator, "trainer", None), "scaler"):
            trainer.scaler = None
        return model
