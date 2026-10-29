# /// script
# requires-python = ">=3.12"
# dependencies = [
#   "deepmd-kit @ git+https://github.com/deepmodeling/deepmd-kit@a78bdfe8d5dc494a9d9a397e4a9e103a84205eaa",
#   "torch>=2.11",
#   "ase>=3.26",
#   "ase-db-backends>=0.11",  # "aselmdb" ASE database type of the PFT release
#   "numpy>=2",
# ]
# ///
r"""Phonon fine-tuning (PFT) of DPA-4.1-Max-OAM to DPA-4.1-Max-OAM-PFT.

This is a SIMPLIFIED, readable reproduction of the training recipe of the submitted
checkpoint. The original trainer is internal. This file documents the recipe and can
be followed to reproduce it; it is not the exact code that produced the checkpoint, and
it has not been run end-to-end in this form.

Recipe (PFT paper, arXiv:2601.07742, Alg. 1, plus a soft-mode anchor term).
One cycle =
  1. phonon step: a batch of PFT structures; for each, ONE random Hessian column j
     (atom, xyz). The model column is the central finite difference of the forces,
     H[:, j] = -(F(x + d e_j) - F(x - d e_j)) / (2 d), d = 0.01 A, with the neighbor
     topology fixed at x. Loss per structure:
       wH * mean|H_model[:, j] - H_dft[:, j]| + wF * mean_a |F_a|
         + wS * mean|sigma| (+ wE * |E - E_dft| / N);
     the PFT structures are equilibrium cells, so the force and stress targets are
     zero; force, stress and energy are the averages of the +d and -d evaluations.
     Batch weights: Hessian and force terms N_i / sum N, stress and energy 1/B.
     In the same step, soft-mode anchors: 64 MPtrj structures where the BASE model
     has very soft or unstable modes. For each, the base model's own Hessian H0 is
     computed once before training; its six lowest non-translational eigenvectors v_k
     (including the unstable ones, i.e. negative eigenvalues) are the protected
     directions. Per phonon step 8 of the 64 anchors are drawn without replacement
     (min(8, 64) = 8); each drawn anchor contributes exactly ONE term, along one
     randomly drawn v_k of its six, with weight w_anchor / 8 = 12.5, so the eight
     terms sum to w_anchor. Loss of one term:
       w * mean|H(theta) v_k - H0 v_k| / (mean|H0 v_k| + 1e-3),
     with H(theta) v_k again a central difference of forces along v_k. No DFT labels
     are used for the anchors.
  2. co-training: K steps on random MPtrj frames,
       cwE * |E - E*| / N + cwF * |F - F*| + cwS * |sigma - sigma*|
     (energy per atom per frame; per-atom force-vector norm summed and divided by the
     batch's total atom count; mean absolute stress component per frame).

Optimizer: AdamW (betas 0.9/0.999, weight decay 1e-3 on parameters with ndim >= 2),
constant lr after a linear warm-up, gradient clipping by global norm, EMA of the
weights (decay min(0.999, (1 + step) / (10 + step))). The model always runs in eval
mode (float32, no bf16 autocast, no TF32) while its parameters are trained. The
submitted checkpoint is the EMA weights after 3 epochs (fixed endpoint).

Hyper-parameters below are the ones of the submitted run. Data:
  * PFT: the public "Preprocessed data for Phonon Fine-tuning (PFT)" release (figshare
    article 31049869; registered in matbench-discovery data/datasets.yml as
    "MDR-MP PBE ω_q"): .aselmdb shards, row data["hessian"] of shape (3N, 3N) in
    eV/A^2. We removed every row whose reduced formula occurs among the 103 PhononDB
    kappa benchmark structures (train 8,510 -> 8,350; val 448 -> 442) and did not use
    rows with more than 11,000 directed neighbor pairs within 6 A (8,350 -> 7,984
    used). Pass the kept row keys '<shard file>:<row id>' as a text file
    (--train-keys).
  * MPtrj: training frames with energy, forces and stress, read here with ase.io
    (e.g. an extxyz export). Frames with more than 300 atoms are not used for
    co-training.
  * Anchors (64): MPtrj training frames selected by a rule that was fixed before any
    run of this recipe and does not use the test sets: candidates are MPtrj training
    frames with 2 to 8 atoms, one frame per reduced formula (smallest cell, then
    lowest energy per atom); the harmonic Gamma-point frequencies of the BASE model
    are computed on each candidate (mass-weighted Hessian by central finite
    differences, d = 0.01 A; rigid translations projected out, the three acoustic
    branches are the three lowest); a candidate qualifies if its lowest
    non-acoustic branch is below 0.5 THz or any frequency is imaginary
    (< -0.01 THz); every reduced formula that occurs among the 103 PhononDB kappa
    structures (70 formulas, asserted) is excluded; the qualifying frames are
    sorted by their lowest non-acoustic frequency, ties by formula, and the first
    64 are kept (5,691 of 25,020 candidate formulas qualified). The 64 frames have
    3 to 8 atoms and all of them have imaginary Gamma modes in the base model. Each
    is used as a 2x2x2 supercell with six candidate directions (the six lowest
    non-translational eigenvectors of the base-model Hessian); each phonon step
    draws 8 of the 64 anchors, one direction each, eight anchor terms per step,
    total weight 100. The Hessian weight 2000 (1000 in an otherwise identical
    earlier run) was chosen after looking at the kappa-103 trend of earlier runs;
    this is disclosed in the model YAML.

Differences from the internal trainer:
  removed: multi-GPU data parallelism (the run used 4 GPUs; gradients were summed
           with global normalizers, which equals the single-GPU update implemented
           here); out-of-memory recovery and retries; resume state; intermediate
           checkpoints and their automatic evaluation; periodic validation; logging
           and monitoring; precision self-checks; internal data formats and paths;
           batched evaluation of co-training frames.
           The internal trainer also skipped co-training frames with more than 8,000
           directed neighbor pairs (about 2.2% of the sampled frames in the submitted
           run); this cap is not reproduced here.
  kept:    the loss terms and weights, the fixed-topology finite-difference Hessian
           column, the anchor construction (base-model Hessian, translations
           projected out, six lowest eigen-directions including negative ones, one
           term per drawn anchor per step), the optimizer / schedule / clipping /
           EMA, the seeded sampling scheme (per-epoch permutation and Hessian column
           seeds from seed + 1000 * epoch, anchors from seed + 99991 * (cycle + 1),
           co-training frames from seed + 7 * (cycle * K + k) + 1; the exact MPtrj
           frames depend on the order of the MPtrj copy) and the epoch-3 EMA
           endpoint.
  differs: gradients come from autograd through both displaced force evaluations in
           one graph; the internal trainer obtained the same gradient with two
           separate backward passes, which only saves memory.

Usage:
  python train_dpa_4_1_max_oam_pft.py --base DPA-4.1-Max-OAM.pt \
      --pft-train PFT/train-aselmdb --train-keys train_keys.txt \
      --mptrj mptrj_train.extxyz --anchors anchors.extxyz \
      --out dpa-4.1-max-oam-pft.pt
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Any

import ase
import ase.db
import ase.io
import numpy as np
import torch

# -------------------------------------------------- hyper-parameters of the run
CFG = dict(
    epochs=3,  # 499 cycles per epoch with 7,984 rows and batch 16 -> 1,497 cycles
    lr=5e-6,
    warmup_steps=50,  # optimizer steps (1 phonon + 4 co-training steps per cycle)
    weight_decay=1e-3,
    clip_norm=1000.0,
    ema_decay=0.999,
    batch=16,  # PFT structures per phonon step
    fd_delta=0.01,  # A, central-difference step for Hessians and anchor directions
    w_energy=0.0,
    w_force=20.0,
    w_stress=5.0,
    w_hessian=2000.0,
    cotrain_steps=4,
    cotrain_batch=16,
    cotrain_max_atoms=300,
    cw_energy=500.0,
    cw_force=200.0,
    cw_stress=50.0,
    anchor_weight=100.0,
    # upper bound on anchors sampled per phonon step; the pool has 64 anchors, so 8
    # are drawn (without replacement), each contributing one term with weight
    # anchor_weight / 8 = 12.5
    anchors_per_step=8,
    anchor_dirs=6,  # candidate directions per anchor (lowest non-translational)
    anchor_supercell=2,
    seed=20261022,
)


# -------------------------------------------------- model access (deepmd-kit PT)
def load_model(path: str, device: torch.device) -> tuple[Any, torch.nn.Module]:
    """Load a deepmd .pt checkpoint; return (wrapper, model) in eval mode."""
    from deepmd.pt.model.model import get_model
    from deepmd.pt.train.wrapper import ModelWrapper

    state = torch.load(path, map_location="cpu", weights_only=True)
    state = state.get("model", state)
    model = get_model(state["_extra_state"]["model_params"]).to(device)
    wrapper = ModelWrapper(model)
    wrapper.load_state_dict(state)
    model.eval()  # deterministic inference path; parameters still receive gradients
    model.atomic_model.descriptor.use_amp = False
    model.enable_tf32 = False
    torch.backends.cuda.matmul.allow_tf32 = False
    torch.backends.cudnn.allow_tf32 = False
    return wrapper, model


class Graph:
    """Neighbor topology of one periodic structure, built once and then kept fixed.

    Edge vectors are rebuilt from the positions as x[src] - x[dst] + shift, so
    autograd reaches the positions (forces) and, with create_graph=True, the
    parameters through the forces. Fixed topology is exact for derivatives at x0.
    """

    def __init__(
        self,
        model: torch.nn.Module,
        pos0: torch.Tensor,
        types: torch.Tensor,
        cell: torch.Tensor,
    ) -> None:
        """Build the neighbor list of ``pos0`` and freeze its topology."""
        n = len(types)
        with torch.no_grad():
            nl = model.build_neighbor_list(
                pos0.view(1, n, 3), types.view(1, n), cell.reshape(1, 9)
            )
        self.nl = nl
        self.src = nl.edge_scatter_index[0]
        self.dst = nl.edge_scatter_index[1]
        self.mask = nl.edge_mask[:, None]
        ev0 = nl.edge_vec.detach().to(pos0.dtype)
        self.shift = torch.where(
            self.mask, ev0 - (pos0[self.src] - pos0[self.dst]), torch.zeros_like(ev0)
        )
        self.ev_pad = torch.where(self.mask, torch.zeros_like(ev0), ev0)
        self.volume = torch.linalg.det(cell.reshape(3, 3)).abs()
        self.charge_spin = model.convert_charge_spin(
            None, nf=1, dtype=pos0.dtype, device=pos0.device
        )


def energy_forces_stress(
    model: torch.nn.Module,
    g: Graph,
    pos: torch.Tensor,
    *,
    create_graph: bool,
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    """Return energy (scalar), forces (N, 3) and stress = virial / volume (3, 3)."""
    if not pos.requires_grad:
        pos = pos.detach().clone().requires_grad_(requires_grad=True)
    edge_vec = torch.where(g.mask, pos[g.src] - pos[g.dst] + g.shift, g.ev_pad)
    nl = g.nl
    out = model.core_compute(
        nl.coord.detach(),
        nl.atype,
        nl.edge_index,
        edge_vec,
        nl.edge_scatter_index,
        nl.edge_mask,
        fparam=None,
        aparam=None,
        charge_spin=g.charge_spin,
        atomic_output_only=True,
    )
    energy = out["energy"].to(pos.dtype).sum()
    g_pos, g_ev = torch.autograd.grad(
        energy, [pos, edge_vec], create_graph=create_graph
    )
    g_ev = torch.where(g.mask, g_ev, torch.zeros_like(g_ev))
    virial = -(g_ev[:, :, None] * edge_vec[:, None, :]).sum(0)  # deepmd convention
    return energy, -g_pos, virial / g.volume


def hessian_times(
    model: torch.nn.Module,
    g: Graph,
    pos0: torch.Tensor,
    direction: torch.Tensor,
    delta: float,
    *,
    create_graph: bool,
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
    """Return H @ direction by a central difference of forces, plus mean E/F/stress."""
    e_p, f_p, s_p = energy_forces_stress(
        model, g, pos0 + delta * direction, create_graph=create_graph
    )
    e_m, f_m, s_m = energy_forces_stress(
        model, g, pos0 - delta * direction, create_graph=create_graph
    )
    hv = -(f_p - f_m).reshape(-1) / (2 * delta)
    return hv, 0.5 * (e_p + e_m), 0.5 * (f_p + f_m), 0.5 * (s_p + s_m)


def to_tensors(
    model: torch.nn.Module, atoms: ase.Atoms, device: torch.device
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    """Return positions, integer types and cell of an ASE Atoms as tensors."""
    tmap = {t: i for i, t in enumerate(model.get_type_map())}
    pos = torch.tensor(atoms.get_positions(), dtype=torch.float64, device=device)
    cell = torch.tensor(
        np.asarray(atoms.get_cell()), dtype=torch.float64, device=device
    )
    types = torch.tensor(
        [tmap[s] for s in atoms.get_chemical_symbols()],
        dtype=torch.int64,
        device=device,
    )
    return pos, types, cell


# -------------------------------------------------- data
class PFTRows:
    """Rows of the PFT .aselmdb shards, addressed as '<shard file>:<row id>'."""

    def __init__(self, root: str) -> None:
        """Remember the shard directory; shards are opened lazily."""
        self.root, self.dbs = root, {}

    def get(self, key: str) -> tuple[ase.Atoms, float, np.ndarray]:
        """Return (Atoms, DFT energy, DFT Hessian) of one '<shard>:<id>' key."""
        shard, rid = key.split(":")
        if shard not in self.dbs:
            self.dbs[shard] = ase.db.connect(
                f"{self.root}/{shard}",
                type="aselmdb",
                readonly=True,
                use_lock_file=False,
            )
        row = self.dbs[shard].get(id=int(rid))
        return (
            row.toatoms(),
            float(row.energy),
            np.asarray(row.data["hessian"], dtype=np.float64),
        )


def build_anchors(
    model: torch.nn.Module, frames: list[ase.Atoms], device: torch.device, cfg: dict
) -> list[dict]:
    """Base-model Hessian targets along the softest non-translational directions."""
    anchors = []
    for frame in frames:
        atoms = frame.repeat(cfg["anchor_supercell"])
        pos, types, cell = to_tensors(model, atoms, device)
        g = Graph(model, pos, types, cell)
        n3 = 3 * len(atoms)
        hess = torch.zeros(n3, n3, dtype=torch.float64, device=device)
        # full base-model Hessian, column by column (supercells of 24 to 64 atoms)
        for j in range(n3):
            e_j = torch.zeros(n3, dtype=torch.float64, device=device)
            e_j[j] = 1.0
            hess[:, j] = hessian_times(
                model, g, pos, e_j.view(-1, 3), cfg["fd_delta"], create_graph=False
            )[0].detach()
        hess = 0.5 * (hess + hess.T)
        trans = torch.zeros(n3, 3, dtype=torch.float64, device=device)
        for c in range(3):  # the 3 rigid translations
            trans[c::3, c] = 1.0 / math.sqrt(len(atoms))
        proj = torch.eye(n3, dtype=torch.float64, device=device) - trans @ trans.T
        eigval, eigvec = torch.linalg.eigh(proj @ hess @ proj)
        # drop the translation modes, then keep the lowest, including negative ones
        keep = torch.nonzero((trans.T @ eigvec).pow(2).sum(0) < 0.5).flatten()
        sel = keep[torch.argsort(eigval[keep])][: cfg["anchor_dirs"]]
        dirs = eigvec[:, sel].T.contiguous()
        anchors.append(
            dict(graph=g, pos=pos, dirs=dirs, target=(hess @ dirs.T).T.contiguous())
        )
    return anchors


# -------------------------------------------------- training
def main() -> None:
    """Run the PFT fine-tuning and save the epoch-3 EMA weights."""
    ap = argparse.ArgumentParser(
        description="PFT fine-tuning of DPA-4.1-Max-OAM (simplified reproduction)"
    )
    ap.add_argument("--base", required=True, help="DPA-4.1-Max-OAM checkpoint (.pt)")
    ap.add_argument(
        "--pft-train", required=True, help="directory of the PFT train .aselmdb shards"
    )
    ap.add_argument(
        "--train-keys",
        required=True,
        help="text file: one '<shard>:<id>' key per used PFT row",
    )
    ap.add_argument(
        "--mptrj",
        required=True,
        help="MPtrj training frames with energy/forces/stress (ase-readable)",
    )
    ap.add_argument(
        "--anchors", required=True, help="the 64 anchor primitive cells (ase-readable)"
    )
    ap.add_argument(
        "--out", required=True, help="output checkpoint (.pt) with the EMA weights"
    )
    args = ap.parse_args()
    cfg = CFG
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    wrapper, model = load_model(args.base, device)
    # tensors frozen by design stay frozen
    params = [p for p in model.parameters() if p.requires_grad]
    opt = torch.optim.AdamW(
        [
            {
                "params": [p for p in params if p.ndim >= 2],
                "weight_decay": cfg["weight_decay"],
            },
            {"params": [p for p in params if p.ndim < 2], "weight_decay": 0.0},
        ],
        lr=cfg["lr"],
        betas=(0.9, 0.999),
        eps=1e-8,
    )
    ema = [p.detach().clone() for p in params]

    pft = PFTRows(args.pft_train)
    keys = [k for k in Path(args.train_keys).read_text().split() if k]
    mptrj = [
        a
        for a in ase.io.iread(args.mptrj, index=":")
        if len(a) <= cfg["cotrain_max_atoms"]
    ]
    anchors = build_anchors(
        model, list(ase.io.iread(args.anchors, index=":")), device, cfg
    )

    cycles_per_epoch = math.ceil(len(keys) / cfg["batch"])
    total_cycles = cycles_per_epoch * cfg["epochs"]
    step = 0

    def optimizer_step() -> None:
        nonlocal step
        torch.nn.utils.clip_grad_norm_(params, cfg["clip_norm"])
        for group in opt.param_groups:
            group["lr"] = cfg["lr"] * min(1.0, (step + 1) / cfg["warmup_steps"])
        opt.step()
        opt.zero_grad(set_to_none=True)
        step += 1
        decay = min(cfg["ema_decay"], (1 + step) / (10 + step))
        with torch.no_grad():
            for e, p in zip(ema, params, strict=True):
                e.mul_(decay).add_(p.detach(), alpha=1 - decay)

    for cycle in range(total_cycles):
        epoch, pos_in_epoch = divmod(cycle, cycles_per_epoch)
        if pos_in_epoch == 0:  # new epoch: new permutation and Hessian column per row
            rng = np.random.default_rng(cfg["seed"] + 1000 * epoch)
            perm = rng.permutation(len(keys))
            col_seed = rng.integers(0, 2**31, size=len(keys))

        # ---- 1. phonon step: Hessian column + force + stress (+ energy) + anchors
        idx = perm[pos_in_epoch * cfg["batch"] : (pos_in_epoch + 1) * cfg["batch"]]
        batch = [(pft.get(keys[i]), int(col_seed[i])) for i in idx]
        sum_n = float(sum(len(atoms) for (atoms, _, _), _ in batch))
        b = len(batch)
        for (atoms, e_dft, hess), cs in batch:
            n = len(atoms)
            pos, types, cell = to_tensors(model, atoms, device)
            g = Graph(model, pos, types, cell)
            j = cs % (3 * n)
            e_j = torch.zeros(3 * n, dtype=torch.float64, device=device)
            e_j[j] = 1.0
            hcol, e, f, s = hessian_times(
                model, g, pos, e_j.view(n, 3), cfg["fd_delta"], create_graph=True
            )
            h_dft = torch.as_tensor(hess[:, j], device=device)
            loss = (
                (cfg["w_hessian"] * n / sum_n) * (hcol - h_dft).abs().mean()
                + (cfg["w_force"] * n / sum_n) * f.norm(dim=1).mean()
                + (cfg["w_stress"] / b) * s.abs().mean()
            )
            if cfg["w_energy"]:
                loss = loss + (cfg["w_energy"] / b) * (e - e_dft).abs() / n
            loss.backward()  # gradients accumulate over the batch
        # each drawn anchor gives one term (one of its six directions), weight / 8
        rng_a = np.random.default_rng(cfg["seed"] + 99991 * (cycle + 1))
        picks = rng_a.choice(
            len(anchors),
            size=min(cfg["anchors_per_step"], len(anchors)),
            replace=False,
        )
        dir_ids = [
            int(rng_a.integers(0, anchors[int(p)]["dirs"].shape[0])) for p in picks
        ]
        for p, k in zip(picks, dir_ids, strict=True):
            a = anchors[int(p)]
            v = a["dirs"][k].view(-1, 3)
            hv = hessian_times(
                model, a["graph"], a["pos"], v, cfg["fd_delta"], create_graph=True
            )[0]
            scale = float(a["target"][k].abs().mean()) + 1e-3
            term = (hv - a["target"][k]).abs().mean() / scale
            ((cfg["anchor_weight"] / len(picks)) * term).backward()
        optimizer_step()

        # ---- 2. co-training steps on MPtrj (energy, forces, stress)
        for kstep in range(cfg["cotrain_steps"]):
            rng_c = np.random.default_rng(
                cfg["seed"] + 7 * (cycle * cfg["cotrain_steps"] + kstep) + 1
            )
            frames = [
                mptrj[i]
                for i in rng_c.choice(
                    len(mptrj), size=cfg["cotrain_batch"], replace=False
                )
            ]
            tot_n = float(sum(len(a) for a in frames))
            for atoms in frames:
                n = len(atoms)
                pos, types, cell = to_tensors(model, atoms, device)
                e, f, s = energy_forces_stress(
                    model, Graph(model, pos, types, cell), pos, create_graph=True
                )
                f_ref = torch.as_tensor(atoms.get_forces(), device=device)
                # ASE stress (tensile positive) = -virial / volume in the deepmd
                # convention used for `s`
                s_ref = -torch.as_tensor(atoms.get_stress(voigt=False), device=device)
                loss = (
                    (cfg["cw_energy"] / cfg["cotrain_batch"])
                    * (e - atoms.get_potential_energy()).abs()
                    / n
                    + (cfg["cw_force"] / tot_n) * (f - f_ref).norm(dim=1).sum()
                    + (cfg["cw_stress"] / cfg["cotrain_batch"])
                    * (s - s_ref).abs().mean()
                )
                loss.backward()
            optimizer_step()

    # ---- save the EMA weights in the deepmd layout (loadable by DP and `dp freeze`)
    with torch.no_grad():
        for p, e in zip(params, ema, strict=True):
            p.copy_(e)
    wrapper.train_infos = dict(
        wrapper.train_infos, pft_cycle=total_cycles, pft_epoch=float(cfg["epochs"])
    )
    torch.save({"model": wrapper.state_dict()}, args.out)


if __name__ == "__main__":
    main()
