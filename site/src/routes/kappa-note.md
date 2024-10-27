> ⚠️ We are working on extending Matbench Discovery to thermal conductivity prediction via modeling of anharmonic phonons.
> Both presentation and results are work in progress.
> More detailed analysis and a new combined metric incorporating both F1 and symmetric relative mean error (SRME) in predicting thermal conductivity κ to follow in the coming weeks.
> Because it tests the 2nd and 3rd order derivatives of the potential energy surface (PES) and higher derivatives expose even subtle discontinuities in the PES, we believe this to be a stricter and more robust metric for measuring both the utility of ML force fields and the physical accuracy of the PES encoded by a model.
> We invite early feedback on this benchmark extension via [GitHub Discussions](https://github.com/janosh/matbench-discovery/discussions).
>
> For details on the modeling task and evaluation method, refer to [arXiv:2408.00755](https://arxiv.org/abs/2408.00755).
> The only difference between the procedure presented by Póta, Ahlawat, Csányi, and Simoncelli, and the results shown here is the relaxation protocol has been updated and unified for all models.
> It it now a sequential cell relax followed by a site relax (change atom positions only). Each relaxation stage has a maximum number of 300 steps, and consists of a single [`FrechetCellFilter`](https://gitlab.com/ase/ase/-/blob/e65782af/ase/filters.py#L495) relaxation with force threshold $10^{-4} \, \text{eV/Å}$. To preserve crystal symmetry, unit-cell angles are not allowed to change. This unified protocol gives the same SRME reported in [arXiv:2408.00755](https://arxiv.org/abs/2408.00755) for all models except M3GNet which improves to a slightly lower error with the new procedure compared with the non-unified relaxation protocol (1.469 -> 1.412).
