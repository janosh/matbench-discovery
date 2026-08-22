# Artifact ensembles

This family contains deterministic combinations of already-published Matbench
Discovery prediction artifacts. These submissions do not claim that a single ASE
calculator reproduces every task. Each model card states the task-level output
operators explicitly.

## DPA4–EqV3 MPtrj Ensemble

All three submitted tasks use fixed weights `(0.5, 0.5)`: row-aligned discovery
formation-energy means, periodic midpoints of the two published relaxed structures,
and aligned means of the two published PhononDB output records. The rules were frozen
before any composite metric was evaluated; generation did not read WBM or PhononDB
DFT targets. The components' already-public leaderboard results were known when they
were selected. This is a task-level direct-output composite, not a claim of one
callable potential, averaged forces, or averaged force constants. The model YAML
credits and links both component model cards; the combination is submitted by
zihenglutech.
