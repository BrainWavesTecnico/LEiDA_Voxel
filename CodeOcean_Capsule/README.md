# LEiDA_Voxel Code Ocean Capsule

A self-contained copy of the voxel-space LEiDA pipeline, starting from
already-extracted leading eigenvectors, structured as a standard Code Ocean
capsule (`code/`, `data/`, `results/`).

This folder is independent of the rest of the repository — everything the
capsule needs lives inside `code/` (copied from the shared function library
one level up, plus `combat/` and `utilities/`), so it can be uploaded to
Code Ocean as-is.

## Structure

- `code/run_LEiDA_Voxel_CodeOcean.m` — entry point. Run this script directly;
  it defaults to reading from `../data/` and writing to `../results/`.
- `data/` — put your demo eigenvector file and `Scores_ADNI` table here (see
  [`data/README.md`](data/README.md)). Not tracked in git.
- `results/` — all outputs (clusters, occupancies, stats, figures) are
  written here automatically.

## Preparing the demo data

On your own machine, with the full cohort's already-extracted eigenvectors
and `Scores_ADNI` table, run [`../Select_Demo_Subsample.m`](../Select_Demo_Subsample.m)
to pick a small, balanced, unique-participant demo subsample (e.g. 30 scans
per condition) and save the two files this capsule expects, then copy them
into `data/`.

See the main repository [`README.md`](../README.md) for the full pipeline
documentation.
