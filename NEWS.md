# fda.clust 0.1.2

## New Features
- Included three real datasets:
  - **ECG200**: Electrical signals for heartbeats (2 classes: normal and myocardial infarction).
  - **ECG5000**: A larger dataset of 5000 heartbeats (4 classes).
  - **growth_ldata**: Longitudinal growth data from the Berkeley Growth Study.

## Changes
- **Removed** the following functions:
  - `rproc2mu`
  - `rproc2clust`
- **New Classification of Functions**:
  - Clustering functions: `fkmeans`, `fdbscan`, `fmeanshift`, `fhclust`
  - Cluster validation: `fclust.measures`
  - Simulation: `rprocKclust`, `rprocKmu`
  - Internal utility functions: `kmeans.assig.groups`, `kmeans.centers.update`, `kmeans.fd.dist`

## Improvements
- **Vignettes** for the **Introduction** and **Simulations**.
- **New datasets**: **ECG200**, **ECG5000**, and **growth_ldata**.
