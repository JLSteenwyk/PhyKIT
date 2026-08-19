# PhyKIT 2.5.0 Release Notes

## Overview

PhyKIT 2.5.0 adds an experimental topology-robust method for estimating
evolutionary-rate covariation between genes whose inferred trees disagree.
The existing published `covarying_evolutionary_rates` (`cover`) calculation is
unchanged.

## Projected covarying rates

- Added `projected_covarying_rates`, with the short alias `pcover`.
- Prune both gene trees and the reference tree to their three-way shared taxon
  set before analysis.
- Convert each gene tree to pairwise patristic distances and project those
  distances onto an identifiable unrooted reference-tree edge basis using
  nonnegative least squares.
- Normalize projected edge lengths by reference edge lengths, filter extreme
  relative rates, standardize retained rates, and report their Pearson
  correlation.
- Combine complementary branches below a degree-two root because pairwise
  distances identify only their summed unrooted edge length.

## Diagnostics and scalability

- Report normalized root mean squared projection error (NRMSE) for each gene.
  An NRMSE of zero indicates an exact additive fit to the reference topology;
  larger values indicate distance signal that the reference basis does not
  explain.
- Support uniform and inverse-reference distance weighting.
- Support deterministic taxon-pair subsampling through `--max-pairs` and
  `--seed`, with explicit rank checks that reject non-identifiable projections.
- Bound distance calculations as well as design-matrix size when pair
  subsampling is active.
- Provide default summary text, per-edge verbose output, structured JSON, and
  optional scatter plots.

## Validation and documentation

- Added exact edge-length recovery tests, discordant-topology tests,
  missing-taxon integration tests, solver fallback tests, plot tests, and a
  synthetic simulation with a known edge-rate correlation.
- Added a canonical command reference, live parser catalog entry, LLM-readable
  documentation, and tutorial guidance for choosing between `cover` and
  experimental `pcover`.
- The new service has 97% unit-only line and branch coverage. Repository-wide
  unit coverage remains above 92%.

## Experimental status

`pcover` is an experimental PhyKIT estimator, not the published CovER
branch-ratio statistic. Interpret its correlation together with both
projection NRMSE values and the retained-edge count. Its ordinary Pearson
p-value is descriptive because projected phylogenetic edges are not
statistically independent.

## Compatibility

- Made JSON Newick serialization from `print_tree` deterministic across
  supported Biopython versions by explicitly formatting branch lengths to five
  decimal places.
- Stabilized generic eigendecomposition validation across supported Python and
  numerical-library versions without changing the production fallback rules.
