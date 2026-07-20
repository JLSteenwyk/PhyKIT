# PhyKIT 2.4.1 Release Notes

## Overview

PhyKIT 2.4.1 corrects how `covarying_evolutionary_rates` handles incomplete
taxon sampling and several branch-level edge cases.

## Correctness fixes

- Prune both gene trees and the reference tree to their shared taxa before
  calculating branch-specific relative rates.
- Restore validation that the pruned gene trees and reference tree have the
  same rooted topology.
- Retain valid zero-length gene branches while excluding the root stem from
  the correlation.
- Preserve branch-label and rate-vector alignment in fallback calculations.
- Avoid premature rounding of branch-specific relative rates.
- Report clear errors when too few comparable branches remain or either rate
  vector is constant.

## Validation and documentation

- Add an independent R validation against `ape` and `stats::cor.test`.
- Expand the command reference and tutorial to describe missing-taxon
  handling and distinguish PhyKIT's CovER ratio method from residual-based
  ERC2 and RERconverge workflows.

The command-line interface and output format remain unchanged.
