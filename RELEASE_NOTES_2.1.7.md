# PhyKIT 2.1.7 Release Notes

## Overview
This release adds a new tree utility for consensus inference with explicit handling of missing taxa across input trees.

## New Functionality
- Added `consensus_tree` (aliases: `consensus`, `ctree`)
- Supports consensus methods:
  - `majority` (default)
  - `strict`
- Supports missing taxa handling:
  - `--missing-taxa error` (default)
  - `--missing-taxa shared` (prune all trees to shared taxa)
- Supports `--json` output with consensus metadata and Newick output

## CLI Entry Points
- `pk_consensus_tree`
- `pk_consensus`
- `pk_ctree`

## Documentation
- Updated usage docs with a new "Consensus tree" section
- Updated changelog for version 2.1.7

## Testing
- Added unit tests for parsing and missing-taxa handling
- Added integration tests for command + aliases
- Updated CLI dispatch coverage and `pk_` help runner coverage
