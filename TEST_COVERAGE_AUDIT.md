# PhyKIT test coverage audit

This audit records behavioral test coverage on `test/coverage-audit`. Raw
coverage percentages are supporting evidence, not the prioritization rule.
Scientific correctness, public CLI contracts, malformed-input handling, and
cross-layer output behavior take precedence.

## Baseline

The baseline was measured from commit `4b2a49e0` on Python 3.11.14 and macOS.

| Selection | Result | Time | Scope | Statements | Branches |
| --- | --- | ---: | --- | ---: | ---: |
| Repository `not integration` | 5,193 passed, 3 failed, 834 deselected | 40m44s | repository line coverage | 88.55% | not collected |
| Repository `integration` | 834 passed, 5,196 deselected | 10m22s | repository line coverage | 40.52% | not collected |
| Core unit (`tests/unit`) | 5,183 passed | 16m19s | `phykit` | 87.76% | 78.34% |
| Core unit + integration | 6,017 passed | 27m41s | `phykit` | 88.79% | 79.54% |

Integration tests add 483 covered statements and 203 covered branches beyond
the unit suite. The three baseline failures are in the optional phytools
threshold-model validator: its parser assumes every R stdout line contains
`=`, and the same 100,000-generation reference run is repeated three times.
They are test-infrastructure failures, not failures in deterministic
ThresholdModel unit or integration tests.

## Inventory

[`tests/coverage_inventory.json`](tests/coverage_inventory.json) maps all 109
canonical commands, 317 public command names, and 113 service modules to direct
unit tests, integration tests, tutorial wheel smoke workflows, parser contract
coverage, and installed module wrappers. Generate or verify it with:

```console
python scripts/build_test_inventory.py
python scripts/build_test_inventory.py --check
```

Every public service has direct unit coverage. Commands without a separate
integration or tutorial workflow are marked `unit_only_documented`; their unit
modules already exercise service `run()` and output paths, so another layer is
recommended only for a specific parser, process, filesystem, packaging, or
serialization contract. This avoids treating filenames as proof of a gap.

## Audit method

- Ran the repository-native unit and integration coverage commands unchanged.
- Added package-scoped line and branch measurements without external validators.
- Compared service modules, command identities, aliases, entry points, unit
  tests, integration tests, and tutorial smoke workflows.
- Reviewed recent performance, character-map, and codon dN/dS history.
- Ran 20 independent fixed-prompt code audits and verified consensus findings
  directly against source and tests.
- Reproduced current defects before assigning critical or high priority.

## Prioritized gaps

| Priority | Area | Existing coverage | Missing behavior and recommended test | Status |
| --- | --- | --- | --- | --- |
| Critical | DNA threading | Unit and CLI workflows covered ordinary sequences and optimized paths, but optimized tests compared related implementations and froze an incorrect internal-gap result. | Independent unit and real-data CLI oracles now cover leading, internal, repeated, terminal, masked, planned, short-DNA, and terminal-stop paths. The planner operates on one codon per aligned amino acid, so `A-C` + `AAACCC` is `AAA---CCC`. | Completed; no-gap fast-path timings remain equivalent to baseline |
| Critical | Installed entry points | Registry/spec/wheel metadata parity is checked, and 121 manually listed tests invoke `python -m phykit`; installed targets are not loaded. | Derive tests from the registry, require every handler and module wrapper to be callable, load all 318 installed wheel scripts, and run `--help` for every canonical `pk_*` executable. | Callable invariant and wheel-target loading completed; canonical installed help smoke remains open |
| High | Concatenation integrity | Valid sequential/process workflows and output files are covered. Existing tests preserve different duplicate behavior between paths. | Reject duplicate taxa and unequal within-locus lengths identically before either path writes output; assert no partial FASTA, partition, occupancy, plot, or temporary artifacts. | Open |
| High | Codon dN/dS ML | NG86/LWL85/YN00 numerical references, argument forwarding, stop/gap behavior, JSON, and CLI errors are covered. | Exercise real ML with F1x4, F3x4, and F61, exact identical-sequence behavior, alternate genetic codes, estimator failure states, and an independent external numerical oracle. | Open |
| High | Shared trait parsing | Comments, filtering, widths, ordinary nonnumeric values, and optimized ordered paths are covered. | Reject duplicate/blank taxa, duplicate/blank trait headers, and `nan`/infinite values with row, taxon, trait, and value context. | Open |
| High | Newick boundaries | Valid fast paths are compared with Bio.Phylo and malformed numbers fall back. | Differential-test unmatched structure, trailing garbage, multiple trees, empty files, duplicate/blank tips, and non-finite lengths; representative CLIs must exit 2 without normal output or traceback. | Open |
| High | Branch-length domains | Missing lengths, defaults, minimum tips, valid zero-length cases, and comparative calculations are covered. | Comparative/covariance methods should reject negative and non-finite lengths while preserving explicitly supported zero lengths. | Open |
| High | Partition subsampling | Valid site/gene/partition modes and output assembly are covered. | Reject start < 1, end < start, end beyond alignment, malformed rows, and unequal sequence lengths before writing either artifact. | Open |
| Medium | JSON finiteness | Built-ins, NumPy conversion, unsupported objects, caching, and broken pipes are covered. | Require RFC-valid output for nested NaN/infinities, using `null` or an explicit command status; verify with a strict decoder. | Open |
| Medium | Character-map validation | Matrix widths, taxa, valid selectors, optimization modes, exact polytomies, JSON, and plotting are covered. | Reject invalid optimization choices, invalid/out-of-range character selectors, blank states, and duplicate/blank headers before computation or plot creation. | Open |
| Medium | Network JSON schema | Valid quartet structures, ties, gamma, and invalid topology strings are covered. | Reject malformed JSON and wrong-shaped/duplicate/nonnumeric/negative quartet records as user errors without result output. | Open |

## Test infrastructure gaps

1. `tests/integration/tree/test_hybridization_integration.py` lacks the
   integration marker, placing seven integration tests in the unit selection.
2. External validation tests have no marker and therefore run under
   `not integration`; they need an explicit `validation` target.
3. The threshold/phytools stdout parser is brittle and repeats the same expensive
   seeded calculations for every assertion.
4. Broken-pipe integration tests use `producer | head` through `os.system`, so
   they assert the consumer's pipeline status rather than PhyKIT's process status.
5. Coverage targets include tests and repository tooling, omit branch coverage,
   and do not provide a package-level threshold.
6. The supported Python matrix runs only on macOS despite the OS-independent
   package classifier; installed-wheel smoke should include Linux and Windows.

## Deliberate non-gaps

- More valid character-map polytomy examples are low value because brute-force
  oracle, child-order, nested multifurcation, ACCTRAN, DELTRAN, service, and plot
  tests already exist.
- Repeating shared trait-parser failures through every consumer would add cost
  without new behavior; one representative CLI rendering test is sufficient.
- Running all 317 aliases through full scientific analyses is unnecessary. Load
  every installed target, run canonical help, and reserve real workflows for
  representative aliases and commands.

## Final reporting

This document will be updated after implementation with completed changes,
remaining gaps, final coverage, supported-Python results, and CI evidence.
