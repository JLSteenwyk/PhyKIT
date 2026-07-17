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
| Critical | Installed entry points | Registry/spec/wheel metadata parity is checked, and 121 manually listed tests invoke `python -m phykit`; installed targets were not loaded. | Registry-derived checks require every handler and module wrapper to be callable, load all 318 entry points from an installed wheel, and run `--help` through all 109 canonical `pk_*` executables. | Completed and enforced in CI |
| High | Concatenation integrity | Valid sequential/process workflows and output files were covered, but the paths handled duplicates differently. | Sequential, process-pool, and deterministic fallback paths now reject duplicate taxa and unequal within-locus lengths with identical context before writing; tests assert no FASTA, partition, occupancy, plot, or temporary artifacts. | Completed |
| High | Codon dN/dS ML | NG86/LWL85/YN00 numerical references, argument forwarding, stop/gap behavior, JSON, and CLI errors are covered. | Real F1x4, F3x4, and F61 estimates now have numerical regressions; identical pairs return exact zeros; alternate-code and estimator-failure paths are covered; F61 runs through the CLI; and optional PAML `codeml` validation supplies an independent F3x4 oracle. Compatibility handling repairs F1x4/F61 on supported Biopython releases whose private ML implementation adds `dict_keys` to a list. | Completed |
| High | Threshold MCMC proposals | Short chains covered bounds, sampling, reproducibility, and ordinary adaptive tuning, but no test reached prolonged high-acceptance adaptation. | A 100,000-generation phytools comparison exposed an effectively nonterminating repeated-reflection loop after proposal variance grew very large. Correlation proposals now use constant-time reflection, reject non-finite values, and cap proposal standard deviation at half the bounded domain width; huge-value and repeated-adaptation regressions cover the failure, and the unchanged external comparison completes. | Completed |
| High | Shared trait parsing | Comments, filtering, widths, ordinary nonnumeric values, and optimized ordered paths were covered. | The shared parser now rejects duplicate/blank taxa, duplicate/blank trait headers, and `nan`/infinite values with logical row, taxon, trait, and raw-value context; a PGLS CLI regression verifies exit code 2 and rendering. | Completed |
| High | Newick boundaries | Valid fast paths are compared with Bio.Phylo and malformed numbers fall back. | The shared loader now rejects empty or multiple documents, unmatched structure, missing terminators, trailing content, invalid/non-finite lengths, and duplicate/blank tips; valid quoted/commented/scientific-notation trees are compared with Bio.Phylo, and `tip_labels` verifies exit code 2 without result output or traceback. | Completed |
| High | Branch-length domains | Missing lengths, defaults, minimum tips, valid zero-length cases, and comparative calculations are covered. | Shared analysis validation now rejects negative and non-finite non-root lengths for required and default-filled branch workflows, names the branch/value/context, preserves zero-length edges, covers standard and fallback trees, and prevents `phylogenetic_signal --json` from emitting a result. | Completed |
| High | Partition subsampling | Valid site/gene/partition modes and output assembly were covered. | Full-row partition parsing and preflight checks now reject start < 1, end < start, end beyond alignment, malformed rows, and unequal sequence lengths before writing either artifact; CLI rendering is covered. | Completed |
| Medium | JSON finiteness | Built-ins, NumPy conversion, unsupported objects, caching, and broken pipes are covered. | The shared encoder now converts nested built-in and NumPy NaN/infinities to `null`, enforces `allow_nan=False`, and is verified with a strict decoder. | Completed |
| Medium | Character-map validation | Matrix widths, taxa, valid selectors, optimization modes, exact polytomies, JSON, and plotting are covered. | Argparse and service callers now reject invalid optimization modes; selector syntax, uniqueness, sign, and matrix bounds are validated; blank states and blank/duplicate headers fail before computation; subprocess tests assert no JSON or plot artifact. | Completed |
| Medium | Network JSON schema | Valid quartet structures, ties, gamma, and invalid topology strings are covered. | `network_signal` now rejects malformed JSON and wrong-shaped roots/rows/fields, duplicate quartet/taxon records, invalid classifications/topologies, and nonnumeric, negative, or non-finite counts with filename/record context; CLI tests assert exit 2 without a result or traceback. | Completed |

## Test infrastructure gaps

1. **Completed:** `tests/integration/tree/test_hybridization_integration.py` is
   explicitly marked, so its seven tests run only in the integration selection.
2. **Completed:** optional external comparisons use a `validation` marker and
   dedicated `make test.validation` target, outside ordinary unit coverage;
   primary test targets use the active `python` environment consistently.
3. **Completed:** threshold/phytools parsing ignores unrelated R output, checks
   required keys, and shares each seeded 100,000-generation run across assertions;
   the complete PAML, phangorn, and phytools suite passes eight tests in 69.82s.
4. **Completed:** broken-pipe integration tests now close an unbuffered
   producer's stdout directly and assert PhyKIT's own return code and stderr;
   this exposed and fixed uncaught completion output in concatenation and
   variable-sites services.
5. **Completed:** coverage targets now measure only `phykit`, collect branch
   coverage for unit and integration selections, and enforce an 80% package
   threshold on the unit suite.
6. **Completed:** the full supported-Python matrix remains on macOS, while a
   dedicated CI matrix builds and installs the wheel on Linux and Windows,
   loads all 318 entry points, and runs all 109 canonical help commands outside
   the source checkout.

## Deliberate non-gaps

- More valid character-map polytomy examples are low value because brute-force
  oracle, child-order, nested multifurcation, ACCTRAN, DELTRAN, service, and plot
  tests already exist.
- Repeating shared trait-parser failures through every consumer would add cost
  without new behavior; one representative CLI rendering test is sufficient.
- Running all 317 aliases through full scientific analyses is unnecessary. Load
  every installed target, run canonical help, and reserve real workflows for
  representative aliases and commands.

## Final results

Final package coverage was measured on Python 3.11.14 and macOS after the
implementation commits. Validation tests that depend on external scientific
software remain outside these percentages.

| Selection | Baseline | Final | Statement coverage | Branch coverage |
| --- | ---: | ---: | ---: | ---: |
| Core unit | 5,183 passed | 5,289 passed | 87.76% -> 87.82% (+0.06 pp) | 78.34% -> 78.52% (+0.18 pp) |
| Core unit + integration | 6,017 passed | 6,140 passed | 88.79% -> 88.86% (+0.08 pp) | 79.54% -> 79.72% (+0.18 pp) |
| Integration alone | not comparably scoped | 851 passed | 61.86% | 47.38% |

The combined suite covers 42,255 of 47,551 statements and 13,593 of
17,050 branches. Relative to baseline, it covers 280 more statements and 154
more branches while production code grew by 274 statements and 154 branches;
missing statements decreased from 5,302 to 5,296 and missing branches remained
at 3,457. The modest percentage movement reflects the added validation code as
well as the 123 added tests; the audit targeted consequential behavior rather
than maximizing a percentage.

The final ordinary suite passed on every supported interpreter:

| Environment | Evidence |
| --- | --- |
| Python 3.10 | Fresh local environment and final macOS CI job passed |
| Python 3.11 | 6,140 local tests passed with combined branch coverage; final macOS CI job passed |
| Python 3.12 | Fresh local environment and final macOS CI job passed |
| Python 3.13 | 6,128 tests passed in a fresh local environment before the final focused regressions; the final macOS CI job passed the complete suite |
| Installed wheel, Linux and Windows | All 318 entry-point targets loaded and all 109 canonical help commands passed in CI |
| Documentation wheel, macOS | All 21 tutorial smoke workflows passed locally and the final docs CI job passed |
| External scientific validation | All 8 available PAML, phangorn, and phytools comparisons passed locally |

All seven jobs passed in the final implementation
[GitHub Actions run](https://github.com/JLSteenwyk/PhyKIT/actions/runs/29545999899):
Python 3.10-3.13, documentation, and installed-wheel smoke tests on Linux and
Windows.

## Remaining gaps

- No unresolved critical or high-priority gap identified by this audit remains.
- `concordance_asr`, `discordance_asymmetry`, `dtt`, `evo_tempo_map`, and
  `version` remain `unit_only_documented`. Their service and output paths have
  direct coverage; add another layer when a concrete parser, process,
  filesystem, packaging, or serialization contract warrants it.
- PAML, phytools, and phangorn comparisons remain opt-in validation tests
  because they require external executables or R packages and can be much
  slower than the ordinary suite. They should be run for relevant scientific
  changes, not made a prerequisite for every commit.
- The generated inventory remains the source for lower-priority command-layer
  breadth. New public commands must update it and should include a direct unit
  module plus integration or tutorial evidence when they cross a system
  boundary.
