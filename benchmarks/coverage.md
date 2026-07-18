# Coverage improvement baseline

This document records the baseline for the effort to raise PhyKIT's combined
coverage to approximately 92%. Coverage gains must come from behavioral tests;
production modules must not be omitted and branch measurement must remain
enabled.

## Baseline

- Date: 2026-07-17
- Commit: `2c6b4d966dbf9f1facac755052872afd5c411bf9`
- Python: 3.11.14
- Unit tests: 5,322 passed, 867 deselected
- Integration tests: 859 passed, 5,330 deselected
- Statements: 42,756 of 48,108 covered (88.87%)
- Branches: 13,755 of 17,244 covered (79.77%)
- Partial branches: 2,299
- Combined statement-plus-branch coverage: 86.47%
- Codecov coverage for the same commit: 84.08%

Codecov merges the unit and integration XML sessions and accounts for partial
branch lines differently from Coverage.py's local aggregate. The branch target
therefore needs local headroom above 92% before the Codecov target can be
considered achieved.

## Progress checkpoint

- Date: 2026-07-17
- Commit: `f75ab4aa`
- Unit tests: 5,520 passed, 867 deselected
- Integration tests: 859 passed, 5,528 deselected
- Statements: 43,706 of 48,110 covered (90.85%)
- Branches: 14,252 of 17,244 covered (82.65%)
- Partial branches: 2,094
- Combined statement-plus-branch coverage: 88.68%

This checkpoint includes focused helper, alignment, and ancestral
reconstruction coverage batches. Relative to the baseline, the suite covers
950 additional statements and 497 additional branches while retaining the
full branch-enabled unit and integration runs.

## Second progress checkpoint

- Date: 2026-07-17
- Commit: `f7de88e6`
- Unit tests: 5,597 passed, 867 deselected
- Integration tests: 859 passed, 5,605 deselected
- Statements: 44,007 of 48,110 covered (91.47%)
- Branches: 14,411 of 17,244 covered (83.57%)
- Partial branches: 2,017
- Combined statement-plus-branch coverage: 89.39%

This checkpoint adds the OUwie and shared VCV utility coverage batches. Since
the baseline, the suite covers 1,251 additional statements and 656 additional
branches, increasing combined coverage by 2.92 percentage points.

## Third progress checkpoint

- Date: 2026-07-17
- Commit: `2a0e5ebb`
- Unit tests: 5,712 passed, 867 deselected
- Integration tests: 859 passed, 5,720 deselected
- Statements: 44,416 of 48,110 covered (92.32%)
- Branches: 14,624 of 17,244 covered (84.81%)
- Partial branches: 1,916
- Combined statement-plus-branch coverage: 90.34%

This checkpoint includes the discrete-model, independent-contrast, polytomy,
stochastic-character-map, and tree-space coverage batches. Since the baseline,
the suite covers 1,660 additional statements and 869 additional branches,
increasing combined coverage by 3.87 percentage points.

## Fourth progress checkpoint

- Date: 2026-07-17
- Commit: `0e5c734f`
- Unit tests: 5,898 passed, 867 deselected
- Integration tests: 859 passed, 5,906 deselected
- Statements: 44,846 of 48,110 covered (93.22%)
- Branches: 14,906 of 17,244 covered (86.44%)
- Partial branches: 1,734
- Combined statement-plus-branch coverage: 91.43%

This checkpoint adds tree-base, internal-branch-statistics, LTT,
threshold-model, quartet-network, concordance-ASR, and rate-heterogeneity
coverage batches. Since the baseline, the suite covers 2,090 additional
statements and 1,151 additional branches, increasing combined coverage by
4.96 percentage points.

## Final checkpoint

- Date: 2026-07-17
- Commit: `dd8bb2e6`
- Unit tests: 5,970 passed, 869 deselected
- Integration tests: 861 passed, 5,978 deselected
- Unit-only combined statement-plus-branch coverage: 90.81%
- Statements: 45,087 of 48,110 covered (93.72%)
- Branches: 15,042 of 17,244 covered (87.23%)
- Partial branches: 1,670
- Combined statement-plus-branch coverage: 92.00508%

The final batches add shared plot configuration, hybridization, OU shift
detection, and tip-to-tip node-distance coverage. Since the baseline, the suite
covers 2,331 additional statements and 1,287 additional branches, increasing
combined coverage by 5.54 percentage points. No production module remains
below 80% coverage.

The exact result has insufficient cross-version headroom for a 92% hard gate.
CI therefore enforces 90% on the unit-only report and 91.5% on the merged
Codecov project report. The local combined target also fails below 91.5%.

## Codecov-accounting checkpoint

- Date: 2026-07-17
- Commit: `980f685b`
- Unit tests: 6,126 passed, 869 deselected
- Integration tests: 861 passed, 6,134 deselected
- Statements: 45,582 of 48,110 covered (94.75%)
- Branches: 15,305 of 17,244 covered (88.76%)
- Partial branches: 1,531
- Combined statement-plus-branch coverage: 93.16492%

The additional batches target partially covered lines because Codecov treats
those lines as uncovered in its merged report. They exercise plotting and
color annotations, numerical fallbacks, malformed input handling, generic
tree traversal, scalar simulation, and user-facing validation behavior. Since
the baseline, the suite covers 2,826 additional statements and 1,550
additional branches, increasing local combined coverage by 6.69 percentage
points. All production modules remain at or above 80% local combined coverage.

## Reproduction

Run the complete unit and integration suites and create an aggregate report:

```shell
make coverage.combined
```

The target prints a two-decimal report, fails below 91.5%, retains raw data in
the ignored `.coverage` file, and writes machine-readable totals to the ignored
`output/combined.coverage.json` file. The existing `make test.coverage` target
continues to generate the separate XML reports uploaded by CI, with a 90%
unit-only floor.
