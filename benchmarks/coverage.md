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

## Reproduction

Run the complete unit and integration suites and create an aggregate report:

```shell
make coverage.combined
```

The target prints a two-decimal report, retains raw data in the ignored
`.coverage` file, and writes machine-readable totals to the ignored
`output/combined.coverage.json` file. The existing `make test.coverage` target
continues to generate the separate XML reports uploaded by CI.
