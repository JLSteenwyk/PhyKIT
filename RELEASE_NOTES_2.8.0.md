# PhyKIT 2.8.0

## Limit threads with `--threads`

Every command, including the `pk_*` entry points, now accepts a global
`--threads N` option. The same limit can be set once with the `PHYKIT_THREADS`
environment variable (GitHub issue #113).

- Caps the worker processes PhyKIT starts and the threads used by numeric
  libraries (NumPy/SciPy BLAS, OpenMP, numba).
- Useful when running many PhyKIT commands in parallel, for example in a
  workflow manager or with `xargs -P`; replaces setting `OMP_NUM_THREADS=1`.
- An explicit `--threads` overrides `OMP_NUM_THREADS` and related variables;
  `PHYKIT_THREADS` alone only sets those that are not already set.
- Without the option or the variable, behavior is unchanged.

```shell
phykit saturation -a alignment.fa -t tree.tre --threads 1
export PHYKIT_THREADS=1
```

## Topology autocorrelation

New command: `topology_autocorrelation` (alias `topo_ac`), also available as
`pk_topology_autocorrelation` and `pk_topo_ac`.

- Genomic distance-binned topology agreement with finite chromosome-specific
  baselines and optional chromosome-stratified marked-block pointwise intervals
  with block-size sensitivity.
- Accepts raw gene trees or classified TSVs, writes separate reference plots,
  and reports explicit missingness diagnostics. Includes simulation
  calibration and a downloadable tutorial.

This is a descriptive analysis: it does not report p-values, clustering ranges,
introgression claims, or breakpoint calls.

## Other changes

- Boolean command-line arguments accept `yes`/`y` and `no`/`n`
  (case-insensitive) in addition to `true`/`t`/`1` and `false`/`f`/`0`.
- Fixed the worker-count calculation in `covarying_evolutionary_rates`
  (performance only; results are unchanged).
- Shared and consolidated trait-file and alignment file-list parsing across
  commands, without changing input contracts, diagnostics, or calculations.
- Removed the unused internal `phykit.helpers.parallel` module.
- Added release-triggered Bioconda recipe synchronization.
