# Topology autocorrelation validation

## Scope and reproducibility

The prespecified design and revision history are in
`docs/plans/2026-09-11-topology-autocorrelation-design.md`. The scientific
estimator uses finite chromosome-specific random-label baselines and a
pair-weighted reference summary. Uncertainty uses centered endpoint marks,
stratified physical-block resampling, and pointwise percentile intervals.
No significance test, p-value, simultaneous band, or clustering-range
estimator is implemented.

Run the stationary calibration from the repository root:

```shell
venv/bin/python benchmarks/calibrate_topology_autocorrelation.py --datasets 1000 --seed 20261203 --scenarios iid_regular iid_irregular_missing iid_imbalanced clustered_regular clustered_irregular_missing --output benchmarks/topology_autocorrelation_calibration.json
```

Each scenario has two chromosomes (6,000 and 4,000 genes), distance edges
0, 201, 501 bp, physical block sizes 5,000 and 10,000 bp, and 499 replicates
per dataset. Irregular spacing is uniform integer 25-175 bp; regular
spacing is 100 bp. Random unresolved labels occur at a rate of 20% in
missing-data scenarios. The stationary clustered process retains its
previous label with probability exp(-distance/250), otherwise drawing
an independent label from that chromosome's probabilities. The imbalanced
case uses probabilities 0.80, 0.15, 0.05; other nonuniform scenarios use
different probabilities on the two chromosomes.

Targets are derived independently conditional on observed positions and
the missing-label mask. In the clustered process, the probability of joint
support for k is p_k^2 + p_k(1-p_k) exp(-distance/250). The target excess
subtracts the expectation of the finite-sample baseline over *all* pairs,
not merely the limiting marginal probability. Direct small-pair enumeration
tests the analytic calculation, separately from the production prefix-sum
implementation. Exhaustive permutations of a small fixed label multiset
also verify exact zero mean excess under uniform reassignment.

The acceptance criterion was set before simulation: each nondegenerate
95% pointwise interval must have empirical coverage 0.90-0.99 and be
available in at least 95% of datasets. JSON records exact binomial Monte
Carlo intervals, available/withheld counts, widths, targets, and estimator
bias for all metrics, bins, and block sizes. These are checks of the
specified scenarios, not a universal guarantee under arbitrary dependence.

## Final stationary calibration

All 80 metric/bin/block-size cells passed the unchanged acceptance criterion.
All intervals were available in all 1,000 datasets per cell. Nominal 95%
intervals had the following empirical coverage ranges; some undercoverage
remains, so these should be described as approximate pointwise intervals,
not exact 95% guarantees.

| Scenario | Coverage range | Maximum absolute estimator bias |
| --- | --- | ---: |
| Independent, regular spacing | 93.2-95.4% | 0.000091 |
| Independent, irregular spacing and missing labels | 92.8-94.2% | 0.000214 |
| Independent, imbalanced frequencies | 93.3-95.2% | 0.000043 |
| Clustered, regular spacing | 92.1-94.7% | 0.000142 |
| Clustered, irregular spacing and missing labels | 92.5-95.0% | 0.000199 |

The run used Python 3.11.14 / NumPy 2.4.2 and took 538 seconds on the
development machine. The statistical helper is unchanged from commit
`3ac463f9`; later commits add the CLI, plotting validation, and reporting.
Full numerical records are in `topology_autocorrelation_calibration.json`.

## Candidate history and confounding

Three earlier experiments are retained rather than overwritten:

- `topology_autocorrelation_calibration_left_marks.json`: one-sided raw
  marks gave dominant-topology coverage 1.00 in the imbalanced case.
- `topology_autocorrelation_calibration_symmetric_marks.json`: equal
  endpoint allocation reduced, but did not eliminate, overcoverage (0.995).
- `topology_autocorrelation_calibration_centered_200.json`: algebraic
  centering passed the imbalanced scenario; one other estimate was
  179/200 = 0.895. The final experiment extends **every** stationary
  scenario to 1,000 datasets with the same seed sequence, including those
  original 200, with no additional algorithm or acceptance-limit changes.

The centered 200-dataset report also contains intentional assumption
violations: changing topology frequencies along chromosomes, and changing
frequencies combined with spatially biased unresolved labels. Labels are
independent conditional on position in these scenarios, yet expected
aggregate excess is about 0.06. This demonstrates confounding by
nonstationarity, not evidence of a biological dependence mechanism.
Coverage in those scenarios is not treated as stationary-model calibration.

## Performance

Run:

```shell
venv/bin/python benchmarks/benchmark_topology_autocorrelation.py --output benchmarks/topology_autocorrelation_performance.json
```

Measured with Python 3.11.14 and NumPy 2.4.2 on macOS. Ten distance bins,
regular 100-bp anchors, and deterministic labels are used. Time includes
tracemalloc overhead and input validation. Peak additional traced allocation
includes validation copies and NumPy buffers, but excludes preconstructed
input records and interpreter/dependency memory; it is not process RSS.

| Genes | Maximum distance (exclusive bp) | Eligible pairs | Seconds | Additional MiB |
| ---: | ---: | ---: | ---: | ---: |
| 1,000 | 1,000 | 8,955 | 0.028 | 0.59 |
| 1,000 | 1,000,000 | 499,500 | 0.025 | 0.47 |
| 10,000 | 1,000 | 89,955 | 0.201 | 4.77 |
| 10,000 | 1,000,000 | 49,995,000 | 0.178 | 4.69 |
| 100,000 | 1,000 | 899,955 | 2.115 | 48.62 |
| 100,000 | 1,000,000 | 949,905,000 | 1.928 | 48.62 |

The estimator uses range searches and topology prefix sums, not pair
enumeration or an n-by-n matrix. For B bins, counting costs O(B n log n)
after sorting, with O(n + KB) working arrays for K physical blocks.
Bootstrap storage is O(RB + KB), plus batches of at most 100 block-weight
draws; work includes multiplying block marks by replicate weights.
The implementation limits bins, block/bin products, and replicate/bin
products to prevent unbounded accidental allocations. A regression test
counts 4,498,500 dense pairs with less than 15 MiB additional traced memory.

## Interpretation limits

Intervals are approximately calibrated for the validated stationary,
short-range processes with sufficient blocks, not for arbitrary genomic
heterogeneity. Block-size sensitivity, sparse-bin withholding, and local
classification diagnostics remain essential. A homogeneous topology gives
a degenerate bootstrap and is explicitly withheld. Rare topologies,
long-range dependence, or nonrandom tree-estimation failures can make
inference unstable even when numerical safeguards pass.

This is descriptive spatial structure of estimated gene-tree labels.
It neither establishes introgression nor identifies recombination
breakpoints or a number of independent loci. See command documentation
for Pollard et al. (2006) and Loh and Stein (2004), the relevant prior work;
the algebraically centered categorical adaptation is separately evaluated
here and is not presented as their original estimator.
