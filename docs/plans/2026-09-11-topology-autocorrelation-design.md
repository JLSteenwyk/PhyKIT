# Topology autocorrelation specification

Status: prespecified implementation and calibration plan; inference is not yet
validated. No version bump or release is part of this work.

## Estimand

Use strict topology_landscape classifications and its midpoint coordinates.
Retain each gene once per reference. Sort chromosomes separately by anchor and
gene ID. Analyze unordered pairs of resolved genes only, including distinct
genes at tied anchors (distance zero). Bins are [lower, upper), including the
last bin; the maximum distance is an exclusive upper bound. No pairs cross
chromosome, reference, or selected-interval boundaries.

For chromosome c with n resolved genes and n_k in topology k, define
q_ck = n_k(n_k-1)/(n(n-1)). This is the exact joint-support expectation under
uniform reassignment of the observed labels to its resolved gene positions.
The aggregate baseline q_c is the sum of q_ck. Chromosomes with n < 2 have
undefined baselines and contribute no pairs.

For distance bin b let M_cb be its resolved pair count and S_cbk the number of
pairs with both genes in topology k. Observed joint support is S_cbk/M_cb;
aggregate agreement is sum_k S_cbk/M_cb. Excess joint support is the observed
joint support minus q_ck, not a conditional probability or a Pearson
correlation. A relative enrichment ratio is secondary and undefined when its
baseline is zero. Aggregate excess is the sum of topology-specific excesses.

Reference-level estimates sum S across chromosomes and divide by sum M.
Their baseline is sum_c M_cb q_ck / sum_c M_cb, never computed from pooled
gene frequencies. Report both chromosome and reference estimates, actual
contributing genes/pairs/chromosomes, and all six classification counts.
Zero-pair estimates are null. All distances refer to the observed gene set,
not an extrapolated continuous genome or chromosome length.

The descriptive baseline assumes label exchangeability only for its
interpretation as a random-label expectation. It does not establish a
biological null: spatially changing frequencies, classification error, and
nonrandom unresolved genes can all cause apparent excess.

## Computation and uncertainty

Compute right-neighbor pair counts with sorted-anchor range searches and
topology prefix sums. Attach each pair contribution to its left endpoint;
ties use gene ID. Retain the original right endpoint even across block
boundaries. Never concatenate resampled genomes or make new pairs.
This permits O(B n log n) computation without an n-by-n distance matrix or
storing all pairs. Aggregate marks into blocks before resampling.

Use a chromosome-stratified, nonoverlapping physical-block marked bootstrap.
For requested block size L, divide each observed span (first through last
selected anchor, inclusive) into floor(span/L) blocks, at least one, using
near-equal integer widths that cover the entire span. This avoids a short
terminal block. Empty blocks remain eligible for sampling. Each block stores
resolved topology counts, pair denominators, and three joint-support counts.
Draw as many blocks with replacement as originally present, independently
within each chromosome; recompute chromosome baselines and pair-weighted
reference estimates for each replicate. The same draws apply to every bin
and topology. Resampling marks, rather than endpoints independently,
preserves local pair dependencies within sampled blocks and never creates
artificial adjacency. Dependence between blocks remains an approximation.

Use percentile, pointwise intervals for excess agreement and excess joint
support. These describe repeated spatial sampling under approximately
stationary, short-range-dependent marked gene processes within chromosomes;
they are not conditional random-label null intervals or simultaneous bands.
No permutation p-values, significance calls, or clustering-range estimates
are implemented initially. Report clustering range as not estimable.

Require user-specified block sizes for inference and evaluate at least two
sizes. Withhold intervals where a contributing chromosome has fewer than
20 occupied blocks, a block size is less than five times the bin upper
distance, fewer than 95% of resamples are valid, or a bootstrap distribution
is degenerate. Report these reasons rather than silently dropping a
chromosome. At least 499 replicates are required for reported intervals.
Flag interval widths changing by more than a factor of two across admissible
block sizes as block-size sensitive. These are safeguards, not proofs of
stationarity or calibration. Users must inspect local class frequencies and
unresolved-gene distributions, which are exported as block diagnostics.

## Prespecified validation

1. Independently enumerate hand fixtures and random small inputs to verify
   bin membership, tied positions, finite baselines, weighted pooling, and
   bootstrap sufficient statistics. Exhaustively permute a small fixed
   label set to verify mean excess equals zero under the stated baseline.
2. Simulate independent labels and stationary persistence/reset labels on
   regular and irregular coordinates, with unequal chromosome frequencies,
   random missingness, and imbalanced topology frequencies. Use known
   conditional pair probabilities (or independent high-precision Monte
   Carlo targets) for the expected finite-sample excess statistic.
3. For 95% pointwise intervals, use at least 200 independent datasets per
   calibration scenario, 499 bootstrap replicates, and at least two block
   sizes. Prespecify acceptable coverage as 0.90-0.99, reporting binomial
   Monte Carlo uncertainty, mean widths, bias, and withholding frequency.
   Evaluate aggregate and all nondegenerate topology-specific estimates.
   Do not tune acceptance limits after observing results. If coverage fails,
   investigate and document changes; rerun using fresh validation seeds.
4. Simulate nonstationary frequencies and spatially biased missingness as
   intentional assumption violations. Demonstrate confounding rather than
   asserting these scenarios validate stationary-process intervals.
5. Test input errors, singleton and empty bins, multiple references,
   reproducibility, CLI parity between raw and classified input, plots,
   documentation examples, and full-suite regression checks. Benchmark
   increasing n and distance horizons with wall time and peak memory.

## Prior work and scope

Pollard et al. (2006), PLoS Genetics 2:e173,
https://doi.org/10.1371/journal.pgen.0020173, investigated spatial clustering
of gene-tree support and topology block lengths. Spatial topology analysis
is established, not a novelty claim here.

Loh and Stein (2004), Statistica Sinica 14:69-101,
https://web.njit.edu/~loh/Papers/Sinica.LohStein.2004.pdf, motivate resampling
precomputed local pair contributions rather than constructing artificial
point configurations. Our categorical, chromosome-stratified estimator and
its finite-frequency baseline are an adaptation requiring separate
validation, not a reproduction of their K-function inference.

Reuse of strict focal classification follows topology_landscape, which
documents its distinction from Twisst (Martin and Van Belleghem, 2017,
https://doi.org/10.1534/genetics.116.194720).

This analysis does not establish introgression, identify recombination
breakpoints, or estimate a number of independent loci. It measures dependence
of estimated gene-tree labels, not necessarily true genealogies. Publish
implementation, calibration evidence, limitations, and runnable tutorial
through incremental commits to main; do not release a new package version.
