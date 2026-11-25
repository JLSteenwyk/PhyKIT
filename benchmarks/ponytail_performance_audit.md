# Ponytail performance audit

## Scope and execution

PhyKIT processes multiple sequence alignments, phylogenetic trees, and
comparative trait data through local command-line jobs. The assumed load is
a researcher running individual commands or combining many independent
jobs in a workflow manager, with a thread cap per job.

Audit instructions: [Ponytail audit skill at ebd8654](https://github.com/DietrichGebert/ponytail/blob/ebd8654c2c627f721eb8f9c0d1bd4090a29f9e10/skills/ponytail-audit/SKILL.md).
Followed the skill directly in this thread. The source under audit is
`cf8c7652`. Repository mapping covered package entry points, dependency
declarations, command families, helpers, scripts, release and CI workflows,
tests, and benchmark inventory. Detailed tracing focused on the five
workloads below. The remaining service implementations, tests, release
scripts, and docs were inventoried and still need detailed review.

## Data paths

- `setup.py` creates the main and command-specific entry points from
  `phykit/cli_registry.py`. `phykit/phykit.py` handles parsing, dispatch,
  and service construction. `phykit/helpers/threads.py` applies the worker
  and numeric-library caps before computation.
- Alignment services use the shared format detector and FASTA reader in
  `phykit/helpers/files.py` or a command-specific streaming parser. Pairwise
  identity turns ASCII sequences into a numerical matrix, computes matching
  site counts, and summarizes or emits pairwise values.
- Tree services use shared Newick ingestion and validation in
  `phykit/services/tree/base.py`. Patristic distances also has a compact
  Newick reader that builds arrays of parent indices and root distances.
  Covarying evolutionary rates maps branches by the tip sets they subtend
  before calculating the correlations.
- Phylogenetic signal maps a trait table to tree tips, builds a
  variance-covariance matrix, and computes the observed and permuted
  statistics. The alignment GWAS path reads the phenotype and alignment,
  calculates site associations, and renders a plot through Matplotlib.

## Ranked findings

Profile times below are instrumented. The uninstrumented times and memory
measurements are in [the baseline](profiling_baseline.md).

1. **Pairwise identity does symmetric work in the default matrix path.**
   In [pairwise_identity.py](../phykit/services/alignment/pairwise_identity.py#L896),
   `(sequence_matrix[start:stop, None, :] == sequence_matrix[None, :, :]).sum(axis=2, dtype=np.int32)`
   compares against every sequence before applying the upper-triangle mask.
   The 1,024-by-1,000 gapped input spent 0.820 of 0.888 profiled seconds in
   the stats matrix routine, including NumPy reductions. Compare each
   row block only with the trailing columns that can contain output pairs.
   This is a small change with a large measured target. Validate denominator,
   gap modes, and condensed pair order before reusing it in verbose output.
2. **Patristic distances repeatedly scans identical ancestor prefixes.**
   In [patristic_distances.py](../phykit/services/tree/patristic_distances.py#L649),
   `for clade_a, clade_b in zip(path_a, tip_paths[j]):` redoes the common
   prefix search for every tip pair. At 1,024 tips the method consumed 0.299
   of 0.932 profiled seconds, with 0.137 seconds in Newick scanning as well.
   Compute common ancestors by subtree groups or by a vectorized batch.
   Preserve arbitrary requested tip order and floating-point operations.
   Even eliminating the method entirely would cap this profile's speedup
   at about 1.47× with the rest held fixed.
3. **GWAS rendering dominates the bundled command.**
   [phylo_gwas.py](../phykit/services/alignment/phylo_gwas.py#L1570) uses
   `fig.savefig(self.output_path, dpi=config.dpi, bbox_inches="tight")`.
   Profiled time concentrates in Agg font measurement, font lookup, and
   image encoding. Rendering work offers more room than changing site
   statistics for this fixture. First compare render passes and image
   bounds before changing layout or DPI, which are user-visible.
4. **Tree parsing is a shared but smaller target in trait analysis.**
   [vcv_utils.py](../phykit/services/tree/vcv_utils.py#L198) builds the
   covariance matrix from descendant indices. On the 1,024-tip signal
   workload this method took 0.074 cumulative seconds and shared Newick
   reading took 0.140 of 2.321 profiled seconds. Four threads improved
   uninstrumented median from 0.501 to 0.433 seconds. Measure native linear
   algebra and larger trees before changing the scientific implementation.
5. **Small covarying-rate jobs spend more time in setup than correlation.**
   On the bundled three-tree case, Biopython import work including
   `CodonTable.list_ambiguous_codons` is visible while
   [`_reference_tipsets_and_order_direct`](../phykit/services/tree/covarying_evolutionary_rates.py#L766)
   accounted for about 0.002 of 0.675 profiled seconds. Profile larger tree
   sets before tuning the comparison loop. Investigate import ownership if
   the target is repeated invocation on small trees.

## Reproduce the profiles

```sh
P=$(mktemp -d)
python benchmarks/profile_scaled_workloads.py --profile-dir "$P"
python -m cProfile -o "$P/gwas.prof" phykit-runner.py --threads 1 phylo_gwas \
  --alignment benchmarks/fixtures/phylo_gwas_multigroup.fa \
  --phenotype benchmarks/fixtures/phylo_gwas_multigroup.tsv \
  --output "$P/gwas.png" --json >/dev/null
python -m cProfile -o "$P/rates.prof" phykit-runner.py --threads 1 \
  covarying_evolutionary_rates tests/sample_files/tree_simple.tre \
  tests/sample_files/tree_simple_1.tre \
  --reference tests/sample_files/tree_simple_2.tre >/dev/null
python benchmarks/profile_function_coverage.py "$P"/*.prof > "$P/functions.tsv"
python -m pstats "$P/pairwise_identity-large.prof"
```

These captures observed 158 of 3,113 named Python function definitions
under `phykit/`. The TSV identifies each definition and its matching
profiles. Rows with an empty profiles column still need an exercised
workload. Coverage counts function definitions, including nested functions.

On Linux, the following checks include native sampling and subprocess
tracking. The signal trace reported native threads and one executable
invocation. A process-pool fallback workload remains to be measured.

```sh
python -c 'from pathlib import Path; import sys; from benchmarks.profile_scaled_workloads import make_inputs; make_inputs(Path(sys.argv[1]), 1024, 1000)' "$P"
py-spy record --subprocesses --native --format speedscope -o "$P/identity.json" -- \
  python phykit-runner.py --threads 4 pairwise_identity "$P/alignment.fa" >/dev/null
strace -f -e trace=clone,clone3,fork,vfork,execve -o "$P/signal-processes.log" \
  python phykit-runner.py --threads 4 phylogenetic_signal -t "$P/tree.nwk" \
  -d "$P/traits.tsv" -p 100 >/dev/null
```
