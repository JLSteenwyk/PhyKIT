# Ponytail performance audit

## Results and priorities

These changes are on the [profiling branch](https://github.com/JLSteenwyk/PhyKIT/tree/profiling/phykit-performance).
The [baseline](profiling_baseline.md) records the machine, dependencies,
input sizes, thread caps and run-to-run spread.

| Rank | Workload or candidate | Function speedup | CLI speedup | Evidence and disposition |
| ---: | --- | ---: | ---: | --- |
| 1 | Patristic distances, 1,024 tips | 4.55× | 1.80× | Adopted. Low implementation effort and exact output across Python versions. |
| 2 | Pairwise identity, 1,024 × 1,000 gapped alignment | 1.84× | 1.57× | Adopted. Small change, lower memory, exact condensed values. |
| 3 | Topology uncertainty, 50,000 genes | 1.72× | 1.16× | Adopted. Moderate time savings and identical reports and intervals. |
| 4 | Polytomy batch, 250 trees | 1.29× | 1.28× with one worker, 1.16× with four | Adopted. Cache isolation required a regression test. |
| 5 | Phylogenetic-signal covariance solves | 5.54× at 1,024 taxa and 100 permutations | Unmeasured | Discarded. Repeatability was good on the synthetic balanced input, but tie and near-singular results changed. |
| 6 | GWAS rendering | Unmeasured | 0.96× for the crop experiment | Discarded. Lower memory came with different pixel dimensions and a longer runtime. |
| 7 | Covarying-rates imports | Unmeasured | Unmeasured | Import attribution is strong. Parser substitutions require a wider fixture matrix and compatibility work. |

The 1.5× target was reached by three measured functions. Two workloads
reached it end to end. Small examples were often dominated by startup,
while larger ones exposed algorithmic or memory costs. Function gains
apply to the stated workloads. Most definitions have no isolated
before-and-after timing.

### Review and rollback

| Commit | Purpose |
| --- | --- |
| `7602229b` | Bundled CLI benchmark and baseline |
| `3da65a26` | Scaled baseline, profiler hooks, and function inventory tool |
| `f08bb31e` | Initial Ponytail audit and unoptimized A/B control |
| `0d9ffbca` | Pairwise identity optimization |
| `9d4701b7` | Patristic distances optimization |
| `c9fb2baa` | Worker-process attribution |
| `c6db9a42` | Polytomy optimization |
| `f2ba240c` | Expanded CLI coverage and rendering and import experiments |
| `c4994b87` | Native signal profiling |
| `fe0a7744`, `c464ac74` | Signal correctness result, reproducible experiment, and topology validation measurements |
| `bfc099af` | Topology uncertainty optimization |

The four performance changes touch independent production modules and
can be reviewed or reverted separately. Their benchmark instructions
use the baseline harness and common control commit. The report
sections accumulate across commits, so reverting a production change
may require removing its summary row here as well.

## Scope and execution

PhyKIT processes multiple sequence alignments, phylogenetic trees, and
comparative trait data through local command-line jobs. The assumed load is
a researcher running individual commands or combining many independent
jobs in a workflow manager, with a thread cap per job.

Audit instructions: [Ponytail audit skill at ebd8654](https://github.com/DietrichGebert/ponytail/blob/ebd8654c2c627f721eb8f9c0d1bd4090a29f9e10/skills/ponytail-audit/SKILL.md).
Followed the skill directly in this thread. The source under audit is
`cf8c7652`. Repository mapping covered package entry points, dependency
declarations, command families, helpers, scripts, release and CI workflows,
tests, and benchmark inventory. The function inventory scanned all 150
Python modules in the runtime package. Detailed runtime attribution
focused on the measured workloads below. Releases and generated docs
were inspected as supporting workflows.

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

The generated command registry lists 114 canonical commands across
13 categories. For broader coverage, the audit mapped tree comparison
and consensus through tree paths and split processing, trait evolution
through shared parsing, covariance and discrete model helpers, and
genomic workflows through coordinate validation, topology classification,
and aggregation. Output paths include text and JSON summaries, TSV
reports, and Matplotlib figures. The packaging and release scripts are
separate from these runtime paths.

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

## Validated prototype: pairwise identity summary

Compare each row block with trailing columns, then retain the original
upper-triangle mask. This keeps summation and condensed row ordering intact
and avoids allocating comparisons that cannot produce output. The change
targets the default stats-only matrix fallback with gaps or nonstandard
symbols.

Two alternating control/prototype rounds on the baseline machine, seven
readings per round after one warmup, with one thread:

| Metric | Control rounds | Prototype rounds | Speedup |
| --- | --- | --- | --- |
| 1,024 × 1,000 CLI median | 0.473, 0.478 s | 0.300, 0.304 s | 1.57× |
| 1,024 × 1,000 CLI median peak RSS | 65.7, 64.7 MiB | 49.9, 52.0 MiB | |
| 1,024 × 1,000 stats matrix function median | 0.339, 0.340 s | 0.185, 0.184 s | 1.84× |
| 32 × 300 CLI median | 0.101, 0.106 s | 0.101, 0.106 s | 1.00× |

The stats matrix function and the large CLI workload both met the 1.5×
target. Across rounds, large CLI samples ranged from 0.471 to 0.496 s in
the control and 0.299 to 0.310 s in the prototype. Function samples
ranged from 0.334 to 0.342 s and 0.184 to 0.186 s, respectively.

Run the following from each checkout, alternating control at `f08bb31e`
and the prototype. The function check constructs the records before
timing and prints the summary for comparison.

```sh
python benchmarks/profile_scaled_workloads.py --command pairwise_identity \
  --threads 1 --repeat 7 --gnu-time /usr/bin/time
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python - <<'PY'
from argparse import Namespace
import random
import statistics
from time import perf_counter
from types import SimpleNamespace

from phykit.services.alignment.pairwise_identity import PairwiseIdentity

rng = random.Random(41)
records = [
    SimpleNamespace(id=f"taxon{i:04d}", seq="".join(rng.choices("ACGT-", k=1000)))
    for i in range(1024)
]
service = PairwiseIdentity(Namespace(alignment="synthetic.fa", verbose=False, exclude_gaps=False))
run = service._calculate_pairwise_identity_stats_matrix
expected = run(records, is_protein=False, exclude_gaps=False)
timings = []
for _ in range(7):
    start = perf_counter()
    result = run(records, is_protein=False, exclude_gaps=False)
    timings.append(perf_counter() - start)
    assert result == expected
print(f"{statistics.median(timings):.6f} [{min(timings):.6f}, {max(timings):.6f}] {expected}")
PY
```

The new regression checks the condensed values against a direct scalar
calculation across three blocks. The targeted suites are
`tests/unit/services/alignment/test_pairwise_identity_service.py` and
`tests/integration/alignment/test_pairwise_identity_integration.py`.

On CPython 3.12.14, 3.13.15, 3.14.7, and 3.15.2, those suites passed 85
tests each. On all four versions, control and prototype CLI JSON matched
byte-for-byte for the large default workload and for gap exclusion and
verbose output on `tests/sample_files/test_alignment_0.fa`.

## Validated prototype: patristic distances on preorder tips

The shallow-tree fast path records the last selected tip in each subtree.
For each output row, it assigns the common ancestor once per contiguous
range as it walks up the tree. It retains the existing per-pair ancestor
search for arbitrary tip order and the binary-lifting path for deep trees.
Distance arithmetic and pair order are identical.

Two alternating control/prototype rounds on the baseline machine, with
seven samples per round after a warmup and one thread:

| Metric | Control rounds | Prototype rounds | Speedup |
| --- | --- | --- | --- |
| 1,024-tip CLI median | 0.371, 0.368 s | 0.209, 0.202 s | 1.80× |
| 1,024-tip CLI median peak RSS | 57.0, 58.1 MiB | 58.1, 58.1 MiB | |
| 1,024-tip path-distance function median | 0.206, 0.216 s | 0.047, 0.046 s | 4.55× |
| 32-tip CLI median | 0.105, 0.101 s | 0.103, 0.100 s | 1.01× |

The function and large CLI both met the 1.5× target. Large CLI samples
ranged from 0.365 to 0.385 s in the control and 0.201 to 0.210 s in the
prototype. Function samples ranged from 0.204 to 0.220 s and 0.045 to
0.052 s, respectively. The full-precision condensed distance lists had
the same SHA-256 digest. The uninstrumented CLI ratio exceeds the bound
estimated from the instrumented profile above because profiling changed
the distribution of execution time.

Run from the control at `f08bb31e` and the prototype:

```sh
python benchmarks/profile_scaled_workloads.py --command patristic_distances \
  --threads 1 --repeat 7 --gnu-time /usr/bin/time
P=$(mktemp -d)
python -c 'from pathlib import Path; import sys; from benchmarks.profile_scaled_workloads import make_inputs; make_inputs(Path(sys.argv[1]), 1024, 1000)' "$P"
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python - "$P/tree.nwk" <<'PY'
import hashlib
import statistics
import sys
from time import perf_counter

from phykit.services.tree.patristic_distances import PatristicDistances

_, tips, parents, _, depths = PatristicDistances._scan_simple_newick_distance_inputs(sys.argv[1])
run = PatristicDistances._pairwise_tip_distance_values_from_paths
expected = run(tips, parents, depths)
timings = []
for _ in range(7):
    start = perf_counter()
    result = run(tips, parents, depths)
    timings.append(perf_counter() - start)
    assert result == expected
print(f"{statistics.median(timings):.6f} [{min(timings):.6f}, {max(timings):.6f}]",
      hashlib.sha256(repr(expected).encode()).hexdigest())
PY
```

All 45 targeted unit and integration tests passed on each of CPython
3.12.14, 3.13.15, 3.14.7, and 3.15.2. The 1,024-tip command returned
byte-identical text and JSON in both summary and verbose modes on every
version. The new test covers omitted tips, a multifurcation, a negative
branch length, arbitrary tip order, and repeated tips.

## Worker-process profile: polytomy test

The bundled tutorial contains 250 listed tree inputs and exercises the
actual pool path. Five uninstrumented runs after a warmup gave these
results. Both settings enter the pool path, so one thread means one pool
worker. The CLI emitted identical JSON at both caps.

| Worker cap | Median (s) | Min (s) | Max (s) | Median peak RSS (MiB) |
| ---: | ---: | ---: | ---: | ---: |
| 1 | 0.81 | 0.80 | 0.84 | 48.5 |
| 4 | 0.34 | 0.34 | 0.35 | 38.6 |

GNU time reports a process high-water mark. Concurrent worker RSS is
not represented by adding these values. The four-worker speedup was
2.38×. A separate py-spy 0.4.1 capture recorded four worker stack
streams with 99 samples across them. In 58 samples, the innermost
PhyKIT frame was `_safe_copy`. This is evidence for a new candidate
between the existing pairwise identity and GWAS candidates.

[`_read_tree_with_cache`](../phykit/services/tree/polytomy_test.py#L185)
copies a newly parsed tree into the cache and then copies it again before
returning it. [`_prepare_tree_for_triplets`](../phykit/services/tree/polytomy_test.py#L280)
copies the returned tree before rerooting it. A worker can keep custody of
the freshly parsed tree or receive one private copy of the cache entry,
then root that private copy. Test repeat paths and worker fallback
behavior before changing which tree is mutated.

```sh
ROOT=$PWD
P=$(mktemp -d)
tar -xzf docs/data/polytomy_tutorial.tar.gz -C "$P"
cd "$P/polytomy_tutorial"
for threads in 1 4; do
  python "$ROOT/phykit-runner.py" --threads "$threads" polytomy_test \
    -t filamentous_fungi_250_trees.txt -g groups_file0.txt --json > "$P/result-$threads.json"
  for sample in 1 2 3 4 5; do
    /usr/bin/time -f "$threads %e %M" -a -o "$P/timings.tsv" \
      python "$ROOT/phykit-runner.py" --threads "$threads" polytomy_test \
      -t filamentous_fungi_250_trees.txt -g groups_file0.txt --json >/dev/null
  done
done
cmp "$P/result-1.json" "$P/result-4.json"
py-spy record --subprocesses --native --format speedscope -o "$P/pool.json" -- \
  python "$ROOT/phykit-runner.py" --threads 4 polytomy_test \
  -t filamentous_fungi_250_trees.txt -g groups_file0.txt --json >/dev/null
```

## Validated prototype: reuse private polytomy trees

The cache retains the newly parsed tree and returns a private copy.
The service reroots its private copy in place, avoiding a second copy.
The preparation helper still copies by default for callers that supply
a tree they need to reuse.

Two alternating control/prototype rounds with the 250-tree tutorial:

| Metric | Control rounds | Prototype rounds | Speedup |
| --- | --- | --- | --- |
| One-worker CLI median | 0.84, 0.82 s | 0.65, 0.65 s | 1.28× |
| One-worker CLI median peak RSS | 48.5, 50.0 MiB | 48.5, 47.0 MiB | |
| Four-worker CLI median | 0.35, 0.36 s | 0.30, 0.31 s | 1.16× |
| Four-worker CLI median peak RSS | 39.5, 39.5 MiB | 37.1, 36.5 MiB | |
| Worker batch function median | 0.927, 0.917 s | 0.711, 0.717 s | 1.29× |

Each round used five timings after a warmup. One-worker CLI timings
ranged from 0.81 to 0.85 s in the control and 0.61 to 0.66 s in the
prototype. Four-worker CLI timings ranged from 0.35 to 0.36 s and 0.30
to 0.31 s. Batch timings ranged from 0.815 to 0.997 s and 0.587 to
0.787 s. The worker batch includes parsing, copies, rerooting, and
triplet calculations. The speedups fell below 1.5×.

Repeat the worker-process commands above for each checkout, alternating
control at `f08bb31e` and prototype. To isolate the worker batch, run
from the extracted tutorial directory with `ROOT` set to the checkout:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH="$ROOT" python - <<'PY'
from argparse import Namespace
import hashlib
import json
from pathlib import Path
import statistics
from time import perf_counter

from phykit.services.tree.polytomy_test import PolytomyTest

paths = Path("filamentous_fungi_250_trees.txt").read_text().splitlines()
args = Namespace(trees="filamentous_fungi_250_trees.txt", groups="groups_file0.txt")
timings = []
for repeat in range(6):
    service = PolytomyTest(args)
    groups, outgroup = service.determine_groups_of_groups(service.read_in_groups())
    start = perf_counter()
    result = service._process_tree_batch(paths, groups, outgroup)
    if repeat:
        timings.append(perf_counter() - start)
    digest = hashlib.sha256(json.dumps(result, sort_keys=True).encode()).hexdigest()
    if repeat == 0:
        expected = digest
    assert digest == expected
print(f"{statistics.median(timings):.6f} [{min(timings):.6f}, {max(timings):.6f}]", expected)
PY
```

The full worker batch digests matched. On each of CPython 3.12.14,
3.13.15, 3.14.7, and 3.15.2, the tutorial JSON matched byte for byte with
one and four workers. All 95 polytomy unit and integration tests passed
on each version. The new regression reroots a tree obtained from the
cache and verifies that the next read still returns the original topology
without reparsing.

## Rendering and startup experiments

**GWAS crop experiment.** The proposed removal of the tight bounding-box
pass reduced median peak RSS from 136.5 to 109.5 MiB but increased median
CLI time from 1.00 s (0.97–1.03 s) to 1.04 s (1.00–1.04 s). Six runs of
each variant were alternated. The JSON reports matched byte for byte.
The output PNG changed from 3,496 × 2,100 to 4,200 × 2,400 pixels because
it retained the uncropped canvas. The experiment was discarded. The
smallest follow-up would compare an equivalent crop produced from one
Agg draw, checking pixel bounds and labels on default, customized, and
partitioned plots.

From the control checkout, create a temporary wrapper that only alters
the plotting call. Both variants load Matplotlib in the same order.

```sh
P=$(mktemp -d)
cat > "$P/render.py" <<'PY'
import runpy
import sys
from matplotlib.figure import Figure

mode, runner, *cli = sys.argv[1:]
savefig = Figure.savefig

def save_variant(self, *args, **kwargs):
    if mode == "no-tight":
        kwargs["bbox_inches"] = None
    return savefig(self, *args, **kwargs)

Figure.savefig = save_variant
sys.argv[:] = [runner, *cli]
runpy.run_path(runner, run_name="__main__")
PY
for sample in 1 2 3 4 5 6; do
  for mode in tight no-tight; do
    PYTHONPATH="$PWD" OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
      /usr/bin/time -f "$mode %e %M" -a -o "$P/timings.tsv" \
      python "$P/render.py" "$mode" "$PWD/phykit-runner.py" --threads 1 phylo_gwas \
      --alignment benchmarks/fixtures/phylo_gwas_multigroup.fa \
      --phenotype benchmarks/fixtures/phylo_gwas_multigroup.tsv \
      --output "$P/$mode.png" --json > "$P/$mode.json"
  done
done
cmp "$P/tight.json" "$P/no-tight.json"
```

For the PNG dimensions, open both outputs with `PIL.Image.open` and
inspect `.size`.

**Covarying rates imports.** Separately running the bundled command with
`python -X importtime` attributed about 104 ms cumulative import time to
`Bio.Phylo`, including 30 ms to its Nexus reader and 59 ms to its PhyloXML
reader. These are import-tree measurements and fluctuate between runs.
Importing only a Biopython Newick submodule still initializes its parent
package. A follow-up would compare the existing compact PhyKIT Newick
reader on ordinary files against Biopython on feature-rich fixtures,
including metadata and rooting, before considering a parser substitution.

## Expanded CLI function coverage

The bundled harness can capture a separate profile for each command.
A pass across its 30 cases observed 355 of 3,114 named definitions on
the optimized checkout at `c6db9a42`. The exact function inventory can
be regenerated as TSV. This records whether a definition was observed
and which CLI profiles observed it.

```sh
P=$(mktemp -d)
env PHYKIT_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 \
  python benchmarks/profile_commands.py --repeat 1 --warmup 0 --profile-dir "$P"
python benchmarks/profile_function_coverage.py "$P"/*.prof > "$P/functions.tsv"
```

## Native profile of phylogenetic signal

A separate one-thread py-spy profile repeated the 1,024-tip CLI command
seven times in one process and excluded output display. Native BLAS
matrix multiplication and triangular solve accounted for 102 of 200
samples. In 81 of those the stack included `_blombergs_k`,
`_inverse_from_spd_matrix`, and SciPy's `cho_solve`. This ranks ahead
of small Newick parser changes for these inputs.

[`_blombergs_k`](../phykit/services/tree/phylogenetic_signal.py#L515)
computes `C_inv = _inverse_from_spd_matrix(C)` before applying it to
trait vectors. A candidate is to factor the covariance once and solve
for the observed and permuted vectors in batches. It can avoid solving
for every column of the identity matrix. Factor reuse might also avoid
separate work in the R² calculation. This affects statistical outputs,
so compare K, permutation p-values and R² against the current path on
near-singular covariance matrices and permutation ties. Start with
100 versus 1,000 permutations because the solve-to-inverse cost balance
depends on the number of right-hand sides.

To reproduce the sample, generate the same tree and traits with
`make_inputs(Path(P), 1024, 1000)` as above, then run:

```sh
cat > "$P/repeat_signal.py" <<'PY'
import contextlib
import os
import sys
from phykit.phykit import Phykit

command = sys.argv[1:]
with open(os.devnull, "w") as sink, contextlib.redirect_stdout(sink):
    for _ in range(7):
        sys.argv[:] = ["phykit", *command]
        Phykit()
PY
env PHYKIT_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH="$PWD" \
  py-spy record --native --rate 100 --format speedscope -o "$P/signal.json" -- \
  python "$P/repeat_signal.py" phylogenetic_signal -t "$P/tree.nwk" \
    -d "$P/traits.tsv" -p 100
```

### Covariance-solve prototype outcome

A scratch prototype used a Cholesky factor and triangular solves to
whiten the observed and permuted trait vectors in batches of 128.
It retained the fixed permutation seed. At 1,024 taxa, with five
runs per variant, the median time in the K calculation fell from
0.138 to 0.025 s for 100 permutations and from 0.212 to 0.087 s for
1,000 permutations. At 128 taxa and 1,000 permutations the prototype
was slightly slower (0.0087 vs 0.0081 s).

The prototype changed permutation p-values from 0.6 to 1.0 for eight
exchangeable tips with repeated traits and from 0.993 to 0.005 for
a near-singular covariance. The latter also materially changed K.
The current production calculation was retained. A separate numerical
investigation with an independent reference is needed to choose a
tie and conditioning policy before deploying a new formula. The
synthetic cases can be reproduced with the following covariance and
trait pairs. The signed differences around ties matter because the
existing count checks `permutations >= observed` exactly.

```python
import numpy as np

repeated = (np.eye(8) + 0.01, np.array([0.0] * 4 + [1.0] * 4))
near_singular = (np.full((8, 8), 1.0) + 1e-9 * np.eye(8), np.arange(8.0))
```

To reproduce function timings and the negative correctness findings:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH="$PWD" \
  python benchmarks/experiments/signal_covariance_solve.py --repeat 5
```

This was a function-level prototype. CLI runtime and peak memory were not
measured for a production integration.

## Large topology-autocorrelation profile

The existing checked benchmark exercises 1,000, 10,000, and 100,000
classified genes with ten distance bins. Across both profiler and
allocation instrumentation, `validate_rows` accounted for 69.73 s
cumulative and `chromosome_marks` for 7.36 s. These elapsed times
include tracing overhead. The 100,000-gene pair totals were 899,955
and 949,905,000 at the two distance cutoffs, with the same 48.6 MiB peak
additional traced allocation. The algorithm aggregates range counts
without allocating a record per pair.

[`uncertainty`](../phykit/helpers/topology_autocorrelation.py#L243)
calls `_, all_marks = analyze(rows, edges, block_size)` for every
block size. Each call normalizes, validates and groups the same gene
rows, then builds point-estimate tables that uncertainty discards.
This creates a candidate to validate and group rows once and compute
only the required marks for each block size. Confirmation needs both
interval and diagnostic output comparisons because group order affects
random seed assignment.

```sh
P=$(mktemp -d)
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH="$PWD" \
  python -m cProfile -o "$P/topology.prof" \
  benchmarks/benchmark_topology_autocorrelation.py --sizes 1000 10000 100000 \
  --maxima 1000 1000000 --output "$P/result.json"
python -m pstats "$P/topology.prof"
```

## Validated prototype: reuse topology uncertainty inputs

The uncertainty routine validates and groups rows once, and computes
only the chromosome marks needed for each block size. It retains the
same group order, mark calculation, and seed derivation.

Two alternating rounds with 50,000 classified genes across two
chromosomes, three block sizes, 499 replicates, one thread and five
samples per variant after warmup:

| Metric | Control rounds | Prototype rounds | Speedup |
| --- | --- | --- | --- |
| CLI median | 2.00, 2.00 s | 1.77, 1.68 s | 1.16× |
| CLI median peak RSS | 111.2, 111.4 MiB | 108.6, 110.3 MiB | |
| Uncertainty function median | 0.644, 0.647 s | 0.405, 0.348 s | 1.72× |

CLI samples ranged from 1.94 to 2.10 s in the control and 1.67 to
1.80 s in the prototype. Function samples ranged from 0.633 to
0.659 s and 0.347 to 0.410 s. The function met the 1.5× target.
The CLI still reads, validates, and writes the full classified
dataset, so the end-to-end improvement is smaller.

Create the same input for both checkouts:

```sh
P=$(mktemp -d)
python - "$P/genes.tsv" <<'PY'
import random
import sys

rng = random.Random(41)
with open(sys.argv[1], "w") as out:
    out.write("reference\tchromosome\tgene_id\tstart\tend\tanchor\tclassification\n")
    for i in range(50_000):
        chrom = f"chr{i // 25_000 + 1}"
        anchor = i % 25_000 * 100
        category = rng.choice(("topology_1", "topology_2", "topology_3"))
        out.write(f"ref\t{chrom}\tg{i}\t{anchor}\t{anchor+1}\t{anchor}\t{category}\n")
PY
```

Run from each checkout, alternating control at `f08bb31e` and the
prototype. Keep the input and output paths fixed for byte comparisons.
Discard sample zero as warmup.

```sh
for sample in 0 1 2 3 4 5; do
  /usr/bin/time -f "$sample %e %M" -a -o "$P/timings.tsv" \
    python phykit-runner.py --threads 1 topology_autocorrelation \
    --classified-genes "$P/genes.tsv" --distance-edges 0 250 500 1000 \
    --block-sizes 10000 20000 40000 --replicates 499 \
    --output-prefix "$P/result" --json > "$P/command.json"
done
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python - "$P/genes.tsv" <<'PY'
import hashlib
import json
from statistics import median
import sys
from time import perf_counter
from phykit.helpers.topology_autocorrelation import read_classified, uncertainty

rows = read_classified(sys.argv[1])
times = []
for i in range(6):
    start = perf_counter()
    result = uncertainty(rows, [0, 250, 500, 1000], [10000, 20000, 40000], 499)
    if i:
        times.append(perf_counter() - start)
    digest = hashlib.sha256(json.dumps(result, sort_keys=True).encode()).hexdigest()
    if i == 0:
        expected = digest
    assert digest == expected
print(f"{median(times):.6f} [{min(times):.6f}, {max(times):.6f}]", expected)
PY
```

The digest of the uncertainty and diagnostics matched. On CPython
3.12.14, 3.13.15, 3.14.7, and 3.15.2, all 78 targeted tests passed.
The 50,000-gene CLI stdout, stderr, and all output file checksums
also matched within each Python version.

## Representative CLI coverage by package area

The bundled and scaled cases plus all 23 portable tutorial workloads
exercised 786 named definitions on the profiling branch. Each tutorial
completed and produced its expected artifacts. This table describes
function observation, rather than a count of functions individually
benchmarked or optimized. Native sampling and child-process checks are
reported above and are separate from this cProfile inventory.

| Package area | Observed | Definitions | Unobserved |
| --- | ---: | ---: | ---: |
| Package entry points and dispatch | 65 | 281 | 216 |
| Shared helpers | 105 | 245 | 140 |
| Base service and package initialization | 0 | 3 | 3 |
| Alignment services | 138 | 646 | 508 |
| Tree and trait services | 478 | 1,940 | 1,462 |
| Total | 786 | 3,115 | 2,329 |

The coverage tool enumerates definitions across all 150 Python modules
in `phykit/`, including nested functions. A blank profile column in
the generated TSV identifies a definition that needs a workload
or a child-process capture. A nonblank column supplies the exact
capture name for investigation. This includes low-cost helpers where
a 1.5× microbenchmark gain may have no detectable effect on total
runtime.

Run this from the optimized checkout with an active benchmark
environment. Each command type writes its own profiles outside Git.

```sh
P=$(mktemp -d)
export PYTHONPATH="$PWD" PHYKIT_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
python benchmarks/profile_commands.py --repeat 1 --warmup 0 --profile-dir "$P/bundled"
python benchmarks/profile_scaled_workloads.py --repeat 1 --threads 1 --profile-dir "$P/scaled"
mkdir "$P/tutorials"
cat > "$P/profile-python" <<'SH'
#!/bin/sh
profile=$(mktemp "$PHY_PROFILE_DIR/${3:-command}.XXXXXX.prof") || exit
exec "$PHY_PROFILE_PYTHON" -m cProfile -o "$profile" "$@"
SH
chmod +x "$P/profile-python"
PHY_PROFILE_DIR="$P/tutorials" PHY_PROFILE_PYTHON="$(command -v python)" \
  python scripts/smoke_tutorials.py --python "$P/profile-python"
python benchmarks/profile_function_coverage.py \
  "$P"/bundled/*.prof "$P"/scaled/*.prof "$P"/tutorials/*.prof > "$P/functions.tsv"
```

## Unit-test profile and final verification

A separate unit-test cProfile capture on the unoptimized control
`f08bb31e` observed 2,832 of 3,113 functions. This test workload
covers many mock, compatibility, and error paths outside the
representative CLI cases. The counts refer to different source
revisions and should not be added. Reproduce the test-path inventory
from the control checkout with a Python 3.12 environment containing
all runtime and test dependencies and installed CLI entry points:

```sh
P=$(mktemp -d)
env PYTHONPATH="$PWD" OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PHYKIT_THREADS=1 \
  python -m cProfile -o "$P/unit.prof" -m pytest -q --basetemp="$P/pytest" \
  -m 'not (integration or validation)'
python benchmarks/profile_function_coverage.py "$P/unit.prof" > "$P/unit-functions.tsv"
```

The forced cap made three worker-count assertions and four UMAP
thread-configuration cases fail in that observation run. It completed
with 6,653 passing tests. The seven cases passed when rerun on the
control with the default PhyKIT thread policy and no profiler. The
profile is used solely for function observation. The test pass
used for correctness verification was separate.

The optimized code at `bfc099af` passed all 7,537 unit and integration
tests on CPython 3.12.14. This run selected all tests except the separate
validation group. The changed modules' targeted suites and CLI output
comparisons were also verified on CPython 3.13.15, 3.14.7, and 3.15.2,
as detailed in each prototype section.

```sh
P=$(mktemp -d)
env PYTHONPATH="$PWD" OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  python -m pytest -q --basetemp="$P/pytest" -m 'not validation'
```

The five independent topology autocorrelation simulation checks also
passed on CPython 3.12.14:

```sh
python -m pytest -q tests/validation/test_topology_autocorrelation_simulation.py
```
