# Relative Rate Test Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add Tajima's relative rate test to PhyKIT, supporting single-alignment and batch modes with multiple testing correction, validated against R's `pegas::rr.test()`.

**Architecture:** A single `RelativeRateTest` service class inheriting from `Tree` base class, accepting `-a` (single alignment) or `-l` (alignment list) plus `-t` (rooted tree). The outgroup is inferred from the rooted tree. All C(n-1, 2) ingroup pairs are tested per alignment. Batch mode aggregates results across genes.

**Tech Stack:** BioPython (alignment/tree I/O), scipy.stats.chi2 (p-values), existing PhyKIT helpers

---

### Task 1: Create test alignment and tree fixtures

**Files:**
- Create: `tests/sample_files/rrt_test.fa`
- Create: `tests/sample_files/rrt_test.tre`
- Create: `tests/sample_files/rrt_test_batch_1.fa`
- Create: `tests/sample_files/rrt_test_batch_2.fa`
- Create: `tests/sample_files/rrt_test_batch_list.txt`

**Step 1: Create a 4-taxon nucleotide alignment with known m1/m2 counts**

File `tests/sample_files/rrt_test.fa`:
```fasta
>A
AACGTACGTAACGTAACGT
>B
AACGTACGTAACGTAACGT
>C
AATGTACGTAATGTAACGT
>O
AACGTACGTAACGTAACGT
```

This alignment has 4 taxa (A, B, C, O) with 19 sites. O is the outgroup.
- A vs B with outgroup O: A matches O everywhere, B matches O everywhere → m1=0, m2=0, chi2=0, p=1.0
- A vs C with outgroup O: sites 3 and 12 have C≠O but A=O → m1=0, m2=2, chi2=2.0, p=0.1573

We need a more interesting test alignment. We'll create the exact alignment in Step 3 after computing reference values with R.

**Step 2: Create a rooted tree**

File `tests/sample_files/rrt_test.tre`:
```
((A:0.1,B:0.1):0.05,(C:0.15,O:0.2):0.05);
```

Wait — for outgroup inference from a rooted tree, the root should separate the outgroup. A better rooted tree:
```
(((A:0.1,B:0.1):0.05,C:0.15):0.1,O:0.2);
```

This tree is rooted such that O is the outgroup (earliest diverging = on the other side of the root from the largest clade {A, B, C}).

**Step 3: Design test alignment with known properties and validate with R**

Create file `tests/sample_files/rrt_test.fa` with this carefully designed alignment:
```fasta
>A
ACGTACGTAAACGTACGTACGT
>B
ACTTACGTAAACGTTCGTACGT
>C
ACGTACCTAAATGTACGTACGT
>O
ACGTACGTAAACGTACGTACGT
```

21 sites. Comparing to outgroup O:
- A: identical to O at all 21 sites → 0 unique changes
- B: differs from O at positions 3(G→T), 15(A→T) → sites where B≠O
  - Position 3: B=T, O=G, A=G, C=G → B differs, A matches O → m2 for (A,B)
  - Position 15: B=T, O=A, A=A, C=A → B differs, A matches O → m2 for (A,B)
- C: differs from O at positions 8(G→C), 12(C→T) → sites where C≠O
  - Position 8: C=C, O=G, A=G, B=G → C differs, A matches O
  - Position 12: C=T, O=C, A=C, B=C → C differs, A matches O

For pair (A, B) with outgroup O: m1=0, m2=2 → chi2 = 4/2 = 2.0, p = 0.1573
For pair (A, C) with outgroup O: m1=0, m2=2 → chi2 = 4/2 = 2.0, p = 0.1573
For pair (B, C) with outgroup O: m1=2, m2=2 → chi2 = 0/4 = 0.0, p = 1.0

These are simple enough to validate by hand and against R.

Actually, let's use a simpler, more carefully designed alignment. We'll generate R reference values during implementation. For the plan, use this:

```fasta
>A
ACGTACGTACGTACGTACGT
>B
ACGTACGTACTTACTTACGT
>C
ACTTACGTACGTACGTACGT
>O
ACGTACGTACGTACGTACGT
```

20 sites. A=O everywhere. B differs from O at positions 11,12,15,16. C differs from O at positions 3,4.

Pair (A,B) outgroup O: m1=0 (A never differs from O uniquely), m2=4 → chi2=16/4=4.0, p≈0.0455
Pair (A,C) outgroup O: m1=0, m2=2 → chi2=4/2=2.0, p≈0.1573
Pair (B,C) outgroup O: m1=4 (B differs from O but C doesn't at pos 11,12,15,16), m2=2 (C differs from O but B doesn't at pos 3,4) → chi2=(4-2)^2/(4+2)=4/6≈0.667, p≈0.4142

**Step 4: Create batch test files**

File `tests/sample_files/rrt_test_batch_1.fa` — same as `rrt_test.fa`.

File `tests/sample_files/rrt_test_batch_2.fa` — a different alignment (could reuse same content or vary it).

File `tests/sample_files/rrt_test_batch_list.txt`:
```
tests/sample_files/rrt_test_batch_1.fa
tests/sample_files/rrt_test_batch_2.fa
```

**Step 5: Commit**

```
git add tests/sample_files/rrt_test*
git commit -m "test: add sample files for relative rate test"
```

---

### Task 2: Write core algorithm unit tests

**Files:**
- Create: `tests/unit/services/tree/test_relative_rate_test.py`

**Step 1: Write failing unit tests for the core Tajima computation**

```python
import pytest
from phykit.services.tree.relative_rate_test import RelativeRateTest


class TestTajimaTest:
    def test_equal_sequences(self):
        """Identical ingroup sequences: m1=m2=0, chi2=0, p=1.0."""
        result = RelativeRateTest._tajima_test(
            seq_x="ACGTACGT",
            seq_y="ACGTACGT",
            seq_out="ACGTACGT",
        )
        assert result["m1"] == 0
        assert result["m2"] == 0
        assert result["chi2"] == 0.0
        assert result["p_value"] == 1.0

    def test_one_lineage_faster(self):
        """One lineage has unique changes, the other doesn't."""
        # x matches outgroup, y differs at 2 sites
        result = RelativeRateTest._tajima_test(
            seq_x="ACGTACGT",
            seq_y="ACTTACTT",
            seq_out="ACGTACGT",
        )
        assert result["m1"] == 0
        assert result["m2"] == 2
        assert result["chi2"] == 2.0
        assert abs(result["p_value"] - 0.1573) < 0.001

    def test_symmetric_changes(self):
        """Equal unique changes on both lineages: chi2=0, p=1."""
        result = RelativeRateTest._tajima_test(
            seq_x="TCGTACGT",
            seq_y="ACGTACTT",
            seq_out="ACGTACGT",
        )
        assert result["m1"] == 1
        assert result["m2"] == 1
        assert result["chi2"] == 0.0
        assert result["p_value"] == 1.0

    def test_highly_unequal_rates(self):
        """Many changes on one lineage: should reject clock."""
        result = RelativeRateTest._tajima_test(
            seq_x="ACGTACGTACGT",
            seq_y="TTTTTTTTACGT",
            seq_out="ACGTACGTACGT",
        )
        assert result["m1"] == 0
        assert result["m2"] == 8
        assert result["chi2"] == 8.0
        assert result["p_value"] < 0.01

    def test_gaps_are_skipped(self):
        """Sites with gaps in any sequence should be excluded."""
        result = RelativeRateTest._tajima_test(
            seq_x="A-GTACGT",
            seq_y="ACTTACGT",
            seq_out="ACGTACGT",
        )
        # Position 1 has gap in x, skip it. Position 2: y=T, out=G, x=G → m2.
        # Only non-gap informative site where y differs: pos 2 (y=T, out=G, x=G)
        assert result["m1"] == 0
        assert result["m2"] == 1

    def test_shared_changes_not_counted(self):
        """Sites where BOTH ingroup differ from outgroup are uninformative."""
        result = RelativeRateTest._tajima_test(
            seq_x="TCGTACGT",
            seq_y="TCGTACGT",
            seq_out="ACGTACGT",
        )
        # Position 0: both x and y differ from outgroup → shared, skip
        assert result["m1"] == 0
        assert result["m2"] == 0
        assert result["chi2"] == 0.0


class TestIdentifyOutgroup:
    def test_simple_rooted_tree(self):
        """Outgroup should be the earliest-diverging taxon."""
        from io import StringIO
        from Bio import Phylo

        tree = Phylo.read(StringIO("(((A:0.1,B:0.1):0.05,C:0.15):0.1,O:0.2);"), "newick")
        outgroup = RelativeRateTest._identify_outgroup(tree)
        assert outgroup == "O"

    def test_two_taxa_at_root(self):
        """With two clades at root, outgroup is the smaller one."""
        from io import StringIO
        from Bio import Phylo

        tree = Phylo.read(StringIO("((A:0.1,B:0.1):0.1,O:0.2);"), "newick")
        outgroup = RelativeRateTest._identify_outgroup(tree)
        assert outgroup == "O"


class TestMultipleTestingCorrection:
    def test_bonferroni(self):
        p_values = [0.01, 0.04, 0.06]
        corrected = RelativeRateTest._bonferroni(p_values)
        assert corrected[0] == 0.03  # 0.01 * 3
        assert corrected[1] == 0.12  # 0.04 * 3
        assert corrected[2] == 0.18  # 0.06 * 3

    def test_bonferroni_capped_at_one(self):
        p_values = [0.5, 0.8]
        corrected = RelativeRateTest._bonferroni(p_values)
        assert corrected[0] == 1.0
        assert corrected[1] == 1.0

    def test_fdr(self):
        p_values = [0.01, 0.04, 0.06]
        corrected = RelativeRateTest._fdr(p_values)
        # BH-FDR: rank 1: 0.01*3/1=0.03, rank 2: 0.04*3/2=0.06, rank 3: 0.06*3/3=0.06
        assert abs(corrected[0] - 0.03) < 1e-10
        assert abs(corrected[1] - 0.06) < 1e-10
        assert abs(corrected[2] - 0.06) < 1e-10
```

**Step 2: Run tests to verify they fail**

Run: `python -m pytest tests/unit/services/tree/test_relative_rate_test.py -v`
Expected: FAIL (ImportError — module doesn't exist yet)

**Step 3: Commit**

```
git add tests/unit/services/tree/test_relative_rate_test.py
git commit -m "test: add unit tests for relative rate test (red)"
```

---

### Task 3: Implement core service class

**Files:**
- Create: `phykit/services/tree/relative_rate_test.py`

**Step 1: Implement the RelativeRateTest service**

```python
import itertools
from pathlib import Path
from typing import Dict, List, Tuple

from Bio import AlignIO, Phylo

from .base import Tree
from ...errors import PhykitUserError
from ...helpers.files import get_alignment_and_format
from ...helpers.json_output import print_json


def _chi2_sf(x, df):
    """Chi-squared survival function (lazy scipy import)."""
    from scipy.stats import chi2
    return float(chi2.sf(x, df))


class RelativeRateTest(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(
            tree_file_path=parsed["tree_file_path"],
            alignment_file_path=parsed.get("alignment_file_path"),
        )
        self.alignment_list = parsed.get("alignment_list")
        self.verbose = parsed["verbose"]
        self.json_output = parsed["json_output"]

    def process_args(self, args) -> Dict:
        return dict(
            tree_file_path=args.tree,
            alignment_file_path=getattr(args, "alignment", None),
            alignment_list=getattr(args, "alignment_list", None),
            verbose=getattr(args, "verbose", False),
            json_output=getattr(args, "json", False),
        )

    # ------------------------------------------------------------------
    # Core Tajima's test
    # ------------------------------------------------------------------

    @staticmethod
    def _tajima_test(seq_x: str, seq_y: str, seq_out: str) -> Dict:
        """Tajima's relative rate test (Equation 4).

        Counts unique substitutions on each ingroup lineage relative
        to the outgroup.  Returns m1, m2, chi-squared, and p-value.
        """
        m1 = 0  # unique changes on lineage x
        m2 = 0  # unique changes on lineage y
        skip_chars = {"-", "?", "N", "X", "n", "x"}

        for cx, cy, co in zip(seq_x, seq_y, seq_out):
            if cx in skip_chars or cy in skip_chars or co in skip_chars:
                continue
            x_diff = cx.upper() != co.upper()
            y_diff = cy.upper() != co.upper()
            if x_diff and not y_diff:
                m1 += 1
            elif y_diff and not x_diff:
                m2 += 1
            # Both differ or neither differs → uninformative

        if m1 + m2 == 0:
            return {"m1": m1, "m2": m2, "chi2": 0.0, "p_value": 1.0}

        chi2 = (m1 - m2) ** 2 / (m1 + m2)
        p_value = _chi2_sf(chi2, df=1)

        return {"m1": m1, "m2": m2, "chi2": chi2, "p_value": p_value}

    # ------------------------------------------------------------------
    # Outgroup identification
    # ------------------------------------------------------------------

    @staticmethod
    def _identify_outgroup(tree) -> str:
        """Identify outgroup as the earliest-diverging taxon in a rooted tree.

        The outgroup is the single taxon on the smaller side of the root split.
        """
        root = tree.root
        children = root.clades
        if len(children) < 2:
            raise PhykitUserError(
                ["Tree does not appear to be rooted (root has < 2 children)."],
                code=2,
            )

        # Find the child clade with fewest terminals
        clade_sizes = []
        for child in children:
            terminals = list(child.get_terminals())
            clade_sizes.append((len(terminals), terminals, child))

        clade_sizes.sort(key=lambda x: x[0])
        smallest_clade = clade_sizes[0]

        if smallest_clade[0] != 1:
            raise PhykitUserError(
                [
                    "Cannot determine a single outgroup taxon from the rooted tree.",
                    "The smallest clade at the root contains multiple taxa.",
                    "Please ensure the tree is rooted with a single outgroup taxon.",
                ],
                code=2,
            )

        return smallest_clade[1][0].name

    # ------------------------------------------------------------------
    # Multiple testing correction
    # ------------------------------------------------------------------

    @staticmethod
    def _bonferroni(p_values: List[float]) -> List[float]:
        n = len(p_values)
        return [min(p * n, 1.0) for p in p_values]

    @staticmethod
    def _fdr(p_values: List[float]) -> List[float]:
        """Benjamini-Hochberg FDR correction."""
        n = len(p_values)
        if n == 0:
            return []
        indexed = sorted(enumerate(p_values), key=lambda x: x[1])
        corrected = [0.0] * n
        prev = 1.0
        for rank_minus_1 in range(n - 1, -1, -1):
            orig_idx, p = indexed[rank_minus_1]
            rank = rank_minus_1 + 1
            adjusted = min(p * n / rank, prev)
            adjusted = min(adjusted, 1.0)
            corrected[orig_idx] = adjusted
            prev = adjusted
        return corrected

    # ------------------------------------------------------------------
    # Single alignment analysis
    # ------------------------------------------------------------------

    def _run_single(self, alignment_path: str, tree) -> List[Dict]:
        """Run Tajima's test for all ingroup pairs in one alignment."""
        alignment, _, _ = get_alignment_and_format(alignment_path)

        # Build sequence dict
        seq_dict = {}
        for record in alignment:
            seq_dict[record.id] = str(record.seq)

        outgroup = self._identify_outgroup(tree)
        if outgroup not in seq_dict:
            raise PhykitUserError(
                [
                    f"Outgroup taxon '{outgroup}' from tree not found in alignment.",
                    f"Alignment taxa: {sorted(seq_dict.keys())}",
                ],
                code=2,
            )

        seq_out = seq_dict[outgroup]
        ingroup = sorted(t for t in seq_dict if t != outgroup)

        if len(ingroup) < 2:
            raise PhykitUserError(
                ["At least 2 ingroup taxa are required (3 total including outgroup)."],
                code=2,
            )

        results = []
        for t1, t2 in itertools.combinations(ingroup, 2):
            test_result = self._tajima_test(seq_dict[t1], seq_dict[t2], seq_out)
            test_result["taxon1"] = t1
            test_result["taxon2"] = t2
            results.append(test_result)

        # Multiple testing correction
        p_values = [r["p_value"] for r in results]
        p_bonf = self._bonferroni(p_values)
        p_fdr = self._fdr(p_values)
        for i, r in enumerate(results):
            r["p_bonf"] = p_bonf[i]
            r["p_fdr"] = p_fdr[i]

        return results

    # ------------------------------------------------------------------
    # Main entry point
    # ------------------------------------------------------------------

    def run(self):
        tree = self.read_tree_file()

        outgroup = self._identify_outgroup(tree)

        if self.alignment_list:
            # Batch mode
            list_path = Path(self.alignment_list)
            if not list_path.exists():
                raise PhykitUserError(
                    [f"{self.alignment_list} not found."], code=2
                )
            aln_paths = [
                line.strip()
                for line in list_path.read_text().splitlines()
                if line.strip() and not line.strip().startswith("#")
            ]
            if not aln_paths:
                raise PhykitUserError(
                    ["No alignment paths found in list file."], code=2
                )

            # Run each alignment
            all_results = {}  # (t1, t2) -> list of per-gene results
            n_alignments = 0
            for aln_path in aln_paths:
                resolved = aln_path
                if not Path(aln_path).is_absolute():
                    resolved = str(list_path.parent / aln_path)
                try:
                    results = self._run_single(resolved, tree)
                    n_alignments += 1
                    for r in results:
                        key = (r["taxon1"], r["taxon2"])
                        all_results.setdefault(key, []).append(r)
                except Exception:
                    continue  # skip problematic alignments

            self._output_batch(outgroup, n_alignments, all_results)

        elif self.alignment_file_path:
            # Single mode
            results = self._run_single(self.alignment_file_path, tree)
            self._output_single(outgroup, results)
        else:
            raise PhykitUserError(
                ["Provide -a (single alignment) or -l (alignment list)."],
                code=2,
            )

    # ------------------------------------------------------------------
    # Output
    # ------------------------------------------------------------------

    def _output_single(self, outgroup: str, results: List[Dict]):
        n_ingroup = len(set(r["taxon1"] for r in results) | set(r["taxon2"] for r in results))

        if self.json_output:
            print_json(dict(
                outgroup=outgroup,
                n_ingroup_taxa=n_ingroup,
                n_tests=len(results),
                results=[
                    dict(
                        taxon1=r["taxon1"], taxon2=r["taxon2"],
                        m1=r["m1"], m2=r["m2"],
                        chi2=round(r["chi2"], 4),
                        p_value=r["p_value"],
                        p_bonferroni=r["p_bonf"],
                        p_fdr=r["p_fdr"],
                    )
                    for r in results
                ],
            ))
        else:
            print(f"Outgroup: {outgroup}")
            print(f"Number of ingroup taxa: {n_ingroup}")
            print(f"Number of pairwise tests: {len(results)}")
            print("---")
            print("taxon1\ttaxon2\tm1\tm2\tchi2\tp_value\tp_bonf\tp_fdr\tsignificant")
            for r in results:
                sig = "*" if r["p_fdr"] < 0.05 else ""
                print(
                    f"{r['taxon1']}\t{r['taxon2']}\t{r['m1']}\t{r['m2']}\t"
                    f"{r['chi2']:.4f}\t{r['p_value']:.4e}\t{r['p_bonf']:.4e}\t"
                    f"{r['p_fdr']:.4e}\t{sig}"
                )

    def _output_batch(self, outgroup: str, n_alignments: int, all_results: Dict):
        if self.json_output:
            pairs = []
            for (t1, t2), gene_results in sorted(all_results.items()):
                n_reject = sum(1 for r in gene_results if r["p_value"] < 0.05)
                chi2_values = [r["chi2"] for r in gene_results]
                chi2_values.sort()
                median_chi2 = chi2_values[len(chi2_values) // 2] if chi2_values else 0.0
                pairs.append(dict(
                    taxon1=t1, taxon2=t2,
                    n_reject=n_reject, n_total=len(gene_results),
                    pct_reject=round(100 * n_reject / len(gene_results), 1),
                    median_chi2=round(median_chi2, 4),
                ))
            print_json(dict(
                outgroup=outgroup,
                n_alignments=n_alignments,
                pairs=pairs,
            ))
        else:
            print(f"Number of alignments: {n_alignments}")
            print(f"Outgroup: {outgroup}")
            print("---")
            print("taxon1\ttaxon2\tn_reject\tn_total\tpct_reject\tmedian_chi2")
            for (t1, t2), gene_results in sorted(all_results.items()):
                n_reject = sum(1 for r in gene_results if r["p_value"] < 0.05)
                chi2_values = sorted(r["chi2"] for r in gene_results)
                median_chi2 = chi2_values[len(chi2_values) // 2] if chi2_values else 0.0
                pct = 100 * n_reject / len(gene_results)
                print(f"{t1}\t{t2}\t{n_reject}\t{len(gene_results)}\t{pct:.1f}\t{median_chi2:.4f}")
```

**Step 2: Run unit tests to verify they pass**

Run: `python -m pytest tests/unit/services/tree/test_relative_rate_test.py -v`
Expected: All tests PASS

**Step 3: Commit**

```
git add phykit/services/tree/relative_rate_test.py
git commit -m "feat: implement RelativeRateTest core service class"
```

---

### Task 4: Register the command in the PhyKIT framework

**Files:**
- Modify: `phykit/service_factories.py`
- Modify: `phykit/cli_registry.py`
- Modify: `phykit/services/tree/__init__.py`
- Modify: `phykit/phykit.py`
- Modify: `setup.py`

**Step 1: Add service factory**

In `phykit/service_factories.py`, add after the existing `QuartetNetwork` entry:
```python
RelativeRateTest = _LazyServiceFactory("phykit.services.tree.relative_rate_test", "RelativeRateTest")
```

**Step 2: Add CLI aliases**

In `phykit/cli_registry.py`, add in the Tree aliases section:
```python
"rrt": "relative_rate_test",
"tajima_rrt": "relative_rate_test",
```

**Step 3: Add to tree exports**

In `phykit/services/tree/__init__.py`, add to `_EXPORTS`:
```python
"RelativeRateTest": "relative_rate_test",
```

**Step 4: Add handler method in phykit.py**

Add the static handler method (near the other tree-based commands). Follow the `saturation` pattern with `-a/--alignment`, `-t/--tree`, plus the new `-l/--alignment-list`:

```python
@staticmethod
def relative_rate_test(argv):
    parser = _new_parser(
        description=textwrap.dedent(
            f"""\
            {help_header}

            Tajima's relative rate test.

            Tests whether two lineages have evolved at equal rates
            since diverging from their common ancestor. The outgroup
            is inferred from the rooted tree as the earliest-diverging
            taxon. All pairwise ingroup combinations are tested with
            Bonferroni and BH-FDR multiple testing correction.

            Single alignment mode:
            phykit relative_rate_test -a <alignment> -t <rooted_tree>

            Batch mode (multiple alignments, one shared tree):
            phykit relative_rate_test -l <alignment_list> -t <rooted_tree>

            Options
            =====================================================
            -a/--alignment              single alignment file

            -l/--alignment-list         file with one alignment path
                                        per line (batch mode)

            -t/--tree                   rooted tree file

            -v/--verbose                print additional detail

            --json                      output results as JSON
            """
        ),
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "-a", "--alignment", type=str, help=SUPPRESS, metavar=""
    )
    group.add_argument(
        "-l", "--alignment-list", type=str, help=SUPPRESS, metavar=""
    )
    parser.add_argument(
        "-t", "--tree", type=str, required=True, help=SUPPRESS, metavar=""
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", required=False, help=SUPPRESS
    )
    _add_json_argument(parser)
    _run_service(parser, argv, RelativeRateTest)
```

**Step 5: Add module-level wrapper function**

Near the other wrapper functions in `phykit.py`:
```python
def relative_rate_test(argv=None):
    Phykit.relative_rate_test(sys.argv[1:])
```

**Step 6: Add help text entry**

In the `Phykit.__init__` banner, add:
```
relative_rate_test (alias: rrt; tajima_rrt)
```

**Step 7: Add entry points in setup.py**

```python
"pk_relative_rate_test = phykit.phykit:relative_rate_test",
"pk_rrt = phykit.phykit:relative_rate_test",
"pk_tajima_rrt = phykit.phykit:relative_rate_test",
```

**Step 8: Run unit tests**

Run: `python -m pytest tests/unit/services/tree/test_relative_rate_test.py -v`
Expected: All PASS

**Step 9: Commit**

```
git add phykit/service_factories.py phykit/cli_registry.py phykit/services/tree/__init__.py phykit/phykit.py setup.py
git commit -m "feat: register relative_rate_test command in PhyKIT framework"
```

---

### Task 5: Write and run integration tests

**Files:**
- Create: `tests/integration/tree/test_relative_rate_test_integration.py`

**Step 1: Write integration tests**

```python
import json
import sys

import pytest
from mock import patch

from phykit.phykit import Phykit


@pytest.mark.integration
class TestRelativeRateTestIntegration:
    @patch("builtins.print")
    def test_single_alignment(self, mocked_print, tmp_path):
        aln = tmp_path / "test.fa"
        aln.write_text(
            ">A\nACGTACGTACGTACGTACGT\n"
            ">B\nACGTACGTACTTACTTACGT\n"
            ">C\nACTTACGTACGTACGTACGT\n"
            ">O\nACGTACGTACGTACGTACGT\n"
        )
        tree = tmp_path / "test.tre"
        tree.write_text("(((A:0.1,B:0.1):0.05,C:0.15):0.1,O:0.2);")

        testargs = ["phykit", "relative_rate_test", "-a", str(aln), "-t", str(tree)]
        with patch.object(sys, "argv", testargs):
            Phykit()

        output = "\n".join(str(call) for call in mocked_print.call_args_list)
        assert "Outgroup: O" in output
        assert "Number of pairwise tests: 3" in output

    @patch("builtins.print")
    def test_rrt_alias(self, mocked_print, tmp_path):
        aln = tmp_path / "test.fa"
        aln.write_text(
            ">A\nACGTACGT\n>B\nACGTACGT\n>O\nACGTACGT\n"
        )
        tree = tmp_path / "test.tre"
        tree.write_text("((A:0.1,B:0.1):0.1,O:0.2);")

        testargs = ["phykit", "rrt", "-a", str(aln), "-t", str(tree)]
        with patch.object(sys, "argv", testargs):
            Phykit()

        output = "\n".join(str(call) for call in mocked_print.call_args_list)
        assert "Outgroup: O" in output

    @patch("builtins.print")
    def test_json_output(self, mocked_print, tmp_path):
        aln = tmp_path / "test.fa"
        aln.write_text(
            ">A\nACGTACGT\n>B\nACTTACGT\n>O\nACGTACGT\n"
        )
        tree = tmp_path / "test.tre"
        tree.write_text("((A:0.1,B:0.1):0.1,O:0.2);")

        testargs = ["phykit", "rrt", "-a", str(aln), "-t", str(tree), "--json"]
        with patch.object(sys, "argv", testargs):
            Phykit()

        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["outgroup"] == "O"
        assert payload["n_tests"] == 1
        assert "chi2" in payload["results"][0]

    @patch("builtins.print")
    def test_batch_mode(self, mocked_print, tmp_path):
        aln_content = ">A\nACGTACGT\n>B\nACTTACGT\n>O\nACGTACGT\n"
        aln1 = tmp_path / "aln1.fa"
        aln1.write_text(aln_content)
        aln2 = tmp_path / "aln2.fa"
        aln2.write_text(aln_content)
        aln_list = tmp_path / "list.txt"
        aln_list.write_text(f"{aln1}\n{aln2}\n")
        tree = tmp_path / "test.tre"
        tree.write_text("((A:0.1,B:0.1):0.1,O:0.2);")

        testargs = ["phykit", "rrt", "-l", str(aln_list), "-t", str(tree)]
        with patch.object(sys, "argv", testargs):
            Phykit()

        output = "\n".join(str(call) for call in mocked_print.call_args_list)
        assert "Number of alignments: 2" in output
```

**Step 2: Run integration tests**

Run: `python -m pytest tests/integration/tree/test_relative_rate_test_integration.py -v`
Expected: All PASS

**Step 3: Commit**

```
git add tests/integration/tree/test_relative_rate_test_integration.py
git commit -m "test: add integration tests for relative_rate_test"
```

---

### Task 6: Validate against R's pegas::rr.test()

**Files:**
- Modify: `tests/unit/services/tree/test_relative_rate_test.py` (add R-validated tests)

**Step 1: Run R validation**

Create a test alignment, run `pegas::rr.test()` in R, and record exact chi2 and p-values.

```r
library(pegas)
library(ape)
# Create test sequences matching our Python test cases
# Write FASTA, read back, run rr.test, record values
```

**Step 2: Add R-validated test cases**

Add a test class `TestMatchesRReference` with hardcoded values from R.

**Step 3: Run all tests**

Run: `python -m pytest tests/unit/services/tree/test_relative_rate_test.py tests/integration/tree/test_relative_rate_test_integration.py -v`
Expected: All PASS

**Step 4: Commit**

```
git add tests/unit/services/tree/test_relative_rate_test.py
git commit -m "test: add R-validated reference test cases for relative rate test"
```

---

### Task 7: Add documentation, changelog, and version bump

**Files:**
- Modify: `docs/usage/index.rst`
- Modify: `docs/change_log/index.rst`
- Modify: `phykit/version.py`

**Step 1: Add usage documentation**

In `docs/usage/index.rst`, add a section for `relative_rate_test` following the pattern of other commands (like `saturation`). Include:
- Function names and aliases
- Description of the algorithm
- CLI usage for both single and batch mode
- Example output

**Step 2: Add changelog entry**

In `docs/change_log/index.rst`, add entry for the new version describing the feature.

**Step 3: Bump version**

In `phykit/version.py`, bump from current version.

**Step 4: Run full test suite**

Run: `python -m pytest tests/unit/services/tree/test_relative_rate_test.py tests/integration/tree/test_relative_rate_test_integration.py -v`
Expected: All PASS

**Step 5: Commit**

```
git add docs/usage/index.rst docs/change_log/index.rst phykit/version.py
git commit -m "docs: add usage docs and changelog for relative_rate_test"
```
