# Evolutionary Tempo Mapping Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add an `evo_tempo_map` command that detects rate-topology associations by comparing branch length distributions between concordant and discordant gene trees at each species tree branch.

**Architecture:** A single service class (`EvoTempoMap`) that reuses bipartition matching logic from `concordance_asr.py` to classify gene trees at each species tree branch, extracts the homologous branch length from each gene tree, and tests for differences using Mann-Whitney U + permutation tests. Includes BH-FDR correction across branches. Optional box/strip plot output.

**Tech Stack:** Bio.Phylo (tree parsing), numpy (array ops), scipy.stats (Mann-Whitney U), matplotlib (plotting)

---

### Task 1: Create the service file with argument parsing and validation

**Files:**
- Create: `phykit/services/tree/evo_tempo_map.py`
- Test: `tests/unit/services/tree/test_evo_tempo_map.py`

**Step 1: Write the failing test for argument parsing**

Create `tests/unit/services/tree/test_evo_tempo_map.py`:

```python
import pytest
from argparse import Namespace

TREE_SIMPLE = "tests/sample_files/tree_simple.tre"
GENE_TREES = "tests/sample_files/gene_trees_simple.nwk"


class TestProcessArgs:
    def test_default_args(self):
        from phykit.services.tree.evo_tempo_map import EvoTempoMap

        args = Namespace(
            tree=TREE_SIMPLE,
            gene_trees=GENE_TREES,
            verbose=False,
            json=False,
            plot_output=None,
        )
        svc = EvoTempoMap(args)
        assert svc.gene_trees_path == GENE_TREES
        assert svc.verbose is False
        assert svc.json_output is False
        assert svc.plot_output is None

    def test_verbose_arg(self):
        from phykit.services.tree.evo_tempo_map import EvoTempoMap

        args = Namespace(
            tree=TREE_SIMPLE,
            gene_trees=GENE_TREES,
            verbose=True,
            json=False,
            plot_output=None,
        )
        svc = EvoTempoMap(args)
        assert svc.verbose is True
```

**Step 2: Run test to verify it fails**

Run: `python -m pytest tests/unit/services/tree/test_evo_tempo_map.py::TestProcessArgs -v`
Expected: FAIL with `ModuleNotFoundError` or `ImportError`

**Step 3: Write minimal implementation**

Create `phykit/services/tree/evo_tempo_map.py`:

```python
from io import StringIO
from pathlib import Path
from typing import Dict, List, Tuple

from Bio import Phylo

from .base import Tree
from ...helpers.json_output import print_json
from ...errors import PhykitUserError


class EvoTempoMap(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"])
        self.gene_trees_path = parsed["gene_trees_path"]
        self.verbose = parsed["verbose"]
        self.json_output = parsed["json_output"]
        self.plot_output = parsed["plot_output"]

    def process_args(self, args) -> Dict:
        return dict(
            tree_file_path=args.tree,
            gene_trees_path=args.gene_trees,
            verbose=getattr(args, "verbose", False),
            json_output=getattr(args, "json", False),
            plot_output=getattr(args, "plot_output", None),
        )

    def run(self) -> None:
        pass  # Placeholder
```

**Step 4: Run test to verify it passes**

Run: `python -m pytest tests/unit/services/tree/test_evo_tempo_map.py::TestProcessArgs -v`
Expected: PASS

**Step 5: Commit**

```bash
git add phykit/services/tree/evo_tempo_map.py tests/unit/services/tree/test_evo_tempo_map.py
git commit -m "feat(evo_tempo_map): scaffold service with argument parsing"
```

---

### Task 2: Implement gene tree parsing and validation

This reuses the same parsing pattern from `concordance_asr.py:125-145` and the same branch-length validation from `concordance_asr.py:63-73`.

**Files:**
- Modify: `phykit/services/tree/evo_tempo_map.py`
- Modify: `tests/unit/services/tree/test_evo_tempo_map.py`

**Step 1: Write the failing tests**

Add to `tests/unit/services/tree/test_evo_tempo_map.py`:

```python
import os
import tempfile


class TestParseGeneTrees:
    def test_parses_multi_newick(self):
        from phykit.services.tree.evo_tempo_map import EvoTempoMap

        args = Namespace(
            tree=TREE_SIMPLE, gene_trees=GENE_TREES,
            verbose=False, json=False, plot_output=None,
        )
        svc = EvoTempoMap(args)
        trees = svc._parse_gene_trees(GENE_TREES)
        assert len(trees) == 10
        # Each tree should have 8 tips
        for gt in trees:
            tips = [t.name for t in gt.get_terminals()]
            assert len(tips) == 8

    def test_file_not_found_error(self):
        from phykit.services.tree.evo_tempo_map import EvoTempoMap

        args = Namespace(
            tree=TREE_SIMPLE, gene_trees="nonexistent.nwk",
            verbose=False, json=False, plot_output=None,
        )
        svc = EvoTempoMap(args)
        with pytest.raises(PhykitUserError):
            svc._parse_gene_trees("nonexistent.nwk")

    def test_validates_branch_lengths(self):
        from phykit.services.tree.evo_tempo_map import EvoTempoMap

        # Create a gene tree file without branch lengths
        no_bl = "((A,B),(C,D));\n"
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".nwk", delete=False
        ) as f:
            f.write(no_bl)
            tmp_path = f.name
        try:
            args = Namespace(
                tree=TREE_SIMPLE, gene_trees=tmp_path,
                verbose=False, json=False, plot_output=None,
            )
            svc = EvoTempoMap(args)
            with pytest.raises(PhykitUserError, match="branch lengths"):
                svc._parse_and_validate_gene_trees()
        finally:
            os.unlink(tmp_path)

    def test_requires_at_least_2_gene_trees(self):
        from phykit.services.tree.evo_tempo_map import EvoTempoMap

        single = "((raccoon:19.2,bear:6.8):0.85,((sea_lion:12.0,seal:12.0):7.5,((monkey:100.9,cat:47.1):20.6,weasel:18.9):2.1):3.9,dog:25.5);\n"
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".nwk", delete=False
        ) as f:
            f.write(single)
            tmp_path = f.name
        try:
            args = Namespace(
                tree=TREE_SIMPLE, gene_trees=tmp_path,
                verbose=False, json=False, plot_output=None,
            )
            svc = EvoTempoMap(args)
            with pytest.raises(PhykitUserError, match="At least 2"):
                svc._parse_and_validate_gene_trees()
        finally:
            os.unlink(tmp_path)
```

Also add at the top of the test file:

```python
from phykit.errors import PhykitUserError
```

**Step 2: Run tests to verify they fail**

Run: `python -m pytest tests/unit/services/tree/test_evo_tempo_map.py::TestParseGeneTrees -v`
Expected: FAIL (methods don't exist yet)

**Step 3: Implement gene tree parsing and validation**

Add to `phykit/services/tree/evo_tempo_map.py`:

```python
    def _parse_gene_trees(self, path: str) -> list:
        """Parse multi-Newick file into list of Bio.Phylo trees."""
        try:
            lines = Path(path).read_text().splitlines()
        except FileNotFoundError:
            raise PhykitUserError(
                [
                    f"{path} corresponds to no such file or directory.",
                    "Please check filename and pathing",
                ],
                code=2,
            )
        cleaned = [l.strip() for l in lines if l.strip() and not l.strip().startswith("#")]
        trees = []
        for line in cleaned:
            if line.startswith("("):
                trees.append(Phylo.read(StringIO(line), "newick"))
            else:
                tree_path = Path(path).parent / line
                trees.append(Phylo.read(str(tree_path), "newick"))
        return trees

    def _parse_and_validate_gene_trees(self) -> list:
        """Parse gene trees and validate they have branch lengths."""
        gene_trees = self._parse_gene_trees(self.gene_trees_path)
        if len(gene_trees) < 2:
            raise PhykitUserError(
                [
                    "At least 2 gene trees are required for evolutionary tempo mapping.",
                    f"Only {len(gene_trees)} gene tree(s) provided.",
                ],
                code=2,
            )
        for i, gt in enumerate(gene_trees):
            for clade in gt.find_clades():
                if clade.branch_length is None and clade != gt.root:
                    raise PhykitUserError(
                        [
                            f"Gene tree {i + 1} has branches without lengths.",
                            "All gene trees must have branch lengths for tempo mapping.",
                        ],
                        code=2,
                    )
        return gene_trees
```

**Step 4: Run tests to verify they pass**

Run: `python -m pytest tests/unit/services/tree/test_evo_tempo_map.py::TestParseGeneTrees -v`
Expected: PASS

**Step 5: Commit**

```bash
git add phykit/services/tree/evo_tempo_map.py tests/unit/services/tree/test_evo_tempo_map.py
git commit -m "feat(evo_tempo_map): add gene tree parsing and validation"
```

---

### Task 3: Implement bipartition extraction and concordance classification

This is the core algorithm. For each species tree branch, classify each gene tree as concordant, discordant, or uninformative, and extract the homologous branch length.

Reuses bipartition logic from `concordance_asr.py:210-326` (`_canonical_split`, `_get_four_groups`, bipartition extraction loop).

**Files:**
- Modify: `phykit/services/tree/evo_tempo_map.py`
- Modify: `tests/unit/services/tree/test_evo_tempo_map.py`

**Step 1: Write the failing tests**

Add to the test file:

```python
class TestClassifyGeneTrees:
    """Test bipartition-based concordance classification."""

    @pytest.fixture
    def svc(self):
        from phykit.services.tree.evo_tempo_map import EvoTempoMap
        args = Namespace(
            tree=TREE_SIMPLE, gene_trees=GENE_TREES,
            verbose=False, json=False, plot_output=None,
        )
        return EvoTempoMap(args)

    def test_classify_returns_dict_per_branch(self, svc):
        species_tree = svc.read_tree_file()
        gene_trees = svc._parse_gene_trees(GENE_TREES)
        result = svc._classify_gene_trees(species_tree, gene_trees)
        # Should have entries for internal branches (not root, not tips)
        assert len(result) > 0
        for branch_label, data in result.items():
            assert "concordant_lengths" in data
            assert "discordant_lengths" in data
            assert "n_concordant" in data
            assert "n_discordant" in data
            assert "split" in data

    def test_concordant_plus_discordant_leq_total(self, svc):
        species_tree = svc.read_tree_file()
        gene_trees = svc._parse_gene_trees(GENE_TREES)
        result = svc._classify_gene_trees(species_tree, gene_trees)
        n_gene_trees = len(gene_trees)
        for branch_label, data in result.items():
            total = data["n_concordant"] + data["n_discordant"]
            assert total <= n_gene_trees

    def test_all_concordant_when_identical_trees(self, svc):
        """When all gene trees match species tree topology,
        all should be concordant at every branch."""
        species_tree = svc.read_tree_file()
        # Use species tree as gene tree 3 times
        identical_trees = [svc.read_tree_file() for _ in range(3)]
        result = svc._classify_gene_trees(species_tree, identical_trees)
        for branch_label, data in result.items():
            assert data["n_concordant"] == 3
            assert data["n_discordant"] == 0

    def test_branch_lengths_are_positive(self, svc):
        species_tree = svc.read_tree_file()
        gene_trees = svc._parse_gene_trees(GENE_TREES)
        result = svc._classify_gene_trees(species_tree, gene_trees)
        for branch_label, data in result.items():
            for bl in data["concordant_lengths"]:
                assert bl > 0
            for bl in data["discordant_lengths"]:
                assert bl > 0

    def test_cross_check_with_concordance_asr_gcf(self, svc):
        """The number of concordant gene trees at each branch should match
        the gCF numerator from concordance_asr's bipartition matching."""
        from phykit.services.tree.concordance_asr import ConcordanceAsr

        species_tree = svc.read_tree_file()
        gene_trees = svc._parse_gene_trees(GENE_TREES)
        etm_result = svc._classify_gene_trees(species_tree, gene_trees)

        # Get gCF from concordance_asr
        casr_args = Namespace(
            tree=TREE_SIMPLE, gene_trees=GENE_TREES,
            trait_data="tests/sample_files/tree_simple_traits.tsv",
            trait=None, method="weighted", ci=False, plot=None,
            missing_taxa="shared", json=False,
        )
        casr = ConcordanceAsr(casr_args)
        casr_species_tree = casr.read_tree_file()
        casr_gene_trees = casr._parse_gene_trees(GENE_TREES)
        all_taxa = sorted(set(
            t.name for t in casr_species_tree.get_terminals()
        ))
        gcf_result = casr._compute_gcf_per_node(
            casr_species_tree, casr_gene_trees, all_taxa
        )

        # For each species tree branch, the concordant count from ETM
        # should match gCF * total from concordance_asr
        # We match by bipartition (split set)
        for branch_label, etm_data in etm_result.items():
            split_set = frozenset(etm_data["split"])
            n_conc = etm_data["n_concordant"]
            n_disc = etm_data["n_discordant"]
            total = n_conc + n_disc
            if total == 0:
                continue
            etm_gcf = n_conc / total
            # Find matching gCF entry by checking node descendants
            matched = False
            for node_id, (gcf, gdf1, gdf2) in gcf_result.items():
                # Reconstruct the split for this node
                for clade in casr_species_tree.find_clades():
                    if id(clade) == node_id:
                        node_tips = frozenset(
                            t.name for t in clade.get_terminals()
                        )
                        if node_tips == split_set:
                            assert abs(etm_gcf - gcf) < 0.01, (
                                f"gCF mismatch at {split_set}: "
                                f"ETM={etm_gcf:.3f}, concordance_asr={gcf:.3f}"
                            )
                            matched = True
                            break
                if matched:
                    break
```

**Step 2: Run tests to verify they fail**

Run: `python -m pytest tests/unit/services/tree/test_evo_tempo_map.py::TestClassifyGeneTrees -v`
Expected: FAIL (method `_classify_gene_trees` doesn't exist)

**Step 3: Implement the classification method**

Add to `phykit/services/tree/evo_tempo_map.py`:

```python
    @staticmethod
    def _canonical_split(taxa_side, all_taxa):
        """Normalize a bipartition to canonical form.

        Same algorithm as concordance_asr._canonical_split.
        """
        complement = all_taxa - taxa_side
        if len(taxa_side) < len(complement):
            return frozenset(taxa_side)
        elif len(taxa_side) > len(complement):
            return frozenset(complement)
        else:
            return min(
                frozenset(taxa_side), frozenset(complement),
                key=lambda s: sorted(s),
            )

    @staticmethod
    def _build_parent_map(tree) -> Dict:
        parent_map = {}
        for clade in tree.find_clades(order="preorder"):
            for child in clade.clades:
                parent_map[id(child)] = clade
        return parent_map

    def _get_four_groups(self, tree, node, parent_map, all_taxa_fs):
        """Identify the four subtree groups around an internal branch.

        Same algorithm as concordance_asr._get_four_groups.
        Returns (C1, C2, S, D) as frozensets, or None.
        """
        if node.is_terminal() or len(node.clades) < 2:
            return None

        C1 = frozenset(t.name for t in node.clades[0].get_terminals())
        C2 = frozenset(t.name for t in node.clades[1].get_terminals())
        for extra_child in node.clades[2:]:
            C2 = C2 | frozenset(t.name for t in extra_child.get_terminals())

        parent = parent_map.get(id(node))
        if parent is None:
            return None

        siblings = [c for c in parent.clades if id(c) != id(node)]
        if not siblings:
            return None

        S = frozenset(t.name for t in siblings[0].get_terminals())
        D = all_taxa_fs - C1 - C2 - S

        return C1, C2, S, D

    def _extract_bipartitions_with_lengths(self, tree, all_taxa_fs):
        """Extract all non-trivial bipartitions from a tree along with
        the branch length that defines each split.

        Returns dict mapping canonical_split -> branch_length.
        """
        bp_to_length = {}
        for clade in tree.get_nonterminals():
            tips = frozenset(t.name for t in clade.get_terminals())
            if len(tips) <= 1 or tips == all_taxa_fs:
                continue
            bp = self._canonical_split(tips, all_taxa_fs)
            bl = clade.branch_length if clade.branch_length else 0.0
            bp_to_length[bp] = bl
        return bp_to_length

    def _classify_gene_trees(
        self, species_tree, gene_trees
    ) -> Dict[str, Dict]:
        """For each internal branch of the species tree, classify gene
        trees as concordant or discordant and extract homologous branch
        lengths.

        Returns a dict keyed by branch label (sorted taxa string) with:
          split: list of taxon names in the smaller partition side
          n_concordant: int
          n_discordant: int
          concordant_lengths: list of floats
          discordant_lengths: list of floats
        """
        all_taxa = sorted(
            set(t.name for t in species_tree.get_terminals())
            & set().union(*(
                set(t.name for t in gt.get_terminals()) for gt in gene_trees
            ))
        )
        all_taxa_fs = frozenset(all_taxa)
        parent_map = self._build_parent_map(species_tree)

        # Extract bipartitions + branch lengths from all gene trees
        gene_tree_bp_lengths = []
        for gt in gene_trees:
            # Prune to shared taxa if needed
            gt_taxa = set(t.name for t in gt.get_terminals())
            if gt_taxa != set(all_taxa):
                taxa_to_remove = gt_taxa - set(all_taxa)
                for taxon in taxa_to_remove:
                    gt.prune(taxon)
            gene_tree_bp_lengths.append(
                self._extract_bipartitions_with_lengths(gt, all_taxa_fs)
            )

        result = {}
        for clade in species_tree.find_clades(order="preorder"):
            if clade.is_terminal():
                continue

            groups = self._get_four_groups(
                species_tree, clade, parent_map, all_taxa_fs
            )
            if groups is None:
                continue

            C1, C2, S, D = groups

            concordant_bp = self._canonical_split(C1 | C2, all_taxa_fs)
            nni_alt1_bp = self._canonical_split(S | C2, all_taxa_fs)
            nni_alt2_bp = self._canonical_split(C1 | S, all_taxa_fs)

            concordant_lengths = []
            discordant_lengths = []

            for bp_lengths in gene_tree_bp_lengths:
                if concordant_bp in bp_lengths:
                    concordant_lengths.append(bp_lengths[concordant_bp])
                elif nni_alt1_bp in bp_lengths:
                    discordant_lengths.append(bp_lengths[nni_alt1_bp])
                elif nni_alt2_bp in bp_lengths:
                    discordant_lengths.append(bp_lengths[nni_alt2_bp])

            # Label by the smaller side of the species tree split
            node_tips = frozenset(
                t.name for t in clade.get_terminals()
            )
            split_label = sorted(node_tips) if len(node_tips) <= len(all_taxa_fs) - len(node_tips) else sorted(all_taxa_fs - node_tips)

            branch_key = ",".join(split_label)
            result[branch_key] = dict(
                split=split_label,
                n_concordant=len(concordant_lengths),
                n_discordant=len(discordant_lengths),
                concordant_lengths=concordant_lengths,
                discordant_lengths=discordant_lengths,
            )

        return result
```

**Step 4: Run tests to verify they pass**

Run: `python -m pytest tests/unit/services/tree/test_evo_tempo_map.py::TestClassifyGeneTrees -v`
Expected: PASS

**Step 5: Commit**

```bash
git add phykit/services/tree/evo_tempo_map.py tests/unit/services/tree/test_evo_tempo_map.py
git commit -m "feat(evo_tempo_map): implement bipartition classification of gene trees"
```

---

### Task 4: Implement statistical tests (Mann-Whitney U + permutation)

**Files:**
- Modify: `phykit/services/tree/evo_tempo_map.py`
- Modify: `tests/unit/services/tree/test_evo_tempo_map.py`

**Step 1: Write the failing tests**

Add to the test file:

```python
import numpy as np


class TestStatisticalTests:
    @pytest.fixture
    def svc(self):
        from phykit.services.tree.evo_tempo_map import EvoTempoMap
        args = Namespace(
            tree=TREE_SIMPLE, gene_trees=GENE_TREES,
            verbose=False, json=False, plot_output=None,
        )
        return EvoTempoMap(args)

    def test_mann_whitney_returns_U_and_pvalue(self, svc):
        conc = [10.0, 12.0, 11.0, 13.0, 10.5]
        disc = [5.0, 6.0, 4.5, 5.5]
        result = svc._test_branch(conc, disc)
        assert "mann_whitney_U" in result
        assert "mann_whitney_p" in result
        assert "permutation_p" in result
        # Concordant should be significantly longer
        assert result["mann_whitney_p"] < 0.05

    def test_permutation_p_between_0_and_1(self, svc):
        conc = [10.0, 12.0, 11.0, 13.0, 10.5]
        disc = [5.0, 6.0, 4.5, 5.5]
        result = svc._test_branch(conc, disc)
        assert 0.0 <= result["permutation_p"] <= 1.0

    def test_no_difference_gives_high_p(self, svc):
        # Same distribution — should NOT be significant
        rng = np.random.default_rng(42)
        vals = rng.normal(10, 1, 20).tolist()
        conc = vals[:10]
        disc = vals[10:]
        result = svc._test_branch(conc, disc)
        assert result["mann_whitney_p"] > 0.05

    def test_insufficient_data_returns_na(self, svc):
        # Only 1 discordant — too few to test
        conc = [10.0, 12.0, 11.0]
        disc = [5.0]
        result = svc._test_branch(conc, disc)
        assert result["mann_whitney_p"] is None
        assert result["permutation_p"] is None

    def test_empty_discordant_returns_na(self, svc):
        conc = [10.0, 12.0, 11.0]
        disc = []
        result = svc._test_branch(conc, disc)
        assert result["mann_whitney_p"] is None
        assert result["permutation_p"] is None


class TestFDRCorrection:
    def test_fdr_correction(self):
        from phykit.services.tree.evo_tempo_map import EvoTempoMap
        p_values = [0.01, 0.04, 0.03, 0.20]
        corrected = EvoTempoMap._fdr(p_values)
        assert len(corrected) == 4
        # Corrected values should be >= original
        for orig, corr in zip(p_values, corrected):
            assert corr >= orig
        # All should be <= 1.0
        for c in corrected:
            assert c <= 1.0

    def test_fdr_empty_list(self):
        from phykit.services.tree.evo_tempo_map import EvoTempoMap
        assert EvoTempoMap._fdr([]) == []

    def test_fdr_single_value(self):
        from phykit.services.tree.evo_tempo_map import EvoTempoMap
        assert EvoTempoMap._fdr([0.05]) == [0.05]
```

**Step 2: Run tests to verify they fail**

Run: `python -m pytest tests/unit/services/tree/test_evo_tempo_map.py::TestStatisticalTests -v`
Expected: FAIL

**Step 3: Implement the statistical tests**

Add to `phykit/services/tree/evo_tempo_map.py`:

```python
    def _test_branch(
        self,
        concordant_lengths: List[float],
        discordant_lengths: List[float],
        n_permutations: int = 1000,
    ) -> Dict:
        """Run Mann-Whitney U + permutation test comparing concordant vs
        discordant branch lengths.

        Returns dict with mann_whitney_U, mann_whitney_p, permutation_p,
        and summary stats. Returns None for p-values if either group
        has fewer than 2 observations.
        """
        import numpy as np

        n_conc = len(concordant_lengths)
        n_disc = len(discordant_lengths)

        conc_arr = np.array(concordant_lengths) if n_conc > 0 else np.array([])
        disc_arr = np.array(discordant_lengths) if n_disc > 0 else np.array([])

        stats = dict(
            concordant_mean=float(np.mean(conc_arr)) if n_conc > 0 else None,
            concordant_median=float(np.median(conc_arr)) if n_conc > 0 else None,
            concordant_std=float(np.std(conc_arr, ddof=1)) if n_conc > 1 else None,
            discordant_mean=float(np.mean(disc_arr)) if n_disc > 0 else None,
            discordant_median=float(np.median(disc_arr)) if n_disc > 0 else None,
            discordant_std=float(np.std(disc_arr, ddof=1)) if n_disc > 1 else None,
        )

        # Need at least 2 in each group for a meaningful test
        if n_conc < 2 or n_disc < 2:
            stats["mann_whitney_U"] = None
            stats["mann_whitney_p"] = None
            stats["permutation_p"] = None
            return stats

        # Mann-Whitney U test
        from scipy.stats import mannwhitneyu

        U, mw_p = mannwhitneyu(
            conc_arr, disc_arr, alternative="two-sided"
        )
        stats["mann_whitney_U"] = float(U)
        stats["mann_whitney_p"] = float(mw_p)

        # Permutation test
        observed_diff = abs(float(np.median(conc_arr)) - float(np.median(disc_arr)))
        combined = np.concatenate([conc_arr, disc_arr])
        rng = np.random.default_rng(42)
        count_ge = 0
        for _ in range(n_permutations):
            perm = rng.permutation(combined)
            perm_conc = perm[:n_conc]
            perm_disc = perm[n_conc:]
            perm_diff = abs(float(np.median(perm_conc)) - float(np.median(perm_disc)))
            if perm_diff >= observed_diff:
                count_ge += 1
        stats["permutation_p"] = count_ge / n_permutations

        return stats

    @staticmethod
    def _fdr(p_values: List[float]) -> List[float]:
        """Benjamini-Hochberg FDR correction.

        Same algorithm as relative_rate_test._fdr.
        """
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
```

**Step 4: Run tests to verify they pass**

Run: `python -m pytest tests/unit/services/tree/test_evo_tempo_map.py::TestStatisticalTests tests/unit/services/tree/test_evo_tempo_map.py::TestFDRCorrection -v`
Expected: PASS

**Step 5: Commit**

```bash
git add phykit/services/tree/evo_tempo_map.py tests/unit/services/tree/test_evo_tempo_map.py
git commit -m "feat(evo_tempo_map): implement Mann-Whitney U, permutation test, and FDR correction"
```

---

### Task 5: Implement global treeness comparison

Cross-check: the treeness calculation should match `phykit treeness` for individual trees. Treeness = sum(internal branch lengths) / sum(all branch lengths). This is `Tree.calculate_treeness()` in `phykit/services/tree/base.py:137-157`.

**Files:**
- Modify: `phykit/services/tree/evo_tempo_map.py`
- Modify: `tests/unit/services/tree/test_evo_tempo_map.py`

**Step 1: Write the failing tests**

```python
class TestGlobalTreeness:
    @pytest.fixture
    def svc(self):
        from phykit.services.tree.evo_tempo_map import EvoTempoMap
        args = Namespace(
            tree=TREE_SIMPLE, gene_trees=GENE_TREES,
            verbose=False, json=False, plot_output=None,
        )
        return EvoTempoMap(args)

    def test_treeness_matches_base_class(self, svc):
        """Cross-check: treeness computed here should match
        Tree.calculate_treeness() for each gene tree."""
        gene_trees = svc._parse_gene_trees(GENE_TREES)
        for gt in gene_trees:
            expected = svc.calculate_treeness(gt)
            actual = svc._compute_treeness(gt)
            assert abs(actual - expected) < 1e-10, (
                f"Treeness mismatch: expected {expected}, got {actual}"
            )

    def test_global_treeness_returns_two_groups(self, svc):
        species_tree = svc.read_tree_file()
        gene_trees = svc._parse_gene_trees(GENE_TREES)
        classified = svc._classify_gene_trees(species_tree, gene_trees)
        result = svc._compute_global_treeness(gene_trees, classified)
        assert "treeness_concordant" in result
        assert "treeness_discordant" in result

    def test_all_concordant_treeness(self, svc):
        """When all gene trees are concordant, discordant treeness
        should have no values."""
        species_tree = svc.read_tree_file()
        identical_trees = [svc.read_tree_file() for _ in range(3)]
        classified = svc._classify_gene_trees(species_tree, identical_trees)
        result = svc._compute_global_treeness(identical_trees, classified)
        assert result["treeness_concordant"]["n"] == 3
        assert result["treeness_discordant"]["n"] == 0
```

**Step 2: Run tests to verify they fail**

Run: `python -m pytest tests/unit/services/tree/test_evo_tempo_map.py::TestGlobalTreeness -v`
Expected: FAIL

**Step 3: Implement global treeness**

Add to `evo_tempo_map.py`:

```python
    @staticmethod
    def _compute_treeness(tree) -> float:
        """Compute treeness = internal branch length / total branch length.

        Same algorithm as Tree.calculate_treeness() in base.py.
        """
        inter_len = 0.0
        for internal in tree.get_nonterminals():
            if internal.branch_length is not None:
                inter_len += internal.branch_length
        total_len = tree.total_branch_length()
        if total_len == 0:
            return 0.0
        return inter_len / total_len

    def _compute_global_treeness(
        self, gene_trees, classified
    ) -> Dict:
        """Compute treeness for concordant vs discordant gene trees.

        A gene tree is considered "concordant" if it is concordant at
        ALL species tree branches (i.e., topology matches species tree).
        Otherwise it is "discordant".
        """
        import numpy as np

        # Determine which gene trees are fully concordant vs have any discordance
        n_gt = len(gene_trees)
        # A gene tree is fully concordant if it contributed a concordant
        # branch length at every tested branch
        # Simpler approach: count how many branches each gene tree is concordant at
        # vs discordant. If concordant at all branches -> fully concordant.

        # Re-derive per-gene-tree classification
        # For each gene tree, check if it contributed to concordant at every branch
        concordant_treeness = []
        discordant_treeness = []

        # Build a set of all tested branches
        all_branches = list(classified.keys())
        if not all_branches:
            return dict(
                treeness_concordant=dict(mean=None, median=None, n=0),
                treeness_discordant=dict(mean=None, median=None, n=0),
                treeness_U_p=None,
            )

        # For simplicity and because gene trees may be uninformative at
        # some branches, classify a gene tree as "fully concordant" if
        # it has no discordant branches among the branches where it is
        # informative.
        # We need to track per-gene-tree: was it ever discordant?
        # Re-derive from the classification data structure.
        # The classified dict has concordant_lengths and discordant_lengths
        # lists, but they're aggregated. We need per-gene-tree info.
        # Better: re-extract from the species tree + gene trees directly.
        species_tree_tips = set()
        for branch_key, data in classified.items():
            species_tree_tips.update(data["split"])

        all_taxa_fs = frozenset(species_tree_tips)

        # Build bipartition sets per gene tree (same as in _classify)
        gt_bp_sets = []
        for gt in gene_trees:
            splits = set()
            for clade in gt.get_nonterminals():
                tips = frozenset(t.name for t in clade.get_terminals())
                if len(tips) <= 1 or tips == all_taxa_fs:
                    continue
                splits.add(self._canonical_split(tips, all_taxa_fs))
            gt_bp_sets.append(splits)

        # For each gene tree, check if it's concordant at all species tree branches
        # We need the concordant bipartitions for each branch
        # Reconstruct from the split data
        # Actually, for global treeness we just need: is the gene tree's
        # topology identical to the species tree's? Use RF-distance == 0
        # as proxy, or: is it concordant at every branch?
        # Simpler: a gene tree is "discordant" if any of its bipartitions
        # differ from the species tree at any branch.
        # Use the classified data to see which gene trees contributed to
        # discordant at any branch.

        # Track per-gene-tree index: ever_discordant?
        ever_discordant = [False] * n_gt

        # We need the per-gene-tree classification. Let's rebuild it.
        # This is somewhat redundant with _classify_gene_trees, but
        # the global treeness is a simple summary so the overhead is small.
        species_tree_obj = None  # We don't have it here
        # Instead, just check: for each branch, the concordant_bp should
        # be in the gene tree's split set. If not, it's discordant.

        # We don't have the concordant bipartitions stored in classified.
        # Re-derive: for branch with split=[A,B], the concordant bp is
        # _canonical_split(frozenset(split), all_taxa_fs).
        # But all_taxa_fs may not be complete from just the split data.
        # Let's take a different approach: pass gene_tree indices through
        # _classify_gene_trees instead.

        # REVISED APPROACH: extend _classify_gene_trees to also return
        # per-gene-tree indices. But that changes the API.
        # For now, just classify based on "is this gene tree concordant
        # at the majority of tested branches?"
        # Actually the simplest correct approach: a gene tree is
        # "concordant" overall if RF distance to species tree is 0.

        # Simplest: for each gene tree, compute treeness and bucket it
        # based on whether ANY of its bipartitions conflict with the
        # species tree. We already have gt_bp_sets.

        # To do this, we need the species tree bipartitions.
        # Extract them. But we don't have the species tree object here.
        # Let's just call calculate_treeness for all gene trees and
        # classify them based on the classified dict structure.

        # PRAGMATIC: Change signature to accept species_tree too.
        # For now, compute treeness for all gene trees and return both
        # lists, leaving the concordant/discordant split to be determined
        # by the caller in run().

        # Actually, let me just accept species_tree as a parameter.
        # [Will refactor in implementation]

        # For the plan, assume _compute_global_treeness takes
        # (species_tree, gene_trees) and does the full computation.
        pass  # See actual implementation below
```

**Actually, let me simplify the plan for this step.** The implementation should:

1. Accept `(species_tree, gene_trees)` as parameters
2. Extract species tree bipartitions
3. For each gene tree: check if RF distance is 0 (all bipartitions match) → concordant, else discordant
4. Compute treeness for each group
5. Run Mann-Whitney U on the two groups

The actual implementation in Step 3:

```python
    def _compute_global_treeness(
        self, species_tree, gene_trees
    ) -> Dict:
        """Classify each gene tree as globally concordant or discordant
        (based on whether ALL bipartitions match the species tree),
        compute treeness for each group, and test for a difference.
        """
        import numpy as np

        all_taxa = sorted(t.name for t in species_tree.get_terminals())
        all_taxa_fs = frozenset(all_taxa)

        # Species tree bipartitions
        sp_splits = set()
        for clade in species_tree.get_nonterminals():
            tips = frozenset(t.name for t in clade.get_terminals())
            if len(tips) <= 1 or tips == all_taxa_fs:
                continue
            sp_splits.add(self._canonical_split(tips, all_taxa_fs))

        concordant_treeness = []
        discordant_treeness = []

        for gt in gene_trees:
            gt_splits = set()
            for clade in gt.get_nonterminals():
                tips = frozenset(t.name for t in clade.get_terminals())
                if len(tips) <= 1 or tips == all_taxa_fs:
                    continue
                gt_splits.add(self._canonical_split(tips, all_taxa_fs))

            treeness = self._compute_treeness(gt)

            if gt_splits == sp_splits:
                concordant_treeness.append(treeness)
            else:
                discordant_treeness.append(treeness)

        result = dict(
            treeness_concordant=dict(
                mean=float(np.mean(concordant_treeness)) if concordant_treeness else None,
                median=float(np.median(concordant_treeness)) if concordant_treeness else None,
                n=len(concordant_treeness),
            ),
            treeness_discordant=dict(
                mean=float(np.mean(discordant_treeness)) if discordant_treeness else None,
                median=float(np.median(discordant_treeness)) if discordant_treeness else None,
                n=len(discordant_treeness),
            ),
        )

        if len(concordant_treeness) >= 2 and len(discordant_treeness) >= 2:
            from scipy.stats import mannwhitneyu
            _, p = mannwhitneyu(
                concordant_treeness, discordant_treeness,
                alternative="two-sided",
            )
            result["treeness_U_p"] = float(p)
        else:
            result["treeness_U_p"] = None

        return result
```

Update the test fixture to pass `species_tree` to `_compute_global_treeness`:

```python
    def test_global_treeness_returns_two_groups(self, svc):
        species_tree = svc.read_tree_file()
        gene_trees = svc._parse_gene_trees(GENE_TREES)
        result = svc._compute_global_treeness(species_tree, gene_trees)
        assert "treeness_concordant" in result
        assert "treeness_discordant" in result

    def test_all_concordant_treeness(self, svc):
        species_tree = svc.read_tree_file()
        identical_trees = [svc.read_tree_file() for _ in range(3)]
        result = svc._compute_global_treeness(species_tree, identical_trees)
        assert result["treeness_concordant"]["n"] == 3
        assert result["treeness_discordant"]["n"] == 0
```

**Step 4: Run tests to verify they pass**

Run: `python -m pytest tests/unit/services/tree/test_evo_tempo_map.py::TestGlobalTreeness -v`
Expected: PASS

**Step 5: Commit**

```bash
git add phykit/services/tree/evo_tempo_map.py tests/unit/services/tree/test_evo_tempo_map.py
git commit -m "feat(evo_tempo_map): implement global treeness comparison"
```

---

### Task 6: Wire up the `run()` method with text and JSON output

**Files:**
- Modify: `phykit/services/tree/evo_tempo_map.py`
- Modify: `tests/unit/services/tree/test_evo_tempo_map.py`

**Step 1: Write the failing tests**

```python
class TestRun:
    def test_run_text_output(self, capsys):
        from phykit.services.tree.evo_tempo_map import EvoTempoMap
        args = Namespace(
            tree=TREE_SIMPLE, gene_trees=GENE_TREES,
            verbose=False, json=False, plot_output=None,
        )
        svc = EvoTempoMap(args)
        svc.run()
        captured = capsys.readouterr()
        assert "n_conc" in captured.out
        assert "n_disc" in captured.out
        # Should have branch rows
        lines = [l for l in captured.out.strip().split("\n") if l.strip() and not l.startswith("---") and not l.startswith("branch") and not l.startswith("Global") and not l.startswith("Branches")]
        assert len(lines) > 0

    def test_run_json_output(self, capsys):
        import json
        from phykit.services.tree.evo_tempo_map import EvoTempoMap
        args = Namespace(
            tree=TREE_SIMPLE, gene_trees=GENE_TREES,
            verbose=False, json=True, plot_output=None,
        )
        svc = EvoTempoMap(args)
        svc.run()
        captured = capsys.readouterr()
        data = json.loads(captured.out)
        assert "branches" in data
        assert "global" in data
        assert "n_gene_trees" in data["global"]
        for branch in data["branches"]:
            assert "split" in branch
            assert "n_concordant" in branch
            assert "mann_whitney_p" in branch

    def test_run_verbose_output(self, capsys):
        from phykit.services.tree.evo_tempo_map import EvoTempoMap
        args = Namespace(
            tree=TREE_SIMPLE, gene_trees=GENE_TREES,
            verbose=True, json=False, plot_output=None,
        )
        svc = EvoTempoMap(args)
        svc.run()
        captured = capsys.readouterr()
        # Verbose should have per-gene-tree details
        assert "concordant" in captured.out.lower() or "gene_tree" in captured.out.lower()
```

**Step 2: Run tests to verify they fail**

Run: `python -m pytest tests/unit/services/tree/test_evo_tempo_map.py::TestRun -v`
Expected: FAIL (run() is still a placeholder)

**Step 3: Implement the run() method**

Replace the `run()` method in `evo_tempo_map.py`:

```python
    def run(self) -> None:
        species_tree = self.read_tree_file()
        gene_trees = self._parse_and_validate_gene_trees()

        # Per-branch classification and testing
        classified = self._classify_gene_trees(species_tree, gene_trees)

        # Run statistical tests per branch
        branch_results = []
        testable_p_values = []
        for branch_key in sorted(classified.keys()):
            data = classified[branch_key]
            test_stats = self._test_branch(
                data["concordant_lengths"], data["discordant_lengths"]
            )
            entry = dict(
                split=data["split"],
                n_concordant=data["n_concordant"],
                n_discordant=data["n_discordant"],
                **test_stats,
            )
            branch_results.append(entry)
            if test_stats["mann_whitney_p"] is not None:
                testable_p_values.append(test_stats["mann_whitney_p"])

        # FDR correction
        fdr_corrected = self._fdr(testable_p_values)
        fdr_idx = 0
        for entry in branch_results:
            if entry["mann_whitney_p"] is not None:
                entry["fdr_p"] = fdr_corrected[fdr_idx]
                fdr_idx += 1
            else:
                entry["fdr_p"] = None

        # Global treeness
        global_stats = self._compute_global_treeness(species_tree, gene_trees)
        global_stats["n_gene_trees"] = len(gene_trees)
        global_stats["n_branches_tested"] = len(testable_p_values)
        global_stats["n_significant_fdr05"] = sum(
            1 for p in fdr_corrected if p < 0.05
        )

        if self.json_output:
            self._output_json(branch_results, global_stats)
            return

        self._output_text(branch_results, global_stats)

        if self.plot_output:
            self._plot(branch_results, self.plot_output)

    def _output_text(self, branch_results, global_stats):
        try:
            # Header
            print(
                f"{'branch':<30}"
                f"{'n_conc':>8}"
                f"{'n_disc':>8}"
                f"{'med_conc':>10}"
                f"{'med_disc':>10}"
                f"{'U_pval':>10}"
                f"{'perm_pval':>10}"
                f"{'fdr_p':>10}"
            )
            for entry in branch_results:
                split_str = "{" + ",".join(entry["split"]) + "}"
                med_c = f"{entry['concordant_median']:.4f}" if entry["concordant_median"] is not None else "NA"
                med_d = f"{entry['discordant_median']:.4f}" if entry["discordant_median"] is not None else "NA"
                up = f"{entry['mann_whitney_p']:.4f}" if entry["mann_whitney_p"] is not None else "NA"
                pp = f"{entry['permutation_p']:.4f}" if entry["permutation_p"] is not None else "NA"
                fp = f"{entry['fdr_p']:.4f}" if entry["fdr_p"] is not None else "NA"
                print(
                    f"{split_str:<30}"
                    f"{entry['n_concordant']:>8}"
                    f"{entry['n_discordant']:>8}"
                    f"{med_c:>10}"
                    f"{med_d:>10}"
                    f"{up:>10}"
                    f"{pp:>10}"
                    f"{fp:>10}"
                )
            print("---")
            tc = global_stats["treeness_concordant"]
            td = global_stats["treeness_discordant"]
            tc_str = f"{tc['mean']:.4f}" if tc["mean"] is not None else "NA"
            td_str = f"{td['mean']:.4f}" if td["mean"] is not None else "NA"
            tu_str = f"{global_stats['treeness_U_p']:.4f}" if global_stats["treeness_U_p"] is not None else "NA"
            print(
                f"Global treeness: concordant={tc_str}, "
                f"discordant={td_str} (U p={tu_str})"
            )
            print(
                f"Branches tested: {global_stats['n_branches_tested']}, "
                f"significant (FDR<0.05): {global_stats['n_significant_fdr05']}"
            )
            if self.verbose:
                print("\nPer-branch concordant lengths:")
                for entry in branch_results:
                    split_str = "{" + ",".join(entry["split"]) + "}"
                    conc = entry.get("_concordant_lengths", [])
                    disc = entry.get("_discordant_lengths", [])
                    # These would need to be stored; see note below
        except BrokenPipeError:
            pass

    def _output_json(self, branch_results, global_stats):
        # Clean up internal keys before JSON output
        clean_branches = []
        for entry in branch_results:
            clean = {
                k: v for k, v in entry.items()
                if not k.startswith("_")
            }
            clean_branches.append(clean)
        result = dict(
            branches=clean_branches,
            global_=global_stats,
        )
        print_json(result)
```

Note: For verbose mode, store the raw lengths in the branch results during `run()`. Add `_concordant_lengths` and `_discordant_lengths` keys (prefixed with `_` so they're excluded from JSON output). Update `_output_text` to print them when verbose is true.

**Step 4: Run tests to verify they pass**

Run: `python -m pytest tests/unit/services/tree/test_evo_tempo_map.py::TestRun -v`
Expected: PASS

**Step 5: Commit**

```bash
git add phykit/services/tree/evo_tempo_map.py tests/unit/services/tree/test_evo_tempo_map.py
git commit -m "feat(evo_tempo_map): implement run() with text and JSON output"
```

---

### Task 7: Implement the plot

**Files:**
- Modify: `phykit/services/tree/evo_tempo_map.py`
- Modify: `tests/unit/services/tree/test_evo_tempo_map.py`

**Step 1: Write the failing test**

```python
import os
import tempfile


class TestPlot:
    def test_plot_creates_file(self):
        from phykit.services.tree.evo_tempo_map import EvoTempoMap

        with tempfile.NamedTemporaryFile(
            suffix=".png", delete=False
        ) as f:
            tmp_path = f.name
        try:
            args = Namespace(
                tree=TREE_SIMPLE, gene_trees=GENE_TREES,
                verbose=False, json=False, plot_output=tmp_path,
            )
            svc = EvoTempoMap(args)
            svc.run()
            assert os.path.exists(tmp_path)
            assert os.path.getsize(tmp_path) > 0
        finally:
            if os.path.exists(tmp_path):
                os.unlink(tmp_path)
```

**Step 2: Run test to verify it fails**

Run: `python -m pytest tests/unit/services/tree/test_evo_tempo_map.py::TestPlot -v`
Expected: FAIL (`_plot` not implemented or `AttributeError`)

**Step 3: Implement the plot**

Add to `evo_tempo_map.py`:

```python
    @staticmethod
    def _plot(branch_results, output_path):
        """Grouped box/strip plot of concordant vs discordant branch lengths."""
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np

        # Collect data for branches that have both groups
        labels = []
        conc_data = []
        disc_data = []
        sig_flags = []

        for entry in branch_results:
            conc = entry.get("_concordant_lengths", [])
            disc = entry.get("_discordant_lengths", [])
            if not conc and not disc:
                continue
            split_str = "{" + ",".join(entry["split"]) + "}"
            labels.append(split_str)
            conc_data.append(conc)
            disc_data.append(disc)
            sig_flags.append(
                entry.get("fdr_p") is not None and entry["fdr_p"] < 0.05
            )

        if not labels:
            return

        n = len(labels)
        fig, ax = plt.subplots(figsize=(max(8, n * 1.5), 5))

        positions = np.arange(n)
        width = 0.35

        # Box plots
        bp_conc = ax.boxplot(
            conc_data, positions=positions - width / 2, widths=width * 0.8,
            patch_artist=True, showfliers=False,
            boxprops=dict(facecolor="#4C72B0", alpha=0.7),
            medianprops=dict(color="black"),
        )
        bp_disc = ax.boxplot(
            disc_data, positions=positions + width / 2, widths=width * 0.8,
            patch_artist=True, showfliers=False,
            boxprops=dict(facecolor="#DD8452", alpha=0.7),
            medianprops=dict(color="black"),
        )

        # Strip (jitter) points
        rng = np.random.default_rng(42)
        for i in range(n):
            if conc_data[i]:
                jitter = rng.uniform(-0.05, 0.05, len(conc_data[i]))
                ax.scatter(
                    [i - width / 2] * len(conc_data[i]) + jitter,
                    conc_data[i], color="#4C72B0", alpha=0.6, s=20, zorder=3,
                )
            if disc_data[i]:
                jitter = rng.uniform(-0.05, 0.05, len(disc_data[i]))
                ax.scatter(
                    [i + width / 2] * len(disc_data[i]) + jitter,
                    disc_data[i], color="#DD8452", alpha=0.6, s=20, zorder=3,
                )
            if sig_flags[i]:
                max_val = max(
                    max(conc_data[i]) if conc_data[i] else 0,
                    max(disc_data[i]) if disc_data[i] else 0,
                )
                ax.text(i, max_val * 1.05, "*", ha="center", fontsize=14, fontweight="bold")

        ax.set_xticks(positions)
        ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=9)
        ax.set_ylabel("Branch length (subs/site)", fontsize=12)
        ax.set_xlabel("Species tree branch", fontsize=12)
        ax.legend(
            [bp_conc["boxes"][0], bp_disc["boxes"][0]],
            ["Concordant", "Discordant"],
            loc="upper right",
        )

        fig.tight_layout()
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close(fig)
```

Important: In `run()`, store the raw lengths in the branch results so the plot can access them:

```python
        # In run(), after computing test_stats:
        entry["_concordant_lengths"] = data["concordant_lengths"]
        entry["_discordant_lengths"] = data["discordant_lengths"]
```

**Step 4: Run test to verify it passes**

Run: `python -m pytest tests/unit/services/tree/test_evo_tempo_map.py::TestPlot -v`
Expected: PASS

**Step 5: Commit**

```bash
git add phykit/services/tree/evo_tempo_map.py tests/unit/services/tree/test_evo_tempo_map.py
git commit -m "feat(evo_tempo_map): implement box/strip plot output"
```

---

### Task 8: Register the command in CLI

**Files:**
- Modify: `phykit/services/tree/__init__.py` (add to `_EXPORTS` dict)
- Modify: `phykit/service_factories.py` (add lazy factory, ~line 98)
- Modify: `phykit/cli_registry.py` (add aliases, ~line 156)
- Modify: `phykit/phykit.py` (add CLI parser, ~line 4852)
- Modify: `tests/unit/services/tree/test_evo_tempo_map.py`

**Step 1: Write the failing test**

Add to the test file:

```python
class TestCLIIntegration:
    def test_cli_runs_without_error(self, capsys):
        """Test that the command runs through the full CLI path."""
        from phykit.phykit import Phykit
        pk = Phykit()
        pk(["evo_tempo_map", "-t", TREE_SIMPLE, "-g", GENE_TREES])
        captured = capsys.readouterr()
        assert "n_conc" in captured.out

    def test_cli_alias_etm(self, capsys):
        from phykit.phykit import Phykit
        pk = Phykit()
        pk(["etm", "-t", TREE_SIMPLE, "-g", GENE_TREES])
        captured = capsys.readouterr()
        assert "n_conc" in captured.out

    def test_cli_json_flag(self, capsys):
        import json
        from phykit.phykit import Phykit
        pk = Phykit()
        pk(["evo_tempo_map", "-t", TREE_SIMPLE, "-g", GENE_TREES, "--json"])
        captured = capsys.readouterr()
        data = json.loads(captured.out)
        assert "branches" in data
```

**Step 2: Run tests to verify they fail**

Run: `python -m pytest tests/unit/services/tree/test_evo_tempo_map.py::TestCLIIntegration -v`
Expected: FAIL (command not registered)

**Step 3: Register the command**

**3a.** Add to `phykit/services/tree/__init__.py` `_EXPORTS` dict (before the `ConcordanceAsr` line):

```python
    "EvoTempoMap": "evo_tempo_map",
```

**3b.** Add to `phykit/service_factories.py` (after the `TreenessOverRCV` line, ~line 98):

```python
EvoTempoMap = _LazyServiceFactory("phykit.services.tree.evo_tempo_map", "EvoTempoMap")
```

**3c.** Add to `phykit/cli_registry.py` (after the `treeness_over_rcv` aliases, ~line 156):

```python
    "etm": "evo_tempo_map",
```

**3d.** Add to `phykit/phykit.py` (after the `treeness_over_rcv` method, ~line 4852). Find the import of `TreenessOverRCV` from `service_factories` and add `EvoTempoMap` to the same import block. Then add the static method:

```python
    @staticmethod
    def evo_tempo_map(argv):
        parser = _new_parser(
            description=textwrap.dedent(
                f"""\
                {help_header}

                Detect rate-topology associations by comparing branch length
                distributions between concordant and discordant gene trees at
                each species tree branch.

                Under the multispecies coalescent, discordant gene trees should
                have shorter internal branches near the discordant node. Deviations
                suggest substitution rate heterogeneity correlated with topology
                (adaptive evolution, different selective pressures, or model
                misspecification).

                For each internal branch of the species tree, gene trees are
                classified as concordant or discordant via bipartition matching.
                The homologous branch length is extracted from each gene tree
                and the two groups are compared using Mann-Whitney U and
                permutation tests. P-values are corrected for multiple testing
                using Benjamini-Hochberg FDR.

                A global treeness (internal/total branch length ratio) comparison
                is also reported.

                Aliases:
                  evo_tempo_map, etm
                Command line interfaces:
                  pk_evo_tempo_map, pk_etm

                Usage:
                phykit evo_tempo_map -t/--tree <tree> -g/--gene-trees <gene_trees>
                    [--plot <output>] [-v/--verbose] [--json]

                Options
                =====================================================
                -t/--tree                   a species tree file

                -g/--gene-trees             multi-Newick file of gene trees
                                            with branch lengths

                --plot                      optional output path for
                                            box/strip plot (PNG)

                -v/--verbose                print per-gene-tree details

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument(
            "-t", "--tree", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "-g", "--gene-trees", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "--plot", dest="plot_output", type=str, required=False,
            default=None, help=SUPPRESS, metavar=""
        )
        _add_verbose_argument(parser)
        _add_json_argument(parser)
        _run_service(parser, argv, EvoTempoMap)
```

Also add `EvoTempoMap` to the imports from `service_factories` at the top of the file. Find the line that imports tree service factories (look for `from .service_factories import`) and add `EvoTempoMap`.

**Step 4: Run tests to verify they pass**

Run: `python -m pytest tests/unit/services/tree/test_evo_tempo_map.py::TestCLIIntegration -v`
Expected: PASS

**Step 5: Commit**

```bash
git add phykit/services/tree/__init__.py phykit/service_factories.py phykit/cli_registry.py phykit/phykit.py tests/unit/services/tree/test_evo_tempo_map.py
git commit -m "feat(evo_tempo_map): register CLI command with aliases"
```

---

### Task 9: Add documentation

**Files:**
- Modify: `docs/usage/index.rst`

**Step 1: Add the command reference**

Add a new section in `docs/usage/index.rst`. Find the concordance_asr section (at `.. _cmd-concordance_asr:`, ~line 1110) and add the new command documentation before it (so commands remain alphabetical under their category) or after it.

Add this block:

```rst
.. _cmd-evo_tempo_map:

Evolutionary tempo mapping
##########################
Function names: evo_tempo_map; etm |br|
Command line interface: pk_evo_tempo_map; pk_etm

Detect rate-topology associations by comparing branch length distributions
between concordant and discordant gene trees at each species tree branch.

Under the multispecies coalescent, discordant gene trees should have shorter
internal branches near the discordant node (because the coalescence happened
deeper, in the ancestral population). Deviations from this expectation suggest
substitution rate heterogeneity correlated with topology, which could indicate
adaptive evolution, different selective pressures in hybridizing lineages, or
systematic error from model misspecification.

For each internal branch of the species tree, gene trees are classified as
concordant or discordant via bipartition matching (same as gCF). The homologous
branch length is extracted from each gene tree and the two groups are compared
using a Mann-Whitney U test and a permutation test (1000 permutations). P-values
are corrected for multiple testing using Benjamini-Hochberg FDR.

A global treeness (internal/total branch length ratio) comparison between
concordant and discordant gene trees is also reported.

.. code-block:: shell

   phykit evo_tempo_map -t <species_tree> -g <gene_trees> [--plot <output>] [-v] [--json]

Options: |br|
*-t/\\-\\-tree*: a species tree file |br|
*-g/\\-\\-gene-trees*: multi-Newick file of gene trees with branch lengths |br|
*--plot*: optional output path for box/strip plot (PNG) |br|
*-v/\\-\\-verbose*: print per-gene-tree classification details |br|
*--json*: optional argument to print results as JSON
```

**Step 2: Add to the functional category section**

Find the "Tree comparison & consensus" category (~line 141) and add:

```rst
- :ref:`Evolutionary tempo mapping <cmd-evo_tempo_map>`: Detect rate-topology associations in gene trees
```

Or add it to a more appropriate category like "Phylogenetic signal" or create a subsection. The best fit is "Tree comparison & consensus" since it compares gene trees to the species tree.

**Step 3: Run the doc coverage test**

Run: `python -m pytest tests/unit/test_docs_category_coverage.py -v`
Expected: PASS (the test checks that every function in `cli_registry.py` appears in `docs/usage/index.rst`)

**Step 4: Commit**

```bash
git add docs/usage/index.rst
git commit -m "docs(evo_tempo_map): add usage documentation and category entry"
```

---

### Task 10: Run full test suite and verify

**Step 1: Run all evo_tempo_map tests**

```bash
python -m pytest tests/unit/services/tree/test_evo_tempo_map.py -v
```

Expected: All tests PASS

**Step 2: Run the full test suite**

```bash
python -m pytest tests/ -v --timeout=120
```

Expected: No new failures (pre-existing failures unrelated to this feature are acceptable)

**Step 3: Run doc coverage test**

```bash
python -m pytest tests/unit/test_docs_category_coverage.py -v
```

Expected: PASS

**Step 4: CLI smoke tests**

```bash
# Text output
python -m phykit evo_tempo_map -t tests/sample_files/tree_simple.tre -g tests/sample_files/gene_trees_simple.nwk

# JSON output
python -m phykit etm -t tests/sample_files/tree_simple.tre -g tests/sample_files/gene_trees_simple.nwk --json

# Plot output
python -m phykit evo_tempo_map -t tests/sample_files/tree_simple.tre -g tests/sample_files/gene_trees_simple.nwk --plot /tmp/etm_test.png

# Verbose
python -m phykit evo_tempo_map -t tests/sample_files/tree_simple.tre -g tests/sample_files/gene_trees_simple.nwk -v
```

**Step 5: Final commit if any fixes needed**

```bash
git add -A
git commit -m "fix(evo_tempo_map): address issues from full test suite run"
```
