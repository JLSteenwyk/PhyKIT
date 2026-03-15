"""
Unit tests for shared quartet concordance utilities.

Tests gCF/gDF computation via bipartition matching and ASTRAL
q1/q2/q3 annotation parsing. Ground truth values computed from
tests/sample_files/tree_simple.tre + gene_trees_simple.nwk.
"""
import pytest

from Bio import Phylo

from phykit.helpers.quartet_utils import (
    canonical_split,
    compute_gcf_per_node,
    parse_astral_annotations,
    _parse_qs_from_label,
)


class TestCanonicalSplit:
    def test_basic(self):
        all_taxa = frozenset(["A", "B", "C", "D"])
        s1 = canonical_split(frozenset(["A", "B"]), all_taxa)
        s2 = canonical_split(frozenset(["C", "D"]), all_taxa)
        assert s1 == s2  # Same split from either side


class TestComputeGcfPerNode:
    def test_ground_truth_values(self):
        """Validate against manually computed ground truth."""
        species_tree = Phylo.read(
            "tests/sample_files/tree_simple.tre", "newick"
        )
        gene_trees = list(
            Phylo.parse("tests/sample_files/gene_trees_simple.nwk", "newick")
        )

        result = compute_gcf_per_node(species_tree, gene_trees)

        # Collect by tip composition for stable assertions
        by_tips = {}
        for clade in species_tree.find_clades(order="preorder"):
            if clade.is_terminal() or clade == species_tree.root:
                continue
            cid = id(clade)
            if cid in result:
                tips = frozenset(t.name for t in clade.get_terminals())
                by_tips[tips] = result[cid]

        # bear+raccoon: fully concordant
        br = frozenset(["bear", "raccoon"])
        assert by_tips[br][0] == pytest.approx(1.0)  # gCF
        assert by_tips[br][3] == 7  # concordant count

        # sea_lion+seal: fully concordant
        ss = frozenset(["sea_lion", "seal"])
        assert by_tips[ss][0] == pytest.approx(1.0)

        # monkey+cat: fully concordant
        mc = frozenset(["monkey", "cat"])
        assert by_tips[mc][0] == pytest.approx(1.0)

        # monkey+cat+weasel: mixed
        mcw = frozenset(["monkey", "cat", "weasel"])
        assert by_tips[mcw][0] == pytest.approx(9 / 19, abs=0.01)
        assert by_tips[mcw][1] == pytest.approx(10 / 19, abs=0.01)

    def test_empty_gene_trees(self):
        species_tree = Phylo.read(
            "tests/sample_files/tree_simple.tre", "newick"
        )
        result = compute_gcf_per_node(species_tree, [])
        # All nodes should default to gCF=1.0 (no conflicting data)
        for cid, vals in result.items():
            assert vals[0] == 1.0


class TestParseAstralAnnotations:
    def test_bracket_format(self):
        from io import StringIO
        # Simulate ASTRAL -t 2 output with q1/q2/q3 in node names
        newick = "((A:1,B:1)'[q1=0.5;q2=0.3;q3=0.2]':0.5,(C:1,D:1)'[q1=0.8;q2=0.1;q3=0.1]':0.3);"
        tree = Phylo.read(StringIO(newick), "newick")
        result = parse_astral_annotations(tree)
        assert len(result) >= 1  # At least one annotated node

    def test_no_annotations_returns_empty(self):
        from io import StringIO
        tree = Phylo.read(StringIO("((A:1,B:1):0.5,(C:1,D:1):0.3);"), "newick")
        result = parse_astral_annotations(tree)
        assert len(result) == 0


class TestParseQsFromLabel:
    def test_bracket_format(self):
        assert _parse_qs_from_label("[q1=0.5;q2=0.3;q3=0.2]") == (0.5, 0.3, 0.2)

    def test_quoted_bracket_format(self):
        assert _parse_qs_from_label("'[q1=0.5;q2=0.3;q3=0.2;f1=50]'") == (0.5, 0.3, 0.2)

    def test_bare_format(self):
        assert _parse_qs_from_label("q1=0.5;q2=0.3;q3=0.2") == (0.5, 0.3, 0.2)

    def test_missing_q3(self):
        assert _parse_qs_from_label("q1=0.5;q2=0.3") is None

    def test_empty(self):
        assert _parse_qs_from_label("") is None

    def test_garbage(self):
        assert _parse_qs_from_label("support=95") is None
