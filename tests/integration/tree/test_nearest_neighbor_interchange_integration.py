import pytest
import sys
import json
from io import StringIO
from mock import patch
from pathlib import Path

from Bio import Phylo

from phykit.phykit import Phykit

here = Path(__file__)


def _bipartitions(tree):
    """Return frozenset of non-trivial bipartitions (sorted tip-name tuples)."""
    all_tips = frozenset(t.name for t in tree.get_terminals())
    biparts = set()
    for clade in tree.get_nonterminals():
        if clade is tree.root:
            continue
        a = frozenset(t.name for t in clade.get_terminals())
        b = all_tips - a
        if not a or not b:
            continue
        biparts.add(frozenset((a, b)))
    return frozenset(biparts)


def _read_trees(path):
    return list(Phylo.parse(path, "newick"))


@pytest.mark.integration
class TestNearestNeighborInterchange(object):
    @patch("builtins.print")
    def test_nni(self, mocked_print):
        testargs = [
            "phykit",
            "nearest_neighbor_interchange",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tre_rooted.tree",
        ]

        with patch.object(sys, "argv", testargs):
            Phykit()

        with open(f"{here.parent.parent}/expected/small_Aspergillus_tre_rooted.nni_moves.tre", "r") as expected_tree:
            expected_tree_content = expected_tree.read()

        with open(f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tre_rooted.tree.nnis", "r") as out_tree:
            out_tree_content = out_tree.read()

        assert expected_tree_content == out_tree_content

    @patch("builtins.print")
    def test_nni_custom_out(self, mocked_print):
        testargs = [
            "phykit",
            "nearest_neighbor_interchange",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tre_rooted.tree",
            "-o",
            f"{here.parent.parent.parent}/sample_files/nni_custom_out.tre"
        ]

        with patch.object(sys, "argv", testargs):
            Phykit()

        with open(f"{here.parent.parent}/expected/small_Aspergillus_tre_rooted.nni_moves.tre", "r") as expected_tree:
            expected_tree_content = expected_tree.read()

        with open(f"{here.parent.parent.parent}/sample_files/nni_custom_out.tre", "r") as out_tree:
            out_tree_content = out_tree.read()

        assert expected_tree_content == out_tree_content

    @patch("builtins.print")
    def test_nni_alias(self, mocked_print):
        testargs = [
            "phykit",
            "nni",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tre_rooted.tree",
        ]

        with patch.object(sys, "argv", testargs):
            Phykit()

        with open(f"{here.parent.parent}/expected/small_Aspergillus_tre_rooted.nni_moves.tre", "r") as expected_tree:
            expected_tree_content = expected_tree.read()

        with open(f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tre_rooted.tree.nnis", "r") as out_tree:
            out_tree_content = out_tree.read()

        assert expected_tree_content == out_tree_content

    @patch("builtins.print")
    def test_nni_wrong_path_tre(self, mocked_print):
        testargs = [
            "phykit",
            "nni",
            "bad_path",
        ]

        with patch.object(sys, "argv", testargs):
            with pytest.raises(SystemExit) as pytest_wrapped_e:
                Phykit()

        assert pytest_wrapped_e.type is SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_nni_json(self, mocked_print):
        testargs = [
            "phykit",
            "nearest_neighbor_interchange",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tre_rooted.tree",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["nni_neighbors"] == 14
        assert payload["total_trees"] == 15
        assert payload["output_file"].endswith(".nnis")

    @patch("builtins.print")
    def test_nni_targeted_single_branch(self, mocked_print, tmp_path):
        tree_path = tmp_path / "in.tre"
        tree_path.write_text("(((A,B),(C,D)),(E,F));\n")
        out_path = tmp_path / "out.tre"

        testargs = [
            "phykit",
            "nni",
            str(tree_path),
            "--branch",
            "A,B",
            "--no-input-tree",
            "-o",
            str(out_path),
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        trees = _read_trees(out_path)
        assert len(trees) == 2

        # Original (A,B) bipartition splits {A,B} | {C,D,E,F}. Both NNIs
        # around that branch must NOT contain that bipartition any more.
        original_biparts = _bipartitions(
            Phylo.read(StringIO("(((A,B),(C,D)),(E,F));"), "newick")
        )
        ab_split = frozenset(
            (frozenset({"A", "B"}), frozenset({"C", "D", "E", "F"}))
        )
        assert ab_split in original_biparts
        for tr in trees:
            tips = sorted(t.name for t in tr.get_terminals())
            assert tips == ["A", "B", "C", "D", "E", "F"]
            assert ab_split not in _bipartitions(tr)

        # The two NNIs must be distinct from each other.
        biparts_a = _bipartitions(trees[0])
        biparts_b = _bipartitions(trees[1])
        assert biparts_a != biparts_b

    @patch("builtins.print")
    def test_nni_targeted_includes_input_tree_by_default(
        self, mocked_print, tmp_path
    ):
        tree_path = tmp_path / "in.tre"
        tree_path.write_text("(((A,B),(C,D)),(E,F));\n")
        out_path = tmp_path / "out.tre"

        testargs = [
            "phykit",
            "nni",
            str(tree_path),
            "--branch",
            "A,B",
            "-o",
            str(out_path),
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        trees = _read_trees(out_path)
        assert len(trees) == 3  # input + 2 NNIs

    @patch("builtins.print")
    def test_nni_targeted_multiple_branches_via_file(
        self, mocked_print, tmp_path
    ):
        tree_path = tmp_path / "in.tre"
        tree_path.write_text("(((A,B),(C,D)),(E,F));\n")
        branches_path = tmp_path / "branches.tsv"
        branches_path.write_text(
            "# pairs defining branches of interest\n"
            "ab_clade\tA\tB\n"
            "cd_clade,C,D\n"
        )
        out_path = tmp_path / "out.tre"

        testargs = [
            "phykit",
            "nni",
            str(tree_path),
            "--branches",
            str(branches_path),
            "--no-input-tree",
            "-o",
            str(out_path),
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        trees = _read_trees(out_path)
        assert len(trees) == 4  # 2 NNIs per branch * 2 branches

        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["nni_neighbors"] == 4
        assert payload["total_trees"] == 4
        labels = [b["label"] for b in payload["branches"]]
        assert labels == ["ab_clade", "cd_clade"]
        for b in payload["branches"]:
            assert b["n_nnis"] == 2

    @patch("builtins.print")
    def test_nni_targeted_missing_taxon_errors(self, mocked_print, tmp_path):
        tree_path = tmp_path / "in.tre"
        tree_path.write_text("(((A,B),(C,D)),(E,F));\n")
        testargs = [
            "phykit",
            "nni",
            str(tree_path),
            "--branch",
            "A,ZZZ",
        ]
        with patch.object(sys, "argv", testargs):
            with pytest.raises(SystemExit) as exc:
                Phykit()
        assert exc.value.code == 2

    @patch("builtins.print")
    def test_nni_targeted_root_mrca_errors(self, mocked_print, tmp_path):
        tree_path = tmp_path / "in.tre"
        tree_path.write_text("(((A,B),(C,D)),(E,F));\n")
        testargs = [
            "phykit",
            "nni",
            str(tree_path),
            "--branch",
            "A,E",
        ]
        with patch.object(sys, "argv", testargs):
            with pytest.raises(SystemExit) as exc:
                Phykit()
        assert exc.value.code == 2

    @patch("builtins.print")
    def test_nni_targeted_bad_branch_format_errors(self, mocked_print, tmp_path):
        tree_path = tmp_path / "in.tre"
        tree_path.write_text("(((A,B),(C,D)),(E,F));\n")
        testargs = [
            "phykit",
            "nni",
            str(tree_path),
            "--branch",
            "A",
        ]
        with patch.object(sys, "argv", testargs):
            with pytest.raises(SystemExit) as exc:
                Phykit()
        assert exc.value.code == 2
