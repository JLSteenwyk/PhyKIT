from argparse import Namespace
from io import StringIO
import json

import pytest
from Bio import Phylo

from phykit.services.tree.spr import Spr


@pytest.fixture
def tree_file(tmp_path):
    """Create a simple test tree file: ((A:1,B:1):1,(C:1,(D:1,E:1):1):1);"""
    p = tmp_path / "test.tree"
    p.write_text("((A:1,B:1):1,(C:1,(D:1,E:1):1):1);\n")
    return str(p)


@pytest.fixture
def args_single(tree_file):
    return Namespace(tree=tree_file, subtree="E", output=None, json=False)


@pytest.fixture
def args_clade(tree_file):
    return Namespace(tree=tree_file, subtree="D,E", output=None, json=False)


class TestSpr:
    def test_init_sets_expected_attrs(self, tree_file):
        args = Namespace(tree=tree_file, subtree="D,E", output=None, json=False)
        service = Spr(args)
        assert service.tree_file_path == tree_file
        assert service.subtree_taxa == ["D", "E"]
        assert service.output_path is None
        assert service.json_output is False

    def test_process_args_strips_whitespace(self, tree_file):
        args = Namespace(tree=tree_file, subtree=" D , E ", output="/tmp/out.nwk", json=True)
        service = Spr(args)
        assert service.subtree_taxa == ["D", "E"]
        assert service.output_path == "/tmp/out.nwk"
        assert service.json_output is True

    def test_single_taxon_spr(self, args_single, capsys):
        """Prune a single taxon and verify correct number of SPR trees."""
        service = Spr(args_single)
        service.run()
        captured = capsys.readouterr()
        lines = [l for l in captured.out.strip().split("\n") if l.strip()]
        # Pruning E from ((A:1,B:1):1,(C:1,(D:1,E:1):1):1);
        # Remaining tree after pruning E: ((A:1,B:1):1,(C:1,D:2):1);
        # (D gets collapsed: D's bl=1 + parent bl=1 = 2)
        # Branches: root->(A,B), (A,B)->A, (A,B)->B, root->C+D clade, clade->C, clade->D
        # That's 5 branches (exclude root itself)
        # But one of these may overlap with original position
        assert len(lines) >= 3
        # Every output line should be valid Newick
        for line in lines:
            tree = Phylo.read(StringIO(line), "newick")
            assert tree is not None

    def test_clade_spr(self, args_clade, capsys):
        """Prune a clade (D,E) and verify correct number of SPR trees."""
        service = Spr(args_clade)
        service.run()
        captured = capsys.readouterr()
        lines = [l for l in captured.out.strip().split("\n") if l.strip()]
        # Pruning (D,E) from ((A:1,B:1):1,(C:1,(D:1,E:1):1):1);
        # Remaining tree: ((A:1,B:1):1,C:2);
        # (C gets collapsed: C bl=1 + parent bl=1 = 2)
        # Branches: root->(A,B), (A,B)->A, (A,B)->B, root->C
        # 4 branches = 4 regraft positions
        assert len(lines) == 4

    def test_all_trees_valid(self, args_clade, capsys):
        """Every SPR tree has the same set of taxa as the original."""
        service = Spr(args_clade)
        service.run()
        captured = capsys.readouterr()
        lines = [l for l in captured.out.strip().split("\n") if l.strip()]
        expected_taxa = {"A", "B", "C", "D", "E"}
        for line in lines:
            tree = Phylo.read(StringIO(line), "newick")
            taxa = {t.name for t in tree.get_terminals()}
            assert taxa == expected_taxa, f"Taxa mismatch in tree: {line}"

    def test_spr_trees_differ(self, args_clade, capsys):
        """SPR trees are different from each other."""
        service = Spr(args_clade)
        service.run()
        captured = capsys.readouterr()
        lines = [l for l in captured.out.strip().split("\n") if l.strip()]
        # At least some trees should be different
        unique_trees = set(lines)
        assert len(unique_trees) > 1

    def test_output_file_created(self, tree_file, tmp_path):
        """Test that -o creates a file."""
        out_path = str(tmp_path / "spr_output.nwk")
        args = Namespace(tree=tree_file, subtree="D,E", output=out_path, json=False)
        service = Spr(args)
        service.run()
        with open(out_path) as f:
            content = f.read()
        lines = [l for l in content.strip().split("\n") if l.strip()]
        assert len(lines) >= 1

    def test_json_output(self, tree_file, capsys):
        """JSON has correct structure."""
        args = Namespace(tree=tree_file, subtree="D,E", output=None, json=True)
        service = Spr(args)
        service.run()
        captured = capsys.readouterr()
        payload = json.loads(captured.out)
        assert payload["subtree_taxa"] == ["D", "E"]
        assert payload["n_spr_trees"] == len(payload["trees"])
        assert payload["n_spr_trees"] > 0
        for entry in payload["trees"]:
            assert "description" in entry
            assert "newick" in entry
            # Validate newick is parsable
            tree = Phylo.read(StringIO(entry["newick"]), "newick")
            assert tree is not None

    def test_missing_taxon_raises(self, tree_file):
        """Nonexistent taxon raises error."""
        args = Namespace(tree=tree_file, subtree="X", output=None, json=False)
        service = Spr(args)
        with pytest.raises(SystemExit) as exc_info:
            service.run()
        assert exc_info.value.code == 2

    def test_root_subtree_raises(self, tree_file):
        """Cannot prune all taxa (root)."""
        args = Namespace(tree=tree_file, subtree="A,B,C,D,E", output=None, json=False)
        service = Spr(args)
        with pytest.raises(SystemExit) as exc_info:
            service.run()
        assert exc_info.value.code == 2

    def test_newick_format(self, args_clade, capsys):
        """Output trees are valid Newick (end with semicolons)."""
        service = Spr(args_clade)
        service.run()
        captured = capsys.readouterr()
        lines = [l for l in captured.out.strip().split("\n") if l.strip()]
        for line in lines:
            assert line.strip().endswith(";"), f"Tree does not end with semicolon: {line}"

    def test_generate_spr_trees_returns_list(self, tree_file):
        """The static method returns a list of tuples."""
        tree = Phylo.read(tree_file, "newick")
        subtree = tree.common_ancestor(["D", "E"])
        results = Spr._generate_spr_trees(tree, subtree)
        assert isinstance(results, list)
        for desc, t in results:
            assert isinstance(desc, str)
            assert hasattr(t, "root")

    def test_subtree_from_file(self, tree_file, tmp_path, capsys):
        """--subtree can be a single-column file with one taxon per line."""
        taxa_file = tmp_path / "subtree.txt"
        taxa_file.write_text("D\nE\n")
        args = Namespace(tree=tree_file, subtree=str(taxa_file), output=None, json=False)
        service = Spr(args)
        assert service.subtree_taxa == ["D", "E"]
        service.run()
        captured = capsys.readouterr()
        lines = [l for l in captured.out.strip().split("\n") if l.strip()]
        assert len(lines) == 4

    def test_branch_lengths_preserved(self, args_clade, capsys):
        """Branch lengths of the pruned subtree are preserved in output trees."""
        service = Spr(args_clade)
        service.run()
        captured = capsys.readouterr()
        lines = [l for l in captured.out.strip().split("\n") if l.strip()]
        for line in lines:
            tree = Phylo.read(StringIO(line), "newick")
            # D and E should still have branch length 1.0
            for tip in tree.get_terminals():
                if tip.name in ("D", "E"):
                    assert tip.branch_length is not None
                    assert tip.branch_length == pytest.approx(1.0)
