import json
import os
import pytest
from argparse import Namespace
from unittest.mock import patch

from phykit.services.alignment.taxon_groups import TaxonGroups
from phykit.errors import PhykitUserError


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def _write_tree(path, newick):
    """Write a Newick tree string to a file."""
    path.write_text(newick + "\n")


def _write_fasta(path, seqs):
    """Write a dict of {name: seq} as FASTA."""
    with open(path, "w") as f:
        for name, seq in seqs.items():
            f.write(f">{name}\n{seq}\n")


def _write_list(path, files):
    """Write a list of file paths (one per line)."""
    path.write_text("\n".join(str(f) for f in files) + "\n")


def _make_args(list_path, fmt="trees", json_output=False):
    return Namespace(list=str(list_path), format=fmt, json=json_output)


# ---------------------------------------------------------------------------
# tests — tree mode
# ---------------------------------------------------------------------------

class TestTreeMode:
    def test_identical_taxa_one_group(self, tmp_path):
        """Three trees with the same taxa should form one group."""
        for i in range(3):
            _write_tree(tmp_path / f"tree{i}.nwk", "((A,B),(C,D));")
        list_file = tmp_path / "trees.txt"
        _write_list(list_file, [tmp_path / f"tree{i}.nwk" for i in range(3)])

        args = _make_args(list_file)
        svc = TaxonGroups(args)

        with patch("builtins.print") as mock_print:
            svc.run()

        output = "\n".join(str(c) for c in mock_print.call_args_list)
        assert "Total groups: 1" in output
        assert "Total files: 3" in output

    def test_different_taxa_multiple_groups(self, tmp_path):
        """Trees with different taxa should form separate groups."""
        _write_tree(tmp_path / "tree0.nwk", "((A,B),(C,D));")
        _write_tree(tmp_path / "tree1.nwk", "((A,B),(C,D));")
        _write_tree(tmp_path / "tree2.nwk", "((X,Y),(Z,W));")

        list_file = tmp_path / "trees.txt"
        _write_list(list_file, [
            tmp_path / "tree0.nwk",
            tmp_path / "tree1.nwk",
            tmp_path / "tree2.nwk",
        ])

        args = _make_args(list_file)
        svc = TaxonGroups(args)

        with patch("builtins.print") as mock_print:
            svc.run()

        output = "\n".join(str(c) for c in mock_print.call_args_list)
        assert "Total groups: 2" in output
        assert "Total files: 3" in output

    def test_singleton_group(self, tmp_path):
        """A tree with unique taxa should be its own group."""
        _write_tree(tmp_path / "tree0.nwk", "((A,B),(C,D));")
        _write_tree(tmp_path / "tree1.nwk", "((E,F),(G,H));")

        list_file = tmp_path / "trees.txt"
        _write_list(list_file, [
            tmp_path / "tree0.nwk",
            tmp_path / "tree1.nwk",
        ])

        args = _make_args(list_file)
        svc = TaxonGroups(args)

        with patch("builtins.print") as mock_print:
            svc.run()

        output = "\n".join(str(c) for c in mock_print.call_args_list)
        assert "Total groups: 2" in output

    def test_sorted_by_size(self, tmp_path):
        """The largest group should be reported first."""
        # 3 trees with taxa {A,B,C,D}, 1 tree with {X,Y}
        for i in range(3):
            _write_tree(tmp_path / f"tree{i}.nwk", "((A,B),(C,D));")
        _write_tree(tmp_path / "tree3.nwk", "(X,Y);")

        list_file = tmp_path / "trees.txt"
        _write_list(list_file, [tmp_path / f"tree{i}.nwk" for i in range(4)])

        args = _make_args(list_file)
        svc = TaxonGroups(args)

        with patch("builtins.print") as mock_print:
            svc.run()

        output = "\n".join(str(c) for c in mock_print.call_args_list)
        # Group 1 should have 3 files (the larger group)
        assert "Group 1 (4 taxa, 3 files):" in output
        assert "Group 2 (2 taxa, 1 files):" in output

    def test_taxa_reported_in_output(self, tmp_path):
        """Taxa names should appear in the text output."""
        _write_tree(tmp_path / "tree0.nwk", "((Alpha,Beta),(Gamma,Delta));")

        list_file = tmp_path / "trees.txt"
        _write_list(list_file, [tmp_path / "tree0.nwk"])

        args = _make_args(list_file)
        svc = TaxonGroups(args)

        with patch("builtins.print") as mock_print:
            svc.run()

        output = "\n".join(str(c) for c in mock_print.call_args_list)
        assert "Alpha" in output
        assert "Beta" in output
        assert "Gamma" in output
        assert "Delta" in output


# ---------------------------------------------------------------------------
# tests — fasta mode
# ---------------------------------------------------------------------------

class TestFastaMode:
    def test_fasta_mode(self, tmp_path):
        """FASTA files should be grouped by their sequence IDs."""
        _write_fasta(tmp_path / "aln0.fa", {"sp1": "ATG", "sp2": "ATG", "sp3": "ATG"})
        _write_fasta(tmp_path / "aln1.fa", {"sp1": "GGG", "sp2": "GGG", "sp3": "GGG"})
        _write_fasta(tmp_path / "aln2.fa", {"sp1": "CCC", "sp4": "CCC"})

        list_file = tmp_path / "alns.txt"
        _write_list(list_file, [
            tmp_path / "aln0.fa",
            tmp_path / "aln1.fa",
            tmp_path / "aln2.fa",
        ])

        args = _make_args(list_file, fmt="fasta")
        svc = TaxonGroups(args)

        with patch("builtins.print") as mock_print:
            svc.run()

        output = "\n".join(str(c) for c in mock_print.call_args_list)
        assert "Total groups: 2" in output
        assert "Total files: 3" in output


# ---------------------------------------------------------------------------
# tests — JSON output
# ---------------------------------------------------------------------------

class TestJsonOutput:
    def test_json_output(self, tmp_path):
        """JSON output should have the correct structure."""
        for i in range(2):
            _write_tree(tmp_path / f"tree{i}.nwk", "((A,B),(C,D));")
        _write_tree(tmp_path / "tree2.nwk", "(X,Y);")

        list_file = tmp_path / "trees.txt"
        _write_list(list_file, [tmp_path / f"tree{i}.nwk" for i in range(3)])

        args = _make_args(list_file, json_output=True)
        svc = TaxonGroups(args)

        with patch("builtins.print") as mock_print:
            svc.run()

        payload = json.loads(mock_print.call_args.args[0])
        assert payload["total_files"] == 3
        assert payload["total_groups"] == 2
        assert len(payload["groups"]) == 2

        # Largest group first
        g1 = payload["groups"][0]
        assert g1["group"] == 1
        assert g1["n_files"] == 2
        assert g1["n_taxa"] == 4
        assert sorted(g1["taxa"]) == ["A", "B", "C", "D"]

        g2 = payload["groups"][1]
        assert g2["group"] == 2
        assert g2["n_files"] == 1
        assert g2["n_taxa"] == 2
        assert sorted(g2["taxa"]) == ["X", "Y"]


# ---------------------------------------------------------------------------
# tests — error handling
# ---------------------------------------------------------------------------

class TestErrorHandling:
    def test_empty_list_raises(self, tmp_path):
        """An empty list file should raise an error."""
        list_file = tmp_path / "empty.txt"
        list_file.write_text("\n")

        args = _make_args(list_file)
        svc = TaxonGroups(args)

        with pytest.raises(PhykitUserError):
            svc.run()

    def test_missing_file_raises(self, tmp_path):
        """A nonexistent file in the list should raise an error."""
        list_file = tmp_path / "trees.txt"
        list_file.write_text(str(tmp_path / "nonexistent.nwk") + "\n")

        args = _make_args(list_file)
        svc = TaxonGroups(args)

        with pytest.raises(PhykitUserError):
            svc.run()

    def test_missing_list_file_raises(self, tmp_path):
        """A nonexistent list file should raise an error."""
        args = _make_args(tmp_path / "no_such_file.txt")
        svc = TaxonGroups(args)

        with pytest.raises(PhykitUserError):
            svc.run()

    def test_comments_skipped(self, tmp_path):
        """Lines starting with # should be ignored."""
        _write_tree(tmp_path / "tree0.nwk", "((A,B),(C,D));")

        list_file = tmp_path / "trees.txt"
        list_file.write_text(
            "# this is a comment\n"
            + str(tmp_path / "tree0.nwk") + "\n"
            + "# another comment\n"
        )

        args = _make_args(list_file)
        svc = TaxonGroups(args)

        with patch("builtins.print") as mock_print:
            svc.run()

        output = "\n".join(str(c) for c in mock_print.call_args_list)
        assert "Total files: 1" in output
        assert "Total groups: 1" in output
