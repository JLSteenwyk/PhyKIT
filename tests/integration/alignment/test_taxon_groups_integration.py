import json
import pytest
import sys
from mock import patch
from pathlib import Path

from phykit.phykit import Phykit


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def _write_tree(path, newick):
    path.write_text(newick + "\n")


def _write_fasta(path, seqs):
    with open(path, "w") as f:
        for name, seq in seqs.items():
            f.write(f">{name}\n{seq}\n")


def _write_list(path, files):
    path.write_text("\n".join(str(f) for f in files) + "\n")


# ---------------------------------------------------------------------------
# tests
# ---------------------------------------------------------------------------

@pytest.mark.integration
class TestTaxonGroupsIntegration:

    @patch("builtins.print")
    def test_basic_invocation(self, mocked_print, tmp_path):
        """End-to-end test: three trees, two groups."""
        for i in range(2):
            _write_tree(tmp_path / f"tree{i}.nwk", "((A,B),(C,D));")
        _write_tree(tmp_path / "tree2.nwk", "(X,Y);")

        list_file = tmp_path / "trees.txt"
        _write_list(list_file, [tmp_path / f"tree{i}.nwk" for i in range(3)])

        testargs = [
            "phykit",
            "taxon_groups",
            "-l", str(list_file),
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        output = "\n".join(str(c) for c in mocked_print.call_args_list)
        assert "Total files: 3" in output
        assert "Total groups: 2" in output

    @patch("builtins.print")
    def test_alias_tgroups(self, mocked_print, tmp_path):
        """The 'tgroups' alias should work."""
        _write_tree(tmp_path / "tree0.nwk", "((A,B),(C,D));")
        list_file = tmp_path / "trees.txt"
        _write_list(list_file, [tmp_path / "tree0.nwk"])

        testargs = [
            "phykit",
            "tgroups",
            "-l", str(list_file),
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        output = "\n".join(str(c) for c in mocked_print.call_args_list)
        assert "Total files: 1" in output
        assert "Total groups: 1" in output

    @patch("builtins.print")
    def test_alias_shared_taxa(self, mocked_print, tmp_path):
        """The 'shared_taxa' alias should work."""
        _write_tree(tmp_path / "tree0.nwk", "((A,B),(C,D));")
        list_file = tmp_path / "trees.txt"
        _write_list(list_file, [tmp_path / "tree0.nwk"])

        testargs = [
            "phykit",
            "shared_taxa",
            "-l", str(list_file),
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        output = "\n".join(str(c) for c in mocked_print.call_args_list)
        assert "Total files: 1" in output

    @patch("builtins.print")
    def test_json_output(self, mocked_print, tmp_path):
        """JSON output should have correct structure."""
        for i in range(2):
            _write_tree(tmp_path / f"tree{i}.nwk", "((A,B),(C,D));")
        _write_tree(tmp_path / "tree2.nwk", "(X,Y);")

        list_file = tmp_path / "trees.txt"
        _write_list(list_file, [tmp_path / f"tree{i}.nwk" for i in range(3)])

        testargs = [
            "phykit",
            "taxon_groups",
            "-l", str(list_file),
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["total_files"] == 3
        assert payload["total_groups"] == 2
        assert len(payload["groups"]) == 2

    @patch("builtins.print")
    def test_fasta_format(self, mocked_print, tmp_path):
        """The -f fasta option should work."""
        _write_fasta(tmp_path / "aln0.fa", {"sp1": "ATG", "sp2": "ATG"})
        _write_fasta(tmp_path / "aln1.fa", {"sp1": "GGG", "sp2": "GGG"})
        _write_fasta(tmp_path / "aln2.fa", {"sp3": "CCC", "sp4": "CCC"})

        list_file = tmp_path / "alns.txt"
        _write_list(list_file, [
            tmp_path / "aln0.fa",
            tmp_path / "aln1.fa",
            tmp_path / "aln2.fa",
        ])

        testargs = [
            "phykit",
            "taxon_groups",
            "-l", str(list_file),
            "-f", "fasta",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        output = "\n".join(str(c) for c in mocked_print.call_args_list)
        assert "Total files: 3" in output
        assert "Total groups: 2" in output
