import json
import os
import pytest
import subprocess
import sys
from argparse import Namespace
from unittest.mock import patch

from Bio.Phylo.BaseTree import TreeMixin

from phykit.services.alignment.taxon_groups import TaxonGroups
from phykit.errors import PhykitUserError


def test_module_import_does_not_import_tree_or_biopython_helpers():
    code = """
import sys
import phykit.services.alignment.taxon_groups
assert "hashlib" not in sys.modules
assert "phykit.services.tree.base" not in sys.modules
assert "Bio.SeqIO.FastaIO" not in sys.modules
assert "Bio.AlignIO" not in sys.modules
assert "Bio.Phylo" not in sys.modules
assert "json" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


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

    def test_text_output_is_single_compatible_report(self, tmp_path):
        _write_tree(tmp_path / "tree0.nwk", "((A,B),(C,D));")

        list_file = tmp_path / "trees.txt"
        _write_list(list_file, [tmp_path / "tree0.nwk"])

        args = _make_args(list_file)
        svc = TaxonGroups(args)

        with patch("builtins.print") as mock_print:
            svc.run()

        mock_print.assert_called_once_with(
            "\n".join(
                [
                    "Total files: 1",
                    "Total groups: 1",
                    "",
                    "Group 1 (4 taxa, 1 files):",
                    f"  Files: {tmp_path / 'tree0.nwk'}",
                    "  Taxa: A, B, C, D",
                    "",
                ]
            )
        )

    def test_run_groups_repeated_taxon_sets_from_extraction(
        self, tmp_path, monkeypatch
    ):
        svc = TaxonGroups(_make_args(tmp_path / "unused.txt"))
        monkeypatch.setattr(
            svc,
            "_read_file_list",
            lambda _path: ["tree0.nwk", "tree1.nwk", "tree2.nwk"],
        )
        taxon_sets = {
            "tree0.nwk": ["A", "B"],
            "tree1.nwk": ["A", "B"],
            "tree2.nwk": ["C", "D"],
        }
        monkeypatch.setattr(svc, "_extract_taxa", taxon_sets.__getitem__)

        with patch("builtins.print") as mock_print:
            svc.run()

        output = mock_print.call_args.args[0]
        assert "Total files: 3" in output
        assert "Total groups: 2" in output
        assert "Group 1 (2 taxa, 2 files):" in output
        assert "  Files: tree0.nwk, tree1.nwk" in output
        assert "Group 2 (2 taxa, 1 files):" in output

    def test_equal_size_groups_keep_first_seen_order(self, tmp_path, monkeypatch):
        svc = TaxonGroups(_make_args(tmp_path / "unused.txt"))
        monkeypatch.setattr(
            svc,
            "_read_file_list",
            lambda _path: ["first.nwk", "second.nwk", "third.nwk"],
        )
        taxon_sets = {
            "first.nwk": ["A", "B"],
            "second.nwk": ["C", "D"],
            "third.nwk": ["E", "F"],
        }
        monkeypatch.setattr(svc, "_extract_taxa", taxon_sets.__getitem__)

        with patch("builtins.print") as mock_print:
            svc.run()

        output = mock_print.call_args.args[0]
        assert output.index("  Files: first.nwk") < output.index("  Files: second.nwk")
        assert output.index("  Files: second.nwk") < output.index("  Files: third.nwk")

    def test_extract_tree_taxa_uses_direct_terminal_name_traversal(
        self, tmp_path, monkeypatch
    ):
        tree_path = tmp_path / "tree.nwk"
        _write_tree(tree_path, "((A,B),(C,D));")
        svc = TaxonGroups(_make_args(tmp_path / "unused.txt"))

        def fail_get_terminals(*_args, **_kwargs):
            raise AssertionError("generic terminal traversal should not be used")

        monkeypatch.setattr(TreeMixin, "get_terminals", fail_get_terminals)

        assert svc._extract_taxa(str(tree_path)) == ["A", "B", "C", "D"]

    def test_extract_taxa_uses_os_exists_instead_of_path_per_file(
        self, tmp_path, monkeypatch, mocker
    ):
        svc = TaxonGroups(_make_args(tmp_path / "unused.txt", fmt="fasta"))
        mocked_exists = mocker.patch(
            "phykit.services.alignment.taxon_groups._path_exists",
            return_value=True,
        )
        mocker.patch(
            "phykit.services.alignment.taxon_groups.read_fasta_first_tokens",
            return_value=["A", "B"],
        )
        monkeypatch.setattr(
            "phykit.services.alignment.taxon_groups.Path",
            lambda *_args, **_kwargs: (_ for _ in ()).throw(
                AssertionError("extract taxa should use os.path.exists")
            ),
        )

        assert svc._extract_taxa("alignment.fa") == ["A", "B"]
        mocked_exists.assert_called_once_with("alignment.fa")


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

    def test_fasta_headers_use_first_whitespace_delimited_token(self, tmp_path):
        """FASTA mode should match SeqIO's record.id header semantics."""
        fasta = tmp_path / "aln.fa"
        fasta.write_text(
            ">sp1 description text\nAT G\nC\n>sp2 another description\nGG G\n"
        )

        args = _make_args(tmp_path / "unused.txt", fmt="fasta")
        svc = TaxonGroups(args)

        assert svc._extract_taxa(str(fasta)) == ["sp1", "sp2"]


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
    def test_read_file_list_skips_comments_and_resolves_relative_paths(self, tmp_path):
        """List parsing should skip comments/blanks and resolve relative paths."""
        list_file = tmp_path / "trees.txt"
        list_file.write_text(
            "# comment\n\n"
            ".\n"
            "relative.nwk\n"
            "nested//relative.nwk\n"
            "./dot-relative.nwk\n"
            "/subdir//absolute.nwk"
        )

        args = _make_args(list_file)
        svc = TaxonGroups(args)

        assert svc._read_file_list(str(list_file)) == [
            str(tmp_path),
            str(tmp_path / "relative.nwk"),
            str(tmp_path / "nested/relative.nwk"),
            str(tmp_path / "dot-relative.nwk"),
            "/subdir/absolute.nwk",
        ]

    def test_read_file_list_avoids_per_line_path_objects(self, tmp_path, mocker):
        """Large file lists should not construct Path objects for every row."""
        list_file = tmp_path / "trees.txt"
        list_file.write_text("one.nwk\ntwo.nwk\n/subdir/absolute.nwk\n")
        path_calls = 0
        original_path = __import__(
            "phykit.services.alignment.taxon_groups",
            fromlist=["Path"],
        ).Path

        def counting_path(*args, **kwargs):
            nonlocal path_calls
            path_calls += 1
            return original_path(*args, **kwargs)

        mocker.patch("phykit.services.alignment.taxon_groups.Path", side_effect=counting_path)

        args = _make_args(list_file)
        svc = TaxonGroups(args)

        assert svc._read_file_list(str(list_file)) == [
            str(tmp_path / "one.nwk"),
            str(tmp_path / "two.nwk"),
            "/subdir/absolute.nwk",
        ]
        assert path_calls == 1

    def test_read_file_list_skips_normalization_for_simple_paths(self, tmp_path, mocker):
        """Common simple paths should avoid per-row normalization overhead."""
        list_file = tmp_path / "trees.txt"
        list_file.write_text("one.nwk\ntwo.nwk\n/subdir/absolute.nwk\n")
        mocker.patch(
            "phykit.services.alignment.taxon_groups._normalize_list_path",
            side_effect=AssertionError("simple paths should not be normalized"),
        )

        args = _make_args(list_file)
        svc = TaxonGroups(args)

        assert svc._read_file_list(str(list_file)) == [
            str(tmp_path / "one.nwk"),
            str(tmp_path / "two.nwk"),
            "/subdir/absolute.nwk",
        ]

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
