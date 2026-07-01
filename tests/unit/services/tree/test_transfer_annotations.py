from argparse import Namespace
from io import StringIO
import json
import subprocess
import sys

import pytest
from Bio import Phylo

from phykit.services.tree.transfer_annotations import TransferAnnotations


SOURCE = "tests/sample_files/wastral_annotated.tre"
TARGET = "tests/sample_files/raxml_unannotated.tre"


def _make_args(**overrides):
    defaults = dict(
        source=SOURCE, target=TARGET, output=None, json=False,
    )
    defaults.update(overrides)
    return Namespace(**defaults)


def test_module_import_does_not_import_biophylo_or_numpy():
    code = """
import sys
import phykit.services.tree.transfer_annotations as module

assert callable(module.print_json)
assert hasattr(module.Phylo, "read")
assert hasattr(module.Phylo, "write")
assert "json" not in sys.modules
assert "typing" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
assert "Bio.Phylo" not in sys.modules
assert "numpy" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


class TestTransferAnnotations:
    def test_run_reads_source_tree_unmodified(self, mocker, tmp_path):
        out = str(tmp_path / "out.tre")
        svc = TransferAnnotations(_make_args(output=out))
        target = object()
        source = object()
        tips = ["A", "B"]
        annotations = {frozenset({"A"}): "q1=1"}

        read_target = mocker.patch.object(svc, "read_tree_file", return_value=target)
        read_source = mocker.patch.object(
            svc, "_read_tree_with_error", return_value=source
        )
        mocker.patch.object(
            svc, "get_tip_names_from_tree", side_effect=[tips, tips]
        )
        mocker.patch.object(svc, "_extract_annotations", return_value=annotations)
        transfer = mocker.patch.object(svc, "_transfer", return_value=(1, 0))
        mocker.patch.object(svc, "_write_annotated_tree")
        mocker.patch("builtins.print")

        svc.run()

        read_target.assert_called_once_with()
        read_source.assert_called_once_with(SOURCE, "source_path", copy_tree=False)
        transfer.assert_called_once_with(target, annotations, frozenset(tips))

    def test_run_uses_fast_terminal_names_for_taxa_validation(
        self, monkeypatch, tmp_path, capsys
    ):
        svc = TransferAnnotations(_make_args(output=str(tmp_path / "out.tre")))
        target = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        source = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        source.root.clades[0].comment = "q1=1"

        def fail_get_terminals(*args, **kwargs):
            raise AssertionError("parsed trees should use fast terminal names")

        monkeypatch.setattr(target, "get_terminals", fail_get_terminals)
        monkeypatch.setattr(source, "get_terminals", fail_get_terminals)
        monkeypatch.setattr(svc, "read_tree_file", lambda: target)

        def read_tree_with_error(path, attr_name, copy_tree=True):
            assert path == SOURCE
            assert attr_name == "source_path"
            assert copy_tree is False
            return source

        monkeypatch.setattr(svc, "_read_tree_with_error", read_tree_with_error)

        svc.run()

        captured = capsys.readouterr()
        assert "Annotations transferred: 1" in captured.out

    def test_cached_bipartitions_match_uncached(self):
        svc = TransferAnnotations(_make_args())
        tree = Phylo.read(TARGET, "newick")
        all_taxa = frozenset(t.name for t in tree.get_terminals())
        clade_taxa = svc._collect_clade_taxa(tree)

        for clade in tree.find_clades(order="preorder"):
            uncached = svc._get_bipartition(clade, all_taxa)
            cached = svc._get_bipartition(clade, all_taxa, clade_taxa)
            assert cached == uncached

    def test_collect_clade_taxa_handles_unary_internal_nodes(self):
        svc = TransferAnnotations(_make_args())
        tree = Phylo.read(StringIO("((A:1):1,(B:1,C:1):1);"), "newick")

        clade_taxa = svc._collect_clade_taxa(tree)
        unary = tree.root.clades[0]
        binary = tree.root.clades[1]

        assert clade_taxa[id(unary)] == frozenset({"A"})
        assert clade_taxa[id(binary)] == frozenset({"B", "C"})
        assert clade_taxa[id(tree.root)] == frozenset({"A", "B", "C"})

    def test_collect_clade_taxa_multifurcation_avoids_temp_set(self, monkeypatch):
        svc = TransferAnnotations(_make_args())
        tree = Phylo.read(StringIO("(A:1,B:1,C:1,D:1);"), "newick")

        def fail_set(*_args, **_kwargs):
            raise AssertionError("multifurcation collection should not build a temp set")

        monkeypatch.setattr("builtins.set", fail_set)

        clade_taxa = svc._collect_clade_taxa(tree)

        assert clade_taxa[id(tree.root)] == frozenset({"A", "B", "C", "D"})

    def test_extract_and_transfer_use_direct_postorder_pass(self, monkeypatch):
        svc = TransferAnnotations(_make_args())
        source = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        target = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        source.root.clades[0].comment = "q1=1"
        all_taxa = frozenset(t.name for t in source.get_terminals())

        def fail_find_clades(*args, **kwargs):
            raise AssertionError("standard trees should use direct postorder")

        monkeypatch.setattr(source, "find_clades", fail_find_clades)
        monkeypatch.setattr(target, "find_clades", fail_find_clades)

        annotations = svc._extract_annotations(source, all_taxa)
        transferred, unmatched = svc._transfer(target, annotations, all_taxa)

        assert transferred == 1
        assert unmatched == 1

    def test_print_text_output_batches_summary(self, mocker):
        svc = TransferAnnotations.__new__(TransferAnnotations)
        printed = mocker.patch("builtins.print")

        svc._print_text_output(12, 3, "target.tre.annotated")

        printed.assert_called_once_with(
            "Annotations transferred: 12\n"
            "Unmatched nodes: 3\n"
            "Output: target.tre.annotated"
        )

    def test_basic_run(self, tmp_path, capsys):
        out = str(tmp_path / "out.tre")
        svc = TransferAnnotations(_make_args(output=out))
        svc.run()
        captured = capsys.readouterr()
        assert "Annotations transferred: 3" in captured.out
        assert "Unmatched nodes: 0" in captured.out

    def test_output_file_created(self, tmp_path):
        out = str(tmp_path / "out.tre")
        svc = TransferAnnotations(_make_args(output=out))
        svc.run()
        assert (tmp_path / "out.tre").exists()

    def test_annotations_present_in_output(self, tmp_path):
        out = str(tmp_path / "out.tre")
        svc = TransferAnnotations(_make_args(output=out))
        svc.run()
        tree = Phylo.read(out, "newick")
        annotated = [
            c for c in tree.find_clades()
            if not c.is_terminal() and c.comment and "q1=" in str(c.comment)
        ]
        assert len(annotated) == 3

    def test_branch_lengths_from_target(self, tmp_path):
        """Output should have target's branch lengths, not source's."""
        out = str(tmp_path / "out.tre")
        svc = TransferAnnotations(_make_args(output=out))
        svc.run()
        tree = Phylo.read(out, "newick")
        # Target has taxA:0.12, source has taxA:0.1
        for tip in tree.get_terminals():
            if tip.name == "taxA":
                assert tip.branch_length == pytest.approx(0.12, abs=0.01)
                break

    def test_all_taxa_preserved(self, tmp_path):
        out = str(tmp_path / "out.tre")
        svc = TransferAnnotations(_make_args(output=out))
        svc.run()
        tree = Phylo.read(out, "newick")
        tips = {t.name for t in tree.get_terminals()}
        assert tips == {"taxA", "taxB", "taxC", "taxD", "taxE", "taxF"}

    def test_json_output(self, tmp_path, capsys):
        out = str(tmp_path / "out.tre")
        svc = TransferAnnotations(_make_args(output=out, json=True))
        svc.run()
        captured = capsys.readouterr()
        payload = json.loads(captured.out)
        assert payload["annotations_transferred"] == 3
        assert payload["unmatched_nodes"] == 0
        assert payload["n_taxa"] == 6

    def test_mismatched_taxa_raises(self, tmp_path):
        """Trees with different taxa should raise error."""
        bad_target = tmp_path / "bad.tre"
        bad_target.write_text("((taxA:1,taxB:1):1,taxZ:1);\n")
        svc = TransferAnnotations(_make_args(target=str(bad_target)))
        with pytest.raises(SystemExit) as exc_info:
            svc.run()
        assert exc_info.value.code == 2

    def test_missing_source_raises(self):
        svc = TransferAnnotations(_make_args(source="nonexistent.tre"))
        with pytest.raises(SystemExit) as exc_info:
            svc.run()
        assert exc_info.value.code == 2

    def test_default_output_path(self, capsys):
        """Without -o, output should be target + '.annotated'."""
        svc = TransferAnnotations(_make_args())
        svc.run()
        captured = capsys.readouterr()
        assert ".annotated" in captured.out

    def test_annotation_format_preserved(self, tmp_path):
        """Annotations should use [key=val;...] format (not [&key=val])."""
        out = str(tmp_path / "out.tre")
        svc = TransferAnnotations(_make_args(output=out))
        svc.run()
        with open(out) as f:
            content = f.read()
        assert "[q1=" in content
        assert "[&" not in content

    def test_write_annotated_tree_formats_in_memory_before_single_file_write(
        self, monkeypatch, tmp_path
    ):
        seen = {}

        def fake_write(tree, handle, fmt):
            seen["tree"] = tree
            seen["handle"] = handle
            seen["format"] = fmt
            handle.write("(A:1)[&q1=1]\\[extra\\];\n")
            return 1

        monkeypatch.setattr(
            "phykit.services.tree.transfer_annotations.Phylo.write",
            fake_write,
        )
        output_path = tmp_path / "out.tre"
        tree = object()

        TransferAnnotations._write_annotated_tree(tree, str(output_path))

        assert seen["tree"] is tree
        assert hasattr(seen["handle"], "getvalue")
        assert seen["format"] == "newick"
        assert output_path.read_text() == "(A:1)[q1=1][extra];\n"
