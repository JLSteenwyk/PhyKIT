from argparse import Namespace
import json

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


class TestTransferAnnotations:
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
