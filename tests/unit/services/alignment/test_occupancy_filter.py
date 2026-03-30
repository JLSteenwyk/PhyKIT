from argparse import Namespace
import json

import pytest

from phykit.services.alignment.occupancy_filter import OccupancyFilter


SAMPLE = "tests/sample_files/occupancy_test"
ALN_LIST = f"{SAMPLE}/aln_list.txt"
TREE_LIST = f"{SAMPLE}/tree_list.txt"


def _make_args(**overrides):
    defaults = dict(
        list=ALN_LIST, format="fasta", threshold=0.5,
        output_dir=None, suffix=".filtered", json=False,
    )
    defaults.update(overrides)
    return Namespace(**defaults)


class TestOccupancyFilter:
    def test_basic_fasta_run(self, tmp_path, capsys):
        """Basic FASTA filtering runs without error."""
        svc = OccupancyFilter(_make_args(output_dir=str(tmp_path)))
        svc.run()
        captured = capsys.readouterr()
        assert "Occupancy filter" in captured.out

    def test_fasta_threshold_2(self, tmp_path, capsys):
        """Threshold=2 removes taxa in only 1 file."""
        svc = OccupancyFilter(_make_args(
            threshold=2, output_dir=str(tmp_path),
        ))
        svc.run()
        captured = capsys.readouterr()
        assert "Taxa kept: 4" in captured.out
        assert "Taxa removed: 2" in captured.out
        assert "dog" in captured.out
        assert "monkey" in captured.out

    def test_fasta_threshold_3(self, tmp_path, capsys):
        """Threshold=3 keeps only raccoon (present in all 3)."""
        svc = OccupancyFilter(_make_args(
            threshold=3, output_dir=str(tmp_path),
        ))
        svc.run()
        captured = capsys.readouterr()
        assert "Taxa kept: 1" in captured.out
        assert "raccoon" in captured.out

    def test_trees_format(self, tmp_path, capsys):
        """Tree format works."""
        svc = OccupancyFilter(_make_args(
            list=TREE_LIST, format="trees",
            threshold=2, output_dir=str(tmp_path),
        ))
        svc.run()
        captured = capsys.readouterr()
        assert "Taxa kept: 4" in captured.out

    def test_output_files_created_fasta(self, tmp_path):
        """Output FASTA files are created."""
        svc = OccupancyFilter(_make_args(
            threshold=2, output_dir=str(tmp_path),
        ))
        svc.run()
        assert (tmp_path / "aln1.filtered.fa").exists()
        assert (tmp_path / "aln2.filtered.fa").exists()
        assert (tmp_path / "aln3.filtered.fa").exists()

    def test_output_files_created_trees(self, tmp_path):
        """Output tree files are created."""
        svc = OccupancyFilter(_make_args(
            list=TREE_LIST, format="trees",
            threshold=2, output_dir=str(tmp_path),
        ))
        svc.run()
        assert (tmp_path / "tree1.filtered.nwk").exists()
        assert (tmp_path / "tree2.filtered.nwk").exists()
        assert (tmp_path / "tree3.filtered.nwk").exists()

    def test_filtered_fasta_has_correct_taxa(self, tmp_path):
        """Filtered FASTA should not contain removed taxa."""
        from Bio import SeqIO
        svc = OccupancyFilter(_make_args(
            threshold=2, output_dir=str(tmp_path),
        ))
        svc.run()
        # aln1 had raccoon, bear, cat, dog → dog removed
        records = list(SeqIO.parse(str(tmp_path / "aln1.filtered.fa"), "fasta"))
        taxa = {r.id for r in records}
        assert "dog" not in taxa
        assert "raccoon" in taxa
        assert "bear" in taxa
        assert "cat" in taxa

    def test_filtered_tree_has_correct_taxa(self, tmp_path):
        """Filtered tree should not contain removed taxa."""
        from Bio import Phylo
        svc = OccupancyFilter(_make_args(
            list=TREE_LIST, format="trees",
            threshold=2, output_dir=str(tmp_path),
        ))
        svc.run()
        tree = Phylo.read(str(tmp_path / "tree1.filtered.nwk"), "newick")
        tips = {t.name for t in tree.get_terminals()}
        assert "dog" not in tips
        assert "raccoon" in tips

    def test_json_output(self, tmp_path, capsys):
        """JSON output has correct structure."""
        svc = OccupancyFilter(_make_args(
            threshold=2, output_dir=str(tmp_path), json=True,
        ))
        svc.run()
        captured = capsys.readouterr()
        payload = json.loads(captured.out)
        assert payload["threshold"] == 2
        assert payload["taxa_kept"] == 4
        assert payload["taxa_removed"] == 2
        assert "dog" in payload["removed_taxa"]
        assert "raccoon" in payload["kept_taxa"]
        assert len(payload["files"]) == 3

    def test_custom_suffix(self, tmp_path):
        """Custom suffix is applied to output filenames."""
        svc = OccupancyFilter(_make_args(
            output_dir=str(tmp_path), suffix=".pruned",
        ))
        svc.run()
        assert (tmp_path / "aln1.pruned.fa").exists()

    def test_threshold_too_high_raises(self, tmp_path):
        """Threshold higher than all occupancies raises error."""
        svc = OccupancyFilter(_make_args(
            threshold=100, output_dir=str(tmp_path),
        ))
        with pytest.raises(SystemExit) as exc_info:
            svc.run()
        assert exc_info.value.code == 2

    def test_missing_list_file_raises(self):
        """Nonexistent list file raises error."""
        svc = OccupancyFilter(_make_args(list="nonexistent.txt"))
        with pytest.raises(SystemExit) as exc_info:
            svc.run()
        assert exc_info.value.code == 2
