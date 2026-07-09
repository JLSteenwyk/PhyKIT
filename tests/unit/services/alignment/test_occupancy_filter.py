from argparse import Namespace
import json
import subprocess
import sys

from Bio.Phylo.BaseTree import TreeMixin
import pytest

from phykit.services.alignment.occupancy_filter import OccupancyFilter
import phykit.services.alignment.occupancy_filter as occupancy_filter_module
from phykit.services.tree.base import Tree


SAMPLE = "tests/sample_files/occupancy_test"
ALN_LIST = f"{SAMPLE}/aln_list.txt"
TREE_LIST = f"{SAMPLE}/tree_list.txt"


def test_module_import_does_not_import_biopython_fasta_parser():
    code = """
import sys
import phykit.services.alignment.occupancy_filter
assert "typing" not in sys.modules
assert "hashlib" not in sys.modules
assert "phykit.services.tree.base" not in sys.modules
assert "Bio.SeqIO.FastaIO" not in sys.modules
assert "Bio.AlignIO" not in sys.modules
assert "Bio.Phylo" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


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

    def test_run_batches_occupancy_counts_with_counter_update(
        self, tmp_path, monkeypatch
    ):
        class RecordingCounter(dict):
            def __init__(self):
                super().__init__()
                self.update_calls = []

            def update(self, iterable):
                taxa = tuple(iterable)
                self.update_calls.append(frozenset(taxa))
                for taxon in taxa:
                    self[taxon] = self.get(taxon, 0) + 1

        occupancy = RecordingCounter()
        svc = OccupancyFilter(_make_args(
            threshold=2,
            output_dir=str(tmp_path),
        ))
        taxa_by_path = {
            "one.fa": ["A", "B"],
            "two.fa": ["A", "C"],
            "three.fa": ["A", "B"],
        }
        printed = {}

        monkeypatch.setattr(
            occupancy_filter_module,
            "Counter",
            lambda: occupancy,
        )
        monkeypatch.setattr(
            svc,
            "_read_file_list",
            lambda _path: list(taxa_by_path),
        )
        monkeypatch.setattr(
            svc,
            "_extract_taxa",
            lambda path: taxa_by_path[path],
        )
        monkeypatch.setattr(
            svc,
            "_filter_fasta",
            lambda _path, _out_path, kept_taxa: (len(kept_taxa), 0),
        )
        monkeypatch.setattr(
            svc,
            "_print_text",
            lambda *args: printed.setdefault("args", args),
        )

        svc.run()

        assert occupancy == {"A": 3, "B": 2, "C": 1}
        assert occupancy.update_calls == [
            frozenset({"A", "B"}),
            frozenset({"A", "C"}),
            frozenset({"A", "B"}),
        ]
        assert printed["args"][5] == {
            "one.fa": [],
            "two.fa": ["C"],
            "three.fa": [],
        }

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

    def test_filter_fasta_streams_counts_and_preserves_formatting(self, tmp_path):
        """Direct FASTA filtering keeps SeqIO formatting and reports counts."""
        input_fasta = tmp_path / "input.fa"
        output_fasta = tmp_path / "output.fa"
        input_fasta.write_text(
            ">keep1 original description\n"
            f"{'A' * 65}\n"
            ">drop1 original description\n"
            "CCCC\n"
            ">keep2 another description\n"
            "GGGG\n"
        )
        svc = OccupancyFilter(_make_args(format="fasta"))

        kept, removed = svc._filter_fasta(
            str(input_fasta),
            str(output_fasta),
            {"keep1", "keep2"},
        )

        assert (kept, removed) == (2, 1)
        assert output_fasta.read_text() == (
            ">keep1 original description\n"
            f"{'A' * 60}\n"
            f"{'A' * 5}\n"
            ">keep2 another description\n"
            "GGGG\n"
        )

    def test_filter_fasta_rewraps_multiline_sequences(self, tmp_path):
        input_fasta = tmp_path / "input.fa"
        output_fasta = tmp_path / "output.fa"
        input_fasta.write_text(
            ">keep1 description\n"
            f"{'A' * 40} \n"
            f"{'C' * 30}\n"
            ">drop1 description\n"
            "GGGG\n"
        )
        svc = OccupancyFilter(_make_args(format="fasta"))

        kept, removed = svc._filter_fasta(
            str(input_fasta),
            str(output_fasta),
            {"keep1"},
        )

        assert (kept, removed) == (1, 1)
        assert output_fasta.read_text() == (
            ">keep1 description\n"
            f"{'A' * 40}{'C' * 20}\n"
            f"{'C' * 10}\n"
        )

    def test_write_wrapped_fasta_sequence_helper(self, tmp_path):
        output_fasta = tmp_path / "wrapped.fa"
        with open(output_fasta, "w") as handle:
            OccupancyFilter._write_wrapped_fasta_sequence(handle, "A" * 65)
            OccupancyFilter._write_wrapped_fasta_sequence(handle, "")

        assert output_fasta.read_text() == f"{'A' * 60}\n{'A' * 5}\n"

    def test_write_wrapped_fasta_sequence_short_sequence_skips_slicing(self, tmp_path):
        class NoSliceStr(str):
            def __getitem__(self, key):
                if isinstance(key, slice):
                    raise AssertionError("short sequences should be written directly")
                return super().__getitem__(key)

        output_fasta = tmp_path / "wrapped_short.fa"
        with open(output_fasta, "w") as handle:
            OccupancyFilter._write_wrapped_fasta_sequence(handle, NoSliceStr("ACGT"))

        assert output_fasta.read_text() == "ACGT\n"

    def test_write_wrapped_fasta_sequence_two_line_sequence(self, tmp_path):
        output_fasta = tmp_path / "wrapped_two_line.fa"
        with open(output_fasta, "w") as handle:
            OccupancyFilter._write_wrapped_fasta_sequence(handle, "A" * 120)
            OccupancyFilter._write_wrapped_fasta_sequence(handle, "C" * 61)

        assert output_fasta.read_text() == (
            f"{'A' * 60}\n"
            f"{'A' * 60}\n"
            f"{'C' * 60}\n"
            "C\n"
        )

    def test_write_wrapped_fasta_sequence_custom_width(self, tmp_path):
        output_fasta = tmp_path / "wrapped_custom.fa"
        with open(output_fasta, "w") as handle:
            OccupancyFilter._write_wrapped_fasta_sequence(handle, "ABCDEFG", width=3)

        assert output_fasta.read_text() == "ABC\nDEF\nG\n"

    def test_write_wrapped_fasta_sequence_batches_large_sequence(
        self, tmp_path, monkeypatch
    ):
        monkeypatch.setattr(
            occupancy_filter_module,
            "_WRAPPED_FASTA_BATCH_MIN_LENGTH",
            1,
        )
        monkeypatch.setattr(
            occupancy_filter_module,
            "_WRAPPED_FASTA_BATCH_CHUNKS",
            1,
        )

        output_fasta = tmp_path / "wrapped_batched.fa"
        with open(output_fasta, "w") as handle:
            OccupancyFilter._write_wrapped_fasta_sequence(
                handle,
                "ABCDEFGHIJKL",
                width=5,
            )

        assert output_fasta.read_text() == "ABCDE\nFGHIJ\nKL\n"

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

    def test_filter_tree_prunes_terminal_objects(self, tmp_path, monkeypatch):
        from Bio import Phylo

        input_tree = tmp_path / "input.nwk"
        output_tree = tmp_path / "output.nwk"
        input_tree.write_text("(A:1,(B:1,C:1):1,(D:1,E:1,F:1):1,G:1);")
        svc = OccupancyFilter(_make_args(format="trees"))

        def fail_get_terminals(*_args, **_kwargs):
            raise AssertionError("generic terminal traversal should not be used")

        def fail_prune(*_args, **_kwargs):
            raise AssertionError("per-tip prune should not be used")

        monkeypatch.setattr(TreeMixin, "get_terminals", fail_get_terminals)
        monkeypatch.setattr(TreeMixin, "prune", fail_prune)

        kept, removed = svc._filter_tree(
            str(input_tree),
            str(output_tree),
            {"A", "B", "D", "F"},
        )

        tree = Phylo.read(str(output_tree), "newick")
        assert (kept, removed) == (4, 3)
        assert set(Tree.calculate_terminal_names_fast(tree)) == {"A", "B", "D", "F"}

    def test_fasta_headers_use_first_whitespace_delimited_token(self, tmp_path):
        """FASTA taxon extraction should match SeqIO's record.id semantics."""
        fasta = tmp_path / "aln.fa"
        fasta.write_text(
            ">sp1 description text\nAT G\nC\n>sp2 another description\nGG G\n"
        )
        svc = OccupancyFilter(_make_args(format="fasta"))

        assert svc._extract_taxa(str(fasta)) == ["sp1", "sp2"]

    def test_extract_tree_taxa_uses_direct_terminal_name_traversal(
        self, tmp_path, monkeypatch
    ):
        tree_path = tmp_path / "tree.nwk"
        tree_path.write_text("((A,B),(C,D));\n")
        svc = OccupancyFilter(_make_args(format="trees"))

        def fail_get_terminals(*_args, **_kwargs):
            raise AssertionError("generic terminal traversal should not be used")

        monkeypatch.setattr(TreeMixin, "get_terminals", fail_get_terminals)

        assert svc._extract_taxa(str(tree_path)) == ["A", "B", "C", "D"]

    def test_extract_tree_taxa_scans_simple_newick_before_biopython(
        self, tmp_path, mocker
    ):
        tree_path = tmp_path / "tree.nwk"
        tree_path.write_text("((A:1,B:2):3,(C:4,D:5):6);\n")
        svc = OccupancyFilter(_make_args(format="trees"))
        mocker.patch(
            "Bio.Phylo.read",
            side_effect=AssertionError("simple Newick should avoid Bio.Phylo"),
        )

        assert svc._extract_taxa(str(tree_path)) == ["A", "B", "C", "D"]

    def test_extract_taxa_uses_os_exists_instead_of_path_per_file(
        self, monkeypatch, mocker
    ):
        svc = OccupancyFilter(_make_args(format="fasta"))
        mocked_exists = mocker.patch(
            "phykit.services.alignment.occupancy_filter._path_exists",
            return_value=True,
        )
        mocker.patch(
            "phykit.services.alignment.occupancy_filter.read_fasta_first_tokens",
            return_value=["A", "B"],
        )
        monkeypatch.setattr(
            "phykit.services.alignment.occupancy_filter.Path",
            lambda *_args, **_kwargs: (_ for _ in ()).throw(
                AssertionError("extract taxa should use os.path.exists")
            ),
        )

        assert svc._extract_taxa("alignment.fa") == ["A", "B"]
        mocked_exists.assert_called_once_with("alignment.fa")

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

    def test_print_text_outputs_single_compatible_report(self, mocker):
        svc = OccupancyFilter(_make_args(threshold=0.5))
        mocked_print = mocker.patch("builtins.print")

        svc._print_text(
            file_paths=["in1.fa", "in2.fa"],
            output_paths=["out1.fa", "out2.fa"],
            occupancy={"taxon_a": 2, "taxon_b": 1},
            kept_taxa={"taxon_a"},
            removed_taxa={"taxon_b"},
            taxa_removed_per_file={"in1.fa": ["taxon_b"], "in2.fa": []},
            total_files=2,
            files_modified=1,
            min_count=1,
        )

        mocked_print.assert_called_once_with(
            "\n".join(
                [
                    "Occupancy filter (threshold: 0.5 = 1/2 files)",
                    "Input files: 2",
                    "Total taxa: 2",
                    "Taxa kept: 1",
                    "Taxa removed: 1",
                    "Files modified: 1",
                    "",
                    "Removed taxa (occupancy):",
                    "  taxon_b: 1/2",
                    "",
                    "Occupancy per taxon:",
                    "  taxon_a: 2/2 (kept)",
                    "  taxon_b: 1/2 (REMOVED)",
                    "",
                    "Output files:",
                    "  out1.fa (1 taxa removed)",
                    "  out2.fa",
                ]
            )
        )

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

    def test_read_file_list_skips_comments_and_resolves_relative_paths(self, tmp_path):
        """List parsing skips comments/blanks and resolves relative paths."""
        list_file = tmp_path / "files.txt"
        list_file.write_text(
            "# comment\n\n"
            ".\n"
            "relative.fa\n"
            "nested//relative.fa\n"
            "./dot-relative.fa\n"
            "/subdir//absolute.fa"
        )
        svc = OccupancyFilter(_make_args(list=str(list_file)))

        assert svc._read_file_list(str(list_file)) == [
            str(tmp_path),
            str(tmp_path / "relative.fa"),
            str(tmp_path / "nested/relative.fa"),
            str(tmp_path / "dot-relative.fa"),
            "/subdir/absolute.fa",
        ]

    def test_read_file_list_avoids_per_line_path_objects(self, tmp_path, mocker):
        """Large file lists should not construct Path objects for every row."""
        list_file = tmp_path / "files.txt"
        list_file.write_text("one.fa\ntwo.fa\n/subdir/absolute.fa\n")
        path_calls = 0
        original_path = __import__(
            "phykit.services.alignment.occupancy_filter",
            fromlist=["Path"],
        ).Path

        def counting_path(*args, **kwargs):
            nonlocal path_calls
            path_calls += 1
            return original_path(*args, **kwargs)

        mocker.patch("phykit.services.alignment.occupancy_filter.Path", side_effect=counting_path)
        svc = OccupancyFilter(_make_args(list=str(list_file)))

        assert svc._read_file_list(str(list_file)) == [
            str(tmp_path / "one.fa"),
            str(tmp_path / "two.fa"),
            "/subdir/absolute.fa",
        ]
        assert path_calls == 1

    def test_read_file_list_skips_normalization_for_simple_paths(self, tmp_path, mocker):
        """Common simple paths avoid per-row normalization overhead."""
        list_file = tmp_path / "files.txt"
        list_file.write_text("one.fa\ntwo.fa\n/subdir/absolute.fa\n")
        mocker.patch(
            "phykit.services.alignment.occupancy_filter._normalize_list_path",
            side_effect=AssertionError("simple paths should not be normalized"),
        )
        svc = OccupancyFilter(_make_args(list=str(list_file)))

        assert svc._read_file_list(str(list_file)) == [
            str(tmp_path / "one.fa"),
            str(tmp_path / "two.fa"),
            "/subdir/absolute.fa",
        ]

    def test_missing_list_file_raises(self):
        """Nonexistent list file raises error."""
        svc = OccupancyFilter(_make_args(list="nonexistent.txt"))
        with pytest.raises(SystemExit) as exc_info:
            svc.run()
        assert exc_info.value.code == 2
