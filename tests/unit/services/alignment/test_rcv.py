import pytest
import subprocess
import sys
from argparse import Namespace
from math import isclose
from types import SimpleNamespace

from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from phykit.services.alignment.rcv import RelativeCompositionVariability
import phykit.services.alignment.rcv as rcv_module
import phykit.services.alignment.base as alignment_base_module


@pytest.fixture
def args():
    kwargs = dict(alignment="/some/path/to/file.fa")
    return Namespace(**kwargs)


class TestRelativeCompositionVariability(object):
    def test_alignment_base_import_does_not_import_numpy_or_biopython_align(self):
        code = """
import sys
import phykit.services.alignment.base as module
assert hasattr(module.np, "__getattr__")
assert callable(module.get_alignment_and_format_helper)
assert "typing" not in sys.modules
assert "numpy" not in sys.modules
assert "Bio.Align" not in sys.modules
assert "Bio.AlignIO" not in sys.modules
assert "phykit.helpers.files" not in sys.modules
"""
        subprocess.run([sys.executable, "-c", code], check=True)

    def test_alignment_json_output_modules_do_not_import_heavy_helpers(self):
        code = """
import importlib
import sys

modules = [
    "phykit.services.alignment.taxon_groups",
    "phykit.services.alignment.dna_threader",
    "phykit.services.alignment.occupancy_filter",
    "phykit.services.alignment.rcv",
]

for module_name in modules:
    module = importlib.import_module(module_name)
    assert callable(module.print_json)

assert "json" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
assert "Bio.SeqIO" not in sys.modules
assert "Bio.SeqIO.FastaIO" not in sys.modules
assert "Bio.Phylo" not in sys.modules
assert "numpy" not in sys.modules
"""
        subprocess.run([sys.executable, "-c", code], check=True)

    def test_init_sets_alignment_file_path(self, args):
        rcv = RelativeCompositionVariability(args)
        assert rcv.alignment_file_path == args.alignment
        assert rcv.output_file_path is None

    def test_relative_composition_variability(self, mocker, alignment_simple, args):
        mocker.patch(
            "phykit.services.alignment.rcv.RelativeCompositionVariability.get_alignment_and_format",
            return_value=(alignment_simple, "fa", True),
        )
        rcv = RelativeCompositionVariability(args)
        relative_composition_variability = rcv.calculate_rcv()
        assert isinstance(relative_composition_variability, float)
        assert isclose(relative_composition_variability, 0.292, rel_tol=0.001)

    def test_relative_composition_variability_handles_lowercase_ambiguity(self, mocker, args):
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("AnXT"), id="t1"),
                SeqRecord(Seq("ACxT"), id="t2"),
                SeqRecord(Seq("AG-T"), id="t3"),
            ]
        )
        rcv = RelativeCompositionVariability(args)

        mocked_get_alignment = mocker.patch.object(
            RelativeCompositionVariability,
            "get_alignment_and_format",
            return_value=(alignment, "fasta", False),
        )
        assert isclose(rcv.calculate_rcv(), 0.3333, rel_tol=0.001)

        mocked_get_alignment.return_value = (alignment, "fasta", True)
        assert isclose(rcv.calculate_rcv(), 0.4444, rel_tol=0.001)

    def test_relative_composition_variability_ascii_path_uses_lookup(self, mocker, args):
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("AnXT"), id="t1"),
                SeqRecord(Seq("ACxT"), id="t2"),
                SeqRecord(Seq("AG-T"), id="t3"),
            ]
        )
        rcv = RelativeCompositionVariability(args)
        mocker.patch.object(
            RelativeCompositionVariability,
            "get_alignment_and_format",
            return_value=(alignment, "fasta", False),
        )
        mocker.patch(
            "phykit.services.alignment.base.np.isin",
            side_effect=AssertionError("ASCII path should use the lookup table"),
        )

        assert isclose(rcv.calculate_rcv(), 0.3333, rel_tol=0.001)

    def test_relative_composition_variability_protein_ascii_uses_bincount(
        self, mocker, args
    ):
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("ACDEFGHIKLMNPQRSTVWY"), id="t1"),
                SeqRecord(Seq("ACDEFGHIKLMNPQRSTVWA"), id="t2"),
            ]
        )
        rcv = RelativeCompositionVariability(args)
        mocker.patch.object(
            RelativeCompositionVariability,
            "get_alignment_and_format",
            return_value=(alignment, "fasta", True),
        )
        bincount_spy = mocker.spy(alignment_base_module.np, "bincount")

        assert rcv.calculate_rcv() > 0.0
        assert bincount_spy.call_count == len(alignment)

    def test_relative_composition_variability_no_gap_ascii_uses_full_lengths(
        self, mocker, args
    ):
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("ACGT"), id="t1"),
                SeqRecord(Seq("AGGT"), id="t2"),
            ]
        )
        rcv = RelativeCompositionVariability(args)
        mocker.patch.object(
            RelativeCompositionVariability,
            "get_alignment_and_format",
            return_value=(alignment, "fasta", True),
        )
        full_spy = mocker.spy(alignment_base_module.np, "full")
        mocker.patch(
            "phykit.services.alignment.base.np.isin",
            side_effect=AssertionError("ASCII path should not use Unicode isin"),
        )

        assert rcv.calculate_rcv() > 0.0
        full_spy.assert_called_once_with(
            len(alignment),
            alignment.get_alignment_length(),
            dtype=alignment_base_module.np.float64,
        )

    def test_relative_composition_variability_identical_sequences_skip_matrix(
        self, mocker, args
    ):
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("ACGTACGT"), id="t1"),
                SeqRecord(Seq("acgtacgt"), id="t2"),
                SeqRecord(Seq("ACGTACGT"), id="t3"),
            ]
        )
        rcv = RelativeCompositionVariability(args)
        mocker.patch.object(
            RelativeCompositionVariability,
            "get_alignment_and_format",
            return_value=(alignment, "fasta", False),
        )
        mocker.patch(
            "phykit.services.alignment.base.np.frombuffer",
            side_effect=AssertionError("identical sequences should skip matrix work"),
        )

        assert rcv.calculate_rcv() == 0.0

    def test_identical_sequence_helper_does_not_slice_rows(self):
        class NoSliceList(list):
            def __getitem__(self, key):
                if isinstance(key, slice):
                    raise AssertionError("identical-sequence scan should not slice")
                return super().__getitem__(key)

        sequences = NoSliceList(["ACGT", "ACGT", "ACGT"])

        assert alignment_base_module._all_sequences_identical(sequences) is True

    def test_relative_composition_variability_unicode_fallback(self, mocker, args):
        class DummyAlignment(list):
            def get_alignment_length(self):
                return 2

        alignment = DummyAlignment(
            [
                SimpleNamespace(seq="A\u00d1", id="t1"),
                SimpleNamespace(seq="A\u00d1", id="t2"),
                SimpleNamespace(seq="T\u00d1", id="t3"),
            ]
        )
        rcv = RelativeCompositionVariability(args)
        mocker.patch.object(
            RelativeCompositionVariability,
            "get_alignment_and_format",
            return_value=(alignment, "fasta", False),
        )

        assert isclose(rcv.calculate_rcv(), 0.4444, rel_tol=0.001)

    def test_process_args_defaults_json_false(self):
        parsed = RelativeCompositionVariability(Namespace(alignment="x.fa")).process_args(
            Namespace(alignment="x.fa")
        )
        assert parsed["json_output"] is False

    def test_run_prints_rcv(self, mocker, capsys):
        rcv = RelativeCompositionVariability(Namespace(alignment="x.fa", json=False))
        mocker.patch.object(RelativeCompositionVariability, "calculate_rcv", return_value=0.123456)
        rcv.run()
        out, _ = capsys.readouterr()
        assert out.strip() == "0.1235"

    def test_run_json_output(self, mocker):
        rcv = RelativeCompositionVariability(Namespace(alignment="x.fa", json=True))
        mocker.patch.object(RelativeCompositionVariability, "calculate_rcv", return_value=0.56789)
        mocked_json = mocker.patch.object(rcv_module, "print_json")
        rcv.run()
        mocked_json.assert_called_once_with({"rcv": 0.5679})
