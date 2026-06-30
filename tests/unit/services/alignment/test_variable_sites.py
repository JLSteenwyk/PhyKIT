import pytest
import subprocess
import sys
from argparse import Namespace
from math import isclose
from types import SimpleNamespace

from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from phykit.services.alignment.variable_sites import VariableSites
import phykit.services.alignment.variable_sites as variable_sites_module


@pytest.fixture
def args():
    kwargs = dict(alignment="/some/path/to/file.fa")
    return Namespace(**kwargs)


class TestVariableSites(object):
    def test_alignment_site_count_modules_defer_heavy_imports(self):
        modules = [
            "phykit.services.alignment.variable_sites",
            "phykit.services.alignment.parsimony_informative_sites",
            "phykit.services.alignment.gc_content",
            "phykit.services.alignment.alignment_length_no_gaps",
        ]
        code = f"""
import sys
modules = {modules!r}
for module_name in modules:
    module = __import__(module_name, fromlist=["*"])
    assert hasattr(module.np, "__getattr__")
    assert callable(module.print_json)

assert "json" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
assert "numpy" not in sys.modules
assert "Bio.Align" not in sys.modules
assert "Bio.AlignIO" not in sys.modules
"""
        subprocess.run([sys.executable, "-c", code], check=True)

    def test_variable_sites_import_does_not_import_typing(self):
        code = """
import sys
import phykit.services.alignment.variable_sites

assert "typing" not in sys.modules
"""
        subprocess.run([sys.executable, "-c", code], check=True)

    def test_init_sets_alignment_file_path(self, args):
        vs = VariableSites(args)
        assert vs.alignment_file_path == args.alignment
        assert vs.output_file_path is None

    def test_variable_sites(self, alignment_simple, args):
        vs = VariableSites(args)
        var_sites, aln_len, var_sites_per = vs.calculate_variable_sites(
            alignment_simple
        )
        assert isinstance(var_sites, int)
        assert isinstance(aln_len, int)
        assert isinstance(var_sites_per, float)
        assert var_sites == 4
        assert aln_len == 6
        assert isclose(var_sites_per, 66.66666666666666, rel_tol=0.001)

    def test_variable_sites_ignores_ambiguous_sites_and_uppercases(self, args):
        alignment = MultipleSeqAlignment([
            SeqRecord(Seq("aC-N"), id="t1"),
            SeqRecord(Seq("AGTN"), id="t2"),
            SeqRecord(Seq("AGTN"), id="t3"),
        ])
        vs = VariableSites(args)

        var_sites, aln_len, var_sites_per = vs.calculate_variable_sites(
            alignment,
            is_protein=False,
        )

        assert var_sites == 1
        assert aln_len == 4
        assert isclose(var_sites_per, 25.0, rel_tol=0.001)

    def test_variable_sites_ascii_matrix_path(self, args):
        alignment = MultipleSeqAlignment([
            SeqRecord(Seq("aa"), id="t1"),
            SeqRecord(Seq("at"), id="t2"),
            SeqRecord(Seq("tt"), id="t3"),
        ])
        vs = VariableSites(args)

        var_sites, aln_len, var_sites_per = vs.calculate_variable_sites(
            alignment,
            is_protein=False,
        )

        assert var_sites == 2
        assert aln_len == 2
        assert isclose(var_sites_per, 100.0, rel_tol=0.001)

    def test_variable_sites_ascii_matrix_path_avoids_isin(self, mocker, args):
        alignment = MultipleSeqAlignment([
            SeqRecord(Seq("aa"), id="t1"),
            SeqRecord(Seq("a-"), id="t2"),
            SeqRecord(Seq("tt"), id="t3"),
        ])
        vs = VariableSites(args)
        mocker.patch(
            "phykit.services.alignment.variable_sites.np.isin",
            side_effect=AssertionError("ASCII path should use the lookup table"),
        )

        var_sites, aln_len, var_sites_per = vs.calculate_variable_sites(
            alignment,
            is_protein=False,
        )

        assert var_sites == 2
        assert aln_len == 2
        assert isclose(var_sites_per, 100.0, rel_tol=0.001)

    def test_variable_sites_ascii_path_uses_gap_code_reduction(self, mocker, args):
        alignment = MultipleSeqAlignment([
            SeqRecord(Seq("aa"), id="t1"),
            SeqRecord(Seq("a-"), id="t2"),
            SeqRecord(Seq("tt"), id="t3"),
        ])
        vs = VariableSites(args)
        gap_codes_spy = mocker.spy(variable_sites_module, "_get_gap_codes")

        var_sites, aln_len, var_sites_per = vs.calculate_variable_sites(
            alignment,
            is_protein=False,
        )

        gap_codes_spy.assert_called_once_with(False)
        assert var_sites == 2
        assert aln_len == 2
        assert isclose(var_sites_per, 100.0, rel_tol=0.001)

    def test_variable_sites_ascii_path_uses_minmax_columns(self, mocker, args):
        alignment = MultipleSeqAlignment([
            SeqRecord(Seq("ACD-"), id="t1"),
            SeqRecord(Seq("ACE-"), id="t2"),
            SeqRecord(Seq("ACF-"), id="t3"),
        ])
        vs = VariableSites(args)
        mocker.patch(
            "phykit.services.alignment.variable_sites.np.unique",
            side_effect=AssertionError("ASCII path should not scan each unique symbol"),
        )

        var_sites, aln_len, var_sites_per = vs.calculate_variable_sites(
            alignment,
            is_protein=True,
        )

        assert var_sites == 1
        assert aln_len == 4
        assert isclose(var_sites_per, 25.0, rel_tol=0.001)

    def test_variable_sites_no_gap_ascii_skips_masked_minmax(self, mocker, args):
        alignment = MultipleSeqAlignment([
            SeqRecord(Seq("ACDE"), id="t1"),
            SeqRecord(Seq("ACDF"), id="t2"),
            SeqRecord(Seq("ACDG"), id="t3"),
        ])
        vs = VariableSites(args)
        mocker.patch(
            "phykit.services.alignment.variable_sites.np.where",
            side_effect=AssertionError("no-gap ASCII path should use direct min/max"),
        )

        var_sites, aln_len, var_sites_per = vs.calculate_variable_sites(
            alignment,
            is_protein=True,
        )

        assert var_sites == 1
        assert aln_len == 4
        assert isclose(var_sites_per, 25.0, rel_tol=0.001)

    def test_variable_sites_no_gap_dna_skips_masked_minmax(self, mocker, args):
        alignment = MultipleSeqAlignment([
            SeqRecord(Seq("ACGT"), id="t1"),
            SeqRecord(Seq("ACGA"), id="t2"),
            SeqRecord(Seq("ACGC"), id="t3"),
        ])
        vs = VariableSites(args)
        mocker.patch(
            "phykit.services.alignment.variable_sites.np.where",
            side_effect=AssertionError("no-gap DNA path should use direct min/max"),
        )

        var_sites, aln_len, var_sites_per = vs.calculate_variable_sites(
            alignment,
            is_protein=False,
        )

        assert var_sites == 1
        assert aln_len == 4
        assert isclose(var_sites_per, 25.0, rel_tol=0.001)

    def test_variable_sites_clean_ascii_skips_gap_code_setup(self, mocker, args):
        alignment = MultipleSeqAlignment([
            SeqRecord(Seq("ACGT"), id="t1"),
            SeqRecord(Seq("ACGA"), id="t2"),
            SeqRecord(Seq("ACGC"), id="t3"),
        ])
        vs = VariableSites(args)
        mocker.patch(
            "phykit.services.alignment.variable_sites._get_gap_codes",
            side_effect=AssertionError("clean ASCII path should not build gap codes"),
        )

        var_sites, aln_len, var_sites_per = vs.calculate_variable_sites(
            alignment,
            is_protein=False,
        )

        assert var_sites == 1
        assert aln_len == 4
        assert isclose(var_sites_per, 25.0, rel_tol=0.001)

    def test_variable_sites_identical_sequences_skip_matrix(self, mocker, args):
        alignment = MultipleSeqAlignment([
            SeqRecord(Seq("AcGt"), id="t1"),
            SeqRecord(Seq("aCgT"), id="t2"),
            SeqRecord(Seq("ACGT"), id="t3"),
        ])
        vs = VariableSites(args)
        mocker.patch(
            "phykit.services.alignment.variable_sites.np.frombuffer",
            side_effect=AssertionError("identical alignments should skip matrix setup"),
        )

        var_sites, aln_len, var_sites_per = vs.calculate_variable_sites(
            alignment,
            is_protein=False,
        )

        assert var_sites == 0
        assert aln_len == 4
        assert var_sites_per == 0.0

    def test_identical_sequence_helper_does_not_slice_rows(self):
        class NoSliceList(list):
            def __getitem__(self, key):
                if isinstance(key, slice):
                    raise AssertionError("identical-sequence scan should not slice")
                return super().__getitem__(key)

        sequences = NoSliceList(["ACGT", "ACGT", "ACGT"])

        assert variable_sites_module._all_sequences_identical(sequences) is True

    def test_variable_sites_unicode_fallback(self, args):
        class DummyAlignment(list):
            def get_alignment_length(self):
                return 2

        alignment = DummyAlignment([
            SimpleNamespace(seq="A\u00d1", id="t1"),
            SimpleNamespace(seq="A\u00d1", id="t2"),
            SimpleNamespace(seq="T\u00d1", id="t3"),
        ])
        vs = VariableSites(args)

        var_sites, aln_len, var_sites_per = vs.calculate_variable_sites(
            alignment,
            is_protein=False,
        )

        assert var_sites == 1
        assert aln_len == 2
        assert isclose(var_sites_per, 50.0, rel_tol=0.001)
