import pytest
from argparse import Namespace
from math import isclose
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from phykit.services.alignment.composition_per_taxon import CompositionPerTaxon
import phykit.services.alignment.composition_per_taxon as composition_per_taxon_module


@pytest.fixture
def args():
    kwargs = dict(alignment="/some/path/to/file.fa")
    return Namespace(**kwargs)


class TestCompositionPerTaxon(object):
    def test_init_sets_alignment_file_path(self, args):
        composition = CompositionPerTaxon(args)
        assert composition.alignment_file_path == args.alignment

    def test_composition_per_taxon(self, alignment_simple, args):
        composition = CompositionPerTaxon(args)
        symbols, rows = composition.calculate_composition_per_taxon(
            alignment_simple, is_protein=False
        )

        assert symbols == ["A", "C", "G", "T"]
        by_taxon = {taxon: vals for taxon, vals in rows}

        assert isclose(by_taxon["1"][0], 0.4, rel_tol=0.001)  # A
        assert isclose(by_taxon["1"][1], 0.0, rel_tol=0.001)  # C
        assert isclose(by_taxon["1"][2], 0.2, rel_tol=0.001)  # G
        assert isclose(by_taxon["1"][3], 0.4, rel_tol=0.001)  # T

        assert isclose(by_taxon["4"][0], 0.6, rel_tol=0.001)  # A
        assert isclose(by_taxon["4"][1], 0.0, rel_tol=0.001)  # C
        assert isclose(by_taxon["4"][2], 0.2, rel_tol=0.001)  # G
        assert isclose(by_taxon["4"][3], 0.2, rel_tol=0.001)  # T

    def test_process_args_defaults_json_false(self):
        parsed = CompositionPerTaxon(Namespace(alignment="x.fa")).process_args(
            Namespace(alignment="x.fa")
        )
        assert parsed["json_output"] is False

    def test_calculate_composition_per_taxon_all_invalid(self, args):
        svc = CompositionPerTaxon(args)
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("----"), id="t1"),
                SeqRecord(Seq("NNNN"), id="t2"),
            ]
        )
        symbols, rows = svc.calculate_composition_per_taxon(alignment, is_protein=False)
        assert symbols == []
        assert rows == []

    def test_calculate_composition_per_taxon_protein(self, args):
        svc = CompositionPerTaxon(args)
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("ACDX"), id="t1"),
                SeqRecord(Seq("AC?A"), id="t2"),
            ]
        )
        symbols, rows = svc.calculate_composition_per_taxon(alignment, is_protein=True)
        assert symbols == ["A", "C", "D"]
        assert len(rows) == 2

    def test_run_json_output(self, mocker):
        svc = CompositionPerTaxon(Namespace(alignment="x.fa", json=True))
        mocker.patch.object(
            CompositionPerTaxon,
            "get_alignment_and_format",
            return_value=(object(), "fasta", False),
        )
        mocker.patch.object(
            CompositionPerTaxon,
            "calculate_composition_per_taxon",
            return_value=(["A", "C"], [("t1", [0.5, 0.5])]),
        )
        mocked_json = mocker.patch.object(composition_per_taxon_module, "print_json")
        svc.run()
        payload = mocked_json.call_args.args[0]
        assert payload["symbols"] == ["A", "C"]
        assert payload["rows"] == payload["taxa"]
        assert payload["rows"][0]["composition"] == {"A": 0.5, "C": 0.5}

    def test_run_text_output(self, mocker, capsys):
        svc = CompositionPerTaxon(Namespace(alignment="x.fa", json=False))
        mocker.patch.object(
            CompositionPerTaxon,
            "get_alignment_and_format",
            return_value=(object(), "fasta", False),
        )
        mocker.patch.object(
            CompositionPerTaxon,
            "calculate_composition_per_taxon",
            return_value=(["A", "C"], [("t1", [0.5, 0.5]), ("t2", [1.0, 0.0])]),
        )
        svc.run()
        out, _ = capsys.readouterr()
        assert "t1\tA:0.5;C:0.5" in out
        assert "t2\tA:1.0;C:0.0" in out

    def test_run_returns_early_when_no_symbols(self, mocker, capsys):
        svc = CompositionPerTaxon(Namespace(alignment="x.fa", json=False))
        mocker.patch.object(
            CompositionPerTaxon,
            "get_alignment_and_format",
            return_value=(object(), "fasta", False),
        )
        mocker.patch.object(
            CompositionPerTaxon,
            "calculate_composition_per_taxon",
            return_value=([], []),
        )
        svc.run()
        out, _ = capsys.readouterr()
        assert out == ""
