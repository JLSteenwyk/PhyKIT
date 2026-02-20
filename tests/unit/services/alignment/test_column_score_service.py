from argparse import Namespace

import numpy as np
import pytest
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from phykit.services.alignment.column_score import ColumnScore


@pytest.fixture
def args():
    return Namespace(fasta="/some/path/to/query.fa", reference="/some/path/to/ref.fa")


class TestColumnScore:
    def test_init_sets_expected_attrs(self, args):
        service = ColumnScore(args)
        assert service.fasta == args.fasta
        assert service.reference == args.reference
        assert service.json_output is False

    def test_get_columns_from_alignments(self, args):
        service = ColumnScore(args)
        reference = MultipleSeqAlignment(
            [SeqRecord(Seq("AC"), id="r1"), SeqRecord(Seq("GT"), id="r2")]
        )
        query = MultipleSeqAlignment(
            [SeqRecord(Seq("AT"), id="q1"), SeqRecord(Seq("GT"), id="q2")]
        )
        ref_cols, query_cols = service.get_columns_from_alignments(reference, query)
        assert ref_cols == ["AG", "CT"]
        assert query_cols == ["AG", "TT"]

    def test_calculate_matches_between_ref_and_query_columns(self, args):
        service = ColumnScore(args)
        matches, total = service.calculate_matches_between_ref_and_query_columns(
            ["AA", "CC", "GG"],
            ["AA", "TT", "GG"],
        )
        assert matches == 2
        assert total == 3

    def test_run_prints_score(self, mocker):
        args = Namespace(fasta="/some/path/to/query.fa", reference="/some/path/to/ref.fa", json=False)
        service = ColumnScore(args)
        query = MultipleSeqAlignment(
            [SeqRecord(Seq("AT"), id="q1"), SeqRecord(Seq("GT"), id="q2")]
        )
        reference = MultipleSeqAlignment(
            [SeqRecord(Seq("AC"), id="r1"), SeqRecord(Seq("GT"), id="r2")]
        )
        mocker.patch(
            "phykit.services.alignment.column_score.AlignIO.read",
            side_effect=[query, reference],
        )
        mocked_print = mocker.patch("builtins.print")
        service.run()
        mocked_print.assert_called_once_with(0.5)

    def test_run_json_output(self, mocker):
        args = Namespace(fasta="/some/path/to/query.fa", reference="/some/path/to/ref.fa", json=True)
        service = ColumnScore(args)
        query = MultipleSeqAlignment(
            [SeqRecord(Seq("AT"), id="q1"), SeqRecord(Seq("GT"), id="q2")]
        )
        reference = MultipleSeqAlignment(
            [SeqRecord(Seq("AC"), id="r1"), SeqRecord(Seq("GT"), id="r2")]
        )
        mocker.patch(
            "phykit.services.alignment.column_score.AlignIO.read",
            side_effect=[query, reference],
        )
        mocked_json = mocker.patch("phykit.services.alignment.column_score.print_json")
        service.run()
        mocked_json.assert_called_once_with({"column_score": 0.5})
