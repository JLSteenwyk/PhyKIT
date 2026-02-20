from argparse import Namespace

import numpy as np
import pytest
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from phykit.services.alignment.pairwise_identity import PairwiseIdentity


@pytest.fixture
def args():
    return Namespace(alignment="/some/path/to/file.fa", verbose=False, exclude_gaps=False)


def _simple_alignment():
    return MultipleSeqAlignment(
        [
            SeqRecord(Seq("ACGT"), id="a"),
            SeqRecord(Seq("AC-T"), id="b"),
            SeqRecord(Seq("TCGT"), id="c"),
        ]
    )


class TestPairwiseIdentity:
    def test_init_sets_expected_attrs(self, args):
        service = PairwiseIdentity(args)
        assert service.alignment_file_path == args.alignment
        assert service.verbose is False
        assert service.exclude_gaps is False
        assert service.json_output is False
        assert service.plot is False

    def test_calculate_identity_vectorized(self, args):
        service = PairwiseIdentity(args)
        seq1 = np.array(list("ACGT"), dtype="U1")
        seq2 = np.array(list("AC-T"), dtype="U1")
        identity = service._calculate_identity_vectorized(seq1, seq2)
        assert identity == 0.75

    def test_calculate_identity_vectorized_excluding_gaps(self, args):
        service = PairwiseIdentity(args)
        seq1 = np.array(list("ACGT"), dtype="U1")
        seq2 = np.array(list("AC-T"), dtype="U1")
        gap_mask1 = np.isin(seq1, ["-"])
        gap_mask2 = np.isin(seq2, ["-"])
        identity = service._calculate_identity_vectorized(
            seq1, seq2, (gap_mask1, gap_mask2), exclude_gaps=True
        )
        # denominator is still full length by implementation
        assert identity == 0.75

    def test_process_pair_batch(self, args):
        service = PairwiseIdentity(args)
        alignment_data = [
            {"id": "a", "seq": np.array(list("AAAA"), dtype="U1")},
            {"id": "b", "seq": np.array(list("AAAT"), dtype="U1")},
        ]
        results = service._process_pair_batch(
            alignment_data,
            pair_indices=[(0, 1)],
            exclude_gaps=False,
            gap_chars=set(["-"]),
        )
        assert results[0]["pair_id"] == ["a", "b"]
        assert round(results[0]["identity"], 4) == 0.75

    def test_calculate_pairwise_identities_sequential(self, mocker, args):
        service = PairwiseIdentity(args)
        mocker.patch.object(PairwiseIdentity, "_should_use_multiprocessing", return_value=False)
        pair_ids, identities, stats = service.calculate_pairwise_identities(
            _simple_alignment(), exclude_gaps=False, is_protein=False
        )
        assert len(pair_ids) == 3
        assert "a-b" in identities
        assert round(identities["a-b"], 4) == 0.75
        assert "mean" in stats

    def test_run_json_verbose_and_plot(self, mocker):
        args = Namespace(
            alignment="/some/path/to/file.fa",
            verbose=True,
            exclude_gaps=False,
            json=True,
            plot=True,
            plot_output="out.png",
        )
        service = PairwiseIdentity(args)
        alignment = _simple_alignment()
        mocker.patch.object(
            PairwiseIdentity,
            "get_alignment_and_format",
            return_value=(alignment, "fasta", False),
        )
        mocker.patch.object(
            PairwiseIdentity,
            "calculate_pairwise_identities",
            return_value=(
                [["a", "b"]],
                {"a-b": 0.75},
                {"mean": 0.75},
            ),
        )
        mocked_plot = mocker.patch.object(PairwiseIdentity, "_plot_pairwise_identity_heatmap")
        mocked_json = mocker.patch(
            "phykit.services.alignment.pairwise_identity.print_json"
        )

        service.run()

        mocked_plot.assert_called_once()
        payload = mocked_json.call_args.args[0]
        assert payload["verbose"] is True
        assert payload["pairs"] == payload["rows"]
        assert payload["plot_output"] == "out.png"
        assert payload["rows"][0] == {"taxon_a": "a", "taxon_b": "b", "identity": 0.75}

    def test_run_summary_prints_and_plot_message(self, mocker):
        args = Namespace(
            alignment="/some/path/to/file.fa",
            verbose=False,
            exclude_gaps=False,
            json=False,
            plot=True,
            plot_output="out.png",
        )
        service = PairwiseIdentity(args)
        alignment = _simple_alignment()
        mocker.patch.object(
            PairwiseIdentity,
            "get_alignment_and_format",
            return_value=(alignment, "fasta", False),
        )
        mocker.patch.object(
            PairwiseIdentity,
            "calculate_pairwise_identities",
            return_value=(
                [["a", "b"]],
                {"a-b": 0.75},
                {"mean": 0.75},
            ),
        )
        mocker.patch.object(PairwiseIdentity, "_plot_pairwise_identity_heatmap")
        mocked_summary = mocker.patch(
            "phykit.services.alignment.pairwise_identity.print_summary_statistics"
        )
        mocked_print = mocker.patch("builtins.print")

        service.run()
        mocked_summary.assert_called_once_with({"mean": 0.75})
        mocked_print.assert_called_once_with("Saved pairwise identity heatmap: out.png")
