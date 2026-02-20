from argparse import Namespace

import numpy as np
import pytest
import builtins
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from phykit.services.alignment.pairwise_identity import PairwiseIdentity
import phykit.services.alignment.pairwise_identity as pairwise_identity_module


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

    def test_should_use_multiprocessing_respects_env(self, args, monkeypatch):
        service = PairwiseIdentity(args)
        monkeypatch.setenv("PHYKIT_DISABLE_MP", "1")
        assert service._should_use_multiprocessing(99999) is False
        monkeypatch.delenv("PHYKIT_DISABLE_MP")
        monkeypatch.setenv("PHYKIT_FORCE_MP", "1")
        assert service._should_use_multiprocessing(1) is True
        monkeypatch.delenv("PHYKIT_FORCE_MP")
        assert service._should_use_multiprocessing(service.MP_MIN_PAIRS - 1) is False
        assert service._should_use_multiprocessing(service.MP_MIN_PAIRS) is True

    def test_process_args_defaults(self):
        parsed = PairwiseIdentity(
            Namespace(alignment="x.fa", verbose=False, exclude_gaps=False)
        ).process_args(
            Namespace(alignment="x.fa", verbose=False, exclude_gaps=False)
        )
        assert parsed["json_output"] is False
        assert parsed["plot"] is False
        assert parsed["plot_output"] == "pairwise_identity_heatmap.png"

    def test_print_json_output_non_verbose(self, mocker):
        service = PairwiseIdentity(
            Namespace(alignment="x.fa", verbose=False, exclude_gaps=True, json=True, plot=False)
        )
        mocked_json = mocker.patch("phykit.services.alignment.pairwise_identity.print_json")
        service._print_json_output(pair_ids=[["a", "b"]], pairwise_identities={"a-b": 0.5}, stats={"mean": 0.5})
        payload = mocked_json.call_args.args[0]
        assert payload["verbose"] is False
        assert payload["exclude_gaps"] is True
        assert payload["summary"] == {"mean": 0.5}

    def test_plot_pairwise_identity_heatmap_creates_file(self, tmp_path):
        pytest.importorskip("matplotlib")
        output = tmp_path / "pairwise.png"
        service = PairwiseIdentity(
            Namespace(alignment="x.fa", verbose=False, exclude_gaps=False, plot=True, plot_output=str(output))
        )
        service._plot_pairwise_identity_heatmap(
            taxa=["a", "b", "c"],
            pair_ids=[["a", "b"], ["a", "c"], ["b", "c"]],
            pairwise_identities={"a-b": 0.8, "a-c": 0.5, "b-c": 0.6},
        )
        assert output.exists()

    def test_plot_pairwise_identity_heatmap_importerror(self, monkeypatch, capsys):
        original_import = builtins.__import__

        def fake_import(name, globals=None, locals=None, fromlist=(), level=0):
            if name.startswith("matplotlib"):
                raise ImportError("no matplotlib")
            return original_import(name, globals, locals, fromlist, level)

        monkeypatch.setattr(builtins, "__import__", fake_import)
        service = PairwiseIdentity(
            Namespace(alignment="x.fa", verbose=False, exclude_gaps=False)
        )
        with pytest.raises(SystemExit) as exc:
            service._plot_pairwise_identity_heatmap(["a"], [], {})
        assert exc.value.code == 2
        out, _ = capsys.readouterr()
        assert "matplotlib is required for --plot in pairwise_identity" in out

    def test_plot_pairwise_identity_heatmap_empty_taxa_returns(self):
        pytest.importorskip("matplotlib")
        service = PairwiseIdentity(
            Namespace(alignment="x.fa", verbose=False, exclude_gaps=False)
        )
        service._plot_pairwise_identity_heatmap([], [], {})

    def test_run_verbose_handles_broken_pipe(self, mocker):
        service = PairwiseIdentity(
            Namespace(alignment="x.fa", verbose=True, exclude_gaps=False, json=False, plot=False)
        )
        alignment = _simple_alignment()
        mocker.patch.object(PairwiseIdentity, "get_alignment_and_format", return_value=(alignment, "fasta", False))
        mocker.patch.object(
            PairwiseIdentity,
            "calculate_pairwise_identities",
            return_value=([["a", "b"]], {"a-b": 0.75}, {"mean": 0.75}),
        )
        mocker.patch("builtins.print", side_effect=BrokenPipeError)
        service.run()

    def test_calculate_pairwise_identities_multiprocessing_branch(self, mocker, args, monkeypatch):
        class FakePool:
            def __init__(self, processes):
                self.processes = processes

            def __enter__(self):
                return self

            def __exit__(self, exc_type, exc, tb):
                return False

            def map(self, process_func, chunks):
                return [process_func(chunk) for chunk in chunks]

            def imap(self, process_func, chunks):
                for chunk in chunks:
                    yield process_func(chunk)

        service = PairwiseIdentity(args)
        mocker.patch.object(PairwiseIdentity, "_should_use_multiprocessing", return_value=True)
        mocker.patch("phykit.services.alignment.pairwise_identity.mp.Pool", FakePool)
        mocker.patch("phykit.services.alignment.pairwise_identity.mp.cpu_count", return_value=2)
        monkeypatch.setattr(pairwise_identity_module.sys.stderr, "isatty", lambda: False)

        pair_ids, identities, stats = service.calculate_pairwise_identities(
            _simple_alignment(),
            exclude_gaps=False,
            is_protein=False,
        )
        assert len(pair_ids) == 3
        assert "a-b" in identities
        assert "mean" in stats
