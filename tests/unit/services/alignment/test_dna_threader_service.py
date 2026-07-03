from argparse import Namespace
from io import StringIO
import subprocess
import sys

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from phykit.services.alignment.dna_threader import DNAThreader
import phykit.services.alignment.dna_threader as dna_threader_module


@pytest.fixture
def args():
    return Namespace(
        stop=True,
        protein="/some/path/to/protein.fa",
        nucleotide="/some/path/to/nucleotide.fa",
        clipkit_log_file=None,
    )


class TestDNAThreader:
    def test_module_import_does_not_import_biopython_seqio_or_seq(self):
        code = """
import sys
import phykit.services.alignment.dna_threader as module
assert hasattr(module.SeqIO, "parse")
assert hasattr(module.SeqIO, "to_dict")
assert "typing" not in sys.modules
assert "Bio.SeqIO" not in sys.modules
assert "Bio.Seq" not in sys.modules
assert "Bio.AlignIO" not in sys.modules
"""
        subprocess.run([sys.executable, "-c", code], check=True)

    def test_init_sets_expected_attrs(self, args):
        service = DNAThreader(args)
        assert service.remove_stop_codon is True
        assert service.protein_file_path == args.protein
        assert service.nucleotide_file_path == args.nucleotide
        assert service.clipkit_log_file is None
        assert service.json_output is False

    def test_clipkit_log_data_streams_rows_without_readlines(self, args, monkeypatch):
        args.clipkit_log_file = "clipkit.log"
        service = DNAThreader(args)

        class StreamingOnlyHandle(StringIO):
            def __enter__(self):
                return self

            def __exit__(self, *exc_info):
                self.close()

            def readlines(self, *args, **kwargs):
                raise AssertionError("clipkit log rows should be streamed")

        def fake_open(path):
            assert path == "clipkit.log"
            return StreamingOnlyHandle("1 keep\n2 trim extra\n")

        monkeypatch.setattr("builtins.open", fake_open)

        assert service.clipkit_log_data == [
            ["1", "keep"],
            ["2", "trim", "extra"],
        ]

    def test_create_mask_without_clipkit_log(self, args):
        service = DNAThreader(args)
        assert service.create_mask(6) == [True] * 6

    def test_create_mask_reads_clipkit_log_once(self, args, mocker):
        args.clipkit_log_file = "clipkit.log"
        service = DNAThreader(args)
        mocked_open = mocker.patch(
            "builtins.open",
            mocker.mock_open(read_data=b"1 keep\n2 trim\n"),
        )

        mask = service.create_mask(0)

        assert mask == [True, True, True, False, False, False]
        mocked_open.assert_called_once_with("clipkit.log", "rb")

    def test_create_mask_matches_exact_keep_status_token(self, args, mocker):
        args.clipkit_log_file = "clipkit.log"
        service = DNAThreader(args)
        mocker.patch(
            "builtins.open",
            mocker.mock_open(read_data=b"1 keep\n2 keep extra\n3 keepish\n4 keep"),
        )

        mask = service.create_mask(0)

        assert mask == [
            True, True, True,
            True, True, True,
            False, False, False,
            True, True, True,
        ]

    def test_create_mask_preserves_single_space_status_rule(self, args, mocker):
        args.clipkit_log_file = "clipkit.log"
        service = DNAThreader(args)
        mocker.patch(
            "builtins.open",
            mocker.mock_open(read_data=b"1  keep\n2 keep\n3 keep \n"),
        )

        mask = service.create_mask(0)

        assert mask == [
            False, False, False,
            True, True, True,
            True, True, True,
        ]

    def test_create_mask_returns_plain_boolean_list_from_codon_triplets(
        self, args, mocker
    ):
        args.clipkit_log_file = "clipkit.log"
        service = DNAThreader(args)
        mocker.patch(
            "builtins.open",
            mocker.mock_open(read_data=b"1 keep\n2 trim"),
        )

        mask = service.create_mask(0)

        assert type(mask) is list
        assert mask == [True, True, True, False, False, False]

    def test_normalize_p_seq(self, args):
        service = DNAThreader(args)
        assert service.normalize_p_seq(Seq("A-*")) == "AAA---***"

    def test_normalize_n_seq_handles_gaps(self, args):
        service = DNAThreader(args)
        assert service.normalize_n_seq(Seq("AAACCC"), Seq("A-?")) == "AAA------"

    def test_normalize_n_seq_handles_missing_codons(self, args):
        service = DNAThreader(args)
        assert service.normalize_n_seq(Seq("AAACCC"), Seq("ABC")) == "AAACCC---"

    def test_normalize_n_seq_preserves_partial_codons_and_gap_positions(self, args):
        service = DNAThreader(args)

        assert service.normalize_n_seq(Seq("AAAC"), Seq("A-B")) == "AAA---C"

    def test_thread_sequence_uses_string_mask_without_numpy_path(self, args):
        service = DNAThreader(args)
        keep_mask = [True, False, True, True, True, False, True, True, True]

        threaded = service._thread_sequence(Seq("AC*"), Seq("AAACCCGGGTTTAAA"), keep_mask)

        assert threaded == "AACCGGG"

    def test_thread_sequence_full_length_mask_avoids_mask_slices(self, args):
        class NoSliceMask(list):
            def __getitem__(self, item):
                if isinstance(item, slice):
                    raise AssertionError("mask should not be sliced")
                return super().__getitem__(item)

        service = DNAThreader(args)
        keep_mask = NoSliceMask([True, False, True, True, True, False, True, True, True])

        threaded = service._thread_sequence("AC*", "AAACCCGGGTTTAAA", keep_mask)

        assert threaded == "AACCGGG"

    def test_thread_sequence_terminal_stop_all_keep_uses_existing_prefix_behavior(self, args):
        service = DNAThreader(args)

        threaded = service._thread_sequence("AC*", "AAACCCGGGTTTAAA", [True] * 9)

        assert threaded == "AAACCCGGG"

    def test_thread_sequence_accepts_plain_strings(self, args):
        service = DNAThreader(args)

        assert service._thread_sequence("AC", "AAACCCGGGTTT", [True] * 6) == "AAACCC"

    def test_has_gap_char_detects_threading_gap_symbols(self):
        assert DNAThreader._has_gap_char("ACD") is False
        assert DNAThreader._has_gap_char("A-D") is True
        assert DNAThreader._has_gap_char("A?D") is True
        assert DNAThreader._has_gap_char("A*D") is True
        assert DNAThreader._has_gap_char("AXD") is True
        assert DNAThreader._has_gap_char("AxD") is True

    def test_thread_sequence_masked_without_terminal_stop(self, args):
        service = DNAThreader(args)

        threaded = service._thread_sequence(
            "ACD",
            "AAACCCGGGTTTAAA",
            [True, False, True, True, False, True, True, True, False],
        )

        assert threaded == "AACCGG"

    def test_thread_sequence_no_gap_mask_filters_nucleotide_prefix(self, args):
        service = DNAThreader(args)

        threaded = service._thread_sequence(
            "ACDEF",
            "AAACCCGGGTTTAAA",
            [True, False, True, True, True, False, False, True, True,
             True, False, True, True, True, True],
        )

        assert threaded == "AACCGGTTAAA"

    def test_thread_sequence_no_gap_mask_uses_compress_path(
        self, args, monkeypatch
    ):
        service = DNAThreader(args)
        calls = []

        def tracking_compress(data, selectors):
            calls.append((data, selectors))
            return (item for item, keep in zip(data, selectors) if keep)

        monkeypatch.setattr(
            "phykit.services.alignment.dna_threader.compress",
            tracking_compress,
        )

        threaded = service._thread_sequence(
            "ACDEF",
            "AAACCCGGGTTTAAA",
            [True, False, True, True, True, False, False, True, True,
             True, False, True, True, True, True],
        )

        assert threaded == "AACCGGTTAAA"
        assert calls == [
            (
                "AAACCCGGGTTTAAA",
                [True, False, True, True, True, False, False, True, True,
                 True, False, True, True, True, True],
            )
        ]

    def test_thread_sequence_no_stop_gap_path_matches_legacy_output(self, args):
        service = DNAThreader(args)

        threaded = service._thread_sequence(
            "A-CDEF",
            "AAACCCGGGTTTAAACCC",
            [True, False, True, True, True, False, True, True, True,
             False, True, True, True, False, True, True, True, False],
        )

        assert threaded == "AACCGGG------"

    def test_thread_sequence_uses_precomputed_mask_plan_for_gap_path(self, args):
        class NoIndexMask(list):
            def __getitem__(self, item):
                raise AssertionError("precomputed mask plan should provide keep bits")

        service = DNAThreader(args)
        protein = "A-CDEF"
        nucleotide = "AAACCCGGGTTTAAACCC"
        keep_mask = [
            True, False, True, True, True, False, True, True, True,
            False, True, True, True, False, True, True, True, False,
            True, False, True, True, True, False, True, True, True,
            True, False, True, True, True, False, True, True, True,
            True, False, True, True, True, False, True, True, True,
            True, False, True, True, True, False, True, True, True,
        ]
        expected = service._thread_sequence(protein, nucleotide, keep_mask)
        plan = service._create_thread_mask_plan(keep_mask, len(protein))

        observed = service._thread_sequence(
            protein,
            nucleotide,
            NoIndexMask(keep_mask),
            keep_mask_plan=plan,
        )

        assert observed == expected

    def test_thread_sequence_emit_plan_matches_direct_short_dna_fallback(self, args):
        class NoIndexMask(list):
            def __getitem__(self, item):
                raise AssertionError("precomputed emit plan should provide keep bits")

        service = DNAThreader(args)
        protein = "A-CDEF"
        nucleotide = "AAACCCGGGTTT"
        keep_mask = [
            True, False, True, True, True, False, True, True, True,
            False, True, True, True, False, True, True, True, False,
            True, False, True, True, True, False, True, True, True,
            True, False, True, True, True, False, True, True, True,
            True, False, True, True, True, False, True, True, True,
            True, False, True, True, True, False, True, True, True,
        ]
        expected = service._thread_sequence(protein, nucleotide, keep_mask)
        plan = service._create_thread_mask_plan(keep_mask, len(protein))

        observed = service._thread_sequence(
            protein,
            nucleotide,
            NoIndexMask(keep_mask),
            keep_mask_plan=plan,
        )

        assert len(plan) == 5
        assert observed == expected

    def test_mask_plan_uses_direct_full_segment_emitter(self, args, monkeypatch):
        def fail_itemgetter(*_args, **_kwargs):
            raise AssertionError("full keep groups should not build itemgetters")

        monkeypatch.setattr(dna_threader_module, "itemgetter", fail_itemgetter)

        plan = DNAThreader._create_thread_mask_plan([True] * 9, 3)

        assert plan[4][0][1] is dna_threader_module._FULL_SEGMENT

    def test_thread_sequence_full_trim_plan_group_consumes_nucleotides(self, args):
        service = DNAThreader(args)
        protein = "AC"
        nucleotide = "AAACCCGGGTTTAAACCC"
        keep_mask = [False] * 9 + [True] * 9
        plan = service._create_thread_mask_plan(keep_mask, len(protein))

        observed = service._thread_sequence(
            protein,
            nucleotide,
            keep_mask,
            keep_mask_plan=plan,
        )

        assert observed == "TTTAAACCC"

    def test_thread_sequence_emit_plan_matches_terminal_stop_path(self, args):
        class NoIndexMask(list):
            def __getitem__(self, item):
                raise AssertionError("precomputed emit plan should provide keep bits")

        service = DNAThreader(args)
        protein = "A-CDE*"
        nucleotide = "AAACCCGGGTTTAAA"
        keep_mask = [
            True, False, True, True, True, False, True, True, True,
            False, True, True, True, False, True, True, True, False,
            True, False, True, True, True, False, True, True, True,
            True, False, True, True, True, False, True, True, True,
            True, False, True, True, True, False, True, True, True,
            True, False, True, True, True, False, True, True, True,
        ]
        expected = service._thread_sequence(protein, nucleotide, keep_mask)
        plan = service._create_thread_mask_plan(keep_mask, len(protein))

        observed = service._thread_sequence(
            protein,
            nucleotide,
            NoIndexMask(keep_mask),
            keep_mask_plan=plan,
        )

        assert observed == expected

    def test_uniform_emit_plan_groups_runs_and_preserves_output(self, args):
        class NoIndexMask(list):
            def __getitem__(self, item):
                raise AssertionError("uniform emit plan should provide keep bits")

        service = DNAThreader(args)
        protein = "AA--CCXDD"
        nucleotide = "AAACCCGGGTTTAAACCCGGGTTT"
        keep_mask = [True, False, True, True, True, False, False, True, True] * len(protein)
        expected = service._thread_sequence(protein, nucleotide, keep_mask)
        plan = service._create_thread_mask_plan(keep_mask, len(protein))

        observed = service._thread_sequence(
            protein,
            nucleotide,
            NoIndexMask(keep_mask),
            keep_mask_plan=plan,
        )

        assert plan[4].uniform_emitter is not None
        assert plan[4].uniform_selectors == (
            True, False, True,
            True, True, False,
            False, True, True,
        )
        assert observed == expected

    def test_uniform_emit_plan_preserves_short_dna_padding(self, args):
        service = DNAThreader(args)
        protein = "AACD"
        nucleotide = "AAACCCGGGTT"
        keep_mask = [True, False, True, True, True, False, False, True, True] * len(protein)
        expected = service._thread_sequence(protein, nucleotide, keep_mask)
        plan = service._create_thread_mask_plan(keep_mask, len(protein))

        observed = service._thread_sequence(
            protein,
            nucleotide,
            keep_mask,
            keep_mask_plan=plan,
        )

        assert observed == expected

    def test_uniform_prefix_emit_plan_preserves_partial_tail_group(self, args):
        class NoIndexMask(list):
            def __getitem__(self, item):
                raise AssertionError("uniform prefix plan should provide keep bits")

        service = DNAThreader(args)
        protein = "AA--CCXDDD"
        nucleotide = "AAACCCGGGTTTAAACCCGGGTTTAAACCCGGG"
        keep_mask = [
            True, False, True,
            True, True, False,
            False, True, True,
        ] * len(protein)
        keep_mask = keep_mask[:len(protein) * 3]
        expected = service._thread_sequence(protein, nucleotide, keep_mask)
        plan = service._create_thread_mask_plan(keep_mask, len(protein))

        observed = service._thread_sequence(
            protein,
            nucleotide,
            NoIndexMask(keep_mask),
            keep_mask_plan=plan,
        )

        assert plan[4].uniform_selectors is None
        assert plan[4].uniform_prefix_selectors == (
            True, False, True,
            True, True, False,
            False, True, True,
        )
        assert plan[4].uniform_prefix_length == len(plan[4]) - 1
        assert observed == expected

    def test_thread_sequence_triplet_path_avoids_full_nucleotide_normalization(
        self, args, mocker
    ):
        service = DNAThreader(args)
        normalize_n_seq = mocker.spy(service, "normalize_n_seq")

        threaded = service._thread_sequence("ACDE", "AAACCCGGGTTT", [True] * 12)

        assert threaded == "AAACCCGGGTTT"
        normalize_n_seq.assert_not_called()

    def test_thread_sequence_all_keep_no_gaps_returns_nucleotide_prefix(self, args):
        service = DNAThreader(args)
        nucleotide = "AAACCCGGGTTT"

        threaded = service._thread_sequence("ACDE", nucleotide, [True] * 12)

        assert threaded == nucleotide

    def test_thread_all_sites_kept_groups_runs_and_preserves_padding(self, args):
        service = DNAThreader(args)

        threaded = service._thread_sequence(
            "AA--CCXDD",
            "AAACCCGGGTTTAAACCCGGGTTT",
            [True] * 81,
        )

        assert threaded == (
            "AAACCCGGG"
            "TTTAAACCC"
            "---------"
            "---------"
            "GGGTTT---"
            "---------"
            "---------"
            "---------"
            "---------"
        )

    def test_thread_sequence_all_keep_hint_skips_mask_scan(self, args):
        class NoScanMask(list):
            def __iter__(self):
                raise AssertionError("all-true mask hint should skip scan")

            def __getitem__(self, item):
                if isinstance(item, slice):
                    raise AssertionError("all-true mask hint should skip slicing")
                return super().__getitem__(item)

        service = DNAThreader(args)
        nucleotide = "AAACCCGGGTTT"
        keep_mask = NoScanMask([True] * 12)

        threaded = service._thread_sequence(
            "ACDE",
            nucleotide,
            keep_mask,
            keep_mask_all_true=True,
        )

        assert threaded == nucleotide

    def test_thread_exits_on_empty_proteins(self, mocker, args):
        service = DNAThreader(args)
        mocked_print = mocker.patch("builtins.print")
        with pytest.raises(SystemExit) as excinfo:
            service.thread(iter([]))
        assert excinfo.value.code == 2
        mocked_print.assert_called_once_with("Protein file is empty or incorrectly formatted.")

    def test_thread_basic(self, mocker, args):
        service = DNAThreader(args)
        prot_records = [SeqRecord(Seq("AC"), id="gene1")]
        nucl_records = [SeqRecord(Seq("AAACCC"), id="gene1")]
        mocker.patch("phykit.services.alignment.dna_threader.SeqIO.parse", return_value=iter(nucl_records))

        pal2nal = service.thread(iter(prot_records))
        assert pal2nal["gene1"] == "AAACCC"

    def test_thread_reuses_precomputed_mask_plan(self, mocker, args):
        args.clipkit_log_file = "clipkit.log"
        service = DNAThreader(args)
        prot_records = [SeqRecord(Seq("A-C"), id="gene1")]
        nucl_records = [SeqRecord(Seq("AAACCC"), id="gene1")]
        mocker.patch.object(
            service,
            "create_mask",
            return_value=[
                True, False, True,
                True, True, False,
                True, True, True,
            ],
        )
        mocker.patch(
            "phykit.services.alignment.dna_threader.SeqIO.parse",
            return_value=iter(nucl_records),
        )
        plans = []
        original = DNAThreader._thread_sequence

        def tracking_thread_sequence(self, *call_args, **kwargs):
            plans.append(kwargs.get("keep_mask_plan"))
            return original(self, *call_args, **kwargs)

        mocker.patch.object(
            DNAThreader,
            "_thread_sequence",
            tracking_thread_sequence,
        )

        service.thread(iter(prot_records))

        assert plans
        assert plans[0] is not None

    def test_thread_skips_precomputed_mask_plan_when_all_sites_kept(
        self, mocker, args
    ):
        service = DNAThreader(args)
        prot_records = [SeqRecord(Seq("A-C"), id="gene1")]
        nucl_records = [SeqRecord(Seq("AAACCC"), id="gene1")]
        mocker.patch(
            "phykit.services.alignment.dna_threader.SeqIO.parse",
            return_value=iter(nucl_records),
        )
        create_plan = mocker.patch.object(
            DNAThreader,
            "_create_thread_mask_plan",
            side_effect=AssertionError("all-true mask should not build a plan"),
        )

        pal2nal = service.thread(iter(prot_records))

        assert pal2nal["gene1"] == "AAACCC---"
        create_plan.assert_not_called()

    def test_thread_uses_mask_helper_all_true_flag_without_rescanning_mask(
        self, mocker, args
    ):
        class NoAllMask(list):
            def __iter__(self):
                raise AssertionError("thread should trust the helper all-true flag")

        args.clipkit_log_file = "clipkit.log"
        service = DNAThreader(args)
        prot_records = [SeqRecord(Seq("AC"), id="gene1")]
        nucl_records = [SeqRecord(Seq("AAACCC"), id="gene1")]
        mocker.patch(
            "phykit.services.alignment.dna_threader.SeqIO.parse",
            return_value=iter(nucl_records),
        )
        mocker.patch.object(
            service,
            "_create_mask_and_all_true",
            return_value=(NoAllMask([True] * 6), True),
        )
        create_plan = mocker.patch.object(
            DNAThreader,
            "_create_thread_mask_plan",
            side_effect=AssertionError("all-true mask should not build a plan"),
        )

        pal2nal = service.thread(iter(prot_records))

        assert pal2nal["gene1"] == "AAACCC"
        create_plan.assert_not_called()

    def test_run_text_output(self, mocker, args, capsys):
        service = DNAThreader(args)
        mocker.patch("phykit.services.alignment.dna_threader.SeqIO.parse", return_value=iter([]))
        mocker.patch.object(
            DNAThreader,
            "thread",
            return_value={"gene1": "AAACCC", "gene2": "GGGTTT"},
        )

        service.run()

        out, _ = capsys.readouterr()
        assert out == ">gene1\nAAACCC\n>gene2\nGGGTTT\n"

    def test_run_json_output(self, mocker):
        args = Namespace(
            stop=True,
            protein="/some/path/to/protein.fa",
            nucleotide="/some/path/to/nucleotide.fa",
            clipkit_log_file=None,
            json=True,
        )
        service = DNAThreader(args)
        mocker.patch("phykit.services.alignment.dna_threader.SeqIO.parse", return_value=iter([]))
        mocker.patch.object(DNAThreader, "thread", return_value={"gene1": "AAACCC"})
        mocked_json = mocker.patch("phykit.services.alignment.dna_threader.print_json")

        service.run()

        payload = mocked_json.call_args.args[0]
        assert payload["remove_stop_codon"] is True
        assert payload["rows"] == payload["taxa"]
        assert payload["rows"][0] == {"taxon": "gene1", "sequence": "AAACCC"}
