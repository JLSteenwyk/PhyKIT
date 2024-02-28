import pytest
import textwrap

from mock import patch, call
from pathlib import Path
import sys

from phykit.phykit import Phykit

here = Path(__file__)


@pytest.mark.integration
class TestAlignmentRecoding(object):
    @patch("builtins.print")
    def test_alignment_recoding_invalid_input(self, mocked_print):  # noqa
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Phykit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_alignment_recoding_ry_nucleotide(self, mocked_print):
        expected_result_0 = textwrap.dedent(
            """>1\nR-RYRY"""
        )
        expected_result_1 = textwrap.dedent(
            """>2\nR-R-RY"""
        )
        expected_result_2 = textwrap.dedent(
            """>3\nR-R-YR"""
        )
        expected_result_3 = textwrap.dedent(
            """>4\nRRR-YR"""
        )
        expected_result_4 = textwrap.dedent(
            """>5\nRYR-Y-"""
        )

        testargs = [
            "phykit",
            "alignment_recoding",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
            "-c",
            "RY-nucleotide",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [
            call(expected_result_0),
            call(expected_result_1),
            call(expected_result_2),
            call(expected_result_3),
            call(expected_result_4),
        ]

    @patch("builtins.print")
    def test_alignment_recoding_dayhoff_6(self, mocked_print):
        expected_result_0 = textwrap.dedent(
            """>200_S38\n301330131053113001440030430042110000000--130100003023021000003000301230000030010-000100---0---------------11000---00-02031011013000010000400212130213332111314330330030001011102311310132131212002321322320231133003030320104011---"""
        )
        expected_result_1 = textwrap.dedent(
            """>203_S40\n301330131053113001440030430042110000000--130100003023021000003000301230000030010-000100---0---------------11000---00-02031011013000010000400212130213332111314330330030001011102311310132131212002321322320231133003030320104011---"""
        )

        testargs = [
            "phykit",
            "alignment_recoding",
            f"{here.parent.parent.parent}/sample_files/EOG091N44MS.fa.mafft",
            "-c",
            "Dayhoff-6",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [
            call(expected_result_0),
            call(expected_result_1),
        ]

    @patch("builtins.print")
    def test_alignment_recoding_dayhoff_9(self, mocked_print):
        expected_result_0 = textwrap.dedent(
            """>200_S38\n130113010371001330223315213320003633666--016033631341340333661633130415533313503-366036---6---------------00653---36-54310500601663603633233404013401114000102113116515330300034100103014010404334140144143410011531315105052300---"""
        )
        expected_result_1 = textwrap.dedent(
            """>203_S40\n130113010371001330223315213320003633666--016033631341340333661633130415533313503-366036---6---------------00653---36-54310500601663603633233404013401114000102113116515330300034100103014010404334140144143410011531315105052300---"""
        )

        testargs = [
            "phykit",
            "alignment_recoding",
            f"{here.parent.parent.parent}/sample_files/EOG091N44MS.fa.mafft",
            "-c",
            "Dayhoff-9",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [
            call(expected_result_0),
            call(expected_result_1),
        ]

    @patch("builtins.print")
    def test_alignment_recoding_dayhoff_12(self, mocked_print):
        expected_result_0 = textwrap.dedent(
            """>200_S38\n1501180108B1001580225814218823095688666--916055651531839777661655159314455515495-766056---6---------------00645---56-43510400691665607678257303015301113000102117116414570500083100108013010303553130133138310011451514134042709---"""
        )
        expected_result_1 = textwrap.dedent(
            """>203_S40\n1501180108B1001580225814218823095688666--916055651531839777661655159314455515495-766056---6---------------00645---56-43510400691665607678257303015301113000102117116414570500083100108013010303553130133138310011451514134042709---"""
        )

        testargs = [
            "phykit",
            "alignment_recoding",
            f"{here.parent.parent.parent}/sample_files/EOG091N44MS.fa.mafft",
            "-c",
            "Dayhoff-12",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [
            call(expected_result_0),
            call(expected_result_1),
        ]

    @patch("builtins.print")
    def test_alignment_recoding_dayhoff_15(self, mocked_print):
        expected_result_0 = textwrap.dedent(
            """>200_S38\n1502180108E100158033581431883B095688666--9260556515A28A9777662655259A24455525495-766056---6---------------00645---56-4C520400691665607678357C0C015C0122A00020312721642457050008C20010801C020A0C55A2C01AA1C8C100214525242B4043709---"""
        )
        expected_result_1 = textwrap.dedent(
            """>203_S40\n1502180108E100158033581431883B095688666--9260556515A28A9777662655259A24455525495-766056---6---------------00645---56-4C520400691665607678357C0C015C0122A00020312721642457050008C20010801C020A0C55A2C01AA1C8C100214525242B4043709---"""
        )

        testargs = [
            "phykit",
            "alignment_recoding",
            f"{here.parent.parent.parent}/sample_files/EOG091N44MS.fa.mafft",
            "-c",
            "Dayhoff-15",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [
            call(expected_result_0),
            call(expected_result_1),
        ]

    @patch("builtins.print")
    def test_alignment_recoding_dayhoff_18(self, mocked_print):
        expected_result_0 = textwrap.dedent(
            """>200_S38\n159218B1B8H19B158B00581401880E9C5688666--C369556515D28DC77766365525CD344555354C5-766B56---6---------------BB645---56-4F53A4A96C1665697678057FBFA15F9122DABB2A012731642457A5ABA8F2B91A8A1F93ADAF55D3FA1DD1F8F1A9314535342E49407BC---"""
        )
        expected_result_1 = textwrap.dedent(
            """>203_S40\n159218B1B8H19B158B00581401880E9C5688666--C369556515D28DC77766365525CD344555354C5-766B56---6---------------BB645---56-4F53A4A96C1665697678057FBFA15F9122DABB2A012731642457A5ABA8F2B91A8A1F93ADAF55D3FA1DD1F8F1A9314535342E49407BC---"""
        )

        testargs = [
            "phykit",
            "alignment_recoding",
            f"{here.parent.parent.parent}/sample_files/EOG091N44MS.fa.mafft",
            "-c",
            "Dayhoff-18",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [
            call(expected_result_0),
            call(expected_result_1),
        ]

    @patch("builtins.print")
    def test_alignment_recoding_sandr_6(self, mocked_print):
        expected_result_0 = textwrap.dedent(
            """>200_S38\n301330232043123002550031530055110000000--130100003023021000003000301231100030110-000200---0---------------22010---00-12031111013000010000500222130213332122315330330131001012102321310132131212002321322320231133103031351115021---"""
        )
        expected_result_1 = textwrap.dedent(
            """>203_S40\n301330232043123002550031530055110000000--130100003023021000003000301231100030110-000200---0---------------22010---00-12031111013000010000500222130213332122315330330131001012102321310132131212002321322320231133103031351115021---"""
        )

        testargs = [
            "phykit",
            "alignment_recoding",
            f"{here.parent.parent.parent}/sample_files/EOG091N44MS.fa.mafft",
            "-c",
            "SandR-6",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [
            call(expected_result_0),
            call(expected_result_1),
        ]

    @patch("builtins.print")
    def test_alignment_recoding_KGB_6(self, mocked_print):
        expected_result_0 = textwrap.dedent(
            """>200_S38\n201221121152112011440120421141110011000--150100002012111000005000201150000050010-000100---0---------------11000---00-01051011012000010001400111120112221111214220520020001011111211211121151111001511211211121152005050210104011---"""
        )
        expected_result_1 = textwrap.dedent(
            """>203_S40\n201221121152112011440120421141110011000--150100002012111000005000201150000050010-000100---0---------------11000---00-01051011012000010001400111120112221111214220520020001011111211211121151111001511211211121152005050210104011---"""
        )

        testargs = [
            "phykit",
            "alignment_recoding",
            f"{here.parent.parent.parent}/sample_files/EOG091N44MS.fa.mafft",
            "-c",
            "KGB-6",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [
            call(expected_result_0),
            call(expected_result_1),
        ]

    @patch("builtins.print")
    def test_alignment_recoding_alias0(self, mocked_print):
        expected_result_0 = textwrap.dedent(
            """>1\nR-RYRY"""
        )
        expected_result_1 = textwrap.dedent(
            """>2\nR-R-RY"""
        )
        expected_result_2 = textwrap.dedent(
            """>3\nR-R-YR"""
        )
        expected_result_3 = textwrap.dedent(
            """>4\nRRR-YR"""
        )
        expected_result_4 = textwrap.dedent(
            """>5\nRYR-Y-"""
        )

        testargs = [
            "phykit",
            "aln_recoding",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
            "-c",
            "RY-nucleotide",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [
            call(expected_result_0),
            call(expected_result_1),
            call(expected_result_2),
            call(expected_result_3),
            call(expected_result_4),
        ]

    @patch("builtins.print")
    def test_alignment_recoding_alias1(self, mocked_print):
        expected_result_0 = textwrap.dedent(
            """>1\nR-RYRY"""
        )
        expected_result_1 = textwrap.dedent(
            """>2\nR-R-RY"""
        )
        expected_result_2 = textwrap.dedent(
            """>3\nR-R-YR"""
        )
        expected_result_3 = textwrap.dedent(
            """>4\nRRR-YR"""
        )
        expected_result_4 = textwrap.dedent(
            """>5\nRYR-Y-"""
        )

        testargs = [
            "phykit",
            "recode",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
            "-c",
            "RY-nucleotide",

        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [
            call(expected_result_0),
            call(expected_result_1),
            call(expected_result_2),
            call(expected_result_3),
            call(expected_result_4),
        ]

    @patch("builtins.print")
    def test_alignment_recoding_no_recoding_table(self, mocked_print):
        expected_call = textwrap.dedent(
            """Please specify a recoding table"""
        )

        testargs = [
            "phykit",
            "recode",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
        ]

        with patch.object(sys, "argv", testargs):
            with pytest.raises(SystemExit) as pytest_wrapped_e:
                Phykit()

        assert pytest_wrapped_e.type == SystemExit
        mocked_print.assert_has_calls([
            call(expected_call),
        ])

    @patch("builtins.print")
    def test_alignment_recoding_custom_code(self, mocked_print):
        expected_result_0 = textwrap.dedent(
            """>1\nR-RYRY"""
        )
        expected_result_1 = textwrap.dedent(
            """>2\nR-R-RY"""
        )
        expected_result_2 = textwrap.dedent(
            """>3\nR-R-YR"""
        )
        expected_result_3 = textwrap.dedent(
            """>4\nRRR-YR"""
        )
        expected_result_4 = textwrap.dedent(
            """>5\nRYR-Y-"""
        )

        testargs = [
            "phykit",
            "recode",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
            "-c",
            f"{here.parent.parent.parent.parent}/phykit/recoding_tables/RY-nucleotide.txt",

        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [
            call(expected_result_0),
            call(expected_result_1),
            call(expected_result_2),
            call(expected_result_3),
            call(expected_result_4),
        ]
