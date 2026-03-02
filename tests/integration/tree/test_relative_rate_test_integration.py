import json
import sys

import pytest
from mock import patch

from phykit.phykit import Phykit


@pytest.mark.integration
class TestRelativeRateTestIntegration:
    @patch("builtins.print")
    def test_single_alignment(self, mocked_print, tmp_path):
        aln = tmp_path / "test.fa"
        aln.write_text(
            ">A\nACGTACGTACGTACGTACGT\n"
            ">B\nACGTACGTACTAACTAACGT\n"
            ">C\nACTAACGTACGTACGTACGT\n"
            ">O\nACGTACGTACGTACGTACGT\n"
        )
        tree = tmp_path / "test.tre"
        tree.write_text("(((A:0.1,B:0.1):0.05,C:0.15):0.1,O:0.2);")

        testargs = ["phykit", "relative_rate_test", "-a", str(aln), "-t", str(tree)]
        with patch.object(sys, "argv", testargs):
            Phykit()

        output = "\n".join(str(call) for call in mocked_print.call_args_list)
        assert "Outgroup: O" in output
        assert "Number of pairwise tests: 3" in output

    @patch("builtins.print")
    def test_rrt_alias(self, mocked_print, tmp_path):
        aln = tmp_path / "test.fa"
        aln.write_text(
            ">A\nACGTACGT\n>B\nACGTACGT\n>O\nACGTACGT\n"
        )
        tree = tmp_path / "test.tre"
        tree.write_text("((A:0.1,B:0.1):0.1,O:0.2);")

        testargs = ["phykit", "rrt", "-a", str(aln), "-t", str(tree)]
        with patch.object(sys, "argv", testargs):
            Phykit()

        output = "\n".join(str(call) for call in mocked_print.call_args_list)
        assert "Outgroup: O" in output

    @patch("builtins.print")
    def test_tajima_rrt_alias(self, mocked_print, tmp_path):
        aln = tmp_path / "test.fa"
        aln.write_text(
            ">A\nACGTACGT\n>B\nACGTACGT\n>O\nACGTACGT\n"
        )
        tree = tmp_path / "test.tre"
        tree.write_text("((A:0.1,B:0.1):0.1,O:0.2);")

        testargs = ["phykit", "tajima_rrt", "-a", str(aln), "-t", str(tree)]
        with patch.object(sys, "argv", testargs):
            Phykit()

        output = "\n".join(str(call) for call in mocked_print.call_args_list)
        assert "Outgroup: O" in output

    @patch("builtins.print")
    def test_json_output(self, mocked_print, tmp_path):
        aln = tmp_path / "test.fa"
        aln.write_text(
            ">A\nACGTACGT\n>B\nACTTACGT\n>O\nACGTACGT\n"
        )
        tree = tmp_path / "test.tre"
        tree.write_text("((A:0.1,B:0.1):0.1,O:0.2);")

        testargs = ["phykit", "rrt", "-a", str(aln), "-t", str(tree), "--json"]
        with patch.object(sys, "argv", testargs):
            Phykit()

        # print_json calls print once with the JSON string
        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["outgroup"] == "O"
        assert payload["n_tests"] == 1
        assert "chi2" in payload["results"][0]

    @patch("builtins.print")
    def test_batch_mode(self, mocked_print, tmp_path):
        aln_content = ">A\nACGTACGT\n>B\nACTTACGT\n>O\nACGTACGT\n"
        aln1 = tmp_path / "aln1.fa"
        aln1.write_text(aln_content)
        aln2 = tmp_path / "aln2.fa"
        aln2.write_text(aln_content)
        aln_list = tmp_path / "list.txt"
        aln_list.write_text(f"{aln1}\n{aln2}\n")
        tree = tmp_path / "test.tre"
        tree.write_text("((A:0.1,B:0.1):0.1,O:0.2);")

        testargs = ["phykit", "rrt", "-l", str(aln_list), "-t", str(tree)]
        with patch.object(sys, "argv", testargs):
            Phykit()

        output = "\n".join(str(call) for call in mocked_print.call_args_list)
        assert "Number of alignments: 2" in output
