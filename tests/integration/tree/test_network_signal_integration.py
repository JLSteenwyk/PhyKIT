import json
import shutil
import sys

import pytest
from mock import patch

from phykit.phykit import Phykit


FIVE_TAXON_TREE = "(((A:1,B:1):0.5,C:1.5):0.5,(D:1,E:1):1);"
FIVE_TAXON_TRAITS = "A\t2.5\nB\t2.8\nC\t1.2\nD\t5.1\nE\t4.9\n"

THREE_TAXON_TREE = "((A:1,B:1):0.5,C:1.5);"
THREE_TAXON_TRAITS = "A\t2.5\nB\t2.8\nC\t1.2\n"

QUARTET_JSON_SOURCE = "tests/sample_files/network_signal_quartet.json"


def _output(mocked_print):
    """Collect all print() calls into one string."""
    return "\n".join(str(call) for call in mocked_print.call_args_list)


@pytest.mark.integration
class TestNetworkSignalIntegration:
    @patch("builtins.print")
    def test_explicit_hybrid(self, mocked_print, tmp_path):
        tree = tmp_path / "tree.tre"
        tree.write_text(FIVE_TAXON_TREE)
        traits = tmp_path / "traits.tsv"
        traits.write_text(FIVE_TAXON_TRAITS)

        testargs = [
            "phykit", "network_signal",
            "-t", str(tree), "-d", str(traits),
            "--hybrid", "B:D:0.3",
            "--permutations", "100",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        output = _output(mocked_print)
        assert "Hybrid edge:" in output
        assert "-> D" in output
        assert "Blomberg's K:" in output
        assert "Pagel's lambda:" in output

    @patch("builtins.print")
    def test_netsig_alias(self, mocked_print, tmp_path):
        tree = tmp_path / "tree.tre"
        tree.write_text(THREE_TAXON_TREE)
        traits = tmp_path / "traits.tsv"
        traits.write_text(THREE_TAXON_TRAITS)

        testargs = [
            "phykit", "netsig",
            "-t", str(tree), "-d", str(traits),
            "--hybrid", "A:C:0.2",
            "--permutations", "100",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        output = _output(mocked_print)
        assert "Blomberg's K:" in output

    @patch("builtins.print")
    def test_net_signal_alias(self, mocked_print, tmp_path):
        tree = tmp_path / "tree.tre"
        tree.write_text(THREE_TAXON_TREE)
        traits = tmp_path / "traits.tsv"
        traits.write_text(THREE_TAXON_TRAITS)

        testargs = [
            "phykit", "net_signal",
            "-t", str(tree), "-d", str(traits),
            "--hybrid", "A:C:0.2",
            "--permutations", "100",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        output = _output(mocked_print)
        assert "Blomberg's K:" in output

    @patch("builtins.print")
    def test_json_output(self, mocked_print, tmp_path):
        tree = tmp_path / "tree.tre"
        tree.write_text(THREE_TAXON_TREE)
        traits = tmp_path / "traits.tsv"
        traits.write_text(THREE_TAXON_TRAITS)

        testargs = [
            "phykit", "network_signal",
            "-t", str(tree), "-d", str(traits),
            "--hybrid", "A:C:0.2",
            "--permutations", "100",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        payload = json.loads(mocked_print.call_args.args[0])
        assert "blombergs_k" in payload
        assert "pagels_lambda" in payload

    @patch("builtins.print")
    def test_quartet_json_mode(self, mocked_print, tmp_path):
        tree = tmp_path / "tree.tre"
        tree.write_text(FIVE_TAXON_TREE)
        traits = tmp_path / "traits.tsv"
        traits.write_text(FIVE_TAXON_TRAITS)

        quartet_file = tmp_path / "quartet.json"
        shutil.copy(QUARTET_JSON_SOURCE, str(quartet_file))

        testargs = [
            "phykit", "network_signal",
            "-t", str(tree), "-d", str(traits),
            "--quartet-json", str(quartet_file),
            "--permutations", "100",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        output = _output(mocked_print)
        assert "Hybrid edge:" in output
        assert "Blomberg's K:" in output

    @patch("builtins.print")
    def test_method_k_only(self, mocked_print, tmp_path):
        tree = tmp_path / "tree.tre"
        tree.write_text(THREE_TAXON_TREE)
        traits = tmp_path / "traits.tsv"
        traits.write_text(THREE_TAXON_TRAITS)

        testargs = [
            "phykit", "network_signal",
            "-t", str(tree), "-d", str(traits),
            "--hybrid", "A:C:0.2",
            "--method", "blombergs_k",
            "--permutations", "100",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        output = _output(mocked_print)
        assert "Blomberg's K:" in output
        assert "Pagel's lambda:" not in output

    @patch("builtins.print")
    def test_method_lambda_only(self, mocked_print, tmp_path):
        tree = tmp_path / "tree.tre"
        tree.write_text(THREE_TAXON_TREE)
        traits = tmp_path / "traits.tsv"
        traits.write_text(THREE_TAXON_TRAITS)

        testargs = [
            "phykit", "network_signal",
            "-t", str(tree), "-d", str(traits),
            "--hybrid", "A:C:0.2",
            "--method", "lambda",
            "--permutations", "100",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        output = _output(mocked_print)
        assert "Pagel's lambda:" in output
        assert "Blomberg's K:" not in output
