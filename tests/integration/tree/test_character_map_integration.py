"""
Integration tests for character_map (discrete character mapping via Fitch parsimony).

Test data:
  tree_character_map.tre  — ((A:1,B:1):1,((C:1,D:1):1,E:1):1,F:1);
  character_matrix_simple.tsv — 6 taxa x 8 characters
"""
from mock import patch
from pathlib import Path
import json
import sys

import pytest

from phykit.phykit import Phykit

here = Path(__file__)


@pytest.mark.integration
class TestCharacterMap:
    @patch("builtins.print")
    def test_taxon_mismatch_fails_with_names_from_both_inputs(
        self,
        mocked_print,
        tmp_path,
    ):
        matrix = tmp_path / "mismatched.tsv"
        matrix.write_text(
            "taxon\tchar0\n"
            "A\t0\n"
            "B\t0\n"
            "C\t1\n"
            "off_tree\t1\n"
        )
        output = str(tmp_path / "mismatch.png")
        testargs = [
            "phykit",
            "character_map",
            "-t", f"{here.parent.parent.parent}/sample_files/tree_character_map.tre",
            "-d", str(matrix),
            "-o", output,
        ]

        with patch.object(sys, "argv", testargs), pytest.raises(SystemExit) as exc:
            Phykit()

        assert exc.value.code == 2
        messages = [call.args[0] for call in mocked_print.call_args_list]
        assert "3 taxa in tree but not in character matrix: D, E, F" in messages
        assert "1 taxon in character matrix but not in tree: off_tree" in messages
        assert not Path(output).exists()

    @patch("builtins.print")
    def test_allow_taxon_mismatch_warns_and_uses_shared_taxa(
        self,
        mocked_print,
        tmp_path,
    ):
        matrix = tmp_path / "mismatched.tsv"
        matrix.write_text(
            "taxon\tchar0\n"
            "A\t0\n"
            "B\t0\n"
            "C\t1\n"
            "off_tree\t1\n"
        )
        output = str(tmp_path / "shared.png")
        testargs = [
            "phykit",
            "character_map",
            "-t", f"{here.parent.parent.parent}/sample_files/tree_character_map.tre",
            "-d", str(matrix),
            "-o", output,
            "--allow-taxon-mismatch",
            "--json",
        ]

        with patch.object(sys, "argv", testargs):
            Phykit()

        messages = [call.args[0] for call in mocked_print.call_args_list]
        assert "Warning: 3 taxa in tree but not in character matrix: D, E, F" in messages
        assert "Warning: 1 taxon in character matrix but not in tree: off_tree" in messages
        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["n_characters"] == 1
        assert Path(output).exists()

    @patch("builtins.print")
    def test_basic_invocation(self, mocked_print, tmp_path):
        output = str(tmp_path / "charmap.png")
        testargs = [
            "phykit",
            "character_map",
            "-t", f"{here.parent.parent.parent}/sample_files/tree_character_map.tre",
            "-d", f"{here.parent.parent.parent}/sample_files/character_matrix_simple.tsv",
            "-o", output,
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert Path(output).exists()

    @patch("builtins.print")
    def test_alias_charmap(self, mocked_print, tmp_path):
        output = str(tmp_path / "charmap_alias.png")
        testargs = [
            "phykit",
            "charmap",
            "-t", f"{here.parent.parent.parent}/sample_files/tree_character_map.tre",
            "-d", f"{here.parent.parent.parent}/sample_files/character_matrix_simple.tsv",
            "-o", output,
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert Path(output).exists()

    @patch("builtins.print")
    def test_alias_synapomorphy_map(self, mocked_print, tmp_path):
        output = str(tmp_path / "charmap_synap.png")
        testargs = [
            "phykit",
            "synapomorphy_map",
            "-t", f"{here.parent.parent.parent}/sample_files/tree_character_map.tre",
            "-d", f"{here.parent.parent.parent}/sample_files/character_matrix_simple.tsv",
            "-o", output,
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert Path(output).exists()

    @patch("builtins.print")
    def test_json_ground_truth(self, mocked_print, tmp_path):
        output = str(tmp_path / "charmap.png")
        testargs = [
            "phykit",
            "character_map",
            "-t", f"{here.parent.parent.parent}/sample_files/tree_character_map.tre",
            "-d", f"{here.parent.parent.parent}/sample_files/character_matrix_simple.tsv",
            "-o", output,
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["n_characters"] == 8
        assert payload["optimization"] == "acctran"
        assert "ci" in payload
        assert "ri" in payload

    @patch("builtins.print")
    def test_deltran_option(self, mocked_print, tmp_path):
        output = str(tmp_path / "charmap_deltran.png")
        testargs = [
            "phykit",
            "character_map",
            "-t", f"{here.parent.parent.parent}/sample_files/tree_character_map.tre",
            "-d", f"{here.parent.parent.parent}/sample_files/character_matrix_simple.tsv",
            "-o", output,
            "--optimization", "deltran",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["optimization"] == "deltran"
        assert Path(output).exists()

    @patch("builtins.print")
    def test_ladderize_flag(self, mocked_print, tmp_path):
        output = str(tmp_path / "charmap_ladder.png")
        testargs = [
            "phykit",
            "character_map",
            "-t", f"{here.parent.parent.parent}/sample_files/tree_character_map.tre",
            "-d", f"{here.parent.parent.parent}/sample_files/character_matrix_simple.tsv",
            "-o", output,
            "--ladderize",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert Path(output).exists()

    @patch("builtins.print")
    def test_character_filter(self, mocked_print, tmp_path):
        output = str(tmp_path / "charmap_filter.png")
        testargs = [
            "phykit",
            "character_map",
            "-t", f"{here.parent.parent.parent}/sample_files/tree_character_map.tre",
            "-d", f"{here.parent.parent.parent}/sample_files/character_matrix_simple.tsv",
            "-o", output,
            "--characters", "0,1",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["n_characters"] == 8
        assert Path(output).exists()

    @patch("builtins.print")
    def test_verbose_output(self, mocked_print, tmp_path):
        output = str(tmp_path / "charmap_verbose.png")
        testargs = [
            "phykit",
            "character_map",
            "-t", f"{here.parent.parent.parent}/sample_files/tree_character_map.tre",
            "-d", f"{here.parent.parent.parent}/sample_files/character_matrix_simple.tsv",
            "-o", output,
            "--verbose",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert mocked_print.call_count == 1
        all_output = mocked_print.call_args.args[0]
        assert "Optimization: acctran" in all_output
        assert "Character" in all_output
        assert "CI:" in all_output

    @patch("builtins.print")
    def test_pdf_output(self, mocked_print, tmp_path):
        output = str(tmp_path / "charmap.pdf")
        testargs = [
            "phykit",
            "character_map",
            "-t", f"{here.parent.parent.parent}/sample_files/tree_character_map.tre",
            "-d", f"{here.parent.parent.parent}/sample_files/character_matrix_simple.tsv",
            "-o", output,
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert Path(output).exists()

    @patch("builtins.print")
    def test_character_map_circular(self, mocked_print, tmp_path):
        """--circular flag produces a circular layout plot without error."""
        output = str(tmp_path / "charmap_circular.png")
        testargs = [
            "phykit",
            "character_map",
            "-t", f"{here.parent.parent.parent}/sample_files/tree_character_map.tre",
            "-d", f"{here.parent.parent.parent}/sample_files/character_matrix_simple.tsv",
            "-o", output,
            "--circular",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert Path(output).exists()
        assert Path(output).stat().st_size > 0

    @patch("builtins.print")
    def test_character_filter_preserves_ci_ri(self, mocked_print, tmp_path):
        """CI and RI should be the same with or without --characters filter,
        since the filter only affects the plot, not the computation."""
        output1 = str(tmp_path / "charmap_all.png")
        testargs1 = [
            "phykit",
            "character_map",
            "-t", f"{here.parent.parent.parent}/sample_files/tree_character_map.tre",
            "-d", f"{here.parent.parent.parent}/sample_files/character_matrix_simple.tsv",
            "-o", output1,
            "--json",
        ]
        with patch.object(sys, "argv", testargs1):
            Phykit()
        payload_all = json.loads(mocked_print.call_args.args[0])

        mocked_print.reset_mock()

        output2 = str(tmp_path / "charmap_filtered.png")
        testargs2 = [
            "phykit",
            "character_map",
            "-t", f"{here.parent.parent.parent}/sample_files/tree_character_map.tre",
            "-d", f"{here.parent.parent.parent}/sample_files/character_matrix_simple.tsv",
            "-o", output2,
            "--characters", "0,1",
            "--json",
        ]
        with patch.object(sys, "argv", testargs2):
            Phykit()
        payload_filtered = json.loads(mocked_print.call_args.args[0])

        assert payload_all["ci"] == payload_filtered["ci"]
        assert payload_all["ri"] == payload_filtered["ri"]
