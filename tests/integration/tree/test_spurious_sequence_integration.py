from mock import patch, call
from pathlib import Path
import pytest
import sys
import json

from phykit.phykit import Phykit

here = Path(__file__)


@pytest.mark.integration
class TestSpuriousSequence(object):
    @patch("builtins.print")
    def test_spurious_sequence(self, mocked_print):
        expected_result = "None"
        testargs = [
            "phykit",
            "spurious_sequence",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_spurious_sequence_alias0(self, mocked_print):
        expected_result = "None"
        testargs = [
            "phykit",
            "spurious_seq",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_spurious_sequence_alias1(self, mocked_print):
        expected_result = "None"
        testargs = [
            "phykit",
            "ss",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_spurious_sequence_custom_factor(self, mocked_print):
        expected_result = (
            "monkey\t100.8593\t38.0791\t19.0396\n"
            "cat\t47.1407\t38.0791\t19.0396"
        )
        testargs = [
            "phykit",
            "spurious_sequence",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "-f",
            "2"
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_spurious_sequence_incorrect_file_path(self, mocked_print):
        testargs = [
            "phykit",
            "spurious_sequence",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tr",  # Invalid file path
        ]
        with patch.object(sys, "argv", testargs):
            with pytest.raises(SystemExit) as pytest_wrapped_e:
                Phykit()

        assert pytest_wrapped_e.type is SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_spurious_sequence_json_none(self, mocked_print):
        testargs = [
            "phykit",
            "spurious_sequence",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        payload = json.loads(mocked_print.call_args.args[0])
        assert payload == {"rows": [], "spurious_sequences": []}

    @patch("builtins.print")
    def test_spurious_sequence_json_custom_factor(self, mocked_print):
        testargs = [
            "phykit",
            "spurious_sequence",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "-f",
            "2",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["rows"][0] == payload["spurious_sequences"][0]
        assert payload["spurious_sequences"][0] == {
            "taxon": "monkey",
            "branch_length": 100.8593,
            "threshold": 38.0791,
            "median": 19.0396,
        }

    @patch("builtins.print")
    def test_spurious_sequence_diameter_impact_json(self, mocked_print, tmp_path):
        tree_path = tmp_path / "outlier.tre"
        tree_path.write_text(
            "(outlier:100,"
            + ",".join(f"tip_{index}:1" for index in range(1, 20))
            + ");\n"
        )
        testargs = [
            "phykit",
            "spurious_sequence",
            str(tree_path),
            "--method",
            "diameter-impact",
            "--json",
        ]

        with patch.object(sys, "argv", testargs):
            Phykit()

        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["method"] == "diameter-impact"
        assert payload["scope"] == "per-gene"
        assert payload["tree_count"] == 1
        assert payload["rows"] == payload["spurious_sequences"]
        assert payload["rows"][0]["taxon"] == "outlier"
        assert payload["rows"][0]["diameter_before"] == 101.0
        assert payload["rows"][0]["diameter_after"] == 2.0
        assert payload["rows"][0]["p_value"] < 0.05

    @patch("builtins.print")
    def test_spurious_sequence_diameter_impact_text(self, mocked_print, tmp_path):
        tree_path = tmp_path / "outlier.tre"
        tree_path.write_text(
            "(outlier:100,"
            + ",".join(f"tip_{index}:1" for index in range(1, 20))
            + ");\n"
        )
        testargs = [
            "phykit",
            "spurious_sequence",
            str(tree_path),
            "--method",
            "diameter-impact",
        ]

        with patch.object(sys, "argv", testargs):
            Phykit()

        fields = mocked_print.call_args.args[0].split("\t")
        assert fields[0] == "outlier"
        assert fields[1] == "3.921973"
        assert float(fields[2]) < 0.05
        assert fields[3:] == ["1", "100.0", "101.0", "2.0"]

    @patch("builtins.print")
    def test_spurious_sequence_per_species_tree_list(self, mocked_print, tmp_path):
        tree_names = []
        ordinary_tips = ",".join(
            f"tip_{index}:1" for index in range(1, 20)
        )
        for index in range(5):
            tree_name = f"gene_{index}.tre"
            outlier_length = 100 if index == 0 else 1
            (tmp_path / tree_name).write_text(
                f"(outlier:{outlier_length},{ordinary_tips});\n"
            )
            tree_names.append(tree_name)
        list_path = tmp_path / "trees.txt"
        list_path.write_text("\n".join(tree_names) + "\n")
        testargs = [
            "phykit",
            "spurious_sequence",
            str(list_path),
            "--method",
            "diameter-impact",
            "--tree-list",
            "--per-species",
            "--json",
        ]

        with patch.object(sys, "argv", testargs):
            Phykit()

        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["scope"] == "per-species"
        assert payload["tree_count"] == 5
        assert [row["taxon"] for row in payload["rows"]] == ["outlier"]
        assert payload["rows"][0]["tree"] == str(tmp_path / "gene_0.tre")

    @patch("builtins.print")
    def test_spurious_sequence_diameter_requires_branch_lengths(
        self,
        mocked_print,
        tmp_path,
    ):
        tree_path = tmp_path / "missing_lengths.tre"
        tree_path.write_text("(a:1,b:1,c:1,d);\n")
        testargs = [
            "phykit",
            "spurious_sequence",
            str(tree_path),
            "--method",
            "diameter-impact",
        ]

        with patch.object(sys, "argv", testargs):
            with pytest.raises(SystemExit) as error:
                Phykit()

        assert error.value.code == 2
        assert "must have lengths" in mocked_print.call_args.args[0]
