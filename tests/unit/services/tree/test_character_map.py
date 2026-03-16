"""
Unit tests for character_map (discrete character mapping on a phylogeny).

Tests matrix parsing, output modes, plot creation, and edge cases.
Uses the sample tree and character matrix in tests/sample_files/.
"""
import json
import pytest
from argparse import Namespace
from pathlib import Path

from phykit.services.tree.character_map import CharacterMap


TREE_PATH = "tests/sample_files/tree_character_map.tre"
MATRIX_PATH = "tests/sample_files/character_matrix_simple.tsv"


def _make_args(tmp_path, **overrides):
    """Build a default Namespace for CharacterMap, with optional overrides."""
    defaults = dict(
        tree=TREE_PATH,
        data=MATRIX_PATH,
        output=str(tmp_path / "charmap.png"),
        optimization="acctran",
        phylogram=False,
        characters=None,
        verbose=False,
        json=False,
    )
    defaults.update(overrides)
    return Namespace(**defaults)


class TestParseCharacterMatrix:
    def test_parse_character_matrix(self, tmp_path):
        args = _make_args(tmp_path)
        cm = CharacterMap(args)
        char_names, tip_states = cm._parse_character_matrix(MATRIX_PATH)

        assert char_names == ["char0", "char1", "char2", "char3",
                              "char4", "char5", "char6", "char7"]
        assert len(tip_states) == 6
        assert tip_states["A"] == ["0", "1", "0", "0", "0", "1", "1", "0"]
        assert tip_states["F"] == ["0", "0", "0", "0", "0", "0", "0", "0"]

    def test_missing_data_preserved(self, tmp_path):
        """Question marks and dashes in the matrix should be preserved."""
        matrix_path = tmp_path / "missing.tsv"
        matrix_path.write_text(
            "taxon\tc0\tc1\n"
            "A\t0\t?\n"
            "B\t-\t1\n"
            "C\t1\t0\n"
        )
        args = _make_args(tmp_path, data=str(matrix_path))
        cm = CharacterMap(args)
        char_names, tip_states = cm._parse_character_matrix(str(matrix_path))

        assert tip_states["A"][1] == "?"
        assert tip_states["B"][0] == "-"

    def test_row_length_mismatch_raises(self, tmp_path):
        """Rows with wrong number of columns should raise an error."""
        matrix_path = tmp_path / "bad.tsv"
        matrix_path.write_text(
            "taxon\tc0\tc1\n"
            "A\t0\n"
        )
        args = _make_args(tmp_path, data=str(matrix_path))
        cm = CharacterMap(args)
        with pytest.raises(SystemExit):
            cm._parse_character_matrix(str(matrix_path))


class TestCharacterMapInit:
    def test_init_sets_fields(self, tmp_path):
        args = _make_args(tmp_path)
        cm = CharacterMap(args)
        assert cm.tree_file_path == TREE_PATH
        assert cm.data_path == MATRIX_PATH
        assert cm.optimization == "acctran"
        assert cm.phylogram is False
        assert cm.characters_filter is None
        assert cm.verbose is False
        assert cm.json_output is False

    def test_characters_filter_parsed(self, tmp_path):
        args = _make_args(tmp_path, characters="0,1,3")
        cm = CharacterMap(args)
        assert cm.characters_filter == [0, 1, 3]


class TestCharacterMapPlot:
    def test_creates_png(self, tmp_path):
        args = _make_args(tmp_path)
        cm = CharacterMap(args)
        cm.run()
        out = Path(args.output)
        assert out.exists()
        assert out.stat().st_size > 0

    def test_creates_pdf(self, tmp_path):
        args = _make_args(tmp_path, output=str(tmp_path / "charmap.pdf"))
        cm = CharacterMap(args)
        cm.run()
        assert Path(args.output).exists()
        assert Path(args.output).stat().st_size > 0

    def test_cladogram_default(self, tmp_path):
        """Default mode (no --phylogram) creates a file."""
        args = _make_args(tmp_path)
        cm = CharacterMap(args)
        cm.run()
        assert Path(args.output).exists()

    def test_phylogram_mode(self, tmp_path):
        """--phylogram flag creates a file."""
        args = _make_args(tmp_path, phylogram=True)
        cm = CharacterMap(args)
        cm.run()
        assert Path(args.output).exists()
        assert Path(args.output).stat().st_size > 0

    def test_character_filter(self, tmp_path):
        """--characters '0,1' creates a file (filters plotted changes)."""
        args = _make_args(tmp_path, characters="0,1")
        cm = CharacterMap(args)
        cm.run()
        assert Path(args.output).exists()
        assert Path(args.output).stat().st_size > 0

    def test_deltran_optimization(self, tmp_path):
        """DELTRAN optimization should also produce output."""
        args = _make_args(tmp_path, optimization="deltran")
        cm = CharacterMap(args)
        cm.run()
        assert Path(args.output).exists()


class TestCharacterMapJson:
    def test_json_output(self, mocker, tmp_path):
        args = _make_args(tmp_path, json=True)
        cm = CharacterMap(args)
        mocked_json = mocker.patch(
            "phykit.services.tree.character_map.print_json"
        )
        cm.run()

        payload = mocked_json.call_args.args[0]
        assert "n_characters" in payload
        assert "n_informative" in payload
        assert "tree_length" in payload
        assert "ci" in payload
        assert "ri" in payload
        assert "optimization" in payload
        assert "output_file" in payload
        assert "characters" in payload
        assert payload["n_characters"] == 8
        assert payload["optimization"] == "acctran"
        assert isinstance(payload["characters"], list)
        assert len(payload["characters"]) == 8

        # Each character entry has required keys
        for char_entry in payload["characters"]:
            assert "index" in char_entry
            assert "name" in char_entry
            assert "steps" in char_entry
            assert "ci" in char_entry
            assert "ri" in char_entry
            assert "changes" in char_entry

    def test_json_deltran(self, mocker, tmp_path):
        args = _make_args(tmp_path, json=True, optimization="deltran")
        cm = CharacterMap(args)
        mocked_json = mocker.patch(
            "phykit.services.tree.character_map.print_json"
        )
        cm.run()
        payload = mocked_json.call_args.args[0]
        assert payload["optimization"] == "deltran"


class TestCharacterMapVerbose:
    def test_verbose_output(self, tmp_path, capsys):
        args = _make_args(tmp_path, verbose=True)
        cm = CharacterMap(args)
        cm.run()
        captured = capsys.readouterr()
        assert "Optimization: acctran" in captured.out
        assert "Characters: 8" in captured.out
        assert "Tree length:" in captured.out
        assert "CI:" in captured.out
        assert "RI:" in captured.out

    def test_summary_output(self, tmp_path, capsys):
        args = _make_args(tmp_path)
        cm = CharacterMap(args)
        cm.run()
        captured = capsys.readouterr()
        assert "Optimization: acctran" in captured.out
        assert "Characters: 8" in captured.out
        assert f"Output: {args.output}" in captured.out


class TestCharacterMapEdgeCases:
    def test_too_few_shared_taxa(self, tmp_path):
        """Trees/matrices with fewer than 3 shared taxa should error."""
        tree_path = tmp_path / "small.tre"
        tree_path.write_text("(X:1,Y:1);\n")
        matrix_path = tmp_path / "small.tsv"
        matrix_path.write_text("taxon\tc0\nX\t0\nY\t1\n")
        args = _make_args(
            tmp_path,
            tree=str(tree_path),
            data=str(matrix_path),
        )
        cm = CharacterMap(args)
        with pytest.raises(SystemExit):
            cm.run()

    def test_missing_matrix_file(self, tmp_path):
        """Non-existent matrix file should raise an error."""
        args = _make_args(tmp_path, data=str(tmp_path / "nonexistent.tsv"))
        cm = CharacterMap(args)
        with pytest.raises(SystemExit):
            cm.run()

    def test_ladderize(self, tmp_path):
        """Ladderize option should work without error."""
        args = _make_args(tmp_path, ladderize=True)
        cm = CharacterMap(args)
        cm.run()
        assert Path(args.output).exists()
