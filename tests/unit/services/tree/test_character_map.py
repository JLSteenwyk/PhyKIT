"""
Unit tests for character_map (discrete character mapping on a phylogeny).

Tests matrix parsing, output modes, plot creation, and edge cases.
Uses the sample tree and character matrix in tests/sample_files/.
"""
import builtins
import json
import subprocess
import sys
import pytest
from argparse import Namespace
from collections import Counter
from io import StringIO
from pathlib import Path

from Bio import Phylo
from Bio.Phylo.BaseTree import Clade, TreeMixin

from phykit.helpers.parsimony_utils import build_parent_map, retention_index
import phykit.services.tree.character_map as character_map_module
from phykit.services.tree.character_map import CharacterMap


def test_module_import_does_not_import_numpy():
    code = """
import sys
import phykit.services.tree.character_map as module
assert callable(module.print_json)
assert callable(module.build_parent_map)
assert callable(module.sankoff_downpass)
assert callable(module.sankoff_uppass)
assert "typing" not in sys.modules
assert "pickle" not in sys.modules
assert "json" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
assert "phykit.helpers.plot_config" not in sys.modules
assert "phykit.helpers.circular_layout" not in sys.modules
assert "phykit.helpers.color_annotations" not in sys.modules
assert "phykit.helpers.parsimony_utils" not in sys.modules
assert "numpy" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


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
        allow_taxon_mismatch=False,
        change_marker_size=None,
        change_fontsize=None,
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

    @pytest.mark.parametrize("taxon", ["", "   "])
    def test_empty_taxon_label_raises(self, tmp_path, taxon):
        matrix_path = tmp_path / "empty-taxon.tsv"
        matrix_path.write_text(
            "taxon\tc0\n"
            f"{taxon}\t0\n"
        )

        with pytest.raises(SystemExit) as exc:
            CharacterMap._parse_character_matrix(str(matrix_path))

        assert exc.value.code == 2
        assert "empty taxon label" in exc.value.messages[0]

    def test_duplicate_taxon_label_raises(self, tmp_path):
        matrix_path = tmp_path / "duplicate-taxon.tsv"
        matrix_path.write_text(
            "taxon\tc0\n"
            "A\t0\n"
            "A\t1\n"
        )

        with pytest.raises(SystemExit) as exc:
            CharacterMap._parse_character_matrix(str(matrix_path))

        assert exc.value.code == 2
        assert "duplicate taxon label 'A'" in exc.value.messages[0]

    def test_parse_character_matrix_streams_nonblank_rows(self, monkeypatch, tmp_path):
        class StreamingOnlyFile:
            def __init__(self):
                self.lines = [
                    "\n",
                    "taxon\tc0\tc1\n",
                    "A\t0\t1\n",
                    "\n",
                    "B\t1\t0\n",
                ]
                self.index = 0

            def __enter__(self):
                return self

            def __exit__(self, exc_type, exc, traceback):
                return False

            def __iter__(self):
                return self

            def __next__(self):
                if self.index >= len(self.lines):
                    raise StopIteration
                line = self.lines[self.index]
                self.index += 1
                return line

            def read(self, *_args, **_kwargs):
                raise AssertionError("character matrix parser should stream rows")

            def readlines(self, *_args, **_kwargs):
                raise AssertionError("character matrix parser should stream rows")

        monkeypatch.setattr("builtins.open", lambda path: StreamingOnlyFile())

        char_names, tip_states = CharacterMap._parse_character_matrix("matrix.tsv")

        assert char_names == ["c0", "c1"]
        assert tip_states == {"A": ["0", "1"], "B": ["1", "0"]}


class TestCharacterMapSharedTaxaSetup:
    def test_ordered_all_shared_returns_original_states_without_sets(self, monkeypatch):
        tip_states = {"A": ["0"], "B": ["1"], "C": ["0"]}

        def fail_set(*_args, **_kwargs):
            raise AssertionError("ordered all-shared character data should skip sets")

        monkeypatch.setattr(builtins, "set", fail_set)

        shared_count, tips_to_prune, filtered = (
            CharacterMap._shared_character_taxa_setup(["A", "B", "C"], tip_states)
        )

        assert shared_count == 3
        assert tips_to_prune == []
        assert filtered is tip_states

    def test_partial_overlap_fails_with_names_from_both_inputs(self):
        tip_states = {"A": ["0"], "C": ["1"], "off_tree": ["0"]}

        with pytest.raises(SystemExit) as exc:
            CharacterMap._shared_character_taxa_setup(
                ["A", "B", "C", "D"],
                tip_states,
            )

        assert exc.value.code == 2
        assert exc.value.messages == [
            "Taxon labels differ between the tree and character matrix.",
            "2 taxa in tree but not in character matrix: B, D",
            "1 taxon in character matrix but not in tree: off_tree",
            "Use --allow-taxon-mismatch to analyze only shared taxa.",
        ]

    def test_allow_mismatch_warns_filters_and_preserves_prune_order(self, capsys):
        tip_states = {"A": ["0"], "C": ["1"], "off_tree": ["0"]}

        shared_count, tips_to_prune, filtered = (
            CharacterMap._shared_character_taxa_setup(
                ["A", "B", "C", "D"],
                tip_states,
                allow_taxon_mismatch=True,
            )
        )

        assert capsys.readouterr().err == (
            "Warning: 2 taxa in tree but not in character matrix: B, D\n"
            "Warning: 1 taxon in character matrix but not in tree: off_tree\n"
        )
        assert shared_count == 2
        assert tips_to_prune == ["B", "D"]
        assert filtered == {"A": ["0"], "C": ["1"]}

    def test_matching_taxa_in_different_order_are_accepted(self):
        tip_states = {"C": ["0"], "A": ["1"], "B": ["0"]}

        shared_count, tips_to_prune, filtered = (
            CharacterMap._shared_character_taxa_setup(
                ["A", "B", "C"],
                tip_states,
            )
        )

        assert shared_count == 3
        assert tips_to_prune == []
        assert filtered is tip_states

    def test_taxon_matching_is_case_sensitive(self):
        tip_states = {"a": ["0"], "B": ["1"], "C": ["0"]}

        with pytest.raises(SystemExit) as exc:
            CharacterMap._shared_character_taxa_setup(
                ["A", "B", "C"],
                tip_states,
            )

        assert "1 taxon in tree but not in character matrix: A" in exc.value.messages
        assert "1 taxon in character matrix but not in tree: a" in exc.value.messages

    def test_taxon_matching_does_not_strip_whitespace(self):
        tip_states = {" A": ["0"], "B": ["1"], "C": ["0"]}

        with pytest.raises(SystemExit) as exc:
            CharacterMap._shared_character_taxa_setup(
                ["A", "B", "C"],
                tip_states,
            )

        assert "1 taxon in tree but not in character matrix: A" in exc.value.messages
        assert "1 taxon in character matrix but not in tree:  A" in exc.value.messages

    def test_duplicate_tree_taxon_labels_raise(self):
        with pytest.raises(SystemExit) as exc:
            CharacterMap._shared_character_taxa_setup(
                ["A", "A", "B"],
                {"A": ["0"], "B": ["1"]},
            )

        assert exc.value.messages[0] == "Tree contains duplicate taxon labels: A"

    @pytest.mark.parametrize("empty_name", [None, "", "   "])
    def test_empty_tree_taxon_labels_raise(self, empty_name):
        with pytest.raises(SystemExit) as exc:
            CharacterMap._shared_character_taxa_setup(
                ["A", "B", empty_name],
                {"A": ["0"], "B": ["1"], "C": ["0"]},
            )

        assert "empty or unnamed taxon labels" in exc.value.messages[0]


class TestCharacterMapInit:
    def test_init_sets_fields(self, tmp_path):
        args = _make_args(tmp_path)
        cm = CharacterMap(args)
        assert cm.tree_file_path == TREE_PATH
        assert cm.data_path == MATRIX_PATH
        assert cm.optimization == "acctran"
        assert cm.phylogram is False
        assert cm.characters_filter is None
        assert cm.allow_taxon_mismatch is False
        assert cm.change_marker_size is None
        assert cm.change_fontsize is None
        assert cm.verbose is False
        assert cm.json_output is False

    def test_characters_filter_parsed(self, tmp_path):
        args = _make_args(tmp_path, characters="0,1,3")
        cm = CharacterMap(args)
        assert cm.characters_filter == [0, 1, 3]

    def test_change_plot_controls_are_stored(self, tmp_path):
        args = _make_args(
            tmp_path,
            change_marker_size=125.0,
            change_fontsize=8.5,
        )

        cm = CharacterMap(args)

        assert cm.change_marker_size == 125.0
        assert cm.change_fontsize == 8.5

    @pytest.mark.parametrize(
        ("option", "value", "message"),
        [
            ("change_marker_size", 0, "--change-marker-size"),
            ("change_marker_size", -1, "--change-marker-size"),
            ("change_fontsize", 0, "--change-fontsize"),
            ("change_fontsize", -1, "--change-fontsize"),
        ],
    )
    def test_change_plot_controls_must_be_positive(
        self,
        tmp_path,
        option,
        value,
        message,
    ):
        args = _make_args(tmp_path, **{option: value})

        with pytest.raises(SystemExit) as exc:
            CharacterMap(args)

        assert exc.value.code == 2
        assert message in exc.value.messages[0]


class TestCharacterMapStateSummary:
    def test_summarize_character_states_counts_wildcards_and_informative_sites(self):
        tip_states = {
            "A": ["0", "0", "?", "1"],
            "B": ["0", "1", "-", "1"],
            "C": ["1", "1", "2", "2"],
            "D": ["1", "2", "2", "2"],
            "E": ["?", "2", "2", "-"],
        }

        n_states, states_by_char, n_informative = (
            CharacterMap._summarize_character_states(tip_states)
        )

        assert n_states == [2, 3, 1, 2]
        assert states_by_char == [
            ["0", "0", "1", "1", "?"],
            ["0", "1", "1", "2", "2"],
            ["?", "-", "2", "2", "2"],
            ["1", "1", "2", "2", "-"],
        ]
        assert n_informative == 3

    def test_summarize_character_states_ignores_wildcard_counts(self):
        tip_states = {
            "A": ["?", "0", "0"],
            "B": ["-", "1", "?"],
            "C": ["?", "2", "1"],
            "D": ["-", "3", "-"],
        }

        n_states, states_by_char, n_informative = (
            CharacterMap._summarize_character_states(tip_states)
        )

        assert n_states == [0, 4, 2]
        assert states_by_char == [
            ["?", "-", "?", "-"],
            ["0", "1", "2", "3"],
            ["0", "?", "1", "-"],
        ]
        assert n_informative == 0

    def test_count_aware_retention_index_matches_reference(self):
        tip_states = {
            "A": ["0", "0", "?", "1", "?"],
            "B": ["0", "1", "-", "1", "-"],
            "C": ["1", "1", "2", "2", "?"],
            "D": ["1", "2", "2", "2", "-"],
            "E": ["?", "2", "2", "-", "?"],
        }
        observed_per_char = [1, 3, 0, 2, 0]

        _, states_by_char, _, counts_per_char = (
            CharacterMap._summarize_character_states_with_counts(tip_states)
        )

        assert CharacterMap._retention_index_from_counts(
            counts_per_char,
            observed_per_char,
        ) == retention_index(states_by_char, observed_per_char)

    def test_count_aware_retention_index_scans_counts_once(self):
        class CountingCounter(Counter):
            def __init__(self, *args, **kwargs):
                super().__init__(*args, **kwargs)
                self.values_calls = 0

            def values(self):
                self.values_calls += 1
                return super().values()

        counts_per_char = [
            CountingCounter({"0": 2, "1": 2}),
            CountingCounter({"0": 1, "1": 2, "2": 2}),
            CountingCounter(),
        ]

        ri_per_char, ri_overall = CharacterMap._retention_index_from_counts(
            counts_per_char,
            [1, 3, 0],
        )

        assert ri_per_char == [1.0, 0.0, None]
        assert ri_overall == pytest.approx(0.5)
        assert [counts.values_calls for counts in counts_per_char] == [1, 1, 1]

    def test_counts_only_summary_matches_full_summary_without_state_lists(self):
        tip_states = {
            "A": ["0", "0", "?", "1"],
            "B": ["0", "1", "-", "1"],
            "C": ["1", "1", "2", "2"],
            "D": ["1", "2", "2", "2"],
            "E": ["?", "2", "2", "-"],
        }

        full = CharacterMap._summarize_character_states_with_counts(tip_states)
        counts_only = CharacterMap._summarize_character_states_with_counts(
            tip_states,
            include_states=False,
        )

        assert counts_only[0] == full[0]
        assert counts_only[1] == []
        assert counts_only[2] == full[2]
        assert counts_only[3] == full[3]
        assert counts_only[3][0] == Counter({"0": 2, "1": 2})

    def test_counts_only_ascii_summary_scans_rows_once(self):
        class SinglePassRows:
            def __init__(self, rows):
                self.rows = rows
                self.iterations = 0

            def __iter__(self):
                self.iterations += 1
                if self.iterations > 1:
                    raise AssertionError("rows should be scanned once")
                return iter(self.rows)

        class TipStates:
            def __init__(self, rows):
                self.rows = SinglePassRows(rows)

            def values(self):
                return self.rows

        tip_states = TipStates(
            [
                ["0", "0", "?", "1"],
                ["0", "1", "-", "1"],
                ["1", "1", "2", "2"],
                ["1", "2", "2", "2"],
            ]
        )

        n_states, states_by_char, n_informative, counts_per_char = (
            CharacterMap._summarize_character_ascii_counts_only(tip_states)
        )

        assert tip_states.rows.iterations == 1
        assert n_states == [2, 3, 1, 2]
        assert states_by_char == []
        assert n_informative == 2
        assert counts_per_char[0] == Counter({"0": 2, "1": 2})

    def test_ascii_symbol_counts_by_char_matches_equality_reference(self):
        import numpy as np

        alphabet = b"0123456789ABCDEFGHIJKLMNOP"
        matrix = np.frombuffer(
            alphabet
            + alphabet[::-1]
            + alphabet[4:]
            + alphabet[:4],
            dtype=np.uint8,
        ).reshape(3, len(alphabet))
        symbols = np.unique(matrix)

        observed = CharacterMap._ascii_symbol_counts_by_char(matrix, symbols)
        expected = np.vstack(
            [np.count_nonzero(matrix == symbol, axis=0) for symbol in symbols]
        )

        np.testing.assert_array_equal(observed, expected)

    def test_ascii_symbol_counts_small_alphabet_skips_bincount(self, monkeypatch):
        import numpy as np

        matrix = np.frombuffer(
            b"0101"
            b"1010"
            b"0011",
            dtype=np.uint8,
        ).reshape(3, 4)
        symbols = np.unique(matrix)

        def fail_bincount(*_args, **_kwargs):
            raise AssertionError("small state alphabets should keep equality counts")

        monkeypatch.setattr(np, "bincount", fail_bincount)

        observed = CharacterMap._ascii_symbol_counts_by_char(matrix, symbols)
        expected = np.vstack(
            [np.count_nonzero(matrix == symbol, axis=0) for symbol in symbols]
        )

        np.testing.assert_array_equal(observed, expected)

    def test_ascii_symbol_counts_twelve_state_alphabet_uses_bincount(self, monkeypatch):
        import numpy as np

        alphabet = b"0123456789AB"
        matrix = np.frombuffer(alphabet * 4, dtype=np.uint8).reshape(4, 12)
        symbols = np.unique(matrix)
        bincount_calls = 0
        original_bincount = np.bincount

        def counting_bincount(*args, **kwargs):
            nonlocal bincount_calls
            bincount_calls += 1
            return original_bincount(*args, **kwargs)

        monkeypatch.setattr(np, "bincount", counting_bincount)

        observed = CharacterMap._ascii_symbol_counts_by_char(matrix, symbols)
        expected = np.vstack(
            [np.count_nonzero(matrix == symbol, axis=0) for symbol in symbols]
        )

        assert bincount_calls == 1
        np.testing.assert_array_equal(observed, expected)

    def test_ascii_counts_per_char_matches_reference_and_omits_zeros(self):
        import numpy as np

        symbol_labels = ["0", "1", "2", "3"]
        symbol_counts = np.array(
            [
                [2, 0, 1],
                [0, 4, 1],
                [3, 0, 0],
                [0, 5, 1],
            ]
        )

        observed = CharacterMap._ascii_counts_per_char(
            symbol_counts,
            symbol_labels,
        )

        assert observed == [
            Counter({"0": 2, "2": 3}),
            Counter({"1": 4, "3": 5}),
            Counter({"0": 1, "1": 1, "3": 1}),
        ]
        assert all(count for counts in observed for count in counts.values())

    def test_ascii_counts_per_char_dense_columns_accept_state_counts(self):
        import numpy as np

        symbol_labels = ["0", "1", "2"]
        symbol_counts = np.array(
            [
                [2, 1],
                [4, 3],
                [5, 6],
            ]
        )

        observed = CharacterMap._ascii_counts_per_char(
            symbol_counts,
            symbol_labels,
            [3, 3],
        )

        assert observed == [
            Counter({"0": 2, "1": 4, "2": 5}),
            Counter({"0": 1, "1": 3, "2": 6}),
        ]

    def test_full_summary_reuses_ascii_summary_for_single_character_states(
        self, monkeypatch
    ):
        tip_states = {
            "A": ["0", "1"],
            "B": ["1", "1"],
        }
        expected_counts = (
            [2, 1],
            [],
            0,
            [Counter({"0": 1, "1": 1}), Counter({"1": 2})],
        )
        calls = []

        def fake_ascii_summary(observed_tip_states, include_states=True):
            calls.append((observed_tip_states, include_states))
            return expected_counts

        monkeypatch.setattr(
            CharacterMap,
            "_summarize_character_ascii",
            staticmethod(fake_ascii_summary),
        )

        observed = CharacterMap._summarize_character_states_with_counts(
            tip_states,
        )

        assert calls == [(tip_states, True)]
        assert observed == expected_counts

    def test_counts_only_summary_falls_back_for_multicharacter_states(self):
        tip_states = {
            "A": ["red", "α", "?"],
            "B": ["red", "β", "-"],
            "C": ["blue", "β", "large"],
            "D": ["blue", "α", "large"],
        }

        full = CharacterMap._summarize_character_states_with_counts(tip_states)
        counts_only = CharacterMap._summarize_character_states_with_counts(
            tip_states,
            include_states=False,
        )

        assert counts_only[0] == full[0]
        assert counts_only[1] == []
        assert counts_only[2] == full[2]
        assert counts_only[3] == full[3]


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

    @pytest.mark.parametrize("circular", [False, True])
    def test_plot_character_map_uses_direct_tree_traversal(
        self, monkeypatch, tmp_path, circular
    ):
        args = _make_args(tmp_path, circular=circular)
        cm = CharacterMap(args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        parent_map = build_parent_map(tree)
        node_labels = CharacterMap._assign_node_labels(tree)

        def fail_traversal(*args, **kwargs):
            raise AssertionError("plot setup should use direct traversal")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_traversal)
        monkeypatch.setattr(TreeMixin, "get_terminals", fail_traversal)

        cm._plot_character_map(tree, {}, node_labels, parent_map)

        out = Path(args.output)
        assert out.exists()
        assert out.stat().st_size > 0

    def test_iter_preorder_preserves_order_without_reversed(self):
        class NoReversedList(list):
            def __reversed__(self):
                raise AssertionError("_iter_preorder should push children directly")

        root = Clade(name="root")
        left = Clade(name="left")
        middle = Clade(name="middle")
        right = Clade(name="right")
        left.clades = NoReversedList([Clade(name="left_a"), Clade(name="left_b")])
        middle.clades = NoReversedList([Clade(name="middle_a")])
        right.clades = NoReversedList()
        root.clades = NoReversedList([left, middle, right])

        order = [clade.name for clade in CharacterMap._iter_preorder(root)]

        assert order == [
            "root",
            "left",
            "left_a",
            "left_b",
            "middle",
            "middle_a",
            "right",
        ]

    def test_assign_node_labels_uses_direct_preorder(self, monkeypatch):
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        expected = CharacterMap._assign_node_labels(tree)

        def fail_find_clades(*args, **kwargs):
            raise AssertionError("standard tree labels should use direct traversal")

        def fail_is_terminal(self):
            raise AssertionError("terminal checks should use child lists")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_find_clades)
        monkeypatch.setattr(Clade, "is_terminal", fail_is_terminal)

        assert CharacterMap._assign_node_labels(tree) == expected
        assert list(expected.values()) == [
            "node_0",
            "node_1",
            "A",
            "B",
            "node_2",
            "C",
            "D",
        ]

    def test_plot_character_map_reuses_preorder_for_node_positions(
        self, monkeypatch, tmp_path
    ):
        pytest.importorskip("matplotlib")
        import phykit.helpers.plot_config as plot_config

        args = _make_args(tmp_path)
        cm = CharacterMap(args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        parent_map = build_parent_map(tree)
        node_labels = CharacterMap._assign_node_labels(tree)
        expected_preorder_ids = [
            id(clade) for clade in CharacterMap._iter_preorder(tree.root)
        ]
        original_compute_node_positions = plot_config.compute_node_positions
        calls = []

        def assert_preorder_reused(
            tree_arg, parent_map_arg, cladogram=False, preorder_clades=None
        ):
            assert preorder_clades is not None
            calls.append([id(clade) for clade in preorder_clades])
            return original_compute_node_positions(
                tree_arg,
                parent_map_arg,
                cladogram=cladogram,
                preorder_clades=preorder_clades,
            )

        monkeypatch.setattr(
            plot_config, "compute_node_positions", assert_preorder_reused
        )

        cm._plot_character_map(tree, {}, node_labels, parent_map)

        assert calls == [expected_preorder_ids]
        out = Path(args.output)
        assert out.exists()
        assert out.stat().st_size > 0

    def test_plot_character_map_reuses_clade_lists_for_circular_coords(
        self, monkeypatch, tmp_path
    ):
        pytest.importorskip("matplotlib")
        import phykit.helpers.circular_layout as circular_layout

        args = _make_args(
            tmp_path,
            circular=True,
            ylabel_fontsize=0,
            no_title=True,
        )
        cm = CharacterMap(args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        parent_map = build_parent_map(tree)
        node_labels = CharacterMap._assign_node_labels(tree)
        preorder_clades = list(CharacterMap._iter_preorder(tree.root))
        tips = [clade for clade in preorder_clades if not clade.clades]
        expected_preorder_ids = [id(clade) for clade in preorder_clades]
        expected_terminal_ids = [id(tip) for tip in tips]
        original_compute_circular_coords = circular_layout.compute_circular_coords
        calls = []

        def assert_clade_lists_reused(
            tree_arg,
            node_x_arg,
            parent_map_arg,
            preorder_clades=None,
            terminal_clades=None,
        ):
            assert preorder_clades is not None
            assert terminal_clades is not None
            calls.append(
                (
                    [id(clade) for clade in preorder_clades],
                    [id(clade) for clade in terminal_clades],
                )
            )
            return original_compute_circular_coords(
                tree_arg,
                node_x_arg,
                parent_map_arg,
                preorder_clades=preorder_clades,
                terminal_clades=terminal_clades,
            )

        monkeypatch.setattr(
            circular_layout, "compute_circular_coords", assert_clade_lists_reused
        )

        cm._plot_character_map(tree, {}, node_labels, parent_map)

        assert calls == [(expected_preorder_ids, expected_terminal_ids)]
        out = Path(args.output)
        assert out.exists()
        assert out.stat().st_size > 0

    def test_rectangular_plot_batches_base_branches(self, monkeypatch, tmp_path):
        pytest.importorskip("matplotlib")
        import matplotlib.axes

        args = _make_args(tmp_path)
        cm = CharacterMap(args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        parent_map = build_parent_map(tree)
        node_labels = CharacterMap._assign_node_labels(tree)

        def fail_plot(*args, **kwargs):
            raise AssertionError("base rectangular branches should be batched")

        monkeypatch.setattr(matplotlib.axes.Axes, "plot", fail_plot)

        cm._plot_character_map(tree, {}, node_labels, parent_map)

        out = Path(args.output)
        assert out.exists()
        assert out.stat().st_size > 0

    @pytest.mark.parametrize("circular", [False, True])
    def test_plot_batches_color_clade_overlay(self, monkeypatch, tmp_path, circular):
        pytest.importorskip("matplotlib")
        from matplotlib.axes import Axes
        from matplotlib.collections import LineCollection

        color_file = tmp_path / "colors.tsv"
        color_file.write_text("A,B\tclade\t#ff0000\tAB\n")
        args = _make_args(
            tmp_path,
            circular=circular,
            color_file=str(color_file),
            legend_position="none",
            ylabel_fontsize=0,
            no_title=True,
        )
        cm = CharacterMap(args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        parent_map = build_parent_map(tree)
        node_labels = CharacterMap._assign_node_labels(tree)

        original_add_collection = Axes.add_collection
        line_collections = []

        def fail_plot(*args, **kwargs):
            raise AssertionError("clade overlay branches should use LineCollection")

        def capture_collection(self, collection, *args, **kwargs):
            if isinstance(collection, LineCollection):
                line_collections.append(collection)
            return original_add_collection(self, collection, *args, **kwargs)

        monkeypatch.setattr(Axes, "plot", fail_plot)
        monkeypatch.setattr(Axes, "add_collection", capture_collection)

        cm._plot_character_map(tree, {}, node_labels, parent_map)

        assert len(line_collections) >= 3
        out = Path(args.output)
        assert out.exists()
        assert out.stat().st_size > 0

    @pytest.mark.parametrize("circular", [False, True])
    def test_plot_batches_character_change_markers(
        self, monkeypatch, tmp_path, circular
    ):
        pytest.importorskip("matplotlib")
        from matplotlib.axes import Axes

        args = _make_args(
            tmp_path,
            circular=circular,
            legend_position="none",
            ylabel_fontsize=0,
            no_title=True,
            change_marker_size=123.0,
            change_fontsize=8.5,
        )
        cm = CharacterMap(args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        parent_map = build_parent_map(tree)
        node_labels = CharacterMap._assign_node_labels(tree)
        changed_clades = [
            clade for clade in CharacterMap._iter_preorder(tree.root)
            if id(clade) in parent_map
        ][:3]
        classified = {
            id(changed_clades[0]): [
                (0, "0", "1", "synapomorphy"),
                (1, "1", "0", "convergence"),
            ],
            id(changed_clades[1]): [
                (2, "0", "1", "reversal"),
            ],
            id(changed_clades[2]): [
                (3, "1", "0", "synapomorphy"),
            ],
        }

        original_scatter = Axes.scatter
        original_annotate = Axes.annotate
        scatter_calls = []
        annotation_fontsizes = []

        def capture_scatter(self, x, y, *args, **kwargs):
            scatter_calls.append((x, y, kwargs))
            return original_scatter(self, x, y, *args, **kwargs)

        def capture_annotate(self, *args, **kwargs):
            annotation_fontsizes.append(kwargs.get("fontsize"))
            return original_annotate(self, *args, **kwargs)

        monkeypatch.setattr(Axes, "scatter", capture_scatter)
        monkeypatch.setattr(Axes, "annotate", capture_annotate)

        cm._plot_character_map(tree, classified, node_labels, parent_map)

        assert len(scatter_calls) == 1
        x, y, kwargs = scatter_calls[0]
        assert len(x) == len(y) == 4
        assert len(kwargs["c"]) == 4
        assert kwargs["s"] == 123.0
        assert kwargs["edgecolors"] == "black"
        assert annotation_fontsizes == [8.5] * 8
        out = Path(args.output)
        assert out.exists()
        assert out.stat().st_size > 0


class TestCharacterMapJson:
    def test_run_uses_fast_tip_name_helper_for_setup(self, mocker, tmp_path):
        args = _make_args(tmp_path, json=True)
        cm = CharacterMap(args)
        tip_name_spy = mocker.spy(cm, "get_tip_names_from_tree")
        mocker.patch.object(CharacterMap, "_plot_character_map")
        mocker.patch("phykit.services.tree.character_map.print_json")

        cm.run()

        tip_name_spy.assert_called_once()

    def test_run_uses_unmodified_tree_read(self, mocker, tmp_path):
        args = _make_args(tmp_path, json=True)
        cm = CharacterMap(args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,C:1);"), "newick")
        tip_states = {"A": ["0"], "B": ["1"], "C": ["0"]}

        read_unmodified = mocker.patch.object(
            cm, "read_tree_file_unmodified", return_value=tree
        )
        mocker.patch.object(
            cm,
            "read_tree_file",
            side_effect=AssertionError("run should use read_tree_file_unmodified"),
        )
        mocker.patch.object(cm, "get_tip_names_from_tree", return_value=list(tip_states))
        mocker.patch.object(
            cm, "_parse_character_matrix", return_value=(["c1"], tip_states)
        )
        mocker.patch("phykit.services.tree.character_map.build_parent_map", return_value={})
        mocker.patch("phykit.services.tree.character_map.fitch_downpass", return_value=({}, [0]))
        mocker.patch(
            "phykit.services.tree.character_map.fitch_uppass_acctran",
            return_value={},
        )
        mocker.patch("phykit.services.tree.character_map.detect_changes", return_value={})
        mocker.patch("phykit.services.tree.character_map.classify_changes", return_value={})
        mocker.patch.object(
            cm,
            "_summarize_character_states_with_counts",
            return_value=([2], [], 1, [Counter({"0": 2, "1": 1})]),
        )
        mocker.patch("phykit.services.tree.character_map.consistency_index", return_value=([1.0], 1.0))
        mocker.patch.object(cm, "_retention_index_from_counts", return_value=([1.0], 1.0))
        mocker.patch.object(cm, "_assign_node_labels", return_value={})
        plot_character_map = mocker.patch.object(cm, "_plot_character_map")
        mocker.patch.object(cm, "_print_json")

        cm.run()

        read_unmodified.assert_called_once_with()
        assert plot_character_map.call_args.args[0] is not tree

    def test_run_skips_copy_for_clean_binary_tree(self, mocker, tmp_path):
        args = _make_args(tmp_path, json=True)
        cm = CharacterMap(args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,C:1):0;"), "newick")
        tip_states = {"A": ["0"], "B": ["1"], "C": ["0"]}

        mocker.patch.object(cm, "read_tree_file_unmodified", return_value=tree)
        mocker.patch.object(cm, "get_tip_names_from_tree", return_value=list(tip_states))
        mocker.patch.object(
            cm, "_parse_character_matrix", return_value=(["c1"], tip_states)
        )
        fast_copy = mocker.patch.object(
            cm,
            "_fast_copy",
            side_effect=AssertionError("clean all-shared analysis should not copy tree"),
        )
        mocker.patch(
            "phykit.services.tree.character_map.sankoff_downpass",
            side_effect=AssertionError("binary trees should stay on the Fitch path"),
        )
        mocker.patch("phykit.services.tree.character_map.build_parent_map", return_value={})
        mocker.patch("phykit.services.tree.character_map.fitch_downpass", return_value=({}, [0]))
        mocker.patch(
            "phykit.services.tree.character_map.fitch_uppass_acctran",
            return_value={},
        )
        mocker.patch("phykit.services.tree.character_map.detect_changes", return_value={})
        mocker.patch("phykit.services.tree.character_map.classify_changes", return_value={})
        mocker.patch.object(
            cm,
            "_summarize_character_states_with_counts",
            return_value=([2], [], 1, [Counter({"0": 2, "1": 1})]),
        )
        mocker.patch("phykit.services.tree.character_map.consistency_index", return_value=([1.0], 1.0))
        mocker.patch.object(cm, "_retention_index_from_counts", return_value=([1.0], 1.0))
        mocker.patch.object(cm, "_assign_node_labels", return_value={})
        plot_character_map = mocker.patch.object(cm, "_plot_character_map")
        mocker.patch.object(cm, "_print_json")

        cm.run()

        fast_copy.assert_not_called()
        assert plot_character_map.call_args.args[0] is tree

    def test_run_preserves_polytomy_and_uses_exact_sankoff_path(
        self,
        mocker,
        tmp_path,
    ):
        args = _make_args(tmp_path, json=True)
        cm = CharacterMap(args)
        tree = Phylo.read(StringIO("(A:1,B:1,C:1):0;"), "newick")
        tip_states = {"A": ["1"], "B": ["1"], "C": ["0"]}

        mocker.patch.object(cm, "read_tree_file_unmodified", return_value=tree)
        mocker.patch.object(
            cm,
            "get_tip_names_from_tree",
            return_value=list(tip_states),
        )
        mocker.patch.object(
            cm,
            "_parse_character_matrix",
            return_value=(["c1"], tip_states),
        )
        mocker.patch.object(
            cm,
            "_fast_copy",
            side_effect=AssertionError("a clean polytomy should not be copied"),
        )
        mocker.patch(
            "phykit.services.tree.character_map.fitch_downpass",
            side_effect=AssertionError("polytomies should use exact Sankoff costs"),
        )
        sankoff_spy = mocker.spy(character_map_module, "sankoff_downpass")
        plot_character_map = mocker.patch.object(cm, "_plot_character_map")
        mocked_json = mocker.patch(
            "phykit.services.tree.character_map.print_json"
        )

        cm.run()

        sankoff_spy.assert_called_once()
        plotted_tree = plot_character_map.call_args.args[0]
        assert plotted_tree is tree
        assert len(plotted_tree.root.clades) == 3
        payload = mocked_json.call_args.args[0]
        assert payload["tree_length"] == 1
        assert payload["ci"] == 1.0

    def test_run_copies_before_pruning_missing_tree_tips(self, mocker, tmp_path):
        args = _make_args(tmp_path, json=True, allow_taxon_mismatch=True)
        cm = CharacterMap(args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1):0;"), "newick")
        tree_copy = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1):0;"), "newick")
        pruned_tree = Phylo.read(StringIO("((A:1,B:1):1,C:1):0;"), "newick")
        tip_states = {"A": ["0"], "B": ["1"], "C": ["0"]}

        mocker.patch.object(cm, "read_tree_file_unmodified", return_value=tree)
        mocker.patch.object(
            cm, "get_tip_names_from_tree", return_value=["A", "B", "C", "D"]
        )
        mocker.patch.object(
            cm, "_parse_character_matrix", return_value=(["c1"], tip_states)
        )
        fast_copy = mocker.patch.object(cm, "_fast_copy", return_value=tree_copy)
        prune = mocker.patch.object(
            cm, "prune_tree_using_taxa_list", return_value=pruned_tree
        )
        mocker.patch("phykit.services.tree.character_map.build_parent_map", return_value={})
        mocker.patch("phykit.services.tree.character_map.fitch_downpass", return_value=({}, [0]))
        mocker.patch(
            "phykit.services.tree.character_map.fitch_uppass_acctran",
            return_value={},
        )
        mocker.patch("phykit.services.tree.character_map.detect_changes", return_value={})
        mocker.patch("phykit.services.tree.character_map.classify_changes", return_value={})
        mocker.patch.object(
            cm,
            "_summarize_character_states_with_counts",
            return_value=([2], [], 1, [Counter({"0": 2, "1": 1})]),
        )
        mocker.patch("phykit.services.tree.character_map.consistency_index", return_value=([1.0], 1.0))
        mocker.patch.object(cm, "_retention_index_from_counts", return_value=([1.0], 1.0))
        mocker.patch.object(cm, "_assign_node_labels", return_value={})
        plot_character_map = mocker.patch.object(cm, "_plot_character_map")
        mocker.patch.object(cm, "_print_json")

        cm.run()

        fast_copy.assert_called_once_with(tree)
        prune.assert_called_once_with(tree_copy, ["D"])
        assert plot_character_map.call_args.args[0] is pruned_tree

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

    def test_print_json_groups_changes_with_single_classified_scan(
        self,
        mocker,
        tmp_path,
    ):
        class CountingClassified(dict):
            def __init__(self, *args, **kwargs):
                super().__init__(*args, **kwargs)
                self.items_calls = 0

            def items(self):
                self.items_calls += 1
                return super().items()

        args = _make_args(tmp_path, json=True)
        cm = CharacterMap(args)
        classified = CountingClassified(
            {
                10: [(1, "0", "1", "convergence"), (0, "1", "0", "reversal")],
                20: [(1, "1", "2", "synapomorphy")],
            }
        )
        mocked_json = mocker.patch(
            "phykit.services.tree.character_map.print_json"
        )

        cm._print_json(
            ["char0", "char1", "char2"],
            3,
            2,
            4,
            0.75,
            0.5,
            [1.0, None, 0.5],
            [0.5, 0.25, None],
            [1, 2, 3],
            classified,
            {10: "node_a", 20: "node_b"},
        )

        assert classified.items_calls == 1
        characters = mocked_json.call_args.args[0]["characters"]
        assert characters[0]["changes"] == [
            {
                "branch": "node_a",
                "from": "1",
                "to": "0",
                "type": "reversal",
            }
        ]
        assert characters[1]["changes"] == [
            {
                "branch": "node_a",
                "from": "0",
                "to": "1",
                "type": "convergence",
            },
            {
                "branch": "node_b",
                "from": "1",
                "to": "2",
                "type": "synapomorphy",
            },
        ]
        assert characters[2]["changes"] == []

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

    def test_print_verbose_groups_changes_once_and_batches_output(self, tmp_path, mocker):
        args = _make_args(tmp_path, verbose=True)
        cm = CharacterMap(args)
        printed = mocker.patch("builtins.print")

        cm._print_verbose(
            ["c0", "c1", "c2"],
            3,
            2,
            4,
            0.5,
            None,
            [1.0, None, 0.25],
            [0.5, 0.75, None],
            [1, 2, 1],
            {
                10: [(1, "0", "1", "convergence"), (0, "A", "B", "synapomorphy")],
                20: [(1, "1", "0", "reversal")],
            },
            {10: "node_10"},
        )

        printed.assert_called_once_with(
            f"Optimization: {cm.optimization}\n"
            "Characters: 3\n"
            "Parsimony-informative: 2\n"
            "Tree length: 4\n"
            "CI: 0.5000\n"
            "RI: N/A\n"
            f"Output: {cm.output_path}\n"
            "\n"
            "Character 0 (c0): steps=1, CI=1.0000, RI=0.5000\n"
            "  node_10: A->B (synapomorphy)\n"
            "Character 1 (c1): steps=2, CI=N/A, RI=0.7500\n"
            "  node_10: 0->1 (convergence)\n"
            "  id_20: 1->0 (reversal)\n"
            "Character 2 (c2): steps=1, CI=0.2500, RI=N/A\n"
        )

    def test_print_summary_batches_output(self, tmp_path, mocker):
        args = _make_args(tmp_path)
        cm = CharacterMap(args)
        printed = mocker.patch("builtins.print")

        cm._print_summary(8, 5, 12, 0.625, None)

        printed.assert_called_once_with(
            f"Optimization: {cm.optimization}\n"
            "Characters: 8\n"
            "Parsimony-informative: 5\n"
            "Tree length: 12\n"
            "CI: 0.6250\n"
            "RI: N/A\n"
            f"Output: {cm.output_path}"
        )

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
