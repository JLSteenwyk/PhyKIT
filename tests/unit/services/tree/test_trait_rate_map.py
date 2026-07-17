import builtins
import json
import os
import subprocess
import sys
import tempfile

import pytest
from argparse import Namespace
from io import StringIO
from pathlib import Path

from Bio import Phylo
from Bio.Phylo.BaseTree import Clade, TreeMixin

from phykit.errors import PhykitUserError
from phykit.services.tree.trait_rate_map import TraitRateMap


def test_module_import_does_not_import_numpy():
    code = """
import sys
import phykit.services.tree.trait_rate_map as module
assert callable(module.print_json)
assert "typing" not in sys.modules
assert "pickle" not in sys.modules
assert "json" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
assert "phykit.helpers.plot_config" not in sys.modules
assert "phykit.helpers.circular_layout" not in sys.modules
assert "phykit.helpers.color_annotations" not in sys.modules
assert "numpy" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


here = Path(__file__)
SAMPLE_FILES = here.parent.parent.parent.parent / "sample_files"
TREE_SIMPLE = str(SAMPLE_FILES / "tree_simple.tre")
TRAITS_FILE = str(SAMPLE_FILES / "tree_simple_traits.tsv")
MULTI_TRAITS_FILE = str(SAMPLE_FILES / "tree_simple_multi_traits.tsv")


def _make_args(**kwargs):
    defaults = dict(
        tree=TREE_SIMPLE,
        trait_data=TRAITS_FILE,
        output="/tmp/test_rate_map.png",
        trait=None,
        json=False,
    )
    defaults.update(kwargs)
    return Namespace(**defaults)


class TestProcessArgs:
    def test_defaults(self):
        args = Namespace(
            tree="t.tre",
            trait_data="d.tsv",
            output="out.png",
        )
        svc = TraitRateMap.__new__(TraitRateMap)
        parsed = svc.process_args(args)
        assert parsed["tree_file_path"] == "t.tre"
        assert parsed["trait_data_path"] == "d.tsv"
        assert parsed["output_path"] == "out.png"
        assert parsed["trait_column"] is None
        assert parsed["json_output"] is False

    def test_with_trait_column(self):
        args = Namespace(
            tree="t.tre",
            trait_data="d.tsv",
            output="out.png",
            trait="body_mass",
            json=True,
        )
        svc = TraitRateMap.__new__(TraitRateMap)
        parsed = svc.process_args(args)
        assert parsed["trait_column"] == "body_mass"
        assert parsed["json_output"] is True


class TestSingleTraitParsing:
    def test_missing_file_reports_path(self, tmp_path):
        missing = tmp_path / "missing.tsv"
        svc = TraitRateMap.__new__(TraitRateMap)

        with pytest.raises(PhykitUserError) as exc_info:
            svc._parse_single_trait_data(str(missing), ["A", "B", "C"])

        assert str(missing) in exc_info.value.messages[0]

    def test_comments_and_blanks(self, tmp_path):
        trait_file = tmp_path / "traits.tsv"
        trait_file.write_text(
            "# comment\n\nraccoon\t1.0\nbear\t2.0\nsea_lion\t3.0\n"
            "seal\t4.0\nmonkey\t5.0\ncat\t6.0\nweasel\t7.0\ndog\t8.0\n"
        )
        svc = TraitRateMap(_make_args(trait_data=str(trait_file)))
        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)

        traits = svc._parse_single_trait_data(str(trait_file), tree_tips)

        assert len(traits) == 8
        assert traits["raccoon"] == pytest.approx(1.0)

    def test_whitespace_prefixed_comments(self, tmp_path):
        trait_file = tmp_path / "traits.tsv"
        trait_file.write_text(
            "   # comment\n\nraccoon\t1.0\nbear\t2.0\nsea_lion\t3.0\n"
            "seal\t4.0\nmonkey\t5.0\ncat\t6.0\nweasel\t7.0\ndog\t8.0\n"
        )
        svc = TraitRateMap(_make_args(trait_data=str(trait_file)))
        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)

        traits = svc._parse_single_trait_data(str(trait_file), tree_tips)

        assert len(traits) == 8
        assert traits["raccoon"] == pytest.approx(1.0)

    def test_all_shared_trait_file_emits_no_warnings(self, tmp_path, capsys):
        trait_file = tmp_path / "traits.tsv"
        trait_file.write_text(
            "raccoon\t1.0\nbear\t2.0\nsea_lion\t3.0\nseal\t4.0\n"
            "monkey\t5.0\ncat\t6.0\nweasel\t7.0\ndog\t8.0\n"
        )
        svc = TraitRateMap(_make_args(trait_data=str(trait_file)))
        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)

        traits = svc._parse_single_trait_data(str(trait_file), tree_tips)

        stderr = capsys.readouterr().err
        assert set(traits) == set(tree_tips)
        assert stderr == ""

    def test_ordered_all_shared_trait_file_skips_sets(self, tmp_path, monkeypatch):
        trait_file = tmp_path / "traits.tsv"
        trait_file.write_text("A\t1.0\nB\t2.0\nC\t3.0\n")
        svc = TraitRateMap(_make_args(trait_data=str(trait_file)))

        def fail_set(*args, **kwargs):
            raise AssertionError("ordered exact trait path should not build sets")

        monkeypatch.setattr(builtins, "set", fail_set)
        traits = svc._parse_single_trait_data(str(trait_file), ["A", "B", "C"])

        assert traits == {"A": 1.0, "B": 2.0, "C": 3.0}
        assert builtins.set is fail_set

    def test_extra_columns_error(self, tmp_path):
        trait_file = tmp_path / "traits.tsv"
        trait_file.write_text(
            "raccoon\t1.0\nbear\t2.0\nsea_lion\t3.0\ndog\t4.0\nextra\t1.0\t2.0\n"
        )
        svc = TraitRateMap(_make_args(trait_data=str(trait_file)))
        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)

        with pytest.raises(PhykitUserError) as exc_info:
            svc._parse_single_trait_data(str(trait_file), tree_tips)

        assert (
            "Line 5 in trait file has 3 columns; expected 2."
            in exc_info.value.messages
        )

    def test_non_numeric_trait_error(self, tmp_path):
        trait_file = tmp_path / "traits.tsv"
        trait_file.write_text("raccoon\t1.0\nbear\tbad\nsea_lion\t3.0\ndog\t4.0\n")
        svc = TraitRateMap(_make_args(trait_data=str(trait_file)))
        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)

        with pytest.raises(PhykitUserError) as exc_info:
            svc._parse_single_trait_data(str(trait_file), tree_tips)

        assert (
            "Non-numeric trait value 'bad' for taxon 'bear' on line 2."
            in exc_info.value.messages
        )


class TestTextOutput:
    def test_summarize_branch_rate_stats_single_pass_preserves_ties(self):
        first_min = {"rate": 1.0, "child": "first_min"}
        second_min = {"rate": 1.0, "child": "second_min"}
        first_max = {"rate": 5.0, "child": "first_max"}
        second_max = {"rate": 5.0, "child": "second_max"}
        branch_rates = [first_min, first_max, second_min, second_max]

        mean_rate, min_entry, max_entry = (
            TraitRateMap._summarize_branch_rate_stats(branch_rates)
        )

        assert mean_rate == pytest.approx(3.0)
        assert min_entry is first_min
        assert max_entry is first_max

    def test_summarize_branch_rate_stats_empty(self):
        assert TraitRateMap._summarize_branch_rate_stats([]) == (0.0, None, None)

    def test_print_text_output_batches_summary(self, mocker):
        svc = TraitRateMap.__new__(TraitRateMap)
        svc.output_path = "trait_rate.png"
        printed = mocker.patch("builtins.print")
        min_entry = {"rate": 0.012345, "child": "N1"}
        max_entry = {"rate": 9.876543, "child": "N2"}

        svc._print_text_output(12, "body_mass", 1.234567, min_entry, max_entry)

        printed.assert_called_once_with(
            "Trait Rate Map\n"
            "Taxa: 12\n"
            "Trait: body_mass\n"
            "Mean rate: 1.2346\n"
            "Min rate: 0.0123 (branch to N1)\n"
            "Max rate: 9.8765 (branch to N2)\n"
            "Output: trait_rate.png"
        )

    def test_print_text_output_handles_empty_rates(self, mocker):
        svc = TraitRateMap.__new__(TraitRateMap)
        svc.output_path = "empty.png"
        printed = mocker.patch("builtins.print")

        svc._print_text_output(3, "mass", 0.0, None, None)

        printed.assert_called_once_with(
            "Trait Rate Map\n"
            "Taxa: 3\n"
            "Trait: mass\n"
            "Mean rate: 0.0000\n"
            "Output: empty.png"
        )


class TestTraitDataParsing:
    @pytest.mark.parametrize(
        ("contents", "column", "message"),
        [
            ("# comment\n\n", "mass", "must have a header row"),
            (
                "taxon\tmass\nA\t1\nB\t2\nC\t3\n",
                "height",
                "Column 'height' not found",
            ),
            ("taxon\tmass\nA\n", "mass", "too few columns"),
            ("taxon\tmass\n# no data\n", "mass", "must have a header row"),
        ],
    )
    def test_multi_trait_parser_reports_structural_errors(
        self, tmp_path, contents, column, message
    ):
        trait_file = tmp_path / "traits.tsv"
        trait_file.write_text(contents)
        svc = TraitRateMap.__new__(TraitRateMap)

        with pytest.raises(PhykitUserError) as exc_info:
            svc._parse_multi_trait_data(
                str(trait_file), ["A", "B", "C"], column
            )

        assert message in exc_info.value.messages[0]

    def test_multi_trait_parser_reports_missing_path(self, tmp_path):
        missing = tmp_path / "missing.tsv"
        svc = TraitRateMap.__new__(TraitRateMap)

        with pytest.raises(PhykitUserError) as exc_info:
            svc._parse_multi_trait_data(
                str(missing), ["A", "B", "C"], "mass"
            )

        assert str(missing) in exc_info.value.messages[0]

    def test_intersect_traits_accepts_unordered_exact_taxa(self, capsys):
        svc = TraitRateMap.__new__(TraitRateMap)
        traits = {"C": 3.0, "A": 1.0, "B": 2.0}

        assert svc._intersect_traits(traits, ["A", "B", "C"]) is traits
        assert capsys.readouterr().err == ""

    def test_intersect_traits_warns_for_both_mismatch_directions(self, capsys):
        svc = TraitRateMap.__new__(TraitRateMap)

        shared = svc._intersect_traits(
            {"A": 1.0, "B": 2.0, "C": 3.0, "extra": 4.0},
            ["A", "B", "C", "missing"],
        )

        assert shared == {"A": 1.0, "B": 2.0, "C": 3.0}
        stderr = capsys.readouterr().err
        assert "missing" in stderr
        assert "extra" in stderr

    def test_intersect_traits_rejects_too_few_shared_taxa(self, capsys):
        svc = TraitRateMap.__new__(TraitRateMap)

        with pytest.raises(PhykitUserError) as exc_info:
            svc._intersect_traits(
                {"A": 1.0, "B": 2.0, "extra": 3.0},
                ["A", "B", "missing"],
            )

        assert "Only 2 shared taxa" in exc_info.value.messages[0]
        stderr = capsys.readouterr().err
        assert "missing" in stderr
        assert "extra" in stderr

    def test_multi_trait_parser_streams_without_readlines(self, monkeypatch):
        class StreamingOnlyFile(StringIO):
            def __enter__(self):
                return self

            def __exit__(self, *_args):
                self.close()

            def readlines(self):
                raise AssertionError("multi-trait parser should stream rows")

        data = (
            "# ignored before header\n"
            "\n"
            "taxon\tbody_mass\twing\n"
            "A\t1.0\t10.0\textra\n"
            "# ignored between rows\n"
            "B\t2.0\t20.0\n"
            "C\t3.0\t30.0\n"
        )

        def fake_open(path):
            assert path == "traits.tsv"
            return StreamingOnlyFile(data)

        monkeypatch.setattr("builtins.open", fake_open)
        svc = TraitRateMap.__new__(TraitRateMap)

        traits = svc._parse_multi_trait_data(
            "traits.tsv", ["A", "B", "C"], "body_mass"
        )

        assert traits == {"A": 1.0, "B": 2.0, "C": 3.0}

    def test_multi_trait_parser_skips_whitespace_prefixed_comments(self, tmp_path):
        trait_file = tmp_path / "traits.tsv"
        trait_file.write_text(
            "   # ignored before header\n"
            "taxon\tbody_mass\twing\n"
            "A\t1.0\t10.0\n"
            "  # ignored between rows\n"
            "B\t2.0\t20.0\n"
            "C\t3.0\t30.0\n"
        )
        svc = TraitRateMap.__new__(TraitRateMap)

        traits = svc._parse_multi_trait_data(
            str(trait_file), ["A", "B", "C"], "body_mass"
        )

        assert traits == {"A": 1.0, "B": 2.0, "C": 3.0}

    def test_multi_trait_parser_preserves_logical_line_numbers(self, tmp_path):
        trait_file = tmp_path / "traits.tsv"
        trait_file.write_text(
            "\n".join(
                [
                    "# ignored before header",
                    "taxon\tbody_mass\twing",
                    "A\t1.0\t10.0",
                    "# ignored between rows",
                    "B\tbad\t20.0",
                ]
            )
            + "\n"
        )
        svc = TraitRateMap.__new__(TraitRateMap)

        with pytest.raises(PhykitUserError) as excinfo:
            svc._parse_multi_trait_data(
                str(trait_file), ["A", "B", "C"], "body_mass"
            )

        assert (
            "Non-numeric trait value 'bad' for taxon 'B' on line 3."
            in excinfo.value.messages
        )

    def test_multi_trait_parser_ignores_unused_trailing_columns(self, tmp_path):
        trait_file = tmp_path / "traits.tsv"
        trait_file.write_text(
            "taxon\tbody_mass\twing\tunused\n"
            "A\t1.0\t10.0\tnot_numeric\n"
            "B\t2.0\t20.0\talso_bad\n"
            "C\t3.0\t30.0\tignored\n"
        )
        svc = TraitRateMap.__new__(TraitRateMap)

        traits = svc._parse_multi_trait_data(
            str(trait_file), ["A", "B", "C"], "body_mass"
        )

        assert traits == {"A": 1.0, "B": 2.0, "C": 3.0}


class TestAncestralReconstruction:
    def test_node_values_computed(self):
        try:
            import matplotlib
        except ImportError:
            pytest.skip("matplotlib not installed")

        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
            tmppath = f.name

        try:
            args = _make_args(output=tmppath)
            svc = TraitRateMap(args)
            tree = svc.read_tree_file()
            tree_tips = svc.get_tip_names_from_tree(tree)
            trait_values = svc._parse_single_trait_data(TRAITS_FILE, tree_tips)
            node_values = svc._ancestral_reconstruction(tree, trait_values)

            # Should have values for all nodes (tips + internal)
            all_nodes = list(tree.find_clades())
            for node in all_nodes:
                assert id(node) in node_values, f"Missing value for node {node.name}"
        finally:
            if os.path.exists(tmppath):
                os.unlink(tmppath)

    def test_missing_tip_values_are_ignored(self):
        svc = TraitRateMap.__new__(TraitRateMap)
        tree = Phylo.read(StringIO("(A:1,B:1);"), "newick")

        one_value = svc._ancestral_reconstruction(tree, {"A": 4.0})
        no_values = svc._ancestral_reconstruction(tree, {})

        assert one_value[id(tree.root)] == pytest.approx(4.0)
        assert id(tree.root) not in no_values


class TestBranchRates:
    def test_internal_labels_preserve_existing_names(self):
        svc = TraitRateMap.__new__(TraitRateMap)
        tree = Phylo.read(StringIO("(A:1,B:1)ancestor;"), "newick")

        labels = svc._label_internal_nodes(tree)

        assert labels[id(tree.root)] == "ancestor"

    def test_rate_pipeline_uses_direct_tree_traversal(self, monkeypatch):
        svc = TraitRateMap.__new__(TraitRateMap)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        trait_values = {"A": 1.0, "B": 2.0, "C": 3.0, "D": 4.0}

        def fail_traversal(*args, **kwargs):
            raise AssertionError("optimized path should use direct traversal")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_traversal)
        monkeypatch.setattr(TreeMixin, "get_terminals", fail_traversal)

        node_labels = svc._label_internal_nodes(tree)
        node_values = svc._ancestral_reconstruction(tree, trait_values)
        parent_map = svc._build_parent_map(tree)
        branch_rates = svc._compute_branch_rates(
            tree, node_values, parent_map, node_labels
        )

        assert set(node_labels.values()) == {"N1", "N2", "N3"}
        assert node_values[id(tree.root)] == pytest.approx(2.5)
        assert len(branch_rates) == 6
        assert all(entry["rate"] >= 0 for entry in branch_rates)

    def test_parent_map_handles_mixed_child_counts_without_preorder(
        self, monkeypatch
    ):
        svc = TraitRateMap.__new__(TraitRateMap)
        tree = Phylo.read(
            StringIO("(A:1,(B:1,C:1):1,(D:1,E:1,F:1):1);"),
            "newick",
        )

        def fail_iter_preorder(*args, **kwargs):
            raise AssertionError("parent map should build directly")

        monkeypatch.setattr(svc, "_iter_preorder", fail_iter_preorder)

        parent_map = svc._build_parent_map(tree)

        terminal, binary, trifurcating = tree.root.clades
        assert id(tree.root) not in parent_map
        assert parent_map[id(terminal)] is tree.root
        assert parent_map[id(binary)] is tree.root
        assert parent_map[id(trifurcating)] is tree.root
        assert all(parent_map[id(child)] is binary for child in binary.clades)
        assert all(
            parent_map[id(child)] is trifurcating for child in trifurcating.clades
        )

    def test_iter_preorder_preserves_order_without_reversed(self):
        class NoReversedList(list):
            def __reversed__(self):
                raise AssertionError("_iter_preorder should push children directly")

        svc = TraitRateMap.__new__(TraitRateMap)
        root = Clade(name="root")
        left = Clade(name="left")
        middle = Clade(name="middle")
        right = Clade(name="right")
        left.clades = NoReversedList([Clade(name="left_a"), Clade(name="left_b")])
        middle.clades = NoReversedList([Clade(name="middle_a")])
        right.clades = NoReversedList()
        root.clades = NoReversedList([left, middle, right])

        order = [clade.name for clade in svc._iter_preorder(root)]

        assert order == [
            "root",
            "left",
            "left_a",
            "left_b",
            "middle",
            "middle_a",
            "right",
        ]

    def test_iter_postorder_preserves_biophylo_order(self, monkeypatch):
        tree = Phylo.read(
            StringIO("((A:1,B:1):1,(C:1,(D:1,E:1):1):1);"),
            "newick",
        )
        expected_ids = [
            id(clade) for clade in tree.find_clades(order="postorder")
        ]

        def fail_find_clades(*_args, **_kwargs):
            raise AssertionError("direct postorder helper should not use find_clades")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_find_clades)

        observed_ids = [
            id(clade) for clade in TraitRateMap._iter_postorder(tree.root)
        ]

        assert observed_ids == expected_ids

    def test_branch_rates_label_children_from_direct_child_lists(self, monkeypatch):
        svc = TraitRateMap.__new__(TraitRateMap)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        trait_values = {"A": 1.0, "B": 2.0, "C": 3.0, "D": 4.0}
        preorder_clades = list(svc._iter_preorder(tree.root))
        postorder_clades = list(reversed(preorder_clades))
        node_labels = svc._label_internal_nodes(tree, preorder_clades)
        node_values = svc._ancestral_reconstruction(
            tree, trait_values, postorder_clades
        )
        parent_map = svc._build_parent_map(tree, preorder_clades)

        def fail_is_terminal(self):
            raise AssertionError("branch-rate labels should use child lists")

        monkeypatch.setattr(Clade, "is_terminal", fail_is_terminal)

        branch_rates = svc._compute_branch_rates(
            tree, node_values, parent_map, node_labels, preorder_clades
        )

        assert len(branch_rates) == 6
        assert {entry["child"] for entry in branch_rates} >= {"A", "B", "C", "D"}

    def test_branch_rate_boundaries_and_missing_values(self):
        svc = TraitRateMap.__new__(TraitRateMap)
        valid_tip = Clade(name="A", branch_length=0.0)
        internal = Clade(
            name=None,
            branch_length=-1.0,
            clades=[Clade(name="descendant")],
        )
        root = Clade(name=None, clades=[valid_tip, internal])
        tree = Namespace(root=root)
        orphan = Clade(name="orphan", branch_length=1.0)
        no_parent_value = Clade(name="no-parent-value", branch_length=1.0)
        no_child_value = Clade(name="no-child-value", branch_length=1.0)
        missing_parent = Clade(name="missing-parent")
        parent_map = {
            id(valid_tip): root,
            id(internal): root,
            id(no_parent_value): missing_parent,
            id(no_child_value): root,
        }
        node_values = {
            id(root): 1.0,
            id(valid_tip): 2.0,
            id(internal): 3.0,
            id(no_parent_value): 4.0,
        }

        rates = svc._compute_branch_rates(
            tree,
            node_values,
            parent_map,
            {},
            [
                root,
                orphan,
                no_parent_value,
                no_child_value,
                valid_tip,
                internal,
            ],
        )

        assert [entry["child"] for entry in rates] == ["A", "unknown"]
        assert all(entry["parent"] == "root" for entry in rates)
        assert all(entry["branch_length"] > 0 for entry in rates)

    def test_ancestral_reconstruction_weighted_polytomy(self):
        svc = TraitRateMap.__new__(TraitRateMap)
        tree = Phylo.read(StringIO("(A:1,B:2,C:4);"), "newick")
        trait_values = {"A": 1.0, "B": 5.0, "C": 9.0}

        node_values = svc._ancestral_reconstruction(tree, trait_values)

        expected = (
            1.0 / 1.0 * 1.0
            + 1.0 / 2.0 * 5.0
            + 1.0 / 4.0 * 9.0
        ) / (1.0 + 1.0 / 2.0 + 1.0 / 4.0)
        assert node_values[id(tree.root)] == pytest.approx(expected)

    def test_rates_nonnegative(self):
        try:
            import matplotlib
        except ImportError:
            pytest.skip("matplotlib not installed")

        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
            tmppath = f.name

        try:
            args = _make_args(output=tmppath)
            svc = TraitRateMap(args)
            tree = svc.read_tree_file()
            tree_tips = svc.get_tip_names_from_tree(tree)
            trait_values = svc._parse_single_trait_data(TRAITS_FILE, tree_tips)
            node_values = svc._ancestral_reconstruction(tree, trait_values)
            parent_map = svc._build_parent_map(tree)
            node_labels = svc._label_internal_nodes(tree)
            branch_rates = svc._compute_branch_rates(tree, node_values, parent_map, node_labels)

            for entry in branch_rates:
                assert entry["rate"] >= 0, f"Negative rate for branch {entry['parent']} -> {entry['child']}"
        finally:
            if os.path.exists(tmppath):
                os.unlink(tmppath)

    def test_rate_zero_for_no_change(self):
        """If parent and child have the same value, rate should be 0."""
        try:
            import matplotlib
        except ImportError:
            pytest.skip("matplotlib not installed")

        # Create a trait file where all taxa have the same value
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
            equal_trait_path = f.name
            f.write("raccoon\t1.0\n")
            f.write("bear\t1.0\n")
            f.write("sea_lion\t1.0\n")
            f.write("seal\t1.0\n")
            f.write("monkey\t1.0\n")
            f.write("cat\t1.0\n")
            f.write("weasel\t1.0\n")
            f.write("dog\t1.0\n")

        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
            tmppath = f.name

        try:
            args = _make_args(trait_data=equal_trait_path, output=tmppath)
            svc = TraitRateMap(args)
            tree = svc.read_tree_file()
            tree_tips = svc.get_tip_names_from_tree(tree)
            trait_values = svc._parse_single_trait_data(equal_trait_path, tree_tips)
            node_values = svc._ancestral_reconstruction(tree, trait_values)
            parent_map = svc._build_parent_map(tree)
            node_labels = svc._label_internal_nodes(tree)
            branch_rates = svc._compute_branch_rates(tree, node_values, parent_map, node_labels)

            for entry in branch_rates:
                assert entry["rate"] == pytest.approx(0.0, abs=1e-10), \
                    f"Rate should be 0 for branch {entry['parent']} -> {entry['child']}, got {entry['rate']}"
        finally:
            if os.path.exists(equal_trait_path):
                os.unlink(equal_trait_path)
            if os.path.exists(tmppath):
                os.unlink(tmppath)


class TestPlotRateMap:
    @pytest.mark.parametrize("circular", [False, True])
    def test_plot_rate_map_uses_direct_tree_traversal(self, monkeypatch, circular):
        try:
            import matplotlib
        except ImportError:
            pytest.skip("matplotlib not installed")

        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
            tmppath = f.name

        try:
            args = _make_args(output=tmppath, circular=circular)
            svc = TraitRateMap(args)
            tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
            trait_values = {"A": 1.0, "B": 2.0, "C": 3.0, "D": 4.0}
            preorder_clades = list(svc._iter_preorder(tree.root))
            postorder_clades = list(reversed(preorder_clades))
            node_labels = svc._label_internal_nodes(tree, preorder_clades)
            node_values = svc._ancestral_reconstruction(
                tree, trait_values, postorder_clades
            )
            parent_map = svc._build_parent_map(tree, preorder_clades)
            branch_rates = svc._compute_branch_rates(
                tree, node_values, parent_map, node_labels, preorder_clades
            )

            def fail_traversal(*args, **kwargs):
                raise AssertionError("plot setup should use direct traversal")

            monkeypatch.setattr(TreeMixin, "find_clades", fail_traversal)
            monkeypatch.setattr(TreeMixin, "get_terminals", fail_traversal)

            svc._plot_rate_map(
                tree, node_values, branch_rates, parent_map,
                node_labels, trait_values, "trait", tmppath,
            )

            assert os.path.exists(tmppath)
            assert os.path.getsize(tmppath) > 0
        finally:
            if os.path.exists(tmppath):
                os.unlink(tmppath)

    def test_prepare_rate_map_plot_data_single_pass(self):
        args = _make_args()
        svc = TraitRateMap(args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        preorder_clades = list(svc._iter_preorder(tree.root))
        branch_rates = [
            {"clade_id": id(clade), "rate": float(index)}
            for index, clade in enumerate(preorder_clades)
            if clade is not tree.root
        ]

        rate_by_clade, preorder, tips, all_rates = (
            svc._prepare_rate_map_plot_data(tree, branch_rates)
        )

        assert [id(clade) for clade in preorder] == [
            id(clade) for clade in preorder_clades
        ]
        assert [tip.name for tip in tips] == ["A", "B", "C", "D"]
        assert all_rates == [entry["rate"] for entry in branch_rates]
        assert rate_by_clade == {
            entry["clade_id"]: entry["rate"] for entry in branch_rates
        }

    def test_plot_rate_map_reuses_clade_lists_for_circular_coords(
        self, monkeypatch, tmp_path
    ):
        pytest.importorskip("matplotlib")
        import phykit.helpers.circular_layout as circular_layout

        args = _make_args(
            output=str(tmp_path / "trait_rate_circular_precomputed_coords.png"),
            circular=True,
            ylabel_fontsize=0,
            no_title=True,
        )
        svc = TraitRateMap(args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        trait_values = {"A": 1.0, "B": 2.0, "C": 3.0, "D": 4.0}
        preorder_clades = list(svc._iter_preorder(tree.root))
        postorder_clades = list(reversed(preorder_clades))
        node_labels = svc._label_internal_nodes(tree, preorder_clades)
        node_values = svc._ancestral_reconstruction(
            tree, trait_values, postorder_clades
        )
        parent_map = svc._build_parent_map(tree, preorder_clades)
        branch_rates = svc._compute_branch_rates(
            tree, node_values, parent_map, node_labels, preorder_clades
        )
        _, expected_preorder, expected_tips, _ = svc._prepare_rate_map_plot_data(
            tree, branch_rates
        )
        expected_preorder_ids = [id(clade) for clade in expected_preorder]
        expected_terminal_ids = [id(tip) for tip in expected_tips]
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

        svc._plot_rate_map(
            tree, node_values, branch_rates, parent_map,
            node_labels, trait_values, "trait", args.output,
        )

        assert calls == [(expected_preorder_ids, expected_terminal_ids)]
        assert os.path.exists(args.output)

    def test_rectangular_plot_batches_rate_branches(self, monkeypatch, tmp_path):
        pytest.importorskip("matplotlib")
        import matplotlib.axes
        from matplotlib.collections import LineCollection

        args = _make_args(output=str(tmp_path / "trait_rate_batched.png"))
        svc = TraitRateMap(args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        trait_values = {"A": 1.0, "B": 2.0, "C": 3.0, "D": 4.0}
        preorder_clades = list(svc._iter_preorder(tree.root))
        postorder_clades = list(reversed(preorder_clades))
        node_labels = svc._label_internal_nodes(tree, preorder_clades)
        node_values = svc._ancestral_reconstruction(
            tree, trait_values, postorder_clades
        )
        parent_map = svc._build_parent_map(tree, preorder_clades)
        branch_rates = svc._compute_branch_rates(
            tree, node_values, parent_map, node_labels, preorder_clades
        )
        original_add_collection = matplotlib.axes.Axes.add_collection
        line_collections = []

        def fail_plot(*args, **kwargs):
            raise AssertionError("rectangular rate-map branches should be batched")

        def count_collection(self, collection, *args, **kwargs):
            if isinstance(collection, LineCollection):
                line_collections.append(collection)
            return original_add_collection(self, collection, *args, **kwargs)

        monkeypatch.setattr(matplotlib.axes.Axes, "plot", fail_plot)
        monkeypatch.setattr(
            matplotlib.axes.Axes, "add_collection", count_collection
        )

        svc._plot_rate_map(
            tree, node_values, branch_rates, parent_map,
            node_labels, trait_values, "trait", args.output,
        )

        rate_collections = [
            collection
            for collection in line_collections
            if collection.get_array() is not None
        ]
        assert len(rate_collections) == 1
        assert len(rate_collections[0].get_array()) == len(branch_rates)
        assert os.path.exists(args.output)

    def test_rectangular_plot_uses_scalar_rate_collection(self, monkeypatch, tmp_path):
        pytest.importorskip("matplotlib")
        import matplotlib.axes
        from matplotlib.collections import LineCollection

        args = _make_args(
            output=str(tmp_path / "trait_rate_scalar_collection.png"),
            ylabel_fontsize=0,
            no_title=True,
        )
        svc = TraitRateMap(args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        preorder_clades = list(svc._iter_preorder(tree.root))
        parent_map = svc._build_parent_map(tree, preorder_clades)
        node_values = {id(clade): 1.0 for clade in preorder_clades}
        node_labels = {id(clade): clade.name or "internal" for clade in preorder_clades}
        trait_values = {"A": 1.0, "B": 1.0, "C": 1.0, "D": 1.0}
        non_root_clades = [clade for clade in preorder_clades if clade is not tree.root]
        branch_rates = [
            {
                "clade_id": id(clade),
                "parent": "parent",
                "child": clade.name or "internal",
                "rate": 0.75,
                "change": 0.75,
                "branch_length": 1.0,
            }
            for clade in non_root_clades[1:]
        ]
        original_add_collection = matplotlib.axes.Axes.add_collection
        line_collections = []

        def count_collection(self, collection, *args, **kwargs):
            if isinstance(collection, LineCollection):
                line_collections.append(collection)
            return original_add_collection(self, collection, *args, **kwargs)

        monkeypatch.setattr(
            matplotlib.axes.Axes, "add_collection", count_collection
        )

        svc._plot_rate_map(
            tree, node_values, branch_rates, parent_map,
            node_labels, trait_values, "trait", args.output,
        )

        rate_collections = [
            collection
            for collection in line_collections
            if collection.get_array() is not None
        ]
        assert len(rate_collections) == 1
        assert list(rate_collections[0].get_array()) == [0.75] * len(branch_rates)
        assert any(
            collection.get_array() is None and len(collection.get_segments()) == 1
            for collection in line_collections
        )
        assert os.path.exists(args.output)

    def test_circular_plot_batches_rate_branches(self, monkeypatch, tmp_path):
        pytest.importorskip("matplotlib")
        import matplotlib.axes
        from matplotlib.collections import LineCollection

        args = _make_args(
            output=str(tmp_path / "trait_rate_circular_batched.png"),
            circular=True,
            ylabel_fontsize=0,
            no_title=True,
        )
        svc = TraitRateMap(args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        trait_values = {"A": 1.0, "B": 2.0, "C": 3.0, "D": 4.0}
        preorder_clades = list(svc._iter_preorder(tree.root))
        postorder_clades = list(reversed(preorder_clades))
        node_labels = svc._label_internal_nodes(tree, preorder_clades)
        node_values = svc._ancestral_reconstruction(
            tree, trait_values, postorder_clades
        )
        parent_map = svc._build_parent_map(tree, preorder_clades)
        branch_rates = svc._compute_branch_rates(
            tree, node_values, parent_map, node_labels, preorder_clades
        )
        original_add_collection = matplotlib.axes.Axes.add_collection
        line_collections = []

        def fail_plot(*args, **kwargs):
            raise AssertionError("circular rate-map branches should be batched")

        def count_collection(self, collection, *args, **kwargs):
            if isinstance(collection, LineCollection):
                line_collections.append(collection)
            return original_add_collection(self, collection, *args, **kwargs)

        monkeypatch.setattr(matplotlib.axes.Axes, "plot", fail_plot)
        monkeypatch.setattr(
            matplotlib.axes.Axes, "add_collection", count_collection
        )

        svc._plot_rate_map(
            tree, node_values, branch_rates, parent_map,
            node_labels, trait_values, "trait", args.output,
        )

        rate_collections = [
            collection
            for collection in line_collections
            if collection.get_array() is not None
        ]
        assert len(rate_collections) == 1
        assert len(rate_collections[0].get_array()) == len(branch_rates)
        assert len(line_collections) >= 2
        assert os.path.exists(args.output)

    def test_plot_rate_map_skips_redundant_tight_layout(self, monkeypatch, tmp_path):
        pytest.importorskip("matplotlib")
        from matplotlib.figure import Figure

        args = _make_args(
            output=str(tmp_path / "trait_rate_no_tight_layout.png"),
            ylabel_fontsize=0,
            no_title=True,
        )
        svc = TraitRateMap(args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        trait_values = {"A": 1.0, "B": 2.0, "C": 3.0, "D": 4.0}
        preorder_clades = list(svc._iter_preorder(tree.root))
        postorder_clades = list(reversed(preorder_clades))
        node_labels = svc._label_internal_nodes(tree, preorder_clades)
        node_values = svc._ancestral_reconstruction(
            tree, trait_values, postorder_clades
        )
        parent_map = svc._build_parent_map(tree, preorder_clades)
        branch_rates = svc._compute_branch_rates(
            tree, node_values, parent_map, node_labels, preorder_clades
        )

        def fail_tight_layout(self, *args, **kwargs):
            raise AssertionError("bbox_inches='tight' handles saved bounds")

        monkeypatch.setattr(Figure, "tight_layout", fail_tight_layout)
        svc._plot_rate_map(
            tree, node_values, branch_rates, parent_map,
            node_labels, trait_values, "trait", args.output,
        )

        assert os.path.exists(args.output)
        assert os.path.getsize(args.output) > 0


class TestRun:
    def test_creates_png(self):
        try:
            import matplotlib
        except ImportError:
            pytest.skip("matplotlib not installed")

        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
            tmppath = f.name

        try:
            args = _make_args(output=tmppath)
            svc = TraitRateMap(args)
            svc.run()
            assert os.path.exists(tmppath)
            assert os.path.getsize(tmppath) > 0
        finally:
            if os.path.exists(tmppath):
                os.unlink(tmppath)

    def test_cladogram_mode(self):
        try:
            import matplotlib
        except ImportError:
            pytest.skip("matplotlib not installed")

        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
            tmppath = f.name

        try:
            args = _make_args(output=tmppath, cladogram=True)
            svc = TraitRateMap(args)
            svc.run()
            assert os.path.exists(tmppath)
            assert os.path.getsize(tmppath) > 0
        finally:
            if os.path.exists(tmppath):
                os.unlink(tmppath)

    def test_circular_mode(self):
        try:
            import matplotlib
        except ImportError:
            pytest.skip("matplotlib not installed")

        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
            tmppath = f.name

        try:
            args = _make_args(output=tmppath, circular=True)
            svc = TraitRateMap(args)
            svc.run()
            assert os.path.exists(tmppath)
            assert os.path.getsize(tmppath) > 0
        finally:
            if os.path.exists(tmppath):
                os.unlink(tmppath)

    def test_json_output(self, capsys):
        try:
            import matplotlib
        except ImportError:
            pytest.skip("matplotlib not installed")

        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
            tmppath = f.name

        try:
            args = _make_args(output=tmppath, json=True)
            svc = TraitRateMap(args)
            svc.run()
            captured = capsys.readouterr()
            lines = captured.out.strip().split("\n")
            payload = json.loads(lines[-1])
            assert "n_taxa" in payload
            assert payload["n_taxa"] == 8
            assert "mean_rate" in payload
            assert "min_rate" in payload
            assert "max_rate" in payload
            assert "branches" in payload
            assert len(payload["branches"]) > 0
            assert "output_file" in payload
        finally:
            if os.path.exists(tmppath):
                os.unlink(tmppath)

    def test_run_reuses_initial_tip_names_for_prune_setup(self, mocker):
        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
            tmppath = f.name

        try:
            args = _make_args(output=tmppath, json=True)
            svc = TraitRateMap(args)
            spy = mocker.spy(svc, "get_tip_names_from_tree")
            mocker.patch.object(TraitRateMap, "_plot_rate_map")
            mocker.patch("phykit.services.tree.trait_rate_map.print_json")

            svc.run()

            assert spy.call_count == 1
        finally:
            if os.path.exists(tmppath):
                os.unlink(tmppath)

    def test_run_all_tips_present_uses_read_only_tree_without_copy_or_prune(
        self, mocker
    ):
        args = _make_args()
        svc = TraitRateMap(args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        trait_values = {"A": 1.0, "B": 2.0, "C": 3.0, "D": 4.0}
        branch_rates = [
            {
                "rate": 1.0,
                "child": "A",
                "parent": "N1",
                "change": 1.0,
                "branch_length": 1.0,
            }
        ]

        read_unmodified = mocker.patch.object(
            svc, "read_tree_file_unmodified", return_value=tree
        )
        mocker.patch.object(
            svc,
            "read_tree_file",
            side_effect=AssertionError("run should use read_tree_file_unmodified"),
        )
        validate = mocker.patch.object(svc, "validate_tree")
        mocker.patch.object(
            svc,
            "get_tip_names_from_tree",
            return_value=list(trait_values),
        )
        mocker.patch.object(svc, "_parse_single_trait_data", return_value=trait_values)
        mocker.patch.object(
            svc,
            "_fast_copy",
            side_effect=AssertionError("all-shared tree should not be copied"),
        )
        mocker.patch.object(
            svc,
            "prune_tree_using_taxa_list",
            side_effect=AssertionError("all-shared tree should not be pruned"),
        )
        mocker.patch.object(svc, "_label_internal_nodes", return_value={})
        mocker.patch.object(svc, "_ancestral_reconstruction", return_value={})
        mocker.patch.object(svc, "_build_parent_map", return_value={})
        mocker.patch.object(svc, "_compute_branch_rates", return_value=branch_rates)
        plot_rate_map = mocker.patch.object(svc, "_plot_rate_map")
        mocker.patch("builtins.print")

        svc.run()

        read_unmodified.assert_called_once_with()
        validate.assert_called_once()
        assert plot_rate_map.call_args.args[0] is tree

    def test_run_counts_traits_without_sorting_trait_keys(self, mocker):
        class UnsortableName(str):
            def __lt__(self, other):
                raise AssertionError("trait keys should not be sorted for count")

        args = _make_args()
        svc = TraitRateMap(args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        trait_values = {
            UnsortableName("D"): 4.0,
            UnsortableName("C"): 3.0,
            UnsortableName("B"): 2.0,
            UnsortableName("A"): 1.0,
        }
        branch_rates = [
            {
                "rate": 1.0,
                "child": "A",
                "parent": "N1",
                "change": 1.0,
                "branch_length": 1.0,
            }
        ]

        mocker.patch.object(svc, "read_tree_file_unmodified", return_value=tree)
        mocker.patch.object(svc, "validate_tree")
        mocker.patch.object(svc, "get_tip_names_from_tree", return_value=list("ABCD"))
        mocker.patch.object(svc, "_parse_single_trait_data", return_value=trait_values)
        mocker.patch.object(svc, "_label_internal_nodes", return_value={})
        mocker.patch.object(svc, "_ancestral_reconstruction", return_value={})
        mocker.patch.object(svc, "_build_parent_map", return_value={})
        mocker.patch.object(svc, "_compute_branch_rates", return_value=branch_rates)
        mocker.patch.object(svc, "_plot_rate_map")
        print_text = mocker.patch.object(svc, "_print_text_output")

        svc.run()

        assert print_text.call_args.args[0] == 4

    def test_run_missing_trait_taxa_copies_before_pruning(self, monkeypatch):
        args = _make_args()
        svc = TraitRateMap(args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        original_fast_copy = svc._fast_copy
        copied_trees = []
        captured = {}
        branch_rates = [
            {
                "rate": 1.0,
                "child": "A",
                "parent": "N1",
                "change": 1.0,
                "branch_length": 1.0,
            }
        ]

        def copy_spy(tree_to_copy):
            copied_tree = original_fast_copy(tree_to_copy)
            copied_trees.append(copied_tree)
            return copied_tree

        def plot_rate_map(tree_arg, *_args):
            captured["tree"] = tree_arg

        monkeypatch.setattr(svc, "read_tree_file_unmodified", lambda: tree)
        monkeypatch.setattr(svc, "_fast_copy", copy_spy)
        monkeypatch.setattr(
            svc,
            "_parse_single_trait_data",
            lambda *_args: {"A": 1.0, "B": 2.0, "C": 3.0},
        )
        monkeypatch.setattr(svc, "_label_internal_nodes", lambda *_args: {})
        monkeypatch.setattr(svc, "_ancestral_reconstruction", lambda *_args: {})
        monkeypatch.setattr(svc, "_build_parent_map", lambda *_args: {})
        monkeypatch.setattr(
            svc,
            "_compute_branch_rates",
            lambda *_args: branch_rates,
        )
        monkeypatch.setattr(svc, "_plot_rate_map", plot_rate_map)
        monkeypatch.setattr("builtins.print", lambda *_args, **_kwargs: None)

        svc.run()

        assert len(copied_trees) == 1
        assert captured["tree"] is copied_trees[0]
        assert {tip.name for tip in tree.get_terminals()} == {
            "A",
            "B",
            "C",
            "D",
        }
        assert {tip.name for tip in copied_trees[0].get_terminals()} == {
            "A",
            "B",
            "C",
        }

    def test_run_ladderize_copies_before_mutating_tree(self, monkeypatch):
        args = _make_args(ladderize=True)
        svc = TraitRateMap(args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        original_fast_copy = svc._fast_copy
        copied_trees = []
        captured = {}
        branch_rates = [
            {
                "rate": 1.0,
                "child": "A",
                "parent": "N1",
                "change": 1.0,
                "branch_length": 1.0,
            }
        ]

        def copy_spy(tree_to_copy):
            copied_tree = original_fast_copy(tree_to_copy)
            copied_trees.append(copied_tree)
            return copied_tree

        def plot_rate_map(tree_arg, *_args):
            captured["tree"] = tree_arg

        monkeypatch.setattr(svc, "read_tree_file_unmodified", lambda: tree)
        monkeypatch.setattr(svc, "_fast_copy", copy_spy)
        monkeypatch.setattr(
            svc,
            "_parse_single_trait_data",
            lambda *_args: {"A": 1.0, "B": 2.0, "C": 3.0, "D": 4.0},
        )
        monkeypatch.setattr(svc, "_label_internal_nodes", lambda *_args: {})
        monkeypatch.setattr(svc, "_ancestral_reconstruction", lambda *_args: {})
        monkeypatch.setattr(svc, "_build_parent_map", lambda *_args: {})
        monkeypatch.setattr(
            svc,
            "_compute_branch_rates",
            lambda *_args: branch_rates,
        )
        monkeypatch.setattr(svc, "_plot_rate_map", plot_rate_map)
        monkeypatch.setattr("builtins.print", lambda *_args, **_kwargs: None)

        svc.run()

        assert len(copied_trees) == 1
        assert captured["tree"] is copied_trees[0]
        assert {tip.name for tip in tree.get_terminals()} == {"A", "B", "C", "D"}

    def test_multi_trait_with_column(self):
        try:
            import matplotlib
        except ImportError:
            pytest.skip("matplotlib not installed")

        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
            tmppath = f.name

        try:
            args = _make_args(
                trait_data=MULTI_TRAITS_FILE,
                output=tmppath,
                trait="body_mass",
                json=True,
            )
            svc = TraitRateMap(args)
            svc.run()
            assert os.path.exists(tmppath)
            assert os.path.getsize(tmppath) > 0
        finally:
            if os.path.exists(tmppath):
                os.unlink(tmppath)
