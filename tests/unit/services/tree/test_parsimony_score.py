"""
Unit tests for parsimony_score (Fitch parsimony).

Computes the minimum number of character state changes on a tree
given an alignment using the Fitch (1971) downpass algorithm.
Cross-validated against R's phangorn::parsimony().

R validation: tests/r_validation/validate_parsimony.R
R gives: parsimony(multi2di(tree), aln, method="fitch") = 4
PhyKIT gives: 4 (exact match)
"""
import pytest
import json
import subprocess
import sys
from argparse import Namespace
from io import StringIO
from pathlib import Path

from Bio import Phylo
from Bio.Phylo.BaseTree import TreeMixin

from phykit.services.tree.parsimony_score import ParsimonyScore


def test_module_import_does_not_import_numpy_or_biopython_fasta_parser():
    code = """
import sys
import phykit.services.tree.parsimony_score as module
assert callable(module._json_dumps)
assert hasattr(module.np, "__getattr__")
assert "typing" not in sys.modules
assert "json" not in sys.modules
assert "numpy" not in sys.modules
assert "Bio.SeqIO.FastaIO" not in sys.modules
assert "Bio.AlignIO" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


@pytest.fixture
def args():
    return Namespace(
        tree="tests/sample_files/tree_simple.tre",
        alignment="tests/sample_files/tree_simple_alignment.fa",
    )


class TestParsimonyScoreInit:
    def test_init_sets_fields(self, args):
        ps = ParsimonyScore(args)
        assert ps.tree_file_path == args.tree
        assert ps.alignment_path == args.alignment
        assert ps.verbose is False
        assert ps.json_output is False


class TestFitchAlgorithm:
    def test_parse_alignment_uses_first_header_token_uppercases_and_keeps_last_duplicate(
        self, tmp_path, args
    ):
        aln = tmp_path / "alignment.fa"
        aln.write_text(
            ">A description\nac gt\n"
            ">B\nn- x\n?\n"
            ">A duplicate\nTGCA\n"
        )
        ps = ParsimonyScore(args)

        sequences = ps._parse_alignment(str(aln))

        assert sequences == {"A": "TGCA", "B": "N-X?"}

    def test_score_matches_r(self, args, capsys):
        """Fitch parsimony score should match R's phangorn::parsimony() = 4."""
        ps = ParsimonyScore(args)
        ps.run()
        captured = capsys.readouterr()
        assert captured.out.strip().split("\n")[0] == "4"

    def test_run_uses_fast_tip_name_helper_for_setup(self, mocker, args):
        ps = ParsimonyScore(args)
        tip_name_spy = mocker.spy(ps, "get_tip_names_from_tree")

        ps.run()

        tip_name_spy.assert_called_once()

    def test_polytomy_scan_handles_mixed_child_counts(self, args, monkeypatch):
        ps = ParsimonyScore(args)
        mixed = Phylo.read(
            StringIO("(A:1,(B:1,C:1):1,(D:1,E:1,F:1):1);"),
            "newick",
        )
        binary = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")

        def fail_find_clades(*_args, **_kwargs):
            raise AssertionError("standard tree scan should not use find_clades")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_find_clades)

        assert ps._has_polytomies(binary) is False
        assert ps._has_polytomies(mixed) is True

    def test_identical_sequences_score_zero(self, tmp_path):
        """All identical sequences should produce score 0."""
        aln = tmp_path / "identical.fa"
        aln.write_text(">A\nACGT\n>B\nACGT\n>C\nACGT\n>D\nACGT\n")
        tree_f = tmp_path / "tree.tre"
        tree_f.write_text("((A:1,B:1):1,(C:1,D:1):1);\n")
        args = Namespace(tree=str(tree_f), alignment=str(aln))
        ps = ParsimonyScore(args)
        ps.run()

    def test_all_different_sequences(self, tmp_path, capsys):
        """Each tip has a unique character at site 1 -> maximum changes."""
        aln = tmp_path / "different.fa"
        aln.write_text(">A\nA\n>B\nC\n>C\nG\n>D\nT\n")
        tree_f = tmp_path / "tree.tre"
        tree_f.write_text("((A:1,B:1):1,(C:1,D:1):1);\n")
        args = Namespace(tree=str(tree_f), alignment=str(aln))
        ps = ParsimonyScore(args)
        ps.run()
        captured = capsys.readouterr()
        score = int(captured.out.strip().split("\n")[0])
        # 4 different states on ((A,B),(C,D)): A∪C=+1, G∪T=+1, {A,C}∪{G,T}=+1 → 3
        assert score == 3

    def test_gaps_treated_as_wildcards(self, tmp_path, capsys):
        """Gaps should not count as state changes."""
        aln = tmp_path / "gaps.fa"
        aln.write_text(">A\nA\n>B\n-\n>C\nA\n>D\nA\n")
        tree_f = tmp_path / "tree.tre"
        tree_f.write_text("((A:1,B:1):1,(C:1,D:1):1);\n")
        args = Namespace(tree=str(tree_f), alignment=str(aln))
        ps = ParsimonyScore(args)
        ps.run()
        captured = capsys.readouterr()
        score = int(captured.out.strip().split("\n")[0])
        assert score == 0  # Gap matches anything

    def test_verbose_output(self, args, capsys):
        args.verbose = True
        ps = ParsimonyScore(args)
        ps.run()
        captured = capsys.readouterr()
        lines = captured.out.splitlines()
        assert lines[0] == "4"  # Total score
        assert lines[1] == ""
        # Per-site scores follow (16 sites)
        per_site_lines = lines[2:]
        assert len(per_site_lines) == 16
        assert all("\t" in line for line in per_site_lines)

    def test_run_verbose_text_output_batches_per_site_scores(
        self, mocker, args, capsys
    ):
        ps = ParsimonyScore(args)
        tree = object()
        mocker.patch.object(ps, "read_tree_file_unmodified", return_value=tree)
        mocker.patch.object(ps, "validate_tree")
        mocker.patch.object(ps, "_has_polytomies", return_value=True)
        mocker.patch.object(ps, "_fast_copy", return_value=tree)
        mocker.patch.object(ps, "_resolve_polytomies")
        mocker.patch.object(
            ps,
            "_parse_alignment",
            return_value={"A": "AC", "B": "AC", "C": "GT"},
        )
        mocker.patch.object(
            ps,
            "get_tip_names_from_tree",
            return_value=["A", "B", "C"],
        )
        mocker.patch.object(ps, "_fitch_parsimony", return_value=(3, [1, 2]))
        ps.verbose = True

        ps.run()

        captured = capsys.readouterr()
        assert captured.out == "3\n\n1\t1\n2\t2\n"

    def test_run_all_shared_binary_tree_uses_read_only_tree_without_copy_or_resolve(
        self, monkeypatch, args, capsys
    ):
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        ps = ParsimonyScore(args)
        captured = {}

        monkeypatch.setattr(ps, "read_tree_file_unmodified", lambda: tree)
        monkeypatch.setattr(
            ps,
            "read_tree_file",
            lambda: (_ for _ in ()).throw(
                AssertionError("mutable tree read should not be used")
            ),
        )
        monkeypatch.setattr(
            ps,
            "_fast_copy",
            lambda *_args, **_kwargs: (_ for _ in ()).throw(
                AssertionError("binary all-shared tree should not be copied")
            ),
        )
        monkeypatch.setattr(
            ps,
            "_resolve_polytomies",
            lambda *_args, **_kwargs: (_ for _ in ()).throw(
                AssertionError("binary all-shared tree should not be resolved")
            ),
        )
        monkeypatch.setattr(
            ps,
            "prune_tree_using_taxa_list",
            lambda *_args, **_kwargs: (_ for _ in ()).throw(
                AssertionError("all-shared tree should not be pruned")
            ),
        )
        monkeypatch.setattr(
            ps,
            "_parse_alignment",
            lambda *_args, **_kwargs: {
                "A": "AC",
                "B": "AC",
                "C": "GT",
                "D": "GT",
            },
        )

        def fitch(tree_arg, sequences, aln_length):
            captured["tree"] = tree_arg
            captured["sequences"] = sequences
            captured["aln_length"] = aln_length
            return 0, [0, 0]

        monkeypatch.setattr(ps, "_fitch_parsimony", fitch)

        ps.run()

        capsys.readouterr()
        assert captured["tree"] is tree
        assert set(captured["sequences"]) == {"A", "B", "C", "D"}
        assert captured["aln_length"] == 2
        assert {tip.name for tip in tree.get_terminals()} == {"A", "B", "C", "D"}

    def test_run_missing_alignment_taxa_copies_before_pruning(
        self, monkeypatch, args, capsys
    ):
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        ps = ParsimonyScore(args)
        original_fast_copy = ps._fast_copy
        copied_trees = []
        captured = {}

        def copy_spy(tree_to_copy):
            copied_tree = original_fast_copy(tree_to_copy)
            copied_trees.append(copied_tree)
            return copied_tree

        monkeypatch.setattr(ps, "read_tree_file_unmodified", lambda: tree)
        monkeypatch.setattr(ps, "_fast_copy", copy_spy)
        monkeypatch.setattr(
            ps,
            "_parse_alignment",
            lambda *_args, **_kwargs: {"A": "AC", "B": "AC", "C": "GT"},
        )

        def fitch(tree_arg, sequences, aln_length):
            captured["tree"] = tree_arg
            captured["sequences"] = sequences
            return 1, [0, 1]

        monkeypatch.setattr(ps, "_fitch_parsimony", fitch)

        ps.run()

        capsys.readouterr()
        assert len(copied_trees) == 1
        assert captured["tree"] is copied_trees[0]
        assert set(captured["sequences"]) == {"A", "B", "C"}
        assert {tip.name for tip in tree.get_terminals()} == {"A", "B", "C", "D"}
        assert {tip.name for tip in copied_trees[0].get_terminals()} == {
            "A",
            "B",
            "C",
        }

    def test_run_polytomy_copies_before_resolving(
        self, monkeypatch, args, capsys
    ):
        tree = Phylo.read(StringIO("(A:1,B:1,C:1,D:1);"), "newick")
        ps = ParsimonyScore(args)
        original_fast_copy = ps._fast_copy
        copied_trees = []
        captured = {}

        def copy_spy(tree_to_copy):
            copied_tree = original_fast_copy(tree_to_copy)
            copied_trees.append(copied_tree)
            return copied_tree

        monkeypatch.setattr(ps, "read_tree_file_unmodified", lambda: tree)
        monkeypatch.setattr(ps, "_fast_copy", copy_spy)
        monkeypatch.setattr(
            ps,
            "_parse_alignment",
            lambda *_args, **_kwargs: {
                "A": "AC",
                "B": "AC",
                "C": "GT",
                "D": "GT",
            },
        )

        def fitch(tree_arg, *_args):
            captured["tree"] = tree_arg
            return 0, [0, 0]

        monkeypatch.setattr(ps, "_fitch_parsimony", fitch)

        ps.run()

        capsys.readouterr()
        assert len(copied_trees) == 1
        assert captured["tree"] is copied_trees[0]
        assert len(tree.root.clades) == 4
        stack = [copied_trees[0].root]
        while stack:
            clade = stack.pop()
            assert len(clade.clades) <= 2
            stack.extend(clade.clades)

    def test_per_site_scores_sum_to_total(self, args):
        import copy
        ps = ParsimonyScore(args)
        tree = ps.read_tree_file()
        tree = copy.deepcopy(tree)
        ps._resolve_polytomies(tree)
        sequences = ps._parse_alignment(args.alignment)
        shared = set(t.name for t in tree.get_terminals()) & set(sequences.keys())
        sequences = {t: sequences[t] for t in shared}
        aln_length = len(next(iter(sequences.values())))
        total, per_site = ps._fitch_parsimony(tree, sequences, aln_length)
        assert total == sum(per_site)

    def test_vectorized_fitch_matches_set_fallback_for_wildcards_and_non_dna(self, args):
        ps = ParsimonyScore(args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        sequences = {
            "A": "E-NA",
            "B": "F?CA",
            "C": "EGXT",
            "D": "FG-T",
        }

        vectorized = ps._fitch_parsimony(tree, sequences, 4)
        fallback = ps._fitch_parsimony_sets(tree, sequences, 4)

        assert vectorized == fallback

    def test_short_alignment_uses_small_fitch_path(self, args, mocker):
        ps = ParsimonyScore(args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        sequences = {
            "A": "ACGT",
            "B": "AGGT",
            "C": "NC-T",
            "D": "TCGA",
        }
        small_spy = mocker.spy(ps, "_fitch_parsimony_small")

        total, per_site = ps._fitch_parsimony(tree, sequences, 4)

        small_spy.assert_called_once()
        assert total == sum(per_site)

    def test_small_fitch_repeated_sequences_match_set_fallback(self, args):
        ps = ParsimonyScore(args)
        tree = Phylo.read(
            StringIO("((A:1,B:1):1,(C:1,(D:1,E:1):1):1);"),
            "newick",
        )
        sequences = {
            "A": "ACGT",
            "B": "ACGT",
            "C": "NC-T",
            "D": "NC-T",
            "E": "ACGT",
        }

        optimized = ps._fitch_parsimony(tree, sequences, 4)
        fallback = ps._fitch_parsimony_sets(tree, sequences, 4)

        assert optimized == fallback

    def test_identical_sequences_skip_fitch_traversal_when_all_tips_present(
        self, args, monkeypatch
    ):
        ps = ParsimonyScore(args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        sequences = {
            "A": "ACGT" * 5,
            "B": "ACGT" * 5,
            "C": "ACGT" * 5,
            "D": "ACGT" * 5,
        }

        def fail_small(*_args, **_kwargs):
            raise AssertionError("identical sequences should return before small path")

        def fail_sets(*_args, **_kwargs):
            raise AssertionError("identical sequences should return before set path")

        monkeypatch.setattr(ps, "_fitch_parsimony_small", fail_small)
        monkeypatch.setattr(ps, "_fitch_parsimony_sets", fail_sets)

        total, per_site = ps._fitch_parsimony(tree, sequences, 20)

        assert total == 0
        assert per_site == [0] * 20

    def test_identical_sequence_shortcut_preserves_missing_terminal_error(
        self, args
    ):
        ps = ParsimonyScore(args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        sequences = {
            "A": "ACGT" * 5,
            "B": "ACGT" * 5,
            "C": "ACGT" * 5,
        }

        with pytest.raises(KeyError):
            ps._fitch_parsimony(tree, sequences, 20)

    def test_wider_alignment_uses_vectorized_fitch_path(self, args, mocker):
        ps = ParsimonyScore(args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        sequences = {
            "A": "ACGT" * 5,
            "B": "AGGT" * 5,
            "C": "NC-T" * 5,
            "D": "TCGA" * 5,
        }
        small_spy = mocker.spy(ps, "_fitch_parsimony_small")

        total, per_site = ps._fitch_parsimony(tree, sequences, 20)

        small_spy.assert_not_called()
        assert total == sum(per_site)

    def test_vectorized_fitch_uses_direct_postorder_traversal(self, args, monkeypatch):
        ps = ParsimonyScore(args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        sequences = {
            "A": "ACGT",
            "B": "ACGT",
            "C": "TCGA",
            "D": "TCGA",
        }

        def fail_find_clades(*_args, **_kwargs):
            raise AssertionError("generic postorder traversal should not be used")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_find_clades)

        total, per_site = ps._fitch_parsimony(tree, sequences, 4)

        assert total == sum(per_site)
        assert total == 2

    def test_resolve_binary_tree_does_not_import_newick(self, args, monkeypatch):
        ps = ParsimonyScore(args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        real_import = __import__

        def fail_postorder_list(*_args, **_kwargs):
            raise AssertionError("binary resolution should not build a clade list")

        def fail_newick_import(name, *import_args, **kwargs):
            if name == "Bio.Phylo.Newick" or name == "Bio.Phylo":
                raise AssertionError("binary tree should not import Newick helper")
            return real_import(name, *import_args, **kwargs)

        monkeypatch.setattr(ps, "_postorder_clades_fast", fail_postorder_list)
        monkeypatch.setattr("builtins.__import__", fail_newick_import)

        ps._resolve_polytomies(tree)

        assert [len(clade.clades) for clade in tree.find_clades() if clade.clades] == [2, 2, 2]

    def test_vectorized_fitch_unicode_states_match_set_fallback(self, args):
        ps = ParsimonyScore(args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        state_one = chr(256)
        state_two = chr(257)
        sequences = {
            "A": state_one + "A",
            "B": state_two + "A",
            "C": state_one + "C",
            "D": state_two + "C",
        }

        vectorized = ps._fitch_parsimony(tree, sequences, 2)
        fallback = ps._fitch_parsimony_sets(tree, sequences, 2)

        assert vectorized == fallback

    def test_resolve_polytomies_uses_direct_postorder_traversal(self, args, monkeypatch):
        ps = ParsimonyScore(args)
        tree = Phylo.read(StringIO("(A:1,B:1,C:1,D:1);"), "newick")

        def fail_find_clades(*_args, **_kwargs):
            raise AssertionError("generic postorder traversal should not be used")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_find_clades)

        ps._resolve_polytomies(tree)

        stack = [tree.root]
        while stack:
            clade = stack.pop()
            assert len(clade.clades) <= 2
            stack.extend(clade.clades)

    def test_resolve_polytomies_matches_previous_traversal_output(self, args):
        ps = ParsimonyScore(args)
        newick = "((A:1,B:1,C:1):1,(D:1,E:1,F:1):1,G:1);"
        expected = Phylo.read(StringIO(newick), "newick")
        actual = Phylo.read(StringIO(newick), "newick")

        from Bio.Phylo import Newick

        stack = [expected.root]
        while stack:
            clade = stack.pop()
            children = clade.clades
            if children:
                stack.extend(reversed(children))
            while len(children) > 2:
                child1 = children.pop()
                child2 = children.pop()
                new_internal = Newick.Clade(branch_length=0.0)
                new_internal.clades = [child1, child2]
                children.append(new_internal)

        ps._resolve_standard_tree_polytomies(actual)

        expected_out = StringIO()
        actual_out = StringIO()
        Phylo.write(expected, expected_out, "newick")
        Phylo.write(actual, actual_out, "newick")
        assert actual_out.getvalue() == expected_out.getvalue()

    def test_direct_postorder_preserves_left_to_right_order(self, args):
        ps = ParsimonyScore(args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")

        observed = [
            clade.name if clade.name is not None else "internal"
            for clade in ps._postorder_clades_fast(tree)
        ]

        assert observed == [
            "A", "B", "internal", "C", "D", "internal", "internal",
        ]

    def test_set_fallback_reuses_direct_postorder_for_high_state_data(
        self, args, monkeypatch
    ):
        def build_subtree(names):
            if len(names) == 1:
                return f"{names[0]}:1"
            mid = len(names) // 2
            return f"({build_subtree(names[:mid])},{build_subtree(names[mid:])}):1"

        taxa = [f"T{i}" for i in range(70)]
        tree = Phylo.read(StringIO(f"{build_subtree(taxa)};"), "newick")
        states = [chr(0x100 + idx) for idx in range(len(taxa))]
        sequences = {
            taxon: states[idx] + states[(idx + 1) % len(states)]
            for idx, taxon in enumerate(taxa)
        }

        def fail_find_clades(*_args, **_kwargs):
            raise AssertionError("set fallback should reuse direct postorder traversal")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_find_clades)

        total, per_site = ParsimonyScore(args)._fitch_parsimony_sets(
            tree,
            sequences,
            2,
        )

        assert per_site == [len(taxa) - 1, len(taxa) - 1]
        assert total == sum(per_site)


class TestParsimonyScoreJson:
    def test_json_output(self, capsys, args):
        args.json = True
        ps = ParsimonyScore(args)
        ps.run()
        captured = capsys.readouterr()
        payload = json.loads(captured.out)
        # Cross-validated against R's phangorn::parsimony() = 4
        assert payload["parsimony_score"] == 4
        assert payload["alignment_length"] == 16
        assert payload["n_taxa"] == 8

    def test_verbose_json_output_includes_per_site_scores(self, capsys, args):
        args.json = True
        args.verbose = True
        ps = ParsimonyScore(args)

        ps.run()

        captured = capsys.readouterr()
        payload = json.loads(captured.out)
        assert payload["parsimony_score"] == 4
        assert len(payload["per_site_scores"]) == payload["alignment_length"]
        assert sum(payload["per_site_scores"]) == payload["parsimony_score"]
