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
from argparse import Namespace
from pathlib import Path

from phykit.services.tree.parsimony_score import ParsimonyScore


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
    def test_score_matches_r(self, args, capsys):
        """Fitch parsimony score should match R's phangorn::parsimony() = 4."""
        ps = ParsimonyScore(args)
        ps.run()
        captured = capsys.readouterr()
        assert captured.out.strip().split("\n")[0] == "4"

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
        lines = captured.out.strip().split("\n")
        assert lines[0] == "4"  # Total score
        # Per-site scores follow (16 sites)
        per_site_lines = [l for l in lines[1:] if l.strip()]
        assert len(per_site_lines) == 16

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


class TestParsimonyScoreJson:
    def test_json_output(self, mocker, args):
        args.json = True
        ps = ParsimonyScore(args)
        mocked_json = mocker.patch(
            "phykit.services.tree.parsimony_score.print_json"
        )
        ps.run()
        payload = mocked_json.call_args.args[0]
        # Cross-validated against R's phangorn::parsimony() = 4
        assert payload["parsimony_score"] == 4
        assert payload["alignment_length"] == 16
        assert payload["n_taxa"] == 8
