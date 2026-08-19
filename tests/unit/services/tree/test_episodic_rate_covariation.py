"""Tests for localized evolutionary-rate covariation scans."""

from argparse import Namespace
from io import StringIO
from types import SimpleNamespace

from Bio import Phylo
import numpy as np
import pytest

from phykit.errors import PhykitUserError
from phykit.services.tree.episodic_rate_covariation import (
    CladeScanCandidate,
    EpisodicRateCovariation,
)


def _service(**overrides):
    values = {
        "tree_zero": "gene_zero.tre",
        "tree_one": "gene_one.tre",
        "reference": "reference.tre",
        "verbose": False,
        "weighting": "uniform",
        "max_rate": 5.0,
        "max_pairs": 50_000,
        "seed": 17,
        "json": False,
        "plot": False,
        "permutations": 99,
        "depth_bins": 4,
        "min_edges": 2,
        "alpha": 0.05,
        "output": None,
        "annotated_tree": None,
    }
    values.update(overrides)
    return EpisodicRateCovariation(Namespace(**values))


def _tree(newick):
    return Phylo.read(StringIO(newick), "newick")


def _candidate(label, edge_indices, depth=1):
    return CladeScanCandidate(
        label=label,
        taxa=tuple(label.split(";")),
        edge_indices=tuple(edge_indices),
        root_depth=depth,
    )


@pytest.mark.parametrize(
    ("overrides", "message"),
    [
        ({"permutations": 0}, "--permutations"),
        ({"depth_bins": 0}, "--depth-bins"),
        ({"min_edges": 1}, "--min-edges"),
        ({"alpha": 0.0}, "--alpha"),
        ({"alpha": float("inf")}, "--alpha"),
        ({"alpha": 1.1}, "--alpha"),
    ],
)
def test_invalid_arguments_raise_user_errors(overrides, message):
    with pytest.raises(PhykitUserError) as error:
        _service(**overrides)

    assert message in " ".join(error.value.messages)


def test_clade_scan_basis_excludes_merged_root_edge():
    service = _service()
    tree = _tree(
        "(((A:1,B:1):1,(C:1,D:1):1):1,"
        "((E:1,F:1):1,(G:1,H:1):1):1);"
    )
    taxa = sorted(tip.name for tip in tree.get_terminals())
    edges = service._reference_edge_basis(tree, taxa)

    basis = service._clade_scan_basis(
        tree,
        edges,
        np.arange(len(edges)),
        min_edges=2,
    )

    assert len(edges) == 13
    assert len(basis.rate_positions) == 12
    assert basis.ambiguous_edge_count == 1
    assert len(basis.candidates) == 6
    assert set(basis.edge_labels).isdisjoint({"A;B;C;D", "E;F;G;H"})
    assert {candidate.label for candidate in basis.candidates} == {
        "A;B",
        "A;B;C;D",
        "C;D",
        "E;F",
        "E;F;G;H",
        "G;H",
    }
    assert all(
        len(candidate.edge_indices) >= 2
        and len(basis.rate_positions) - len(candidate.edge_indices) >= 2
        for candidate in basis.candidates
    )


def test_clade_scan_basis_requires_enough_localizable_edges():
    service = _service(min_edges=3)
    tree = _tree("((A:1,B:1):1,(C:1,D:1):1);")
    taxa = sorted(tip.name for tip in tree.get_terminals())
    edges = service._reference_edge_basis(tree, taxa)

    with pytest.raises(PhykitUserError) as error:
        service._clade_scan_basis(
            tree,
            edges,
            np.arange(len(edges)),
            min_edges=3,
        )

    assert "localizable projected edges" in " ".join(error.value.messages)


def test_candidate_membership_matrix_marks_local_edges():
    candidates = (
        _candidate("A;B", (0, 2)),
        _candidate("C;D", (1, 3, 4)),
    )

    observed = EpisodicRateCovariation._candidate_membership_matrix(
        candidates,
        edge_count=5,
    )

    np.testing.assert_array_equal(
        observed,
        [[1.0, 0.0, 1.0, 0.0, 0.0], [0.0, 1.0, 0.0, 1.0, 1.0]],
    )


def test_depth_strata_preserve_ties_and_requested_limit():
    depths = np.asarray([1, 1, 2, 2, 3, 3, 4, 5], dtype=float)

    strata = EpisodicRateCovariation._depth_strata(depths, requested_bins=3)

    assert len(np.unique(strata)) <= 3
    for depth in np.unique(depths):
        assert len(np.unique(strata[depths == depth])) == 1


def test_scan_scores_match_inside_background_contrast():
    products = np.asarray([2.0, 2.0, -1.0, -1.0, -1.0, -1.0])
    membership = np.asarray([[1.0, 1.0, 0.0, 0.0, 0.0, 0.0]])

    scores, local_means, background_means = (
        EpisodicRateCovariation._scan_scores(products, membership)
    )

    assert local_means[0] == pytest.approx(2.0)
    assert background_means[0] == pytest.approx(-1.0)
    assert scores[0] == pytest.approx(4.0 / np.sqrt(4.0 / 3.0))


def test_scan_scores_support_permutation_batches():
    products = np.asarray(
        [[2.0, 2.0, -1.0, -1.0], [-1.0, 2.0, 2.0, -1.0]]
    )
    membership = np.asarray([[1.0, 1.0, 0.0, 0.0]])

    scores = EpisodicRateCovariation._scan_scores(products, membership)[0]

    assert scores.shape == (2, 1)
    assert scores[0, 0] > scores[1, 0]


def test_clade_scan_is_deterministic_and_familywise_corrected():
    candidates = (
        _candidate("A;B", (0, 1, 2)),
        _candidate("C;D", (3, 4, 5)),
        _candidate("E;F", (6, 7, 8)),
    )
    membership = EpisodicRateCovariation._candidate_membership_matrix(
        candidates,
        edge_count=9,
    )
    tree_zero = np.asarray([2.0, 1.5, 1.0, -2.0, -1.5, -1.0, 0.5, -0.5, 0.0])
    tree_one = np.asarray([1.8, 1.4, 0.8, 1.8, 1.4, 0.8, -0.4, 0.4, 0.0])
    strata = np.zeros(9, dtype=np.intp)

    first = _service(seed=23)._scan_clades(
        tree_zero,
        tree_one,
        membership,
        candidates,
        strata,
    )
    second = _service(seed=23)._scan_clades(
        tree_zero,
        tree_one,
        membership,
        candidates,
        strata,
    )

    assert first == second
    assert all(
        row["fwer_p_value"] >= row["unadjusted_p_value"] for row in first
    )
    assert all(
        1.0 / 100.0 <= row["unadjusted_p_value"] <= 1.0 for row in first
    )


def test_clade_scan_rejects_singleton_depth_strata():
    candidates = (_candidate("A;B", (0, 1)),)
    membership = EpisodicRateCovariation._candidate_membership_matrix(
        candidates,
        edge_count=4,
    )

    with pytest.raises(PhykitUserError) as error:
        _service()._scan_clades(
            np.asarray([1.0, 2.0, 3.0, 4.0]),
            np.asarray([4.0, 3.0, 2.0, 1.0]),
            membership,
            candidates,
            np.arange(4),
        )

    assert "exchangeable within depth bins" in " ".join(error.value.messages)


def test_leave_one_out_correlations_match_direct_calculation():
    tree_zero = np.asarray([-1.2, -0.5, 0.2, 0.4, 1.1, 1.7])
    tree_one = np.asarray([-0.8, 0.1, -0.3, 0.7, 1.4, 0.9])

    observed = EpisodicRateCovariation._leave_one_out_correlations(
        tree_zero,
        tree_one,
    )
    expected = [
        np.corrcoef(np.delete(tree_zero, index), np.delete(tree_one, index))[0, 1]
        for index in range(len(tree_zero))
    ]

    np.testing.assert_allclose(observed, expected)


def test_add_edge_diagnostics_marks_excluded_rows_and_largest_influence():
    rows = [
        {"branch": "A", "status": "retained"},
        {"branch": "B", "status": "rate_outlier"},
        {"branch": "C", "status": "retained"},
        {"branch": "D", "status": "retained"},
        {"branch": "E", "status": "retained"},
    ]
    analysis = SimpleNamespace(
        standardized_zero=np.asarray([-1.0, -0.2, 0.3, 1.4]),
        standardized_one=np.asarray([-0.9, 0.4, 0.1, 1.2]),
        retained_edge_indices=np.asarray([0, 2, 3, 4]),
        rows=rows,
        correlation=float(
            np.corrcoef(
                [-1.0, -0.2, 0.3, 1.4],
                [-0.9, 0.4, 0.1, 1.2],
            )[0, 1]
        ),
    )

    influence = _service()._add_edge_diagnostics(analysis)

    assert rows[1]["local_contribution"] is None
    assert rows[0]["local_contribution"] == pytest.approx(0.9)
    assert influence["branch"] in {"A", "C", "D", "E"}
    assert influence["absolute_delta"] == pytest.approx(abs(influence["delta"]))


def test_scan_table_contains_all_diagnostics(tmp_path):
    output = tmp_path / "clades.tsv"
    row = {
        "clade": "A;B",
        "taxon_count": 2,
        "root_depth": 2,
        "edge_count": 3,
        "score": 4.5,
        "local_mean_contribution": 2.0,
        "background_mean_contribution": -1.0,
        "contribution_contrast": 3.0,
        "coupling": "concordant",
        "unadjusted_p_value": 0.01,
        "fwer_p_value": 0.03,
        "significant": True,
    }

    EpisodicRateCovariation._write_scan_table([row], output)

    contents = output.read_text()
    assert "contribution_contrast" in contents
    assert "A;B\t2\t2\t3\t4.5" in contents


def test_annotated_tree_marks_only_significant_clades(tmp_path):
    service = _service()
    output = tmp_path / "annotated.tre"
    tree = _tree("((A:1,B:1):1,(C:1,D:1):1);")
    rows = [
        {
            "taxa": ["A", "B"],
            "score": 4.25,
            "fwer_p_value": 0.01,
            "coupling": "concordant",
            "significant": True,
        },
        {
            "taxa": ["C", "D"],
            "score": 1.0,
            "fwer_p_value": 0.5,
            "coupling": "concordant",
            "significant": False,
        },
    ]

    service._write_annotated_tree(tree, rows, output)

    annotated = Phylo.read(output, "newick")
    comments = {
        tuple(sorted(tip.name for tip in clade.get_terminals())): clade.comment
        for clade in annotated.get_nonterminals()
    }
    assert "erc_score=4.25" in comments[("A", "B")]
    assert comments[("C", "D")] is None
