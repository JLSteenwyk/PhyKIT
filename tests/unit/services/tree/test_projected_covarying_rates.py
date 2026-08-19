"""Tests for reference-edge distance projection."""

from argparse import Namespace
from io import StringIO
import json
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from Bio import Phylo
import numpy as np
import pytest

from phykit.errors import PhykitUserError
from phykit.services.tree.projected_covarying_rates import (
    ProjectedCovaryingRates,
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
        "seed": 0,
        "json": False,
        "plot": False,
    }
    values.update(overrides)
    return ProjectedCovaryingRates(Namespace(**values))


def _tree(newick):
    return Phylo.read(StringIO(newick), "newick")


def test_reference_basis_merges_complementary_root_edges():
    service = _service()
    taxa = ["A", "B", "C", "D"]
    tree = _tree("((A:1,B:2):3,(C:4,D:5):6);")

    edges = service._reference_edge_basis(tree, taxa)

    assert len(edges) == 5
    assert {edge.split: edge.length for edge in edges}[("A", "B")] == 9.0

    pair_i, pair_j = np.triu_indices(len(taxa), k=1)
    design = service._path_design_matrix(taxa, edges, pair_i, pair_j)
    expected = design @ np.asarray([edge.length for edge in edges])
    observed = service._selected_distances(
        tree,
        taxa,
        pair_i,
        pair_j,
        np.arange(len(pair_i)),
    )
    np.testing.assert_allclose(expected, observed)


def test_nnls_projection_exactly_recovers_positive_edge_lengths():
    service = _service()
    taxa = ["A", "B", "C", "D", "E"]
    reference = _tree("((A:1,B:1):1,(C:1,(D:1,E:1):1):1);")
    edges = service._reference_edge_basis(reference, taxa)
    pair_i, pair_j = np.triu_indices(len(taxa), k=1)
    design = service._path_design_matrix(taxa, edges, pair_i, pair_j)
    expected_lengths = np.linspace(0.4, 1.6, len(edges))

    result = service._project_distances(
        design,
        design @ expected_lengths,
        np.ones(len(pair_i)),
    )

    np.testing.assert_allclose(result.edge_lengths, expected_lengths, atol=1e-10)
    assert result.nrmse == pytest.approx(0.0, abs=1e-12)
    assert result.max_absolute_residual == pytest.approx(0.0, abs=1e-12)


def test_discordant_gene_tree_has_nonzero_projection_error():
    service = _service()
    taxa = ["A", "B", "C", "D"]
    reference = _tree("((A:1,B:1):1,(C:1,D:1):1);")
    discordant = _tree("((A:1,C:1):1,(B:1,D:1):1);")
    edges = service._reference_edge_basis(reference, taxa)
    pair_i, pair_j = np.triu_indices(len(taxa), k=1)
    selected = np.arange(len(pair_i))
    design = service._path_design_matrix(taxa, edges, pair_i, pair_j)
    distances = service._selected_distances(
        discordant,
        taxa,
        pair_i,
        pair_j,
        selected,
    )

    result = service._project_distances(
        design,
        distances,
        np.ones(len(distances)),
    )

    assert np.all(result.edge_lengths >= 0.0)
    assert result.nrmse > 0.0


def test_shared_taxa_are_intersection_of_all_three_trees():
    service = _service()
    gene_zero = _tree("((A:1,B:1):1,(C:1,X:1):1);")
    gene_one = _tree("((A:1,B:1):1,(C:1,Y:1):1);")
    reference = _tree("((A:1,B:1):1,(C:1,Z:1):1);")

    assert service._shared_taxa(gene_zero, gene_one, reference) == ["A", "B", "C"]


def test_prune_to_taxa_removes_nonshared_tip():
    service = _service()
    tree = _tree("((A:1,B:1):1,(C:1,X:1):1);")

    pruned = service._prune_to_taxa(tree, ["A", "B", "C"])

    assert sorted(service.get_tip_names_from_tree(pruned)) == ["A", "B", "C"]


def test_shared_taxa_requires_at_least_three_tips():
    service = _service()
    gene_zero = _tree("(A:1,B:1,C:1);")
    gene_one = _tree("(A:1,B:1,D:1);")
    reference = _tree("(A:1,B:1,E:1);")

    with pytest.raises(PhykitUserError) as error:
        service._shared_taxa(gene_zero, gene_one, reference)
    assert "At least 3 shared taxa" in " ".join(error.value.messages)


def test_inverse_reference_weights_downweight_long_distances():
    weights = ProjectedCovaryingRates._pair_weights(
        np.asarray([1.0, 2.0, 4.0]),
        "inverse_reference",
    )

    np.testing.assert_allclose(weights, [1.0, 0.5, 0.25])


def test_pair_subsampling_is_deterministic_and_sorted():
    first = _service(max_pairs=7, seed=23)
    second = _service(max_pairs=7, seed=23)

    first_i, first_j, first_selected, total = first._selected_pair_indices(10, 5)
    second_i, second_j, second_selected, second_total = (
        second._selected_pair_indices(10, 5)
    )

    assert total == second_total == 45
    assert len(first_selected) == 7
    assert np.all(np.diff(first_selected) > 0)
    np.testing.assert_array_equal(first_selected, second_selected)
    np.testing.assert_array_equal(first_i, second_i)
    np.testing.assert_array_equal(first_j, second_j)


def test_rank_deficient_subsample_is_rejected():
    design = np.asarray(
        [
            [1.0, 0.0, 1.0],
            [0.0, 1.0, 1.0],
            [1.0, 1.0, 2.0],
        ]
    )

    with pytest.raises(PhykitUserError) as error:
        ProjectedCovaryingRates._validate_projection_design(
            design,
            edge_count=3,
            was_subsampled=True,
        )

    assert "rank deficient" in " ".join(error.value.messages)


def test_unobserved_reference_edge_is_rejected():
    design = np.asarray([[1.0, 0.0], [1.0, 0.0]])

    with pytest.raises(PhykitUserError) as error:
        ProjectedCovaryingRates._validate_projection_design(
            design,
            edge_count=2,
            was_subsampled=True,
        )

    assert "do not observe every reference edge" in " ".join(
        error.value.messages
    )


def test_full_design_skips_redundant_rank_computation(monkeypatch):
    design = np.eye(3)

    def fail_if_called(*args, **kwargs):
        raise AssertionError("full reference designs have a known identifiable basis")

    monkeypatch.setattr(np.linalg, "matrix_rank", fail_if_called)
    ProjectedCovaryingRates._validate_projection_design(
        design,
        edge_count=3,
        was_subsampled=False,
    )


def test_pair_selection_rejects_too_few_rows_for_edges():
    service = _service(max_pairs=2)

    with pytest.raises(PhykitUserError) as error:
        service._selected_pair_indices(taxon_count=4, edge_count=5)

    assert "At least 5 distance pairs" in " ".join(error.value.messages)


def test_selected_distances_only_calculates_sampled_pairs(monkeypatch):
    service = _service()
    taxa = ["A", "B", "C", "D"]
    tree = _tree("((A:1,B:2):3,(C:4,D:5):6);")
    pair_i, pair_j = np.triu_indices(len(taxa), k=1)
    selected = np.asarray([1, 4])
    monkeypatch.setattr(
        service,
        "calculate_pairwise_tip_distances_fast",
        lambda *args, **kwargs: (_ for _ in ()).throw(
            AssertionError("subsampling should not calculate all distances")
        ),
    )

    observed = service._selected_distances(
        tree,
        taxa,
        pair_i[selected],
        pair_j[selected],
        selected,
    )

    expected = [
        tree.distance(taxa[pair_i[index]], taxa[pair_j[index]])
        for index in selected
    ]
    np.testing.assert_allclose(observed, expected)


def test_selected_distances_falls_back_for_tree_test_double(monkeypatch):
    service = _service()
    tree = SimpleNamespace(distance=lambda first, second: len(first) + len(second))
    monkeypatch.setattr(
        service,
        "calculate_pairwise_tip_distances_fast",
        lambda *args, **kwargs: None,
    )

    observed = service._selected_distances(
        tree,
        ["A", "BB", "CCC"],
        np.asarray([0, 0, 1]),
        np.asarray([1, 2, 2]),
        np.asarray([0, 1, 2]),
    )

    np.testing.assert_allclose(observed, [3.0, 4.0, 5.0])


def test_reference_basis_validation_rejects_mismatched_distances():
    edges = [SimpleNamespace(length=1.0), SimpleNamespace(length=1.0)]

    with pytest.raises(PhykitUserError) as error:
        ProjectedCovaryingRates._validate_reference_basis(
            np.eye(2),
            edges,
            np.asarray([1.0, 2.0]),
        )

    assert "does not reproduce" in " ".join(error.value.messages)


def test_reference_basis_rejects_tree_without_identifiable_edges():
    with pytest.raises(PhykitUserError) as error:
        _service()._reference_edge_basis(_tree("(A:1);"), ["A"])

    assert "fewer than 2 identifiable edges" in " ".join(error.value.messages)


def test_inverse_weights_reject_all_zero_reference_distances():
    with pytest.raises(PhykitUserError) as error:
        ProjectedCovaryingRates._pair_weights(
            np.zeros(3),
            "inverse_reference",
        )

    assert "all zero" in " ".join(error.value.messages)


@pytest.mark.parametrize("distances", [[-1.0, 1.0], [float("nan"), 1.0]])
def test_projection_rejects_invalid_gene_distances(distances):
    with pytest.raises(PhykitUserError) as error:
        ProjectedCovaryingRates._project_distances(
            np.eye(2),
            np.asarray(distances),
            np.ones(2),
        )

    assert "finite and nonnegative" in " ".join(error.value.messages)


def test_projection_falls_back_to_bounded_least_squares():
    with patch("scipy.optimize.nnls", side_effect=RuntimeError("iteration limit")):
        result = ProjectedCovaryingRates._project_distances(
            np.eye(2),
            np.asarray([1.0, 2.0]),
            np.ones(2),
        )

    np.testing.assert_allclose(result.edge_lengths, [1.0, 2.0])
    assert result.nrmse == pytest.approx(0.0)


def test_projection_reports_failure_from_fallback_solver():
    failed = SimpleNamespace(success=False, message="did not converge")
    with (
        patch("scipy.optimize.nnls", side_effect=RuntimeError("iteration limit")),
        patch("scipy.optimize.lsq_linear", return_value=failed),
        pytest.raises(PhykitUserError) as error,
    ):
        ProjectedCovaryingRates._project_distances(
            np.eye(2),
            np.asarray([1.0, 2.0]),
            np.ones(2),
        )

    assert "did not converge" in " ".join(error.value.messages)


def test_zero_distance_projection_has_zero_nrmse():
    result = ProjectedCovaryingRates._project_distances(
        np.eye(2),
        np.zeros(2),
        np.ones(2),
    )

    assert result.nrmse == 0.0


@pytest.mark.parametrize(
    ("overrides", "message"),
    [
        ({"weighting": "bad"}, "--weighting"),
        ({"max_rate": 0}, "--max-rate"),
        ({"max_rate": float("inf")}, "--max-rate"),
        ({"max_pairs": 0}, "--max-pairs"),
    ],
)
def test_invalid_arguments_raise_user_errors(overrides, message):
    with pytest.raises(PhykitUserError) as error:
        _service(**overrides)
    assert message in " ".join(error.value.messages)


def test_rate_rows_filter_zero_reference_edges_and_large_rates():
    service = _service(max_rate=5.0)
    edges = [
        type("Edge", (), {"label": "A", "length": 1.0})(),
        type("Edge", (), {"label": "B", "length": 0.0})(),
        type("Edge", (), {"label": "C", "length": 1.0})(),
    ]

    rows, retained_zero, retained_one = service._rate_rows(
        edges,
        np.asarray([2.0, 1.0, 6.0]),
        np.asarray([3.0, 1.0, 1.0]),
    )

    assert [row["status"] for row in rows] == [
        "retained",
        "zero_reference_length",
        "rate_outlier",
    ]
    assert retained_zero == [2.0]
    assert retained_one == [3.0]


def test_constant_projected_rates_cannot_be_correlated():
    with pytest.raises(PhykitUserError) as error:
        ProjectedCovaryingRates._validate_rate_vectors(
            [1.0, 1.0, 1.0],
            [1.0, 2.0, 3.0],
        )
    assert "no variation" in " ".join(error.value.messages)


def test_rate_vectors_must_have_matching_lengths():
    with pytest.raises(PhykitUserError) as error:
        ProjectedCovaryingRates._validate_rate_vectors(
            [1.0, 2.0, 3.0],
            [1.0, 2.0],
        )
    assert "different lengths" in " ".join(error.value.messages)


def test_rate_vectors_require_three_retained_edges():
    with pytest.raises(PhykitUserError) as error:
        ProjectedCovaryingRates._validate_rate_vectors(
            [1.0, 2.0],
            [2.0, 3.0],
        )
    assert "Fewer than 3" in " ".join(error.value.messages)


def test_run_emits_json_for_discordant_gene_trees(monkeypatch, capsys):
    sample_files = Path(__file__).parents[3] / "sample_files"
    tree_zero = Phylo.read(sample_files / "tree_simple.tre", "newick")
    tree_one = Phylo.read(
        sample_files / "tree_simple_other_topology.tre",
        "newick",
    )
    reference = Phylo.read(sample_files / "tree_simple_2.tre", "newick")
    service = _service(json=True)
    monkeypatch.setattr(service, "read_tree_file_unmodified", lambda: tree_zero)
    monkeypatch.setattr(service, "read_tree1_file_unmodified", lambda: tree_one)
    monkeypatch.setattr(
        service,
        "read_reference_tree_file_unmodified",
        lambda: reference,
    )

    service.run()

    payload = json.loads(capsys.readouterr().out)
    assert payload["shared_taxon_count"] == 8
    assert payload["tree_zero_projection"]["nrmse"] == pytest.approx(
        0.0,
        abs=1e-12,
    )
    assert payload["tree_one_projection"]["nrmse"] > 0.0


def test_verbose_text_output_formats_excluded_values(capsys):
    service = _service(verbose=True)
    payload = {
        "correlation": 0.5,
        "p_value": 0.1,
        "shared_taxon_count": 4,
        "distance_pairs_used": 6,
        "distance_pair_count": 6,
        "reference_edge_count": 5,
        "retained_edge_count": 4,
        "tree_zero_projection": {"nrmse": 0.0},
        "tree_one_projection": {"nrmse": 0.1},
        "branches": [
            {
                "branch": "A",
                "reference_length": 0.0,
                "tree_zero_projected_length": 1.0,
                "tree_one_projected_length": 1.0,
                "tree_zero_relative_rate": None,
                "tree_one_relative_rate": None,
                "tree_zero_zscore": None,
                "tree_one_zscore": None,
                "status": "zero_reference_length",
            }
        ],
    }

    service._print_text(payload)

    output = capsys.readouterr().out
    assert "branch\treference_length\tprojected_zero" in output
    assert "A\t0.000000\t1.000000\t1.000000\tNA\tNA\tNA\tNA" in output


def test_plot_projected_rates_writes_nonempty_image(tmp_path):
    output = tmp_path / "projected_rates.png"
    service = _service(plot=True, plot_output=str(output))

    service._plot_projected_rates(
        [-1.0, -0.2, 0.4, 1.2],
        [-0.8, 0.1, 0.3, 1.0],
        correlation=0.95,
        p_value=0.05,
    )

    assert output.is_file()
    assert output.stat().st_size > 0


def test_text_output_tolerates_closed_pipe(monkeypatch):
    service = _service()
    payload = {
        "correlation": 0.5,
        "p_value": 0.1,
        "shared_taxon_count": 4,
        "distance_pairs_used": 6,
        "distance_pair_count": 6,
        "reference_edge_count": 5,
        "retained_edge_count": 5,
        "tree_zero_projection": {"nrmse": 0.0},
        "tree_one_projection": {"nrmse": 0.1},
    }
    monkeypatch.setattr(
        "builtins.print",
        lambda *args, **kwargs: (_ for _ in ()).throw(BrokenPipeError()),
    )

    service._print_text(payload)
