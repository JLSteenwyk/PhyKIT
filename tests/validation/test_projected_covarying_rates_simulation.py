"""Validate projected rate covariation against known synthetic edge rates."""

from argparse import Namespace
from io import StringIO

from Bio import Phylo
import numpy as np
import pytest

from phykit.services.tree.covarying_evolutionary_rates import _pearsonr, _zscore
from phykit.services.tree.projected_covarying_rates import ProjectedCovaryingRates


pytestmark = pytest.mark.validation


def test_projection_recovers_simulated_edge_rate_correlation():
    service = ProjectedCovaryingRates(
        Namespace(
            tree_zero="gene_zero.tre",
            tree_one="gene_one.tre",
            reference="reference.tre",
            verbose=False,
        )
    )
    taxa = ["A", "B", "C", "D", "E", "F", "G", "H"]
    reference = Phylo.read(
        StringIO(
            "(((A:0.7,B:0.9):0.4,(C:0.8,D:1.1):0.6):0.5,"
            "((E:1.0,F:0.6):0.7,(G:1.2,H:0.8):0.3):0.9);"
        ),
        "newick",
    )
    edges = service._reference_edge_basis(reference, taxa)
    pair_i, pair_j = np.triu_indices(len(taxa), k=1)
    design = service._path_design_matrix(taxa, edges, pair_i, pair_j)
    reference_lengths = np.asarray([edge.length for edge in edges])

    rng = np.random.default_rng(814)
    rates_zero = rng.uniform(0.35, 1.75, len(edges))
    rates_one = 0.7 * rates_zero + 0.3 * rng.uniform(0.35, 1.75, len(edges))
    lengths_zero = reference_lengths * rates_zero
    lengths_one = reference_lengths * rates_one
    weights = np.ones(len(pair_i))

    projected_zero = service._project_distances(
        design,
        design @ lengths_zero,
        weights,
    )
    projected_one = service._project_distances(
        design,
        design @ lengths_one,
        weights,
    )
    recovered_zero = projected_zero.edge_lengths / reference_lengths
    recovered_one = projected_one.edge_lengths / reference_lengths
    recovered_correlation, _ = _pearsonr(
        _zscore(recovered_zero),
        _zscore(recovered_one),
    )
    expected_correlation = float(np.corrcoef(rates_zero, rates_one)[0, 1])

    np.testing.assert_allclose(recovered_zero, rates_zero, atol=1e-10)
    np.testing.assert_allclose(recovered_one, rates_one, atol=1e-10)
    assert projected_zero.nrmse == pytest.approx(0.0, abs=1e-12)
    assert projected_one.nrmse == pytest.approx(0.0, abs=1e-12)
    assert recovered_correlation == pytest.approx(expected_correlation, abs=1e-12)
