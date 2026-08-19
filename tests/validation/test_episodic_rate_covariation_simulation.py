"""Validate episodic rate covariation with planted and null simulations."""

from argparse import Namespace
from io import StringIO

from Bio import Phylo
import numpy as np
import pytest

from phykit.services.tree.episodic_rate_covariation import (
    EpisodicRateCovariation,
)


pytestmark = pytest.mark.validation


REFERENCE_NEWICK = (
    "((((A:1,B:1):1,(C:1,D:1):1):1,"
    "((E:1,F:1):1,(G:1,H:1):1):1):1,"
    "(((I:1,J:1):1,(K:1,L:1):1):1,"
    "((M:1,N:1):1,(O:1,P:1):1):1):1);"
)


def _service(seed):
    return EpisodicRateCovariation(
        Namespace(
            tree_zero="gene_zero.tre",
            tree_one="gene_one.tre",
            reference="reference.tre",
            verbose=False,
            weighting="uniform",
            max_rate=5.0,
            max_pairs=50_000,
            seed=seed,
            json=False,
            plot=False,
            permutations=999,
            depth_bins=4,
            min_edges=3,
            alpha=0.05,
            output=None,
            annotated_tree=None,
        )
    )


def _scaled_tree(service, reference, taxa, edges, rates):
    tree = service._fast_copy(reference)
    all_taxa = frozenset(taxa)
    edge_index = {edge.split: index for index, edge in enumerate(edges)}
    descendant_taxa = {}
    for clade in tree.find_clades(order="postorder"):
        if clade.clades:
            descendants = frozenset().union(
                *(descendant_taxa[id(child)] for child in clade.clades)
            )
        else:
            descendants = frozenset((clade.name,))
        descendant_taxa[id(clade)] = descendants
        if clade is tree.root:
            continue
        split = service._canonical_split(descendants, all_taxa)
        clade.branch_length *= rates[edge_index[split]]
    return tree


def _analyze_rates(service, rates_zero, rates_one):
    reference = Phylo.read(StringIO(REFERENCE_NEWICK), "newick")
    taxa = sorted(tip.name for tip in reference.get_terminals())
    edges = service._reference_edge_basis(reference, taxa)
    tree_zero = _scaled_tree(service, reference, taxa, edges, rates_zero)
    tree_one = _scaled_tree(service, reference, taxa, edges, rates_one)
    service.read_tree_file_unmodified = lambda: tree_zero
    service.read_tree1_file_unmodified = lambda: tree_one
    service.read_reference_tree_file_unmodified = lambda: reference

    analysis = service._analyze_projected_rates()
    basis = service._clade_scan_basis(
        analysis.reference_tree,
        analysis.edges,
        analysis.retained_edge_indices,
        service.min_edges,
    )
    membership = service._candidate_membership_matrix(
        basis.candidates,
        len(basis.rate_positions),
    )
    strata = service._depth_strata(basis.edge_depths, service.depth_bins)
    rows = service._scan_clades(
        analysis.standardized_zero[basis.rate_positions],
        analysis.standardized_one[basis.rate_positions],
        membership,
        basis.candidates,
        strata,
    )
    return analysis, edges, rows


def test_scan_recovers_planted_clade_level_covariation():
    service = _service(seed=91)
    reference = Phylo.read(StringIO(REFERENCE_NEWICK), "newick")
    taxa = sorted(tip.name for tip in reference.get_terminals())
    edges = service._reference_edge_basis(reference, taxa)
    target_taxa = frozenset(("A", "B", "C", "D"))
    target_indices = [
        index
        for index, edge in enumerate(edges)
        if frozenset(edge.split) <= target_taxa
    ]
    rng = np.random.default_rng(481)
    rates_zero = rng.uniform(0.65, 1.35, len(edges))
    rates_one = rng.uniform(0.65, 1.35, len(edges))
    for offset, index in enumerate(target_indices):
        rates_zero[index] = 2.2 + 0.08 * offset
        rates_one[index] = 2.15 + 0.07 * offset

    analysis, _, rows = _analyze_rates(service, rates_zero, rates_one)

    assert analysis.projection_zero.nrmse == pytest.approx(0.0, abs=1e-12)
    assert analysis.projection_one.nrmse == pytest.approx(0.0, abs=1e-12)
    assert rows[0]["clade"] == "A;B;C;D"
    assert rows[0]["coupling"] == "concordant"
    assert rows[0]["fwer_p_value"] <= 0.01
    assert rows[0]["significant"] is True


def test_scan_has_no_familywise_discovery_for_fixed_independent_null():
    service = _service(seed=97)
    edge_count = 29
    rng = np.random.default_rng(221)
    rates_zero = rng.uniform(0.5, 1.5, edge_count)
    rates_one = rng.uniform(0.5, 1.5, edge_count)

    analysis, edges, rows = _analyze_rates(service, rates_zero, rates_one)

    assert len(edges) == edge_count
    assert abs(analysis.correlation) < 0.1
    assert min(row["fwer_p_value"] for row in rows) > 0.05
    assert not any(row["significant"] for row in rows)
