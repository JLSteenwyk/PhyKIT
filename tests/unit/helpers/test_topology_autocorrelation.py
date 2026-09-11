import itertools

import numpy as np
import pytest

from phykit.helpers.topology_autocorrelation import (
    analyze,
    baseline,
    bootstrap_reference,
    chromosome_marks,
    estimate,
    read_classified,
    uncertainty,
    validate_edges,
    validate_rows,
)
from phykit.helpers.topology_landscape import CLASSES


def genes(positions, labels, chromosome="chr1", reference="ref", prefix="g"):
    return [{"gene_id": f"{prefix}{i}", "reference": reference, "chromosome": chromosome,
                 "start": int(x), "end": int(x) + 1, "anchor": int(x), "classification": CLASSES[k]}
            for i, (x, k) in enumerate(zip(positions, labels))]


def brute_force(rows, edges):
    counts = np.zeros((len(edges) - 1, 4), dtype=int)
    participants = [set() for _ in counts]
    for a, b in itertools.combinations(rows, 2):
        if a["classification"] not in CLASSES[:3] or b["classification"] not in CLASSES[:3]:
            continue
        distance = abs(a["anchor"] - b["anchor"])
        for index, (lower, upper) in enumerate(itertools.pairwise(edges)):
            if lower <= distance < upper:
                counts[index, 0] += 1
                participants[index].update((a["gene_id"], b["gene_id"]))
                if a["classification"] == b["classification"]:
                    counts[index, 1 + CLASSES.index(a["classification"])] += 1
    return counts, [len(p) for p in participants]


def test_hand_calculated_baseline_and_bins():
    rows = genes([0, 10, 20, 30], [0, 0, 1, 2])
    result, marks = analyze(rows, [0, 11, 21, 31, 40], block_size=10)
    np.testing.assert_array_equal(marks[0]["pair_marks"].sum(axis=0),
                                  [[3, 1, 0, 0], [2, 0, 0, 0], [1, 0, 0, 0], [0, 0, 0, 0]])
    agreement = result["reference_estimates"][0]
    assert agreement["observed"] == pytest.approx(1 / 3)
    assert agreement["baseline"] == pytest.approx(1 / 6)
    assert agreement["excess"] == pytest.approx(1 / 6)
    assert agreement["enrichment_ratio"] == 2
    assert agreement["contributing_genes"] == 4
    assert result["reference_estimates"][-1]["observed"] is None
    assert result["reference_estimates"][2]["enrichment_ratio"] is None


@pytest.mark.parametrize("seed", range(15))
def test_prefix_counts_match_independent_pair_enumeration(seed):
    rng = np.random.default_rng(seed)
    rows = genes(rng.integers(0, 200, size=35), rng.integers(0, 6, size=35))
    edges = [0, 1, 20, 70, 200]
    expected, participants = brute_force(rows, edges)
    marks = chromosome_marks(rows, edges, block_size=21)
    np.testing.assert_array_equal(marks["pair_marks"].sum(axis=0), expected)
    np.testing.assert_array_equal(marks["participants"], participants)
    assert marks["class_counts"].sum() == len(rows)
    assert np.ptp(np.diff(marks["bounds"])) <= 1


def test_exact_permutation_expectation_is_zero():
    values = []
    for labels in set(itertools.permutations([0, 0, 1, 2])):
        result, _ = analyze(genes([0, 1, 4, 20], labels), [0, 5, 21])
        values.append([r["excess"] for r in result["reference_estimates"]])
    np.testing.assert_allclose(np.mean(values, axis=0), 0, atol=1e-15)


def test_pair_weighted_baseline_not_pooled_frequencies():
    rows = genes([0, 1, 2], [0, 0, 0]) + genes([0, 1], [1, 2], "chr2", prefix="h")
    result, _ = analyze(rows, [0, 10])
    row = result["reference_estimates"][0]
    assert row["pair_count"] == 4
    assert row["contributing_genes"] == 5
    assert row["contributing_chromosomes"] == 2
    assert row["observed"] == row["baseline"] == 0.75
    assert row["excess"] == 0


def test_references_missing_singletons_and_all_unresolved():
    rows = genes([0, 1], [0, 0]) + genes([0, 1], [0, 0], reference="other")
    rows += genes([5], [1], chromosome="singleton", prefix="s")
    rows += genes([5, 9], [3, 4], chromosome="unresolved", prefix="u")
    result, _ = analyze(rows, [0, 10])
    assert [r["pair_count"] for r in result["reference_estimates"][::4]] == [1, 1]
    single = [r for r in result["chromosome_estimates"] if r["chromosome"] == "singleton"]
    assert all(r["baseline"] is None for r in single)
    assert np.isnan(baseline([1, 0, 0])).all()


def test_ties_and_exclusive_last_edge():
    result, _ = analyze(genes([0, 0, 10], [0, 1, 2]), [0, 1, 10])
    assert [r["pair_count"] for r in result["reference_estimates"][::4]] == [1, 0]


@pytest.mark.parametrize("edges", [[], [0], [1, 2], [0, 0], [0, -1], [0, 1.5], [0, True], [0, 2**63]])
def test_invalid_edges(edges):
    with pytest.raises(SystemExit):
        validate_edges(edges)


@pytest.mark.parametrize("field,value", [("anchor", 2), ("start", -1), ("end", 0),
                                        ("end", 2**63), ("anchor", 0.0),
                                        ("classification", "unknown"), ("gene_id", " "),
                                        ("chromosome", None)])
def test_invalid_rows(field, value):
    rows = genes([0], [0])
    rows[0][field] = value
    with pytest.raises(SystemExit):
        validate_rows(rows)


def test_duplicate_and_conflicting_reference_labels():
    with pytest.raises(SystemExit):
        validate_rows(genes([0], [0]) * 2)
    with pytest.raises(SystemExit):
        validate_rows(genes([0], [0]) + genes([0], [1], reference="other"))
    with pytest.raises(SystemExit):
        validate_rows([])


def test_classified_tsv(tmp_path):
    path = tmp_path / "genes.tsv"
    path.write_text("reference\tchromosome\tgene_id\tstart\tend\tanchor\tclassification\n"
                    "ref\tchr1\tg1\t0\t1\t0\ttopology_1\n")
    assert read_classified(path)[0]["anchor"] == 0
    path.write_text("gene_id\nfoo\n")
    with pytest.raises(SystemExit):
        read_classified(path)


def test_block_marks_retain_cross_boundary_pairs_without_new_adjacency():
    m = chromosome_marks(genes([0, 9, 10, 19], [0, 0, 1, 1]), [0, 2, 11], block_size=10)
    np.testing.assert_array_equal(m["bounds"], [0, 10, 20])
    # Pair (9, 10) crosses a boundary but belongs to the original left mark.
    np.testing.assert_array_equal(m["pair_marks"][:, 0, 0], [1, 0])
    duplicate_left = m["pair_marks"][0] * 2
    assert duplicate_left[0, 0] == 2
    obs, exp, excess = estimate(duplicate_left, m["class_counts"][0, :3] * 2)
    assert obs[0, 0] == 0
    assert exp[0, 0] == 1
    assert excess[0, 0] == -1


def test_empty_physical_blocks_are_retained():
    m = chromosome_marks(genes([0, 99], [0, 1]), [0, 100], block_size=10)
    assert len(m["class_counts"]) == 10
    assert np.count_nonzero(m["class_counts"].sum(axis=1)) == 2


def test_bootstrap_reproducible_and_intervals_available():
    rng = np.random.default_rng(27)
    rows = genes(np.arange(5000) * 10, rng.choice(3, 5000, p=[0.5, 0.3, 0.2]))
    a, diagnostics = uncertainty(rows, [0, 11, 21], [500, 1000], replicates=499, seed=82)
    b, _ = uncertainty(list(reversed(rows)), [0, 11, 21], [1000, 500], replicates=499, seed=82)
    assert a == b
    assert all(r["status"] == "available" for r in a)
    assert all(r["lower"] < r["upper"] for r in a)
    assert all(isinstance(r["block_size_sensitive"], bool) for r in a)
    assert sum(r["total_genes"] for r in diagnostics) == len(rows) * 2


def test_uncertainty_withholds_short_sparse_and_degenerate_blocks():
    rows = genes([0, 1, 10, 99], [0, 0, 0, 0])
    result, _ = uncertainty(rows, [0, 20, 100, 200], [10, 20], replicates=499)
    assert all(r["status"] == "withheld" for r in result)
    assert all(r["lower"] is None for r in result)
    assert "no_pairs" in result[-1]["reasons"]
    assert "fewer_than_20_pair_occupied_blocks" in result[0]["reasons"]
    assert "degenerate_bootstrap" in result[0]["reasons"]


def test_bootstrap_original_marks_manual_draw():
    m = chromosome_marks(genes([0, 9, 10, 19], [0, 0, 1, 1]), [0, 2, 11], block_size=10)

    class FixedDraw:
        def multinomial(self, n, p, size):
            assert n == 2
            return np.tile([1, 1], (size, 1))

    values = bootstrap_reference([m], 3, FixedDraw())
    _, _, expected = estimate(m["pair_marks"].sum(axis=0), m["class_counts"][:, :3].sum(axis=0))
    np.testing.assert_allclose(values, np.tile(expected, (3, 1, 1)))


@pytest.mark.parametrize("kwargs", [{"block_sizes": [10]}, {"block_sizes": [10, 10]},
                                   {"block_sizes": [-1, 10]}, {"replicates": 498}, {"seed": -1}])
def test_bad_uncertainty_options(kwargs):
    options = {"block_sizes": [10, 20], "replicates": 499, "seed": 0}
    options.update(kwargs)
    with pytest.raises(SystemExit):
        uncertainty(genes([0, 1], [0, 1]), [0, 10], **options)
