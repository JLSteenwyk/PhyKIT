"""Distance-binned categorical agreement and chromosome-stratified block marks."""

import itertools
from collections import Counter, defaultdict

import numpy as np

from .genomic_coordinates import select_rows, table_rows
from .topology_landscape import CLASSES, fail

METRICS = ("agreement", *CLASSES[:3])


def validate_edges(edges):
    if len(edges) < 2 or edges[0] != 0 or any(
        isinstance(x, bool) or not isinstance(x, (int, np.integer)) for x in edges
    ) or any(b <= a for a, b in itertools.pairwise(edges)):
        fail("Distance edges must be strictly increasing integers starting at zero.")
    if len(edges) > 1001 or edges[-1] > 2**62:
        fail("Use at most 1000 distance bins with maximum distance <= 2**62.")
    return np.asarray(edges, dtype=np.int64)


def validate_rows(rows):
    """Reject duplicate mappings and inconsistent cross-reference labels."""
    seen, labels = set(), {}
    output = []
    for row in rows:
        row = dict(row)
        for key in ("reference", "chromosome", "gene_id"):
            if not isinstance(row.get(key), str) or not row[key].strip() or row[key] == ".":
                fail(f"Invalid or missing {key} in classified genes.")
        key = (row["reference"], row["gene_id"])
        if key in seen:
            fail(f"Duplicate classified gene mapping: {key}.")
        seen.add(key)
        category = row.get("classification")
        if category not in CLASSES:
            fail(f"Unknown topology classification: {category}.")
        if row["gene_id"] in labels and labels[row["gene_id"]] != category:
            fail(f"Conflicting classifications across references: {row['gene_id']}.")
        labels[row["gene_id"]] = category
        for field in ("start", "end", "anchor"):
            try:
                value = row[field]
                if isinstance(value, (bool, float)):
                    raise TypeError
                row[field] = int(value)
            except (KeyError, TypeError, ValueError):
                fail(f"Classified gene {field} must be an integer.")
        if not 0 <= row["start"] < row["end"] <= 2**62:
            fail("Classified coordinates require 0 <= start < end <= 2**62.")
        if row["anchor"] != (row["start"] + row["end"] - 1) // 2:
            fail("Classified gene anchor does not match the gene-span midpoint.")
        output.append(row)
    if not output:
        fail("No classified genes provided.")
    return select_rows(output)


def read_classified(path):
    required = ("reference", "chromosome", "gene_id", "start", "end", "anchor", "classification")
    return validate_rows([row for _, row in table_rows(path, required)])


def baseline(counts):
    """Finite-population joint probabilities, with NaNs for fewer than two genes."""
    counts = np.asarray(counts, dtype=float)
    n = counts.sum(axis=-1)
    denominator = n * (n - 1)
    joint = np.divide(counts * (counts - 1), denominator[..., None],
                      out=np.full_like(counts, np.nan), where=denominator[..., None] > 0)
    return np.concatenate((joint.sum(axis=-1, keepdims=True), joint), axis=-1)


def chromosome_marks(genes, edges, block_size=None):
    """Compute left-endpoint marks without enumerating or storing gene pairs.

    Each bin uses range searches and topology prefix sums. Difference arrays
    identify participating right endpoints, avoiding materialized pair lists.
    """
    edges = validate_edges(edges)
    genes = sorted(genes, key=lambda r: (r["anchor"], r["gene_id"]))
    if not genes:
        fail("Cannot analyze an empty chromosome.")
    start, stop = genes[0]["anchor"], genes[-1]["anchor"] + 1
    span = stop - start
    if block_size is not None and (not isinstance(block_size, int) or block_size <= 0):
        fail("Block sizes must be positive integers.")
    nblocks = max(1, span // block_size) if block_size else 1
    if nblocks > 100_000 or nblocks * (len(edges) - 1) > 2_000_000:
        fail("Too many block/bin combinations; increase block size or reduce bins.")
    # Integer quotient/remainder partitions avoid coordinate multiplication overflow.
    width, remainder = divmod(span, nblocks)
    bounds = start + np.arange(nblocks + 1, dtype=np.int64) * width
    bounds += np.minimum(np.arange(nblocks + 1), remainder)
    all_x = np.asarray([g["anchor"] for g in genes], dtype=np.int64)
    all_block = np.searchsorted(bounds, all_x, side="right") - 1
    classes = np.asarray([CLASSES.index(g["classification"]) for g in genes])
    class_counts = np.zeros((nblocks, 6), dtype=np.int64)
    np.add.at(class_counts, (all_block, classes), 1)
    resolved = classes < 3
    x, y, blocks = all_x[resolved], classes[resolved], all_block[resolved]
    n, nbins = len(x), len(edges) - 1
    pair_marks = np.zeros((nblocks, nbins, 4), dtype=np.int64)
    participants = np.zeros(nbins, dtype=np.int64)
    prefix = np.zeros((n + 1, 3), dtype=np.int64)
    if n:
        prefix[1:] = np.cumsum(np.eye(3, dtype=np.int64)[y], axis=0)
    indices = np.arange(n)
    for b, (lower, upper) in enumerate(itertools.pairwise(edges)):
        left = np.maximum(indices + 1, np.searchsorted(x, x + lower, side="left"))
        right = np.maximum(left, np.searchsorted(x, x + upper, side="left"))
        counts = right - left
        np.add.at(pair_marks[:, b, 0], blocks, counts)
        for k in range(3):
            joint = (prefix[right, k] - prefix[left, k]) * (y == k)
            np.add.at(pair_marks[:, b, k + 1], blocks, joint)
        difference = np.zeros(n + 1, dtype=np.int64)
        np.add.at(difference, left, 1)
        np.add.at(difference, right, -1)
        participants[b] = np.count_nonzero((counts > 0) | (np.cumsum(difference[:-1]) > 0))
    return {
        "reference": genes[0]["reference"], "chromosome": genes[0]["chromosome"],
        "bounds": bounds, "class_counts": class_counts, "pair_marks": pair_marks,
        "participants": participants,
    }


def estimate(pair_counts, gene_counts):
    """Return observed, baseline, and excess arrays with final axis METRICS."""
    pair_counts = np.asarray(pair_counts, dtype=float)
    totals = pair_counts[..., :1]
    successes = np.concatenate((pair_counts[..., 1:].sum(axis=-1, keepdims=True),
                                pair_counts[..., 1:]), axis=-1)
    observed = np.divide(successes, totals, out=np.full_like(successes, np.nan), where=totals > 0)
    expected = np.broadcast_to(np.expand_dims(baseline(gene_counts), -2), observed.shape).copy()
    expected = np.where(totals > 0, expected, np.nan)
    return observed, expected, observed - expected


def _number(value):
    return float(value) if np.isfinite(value) else None


def records(reference, chromosome, pairs, observed, expected, participants, contributors, edges):
    output = []
    for b, (lower, upper) in enumerate(itertools.pairwise(edges)):
        for k, metric in enumerate(METRICS):
            output.append({
                "reference": reference, "chromosome": chromosome, "distance_start": int(lower),
                "distance_end": int(upper), "metric": metric, "pair_count": int(pairs[b]),
                "contributing_genes": int(participants[b]),
                "contributing_chromosomes": int(contributors[b]),
                "observed": _number(observed[b, k]), "baseline": _number(expected[b, k]),
                "excess": _number(observed[b, k] - expected[b, k]),
                "enrichment_ratio": _number(observed[b, k] / expected[b, k])
                if expected[b, k] > 0 else None,
            })
    return output


def analyze(rows, edges, block_size=None):
    """Chromosome and pair-weighted reference estimates, with reusable block marks."""
    rows, edges = validate_rows(rows), validate_edges(edges)
    grouped = defaultdict(list)
    for row in rows:
        grouped[(row["reference"], row["chromosome"])].append(row)
    marks, chromosome_rows, reference_rows, diagnostics = [], [], [], []
    pooled = {}
    for (ref, chrom), genes in sorted(grouped.items()):
        m = chromosome_marks(genes, edges, block_size)
        marks.append(m)
        counts = m["pair_marks"].sum(axis=0)
        obs, exp, _ = estimate(counts, m["class_counts"][:, :3].sum(axis=0))
        pairs = counts[:, 0]
        chromosome_rows.extend(records(ref, chrom, pairs, obs, exp, m["participants"], pairs > 0, edges))
        if ref not in pooled:
            pooled[ref] = [np.zeros_like(counts), np.zeros_like(exp),
                           np.zeros_like(pairs), np.zeros_like(pairs)]
        total, weighted, participants, contributors = pooled[ref]
        total += counts
        weighted += np.nan_to_num(exp) * pairs[:, None]
        participants += m["participants"]
        contributors += pairs > 0
        c = Counter(g["classification"] for g in genes)
        diagnostics.append({"reference": ref, "chromosome": chrom, "total_genes": len(genes),
                                "resolved_genes": sum(c[k] for k in CLASSES[:3]),
                                "classifications": {k: c[k] for k in CLASSES},
                                "first_anchor": genes[0]["anchor"], "last_anchor": genes[-1]["anchor"]})
    for ref, (counts, weighted, participants, contributors) in sorted(pooled.items()):
        pairs = counts[:, 0]
        successes = np.column_stack((counts[:, 1:].sum(axis=1), counts[:, 1:]))
        obs = np.divide(successes, pairs[:, None], out=np.full_like(weighted, np.nan),
                        where=pairs[:, None] > 0)
        exp = np.divide(weighted, pairs[:, None], out=np.full_like(weighted, np.nan),
                        where=pairs[:, None] > 0)
        reference_rows.extend(records(ref, None, pairs, obs, exp, participants, contributors, edges))
    return {"reference_estimates": reference_rows, "chromosome_estimates": chromosome_rows,
            "chromosome_diagnostics": diagnostics}, marks


def bootstrap_reference(marks, replicates, rng):
    """Resample original marks within chromosomes; never assemble new genomes."""
    nbins = marks[0]["pair_marks"].shape[1]
    draws = np.full((replicates, nbins, 4), np.nan)
    for first in range(0, replicates, 100):
        size = min(100, replicates - first)
        totals = np.zeros((size, nbins, 1))
        successes = np.zeros((size, nbins, 4))
        expectations = np.zeros_like(successes)
        valid = np.ones((size, nbins, 1), dtype=bool)
        for m in marks:
            blocks = len(m["class_counts"])
            weights = rng.multinomial(blocks, np.full(blocks, 1 / blocks), size=size)
            counts = weights @ m["class_counts"][:, :3]
            pairs = (weights @ m["pair_marks"].reshape(blocks, -1)).reshape(size, nbins, 4)
            q = baseline(counts)[:, None, :]
            denominator = pairs[:, :, :1]
            valid &= (counts.sum(axis=1)[:, None, None] >= 2) | (denominator == 0)
            totals += denominator
            successes += np.concatenate((pairs[:, :, 1:].sum(axis=2, keepdims=True),
                                         pairs[:, :, 1:]), axis=2)
            expectations += np.nan_to_num(q) * denominator
        np.divide(successes - expectations, totals, out=draws[first:first + size],
                  where=(totals > 0) & valid)
    return draws


def uncertainty(rows, edges, block_sizes, replicates=999, seed=0):
    """Pointwise marked-block intervals and explicit block-size sensitivity."""
    edges = validate_edges(edges)
    if (not isinstance(seed, int) or seed < 0 or not isinstance(replicates, int)
            or replicates < 499 or replicates > 100_000):
        fail("Inference requires a nonnegative integer seed and 499-100000 replicates.")
    if (len(set(block_sizes)) < 2 or any(not isinstance(s, int) or s <= 0 for s in block_sizes)):
        fail("Inference requires at least two distinct positive integer block sizes.")
    if replicates * (len(edges) - 1) > 5_000_000:
        fail("Too many replicate/bin combinations; reduce replicates or bins.")
    output, block_diagnostics = [], []
    for block_size in sorted(set(block_sizes)):
        _, all_marks = analyze(rows, edges, block_size)
        by_reference = defaultdict(list)
        for m in all_marks:
            by_reference[m["reference"]].append(m)
            for b, counts in enumerate(m["class_counts"]):
                block_diagnostics.append({
                    "reference": m["reference"], "chromosome": m["chromosome"],
                    "requested_block_bp": block_size, "block_index": b + 1,
                    "start": int(m["bounds"][b]), "end": int(m["bounds"][b + 1]),
                    "total_genes": int(counts.sum()), "resolved_genes": int(counts[:3].sum()),
                    **{name: int(counts[k]) for k, name in enumerate(CLASSES)},
                })
        for ref_index, (reference, marks) in enumerate(sorted(by_reference.items())):
            rng = np.random.default_rng(np.random.SeedSequence([seed, block_size, ref_index]))
            draws = bootstrap_reference(marks, replicates, rng)
            for b, (lower, upper) in enumerate(itertools.pairwise(edges)):
                contributors = [m for m in marks if m["pair_marks"][:, b, 0].sum() > 0]
                common = []
                if not contributors:
                    common.append("no_pairs")
                if any(np.count_nonzero(m["pair_marks"][:, b, 0]) < 20 for m in contributors):
                    common.append("fewer_than_20_pair_occupied_blocks")
                if any(np.min(np.diff(m["bounds"])) < 5 * int(upper) for m in contributors):
                    common.append("blocks_shorter_than_5x_bin_distance")
                for k, metric in enumerate(METRICS):
                    values = draws[:, b, k]
                    values = values[np.isfinite(values)]
                    reasons = list(common)
                    if len(values) < 0.95 * replicates:
                        reasons.append("too_few_valid_replicates")
                    if len(values) and np.ptp(values) < 1e-12:
                        reasons.append("degenerate_bootstrap")
                    low, high = np.quantile(values, [0.025, 0.975]) if not reasons else (np.nan, np.nan)
                    output.append({
                        "reference": reference, "distance_start": int(lower), "distance_end": int(upper),
                        "metric": metric, "requested_block_bp": block_size,
                        "replicates": replicates, "valid_replicates": len(values), "seed": seed,
                        "confidence_level": 0.95, "interval_type": "pointwise_percentile",
                        "lower": _number(low), "upper": _number(high),
                        "status": "withheld" if reasons else "available",
                        "reasons": reasons, "block_size_sensitive": None,
                    })
    grouped = defaultdict(list)
    for row in output:
        grouped[(row["reference"], row["distance_start"], row["metric"])].append(row)
    for group in grouped.values():
        available = [r for r in group if r["status"] == "available"]
        if len(available) >= 2:
            widths = [r["upper"] - r["lower"] for r in available]
            sensitive = max(widths) > 2 * min(widths)
            for row in group:
                row["block_size_sensitive"] = sensitive
        else:
            for row in group:
                row["reasons"].append("insufficient_admissible_block_sizes_for_sensitivity")
    return output, block_diagnostics
