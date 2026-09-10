"""Validated gene coordinates and nonoverlapping concordance neighborhoods."""

import csv
from collections import Counter, defaultdict
from pathlib import Path

from .topology_landscape import CLASSES, fail


def table_rows(path, required):
    """Read a header-bearing TSV without silently accepting malformed rows."""
    try:
        with open(path, newline="", encoding="utf-8-sig") as handle:
            rows = csv.DictReader(handle, delimiter="\t")
            fields = rows.fieldnames or []
            if len(set(fields)) != len(fields) or not set(required) <= set(fields):
                fail(f"{path}: expected unique TSV columns including {', '.join(required)}.")
            for number, row in enumerate(rows, 2):
                if None in row or any(row.get(key) is None for key in required):
                    fail(f"{path}:{number}: malformed TSV row.")
                if any(not row[key].strip() for key in required):
                    fail(f"{path}:{number}: required fields cannot be empty.")
                yield number, row
    except (OSError, UnicodeError, csv.Error) as exc:
        fail(f"Cannot read {path}: {exc}")


def read_manifest(path):
    genes = {}
    tree_owners = {}
    for number, row in table_rows(path, ("gene_id", "tree")):
        gene = row["gene_id"]
        tree = Path(row["tree"])
        if not tree.is_absolute():
            tree = Path(path).resolve().parent / tree
        tree = tree.resolve()
        if gene in genes:
            fail(f"{path}:{number}: duplicate gene identifier {gene}.")
        if tree in tree_owners:
            fail(f"{path}:{number}: tree file is assigned to multiple genes: {tree}.")
        if not tree.is_file():
            fail(f"{path}:{number}: missing gene-tree file {tree}.")
        genes[gene] = str(tree)
        tree_owners[tree] = gene
    if not genes:
        fail("Gene-tree manifest contains no genes.")
    return genes


def bed_rows(path):
    try:
        with open(path, encoding="utf-8-sig") as handle:
            for number, line in enumerate(handle, 1):
                if not line.strip() or line.startswith(("#", "track ", "browser ")):
                    continue
                fields = line.rstrip("\r\n").split("\t")
                if len(fields) < 4:
                    fail(f"{path}:{number}: BED requires at least four tab-separated columns.")
                yield number, dict(zip(("chromosome", "start", "end", "gene_id"), fields))
    except (OSError, UnicodeError) as exc:
        fail(f"Cannot read {path}: {exc}")


def read_coordinates(sources):
    """Merge repeated same-chromosome records to one gene span per reference."""
    coordinates = {}
    duplicates = Counter()
    references = []
    for reference, fmt, path in sources:
        if not reference.strip():
            fail("Reference genome identifiers cannot be empty.")
        if reference not in references:
            references.append(reference)
        if fmt not in ("bed", "tsv"):
            fail("Coordinate format must be bed or tsv.")
        rows = bed_rows(path) if fmt == "bed" else table_rows(
            path, ("gene_id", "chromosome", "start", "end")
        )
        for number, row in rows:
            gene, chromosome = row["gene_id"], row["chromosome"]
            if not gene.strip() or gene == "." or not chromosome.strip():
                fail(f"{path}:{number}: gene and chromosome identifiers are required.")
            try:
                start, end = int(row["start"]), int(row["end"])
            except ValueError:
                fail(f"{path}:{number}: coordinates must be integers.")
            if start < 0 or end <= start:
                fail(f"{path}:{number}: require 0 <= start < end (half-open coordinates).")
            key = (reference, gene)
            if key in coordinates:
                previous = coordinates[key]
                if chromosome != previous["chromosome"]:
                    fail(f"Ambiguous mapping: {gene} occurs on multiple chromosomes in {reference}.")
                start, end = min(start, previous["start"]), max(end, previous["end"])
                duplicates[reference] += 1
            coordinates[key] = {
                "reference": reference, "gene_id": gene, "chromosome": chromosome,
                "start": start, "end": end, "anchor": (start + end - 1) // 2,
            }
    if not coordinates:
        fail("Coordinate files contain no gene records.")
    return list(coordinates.values()), {
        "references": references, "merged_records": dict(duplicates),
    }


def summarize(rows):
    """Each reference/gene contributes once; fractions are among resolved genes."""
    unique = {}
    for row in rows:
        unique[(row["reference"], row["gene_id"])] = row["classification"]
    counts = Counter(unique.values())
    total = len(unique)
    informative = sum(counts[c] for c in CLASSES[:3])
    result = {"total_genes": total, "informative_genes": informative}
    for category in CLASSES:
        result[category] = counts[category]
    for category in CLASSES[:3]:
        result[f"{category}_proportion"] = counts[category] / informative if informative else None
    return result


def select_rows(rows, chromosomes=None, interval=None):
    if interval is not None:
        start, end = interval
        if start < 0 or end <= start:
            fail("Interval requires 0 <= start < end.")
    selected = [row for row in rows if (
        (not chromosomes or row["chromosome"] in chromosomes)
        and (interval is None or interval[0] <= row["anchor"] < interval[1])
    )]
    return sorted(selected, key=lambda r: (
        r["reference"], r["chromosome"], r["anchor"], r["gene_id"],
    ))


def neighborhoods(rows, window_bp=None, window_genes=None, interval=None):
    """Assign a gene by its span midpoint, once in each nonoverlapping partition.

    Physical windows start at zero, or at the selected interval's start.
    Without chromosome lengths they stop at the last observed gene end.
    Gene-count windows are rank bins; equal anchors can straddle bins.
    """
    if (window_bp is None) == (window_genes is None):
        fail("Specify exactly one of physical-distance or gene-count windows.")
    size = window_bp if window_bp is not None else window_genes
    if size <= 0:
        fail("Window size must be positive.")
    by_chromosome = defaultdict(list)
    for row in select_rows(rows, interval=interval):
        by_chromosome[(row["reference"], row["chromosome"])].append(row)
    output = []
    for (reference, chromosome), genes in sorted(by_chromosome.items()):
        if window_bp is not None:
            start = interval[0] if interval else 0
            stop = interval[1] if interval else max(g["end"] for g in genes)
            n_windows = (stop - start + size - 1) // size
            if n_windows > 1_000_000:
                fail("More than one million windows on a chromosome; increase --window-bp.")
            bins = defaultdict(list)
            for gene in genes:
                bins[(gene["anchor"] - start) // size].append(gene)
            windows = ((start + i * size, min(start + (i + 1) * size, stop), bins[i])
                       for i in range(n_windows))
        else:
            chunks = (genes[i:i + size] for i in range(0, len(genes), size))
            windows = ((chunk[0]["anchor"], chunk[-1]["anchor"] + 1, chunk)
                       for chunk in chunks)
        for index, (window_start, window_end, members) in enumerate(windows, 1):
            output.append(dict(
                reference=reference, chromosome=chromosome, window_index=index,
                start=window_start, end=window_end, window_mode="bp" if window_bp else "genes",
                **summarize(members),
            ))
    return output
