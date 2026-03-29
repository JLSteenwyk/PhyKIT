"""
Shared trait file parsing utilities.

Provides a single implementation for parsing tab-delimited multi-trait
files with a header row, used across phylogenetic regression, signal,
ordination, path analysis, ANOVA, and other comparative methods.
"""
import sys
from typing import Dict, List, Tuple

from ..errors import PhykitUserError


def parse_multi_trait_file(
    path: str,
    tree_tips: List[str],
    min_shared: int = 3,
    min_columns: int = 2,
) -> Tuple[List[str], Dict[str, List[float]]]:
    """Parse a tab-delimited multi-trait file with a header row.

    Format:
        taxon<tab>trait1<tab>trait2<tab>...
        species_A<tab>1.2<tab>3.4<tab>...

    Parameters
    ----------
    path : path to TSV file
    tree_tips : list of tip names from the tree
    min_shared : minimum shared taxa between tree and file (default 3)
    min_columns : minimum columns in header (default 2: taxon + 1 trait)

    Returns
    -------
    (trait_names, traits_dict) where traits_dict maps taxon -> [float values]
    """
    try:
        with open(path) as f:
            lines = f.readlines()
    except FileNotFoundError:
        raise PhykitUserError(
            [
                f"{path} corresponds to no such file or directory.",
                "Please check filename and pathing",
            ],
            code=2,
        )

    # Filter out comments and blank lines
    data_lines = []
    for line in lines:
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        data_lines.append(stripped)

    if len(data_lines) < 2:
        raise PhykitUserError(
            ["Multi-trait file must have a header row and at least one data row."],
            code=2,
        )

    # First data line is the header
    header_parts = data_lines[0].split("\t")
    n_cols = len(header_parts)
    if n_cols < min_columns:
        raise PhykitUserError(
            [
                f"Header must have at least {min_columns} columns "
                f"(taxon + at least {min_columns - 1} trait(s)).",
            ],
            code=2,
        )
    trait_names = header_parts[1:]

    traits = {}
    for line_idx, line in enumerate(data_lines[1:], 2):
        parts = line.split("\t")
        if len(parts) != n_cols:
            raise PhykitUserError(
                [
                    f"Line {line_idx} has {len(parts)} columns; expected {n_cols}.",
                    f"Each line should have: taxon_name<tab>"
                    f"{'<tab>'.join(['trait'] * len(trait_names))}",
                ],
                code=2,
            )
        taxon = parts[0]
        values = []
        for i, val_str in enumerate(parts[1:]):
            try:
                values.append(float(val_str))
            except ValueError:
                raise PhykitUserError(
                    [
                        f"Non-numeric trait value '{val_str}' for taxon '{taxon}' "
                        f"(trait '{trait_names[i]}') on line {line_idx}.",
                    ],
                    code=2,
                )
        traits[taxon] = values

    tree_tip_set = set(tree_tips)
    trait_taxa_set = set(traits.keys())
    shared = tree_tip_set & trait_taxa_set

    tree_only = tree_tip_set - trait_taxa_set
    trait_only = trait_taxa_set - tree_tip_set

    if tree_only:
        print(
            f"Warning: {len(tree_only)} taxa in tree but not in trait file: "
            f"{', '.join(sorted(tree_only))}",
            file=sys.stderr,
        )
    if trait_only:
        print(
            f"Warning: {len(trait_only)} taxa in trait file but not in tree: "
            f"{', '.join(sorted(trait_only))}",
            file=sys.stderr,
        )

    if len(shared) < min_shared:
        raise PhykitUserError(
            [
                f"Only {len(shared)} shared taxa between tree and trait file.",
                f"At least {min_shared} shared taxa are required.",
            ],
            code=2,
        )

    filtered = {taxon: traits[taxon] for taxon in shared}
    return trait_names, filtered
