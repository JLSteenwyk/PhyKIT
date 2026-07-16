"""
Shared trait file parsing utilities.

Provides a single implementation for parsing tab-delimited multi-trait
files with a header row, used across phylogenetic regression, signal,
ordination, path analysis, ANOVA, and other comparative methods.
"""
from math import isfinite
import sys

from ..errors import PhykitUserError


def _split_tsv_line(line: str) -> list[str]:
    return [part.strip() for part in line.rstrip("\r\n").split("\t")]


def _parse_float_values_detailed(
    value_parts: list[str],
    trait_names: list[str],
    taxon: str,
    line_num: int,
) -> list[float]:
    values = []
    for i, val_str in enumerate(value_parts):
        try:
            values.append(float(val_str))
        except ValueError:
            raise PhykitUserError(
                [
                    f"Non-numeric trait value '{val_str}' for taxon '{taxon}' "
                    f"(trait '{trait_names[i]}') on line {line_num}.",
                ],
                code=2,
            )
    return values


def parse_multi_trait_file(
    path: str,
    tree_tips: list[str],
    min_shared: int = 3,
    min_columns: int = 2,
) -> tuple[list[str], dict[str, list[float]]]:
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
            header_parts = None
            for line in f:
                stripped = line.strip()
                if not stripped or stripped[0] == "#":
                    continue
                header_parts = _split_tsv_line(line)
                break

            if header_parts is None:
                raise PhykitUserError(
                    ["Multi-trait file must have a header row and at least one data row."],
                    code=2,
                )

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
            trait_count = len(trait_names)
            trait_columns = {}
            for column, trait_name in enumerate(trait_names, start=2):
                if not trait_name:
                    raise PhykitUserError(
                        [
                            f"Trait header is blank in column {column} on line 1."
                        ],
                        code=2,
                    )
                if trait_name in trait_columns:
                    raise PhykitUserError(
                        [
                            f"Duplicate trait header '{trait_name}' in columns "
                            f"{trait_columns[trait_name]} and {column} on line 1."
                        ],
                        code=2,
                    )
                trait_columns[trait_name] = column

            traits = {}
            taxon_lines = {}
            saw_data = False
            data_line_idx = 2
            to_float = float
            for line in f:
                stripped = line.strip()
                if not stripped or stripped[0] == "#":
                    continue
                saw_data = True
                parts = _split_tsv_line(line)
                if len(parts) != n_cols:
                    raise PhykitUserError(
                        [
                            f"Line {data_line_idx} has {len(parts)} columns; expected {n_cols}.",
                            f"Each line should have: taxon_name<tab>"
                            f"{'<tab>'.join(['trait'] * len(trait_names))}",
                        ],
                        code=2,
                    )
                taxon = parts[0]
                if not taxon:
                    raise PhykitUserError(
                        [f"Taxon name is blank on line {data_line_idx}."],
                        code=2,
                    )
                if taxon in taxon_lines:
                    raise PhykitUserError(
                        [
                            f"Duplicate taxon '{taxon}' on line {data_line_idx}; "
                            f"first seen on line {taxon_lines[taxon]}."
                        ],
                        code=2,
                    )
                try:
                    if trait_count == 1:
                        values = [to_float(parts[1])]
                    elif trait_count == 2:
                        values = [to_float(parts[1]), to_float(parts[2])]
                    elif trait_count == 3:
                        values = [
                            to_float(parts[1]),
                            to_float(parts[2]),
                            to_float(parts[3]),
                        ]
                    elif trait_count == 4:
                        values = [
                            to_float(parts[1]),
                            to_float(parts[2]),
                            to_float(parts[3]),
                            to_float(parts[4]),
                        ]
                    else:
                        values = [to_float(val_str) for val_str in parts[1:]]
                except ValueError:
                    values = _parse_float_values_detailed(
                        parts[1:],
                        trait_names,
                        taxon,
                        data_line_idx,
                    )
                for trait_name, raw_value, value in zip(
                    trait_names, parts[1:], values
                ):
                    if not isfinite(value):
                        raise PhykitUserError(
                            [
                                f"Non-finite trait value '{raw_value}' for taxon "
                                f"'{taxon}' (trait '{trait_name}') on line "
                                f"{data_line_idx}."
                            ],
                            code=2,
                        )
                traits[taxon] = values
                taxon_lines[taxon] = data_line_idx
                data_line_idx += 1

            if not saw_data:
                raise PhykitUserError(
                    ["Multi-trait file must have a header row and at least one data row."],
                    code=2,
                )
    except FileNotFoundError:
        raise PhykitUserError(
            [
                f"{path} corresponds to no such file or directory.",
                "Please check filename and pathing",
            ],
            code=2,
        )

    trait_row_count = len(traits)
    trait_names_in_order = None
    if trait_row_count >= min_shared and len(tree_tips) == trait_row_count:
        if tree_tips[0] == next(iter(traits)):
            if tree_tips[-1] == next(reversed(traits)):
                trait_names_in_order = list(traits)
                if trait_names_in_order == tree_tips:
                    return trait_names, traits

    tree_tip_set = set(tree_tips)
    if (
        len(tree_tip_set) >= min_shared
        and len(tree_tip_set) == trait_row_count
        and tree_tip_set == traits.keys()
    ):
        return trait_names, traits

    trait_taxa_set = (
        set(trait_names_in_order) if trait_names_in_order is not None else set(traits)
    )
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


def trait_matrix_from_rows(traits: dict[str, list[float]], ordered_names: list[str]):
    """Build a numeric matrix from already-parsed trait rows."""
    import numpy as np

    return np.asarray([traits[name] for name in ordered_names], dtype=float)


def trait_column_from_rows(
    traits: dict[str, list[float]],
    ordered_names: list[str],
    column_index: int,
):
    """Build a numeric vector for one parsed trait column."""
    import numpy as np

    return np.fromiter(
        (traits[name][column_index] for name in ordered_names),
        dtype=float,
        count=len(ordered_names),
    )


def response_predictor_arrays(
    traits: dict[str, list[float]],
    ordered_names: list[str],
    response_index: int,
    predictor_indices: list[int],
):
    """Build response and intercept design arrays from selected trait columns."""
    import numpy as np

    n_rows = len(ordered_names)
    if not predictor_indices:
        y = np.array(
            [traits[name][response_index] for name in ordered_names],
            dtype=float,
        )
        return y, np.ones((n_rows, 1), dtype=float)

    predictor_count = len(predictor_indices)
    if predictor_count <= 4:
        y = np.empty(n_rows, dtype=float)
        X = np.ones((n_rows, predictor_count + 1), dtype=float)
        if predictor_count == 1:
            pred0 = predictor_indices[0]
            for row_idx, name in enumerate(ordered_names):
                row = traits[name]
                y[row_idx] = row[response_index]
                X[row_idx, 1] = row[pred0]
        elif predictor_count == 2:
            pred0, pred1 = predictor_indices
            for row_idx, name in enumerate(ordered_names):
                row = traits[name]
                y[row_idx] = row[response_index]
                X[row_idx, 1] = row[pred0]
                X[row_idx, 2] = row[pred1]
        elif predictor_count == 3:
            pred0, pred1, pred2 = predictor_indices
            for row_idx, name in enumerate(ordered_names):
                row = traits[name]
                y[row_idx] = row[response_index]
                X[row_idx, 1] = row[pred0]
                X[row_idx, 2] = row[pred1]
                X[row_idx, 3] = row[pred2]
        else:
            pred0, pred1, pred2, pred3 = predictor_indices
            for row_idx, name in enumerate(ordered_names):
                row = traits[name]
                y[row_idx] = row[response_index]
                X[row_idx, 1] = row[pred0]
                X[row_idx, 2] = row[pred1]
                X[row_idx, 3] = row[pred2]
                X[row_idx, 4] = row[pred3]
        return y, X

    from operator import itemgetter

    getter = itemgetter(response_index, *predictor_indices)
    selected = np.asarray([getter(traits[name]) for name in ordered_names], dtype=float)
    y = selected[:, 0].copy()
    X = np.empty((n_rows, len(predictor_indices) + 1), dtype=float)
    X[:, 0] = 1.0
    X[:, 1:] = selected[:, 1:]
    return y, X


def subset_traits_to_ordered_shared_taxa(
    traits: dict,
    ordered_names: list[str],
    shared_ordered: list[str],
) -> tuple[dict, list[str]]:
    """Subset parsed trait rows when ordered shared taxa exclude some rows."""
    if shared_ordered == ordered_names:
        return traits, ordered_names
    return {taxon: traits[taxon] for taxon in shared_ordered}, shared_ordered
