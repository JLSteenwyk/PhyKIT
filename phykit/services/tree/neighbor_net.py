from __future__ import annotations

import math
import sys
from pathlib import Path

from ..alignment._fasta import read_fasta_first_token_upper_with_taxa
from ...errors import PhykitUserError


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


def _position_extent(positions):
    iterator = iter(positions)
    try:
        min_x, min_y = next(iterator)
    except StopIteration:
        return 0
    max_x = min_x
    max_y = min_y
    for x, y in iterator:
        if x < min_x:
            min_x = x
        elif x > max_x:
            max_x = x
        if y < min_y:
            min_y = y
        elif y > max_y:
            max_y = y
    x_extent = max_x - min_x
    y_extent = max_y - min_y
    return x_extent if x_extent >= y_extent else y_extent


class _LazyNumpy:
    def __getattr__(self, name):
        import numpy as _np

        return getattr(_np, name)


np = _LazyNumpy()

_DISTANCE_SKIP_BYTES = b"-NX?"
_NO_SKIP_DIRECT_MIN_LENGTH = 2048
_NO_SKIP_DIRECT_MAX_TAXA = 700


class NeighborNet:
    """Construct a NeighborNet phylogenetic network from a distance matrix.

    Accepts either a FASTA alignment (from which pairwise distances are
    computed) or a pre-computed CSV distance matrix.  Produces a circular
    split system visualized as a planar splits graph.

    The implementation follows a practical approximation of Bryant & Moulton
    (2004): a Neighbor-Joining tree is built to obtain a circular ordering,
    circular splits compatible with that ordering are enumerated, and split
    weights are estimated via non-negative least squares (NNLS).
    """

    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        self.alignment_path = parsed["alignment_path"]
        self.distance_matrix_path = parsed["distance_matrix_path"]
        self.output_path = parsed["output_path"]
        self.metric = parsed["metric"]
        self.max_splits = parsed["max_splits"]
        self.json_output = parsed["json_output"]
        self.plot_config = parsed["plot_config"]

    def process_args(self, args) -> dict:
        from ...helpers.plot_config import PlotConfig

        aln = getattr(args, "alignment", None)
        dm = getattr(args, "distance_matrix", None)
        if aln is None and dm is None:
            raise PhykitUserError(
                ["Either -a/--alignment or --distance-matrix is required."],
                code=2,
            )
        return dict(
            alignment_path=aln,
            distance_matrix_path=dm,
            output_path=args.output,
            metric=getattr(args, "metric", "p-distance"),
            max_splits=getattr(args, "max_splits", 30),
            json_output=getattr(args, "json", False),
            plot_config=PlotConfig.from_args(args),
        )

    # ------------------------------------------------------------------
    # Distance computation
    # ------------------------------------------------------------------

    def _compute_distances_from_alignment(self) -> tuple[list[str], np.ndarray]:
        """Compute pairwise distance matrix from a FASTA alignment."""
        try:
            taxa, seqs = self._read_alignment(self.alignment_path)
        except FileNotFoundError:
            raise PhykitUserError(
                [
                    f"{self.alignment_path} corresponds to no such file or directory.",
                    "Please check filename and pathing",
                ],
                code=2,
            )

        if len(taxa) == 0:
            raise PhykitUserError(
                ["No sequences found in the alignment file."], code=2
            )

        return taxa, self._compute_distance_matrix_from_sequences(taxa, seqs)

    @staticmethod
    def _read_alignment(path: str) -> tuple[list[str], dict[str, str]]:
        return read_fasta_first_token_upper_with_taxa(path)

    def _compute_distance_matrix_from_sequences(
        self, taxa: list[str], seqs: dict[str, str]
    ) -> np.ndarray:
        lengths = {len(seqs[taxon]) for taxon in taxa}
        if len(lengths) == 1:
            try:
                return self._compute_distance_matrix_from_equal_length_sequences(
                    taxa, seqs
                )
            except UnicodeEncodeError:
                pass
        return self._compute_distance_matrix_from_sequences_legacy(taxa, seqs)

    def _compute_distance_matrix_from_equal_length_sequences(
        self, taxa: list[str], seqs: dict[str, str]
    ) -> np.ndarray:
        n_taxa = len(taxa)
        aln_len = len(seqs[taxa[0]]) if taxa else 0
        seq_matrix = np.frombuffer(
            "".join(seqs[taxon] for taxon in taxa).encode("ascii"),
            dtype=np.uint8,
        ).reshape(n_taxa, aln_len)
        skip = np.frombuffer(_DISTANCE_SKIP_BYTES, dtype=np.uint8)
        valid = ~np.isin(seq_matrix, skip)
        if (
            aln_len >= _NO_SKIP_DIRECT_MIN_LENGTH
            and n_taxa <= _NO_SKIP_DIRECT_MAX_TAXA
            and bool(valid.all())
        ):
            return self._compute_distance_matrix_from_no_skip_matrix(seq_matrix)

        valid_float = valid.astype(np.float64, copy=False)
        compared = valid_float @ valid_float.T

        matches = np.zeros_like(compared, dtype=np.float64)
        valid_symbols = np.unique(seq_matrix[valid])
        for symbol in valid_symbols:
            symbol_mask = seq_matrix == symbol
            symbol_mask &= valid
            symbol_float = symbol_mask.astype(np.float64, copy=False)
            matches += symbol_float @ symbol_float.T

        D = np.zeros(compared.shape, dtype=np.float64)
        np.divide(
            matches,
            compared,
            out=D,
            where=compared > 0,
        )
        if self.metric == "identity":
            np.fill_diagonal(D, 0.0)
            return D

        D = 1.0 - D
        if self.metric == "jc":
            jc = np.full_like(D, 10.0)
            unsaturated = D < 0.75
            jc[unsaturated] = -0.75 * np.log(1 - 4 / 3 * D[unsaturated])
            D = jc

        np.fill_diagonal(D, 0.0)
        return D

    def _compute_distance_matrix_from_no_skip_matrix(self, seq_matrix) -> np.ndarray:
        n_taxa, aln_len = seq_matrix.shape
        target_bytes = 16 * 1024 * 1024
        block_size = max(
            1,
            min(64, target_bytes // max(1, n_taxa * max(1, aln_len))),
        )
        identities = np.empty((n_taxa, n_taxa), dtype=np.float64)
        denominator = float(aln_len)
        for start in range(0, n_taxa, block_size):
            stop = min(n_taxa, start + block_size)
            match_counts = (
                seq_matrix[start:stop, None, :] == seq_matrix[None, :, :]
            ).sum(axis=2, dtype=np.int32)
            identities[start:stop] = match_counts / denominator

        if self.metric == "identity":
            np.fill_diagonal(identities, 0.0)
            return identities

        distances = 1.0 - identities
        if self.metric == "jc":
            jc = np.full_like(distances, 10.0)
            unsaturated = distances < 0.75
            jc[unsaturated] = -0.75 * np.log(
                1 - 4 / 3 * distances[unsaturated]
            )
            distances = jc

        np.fill_diagonal(distances, 0.0)
        return distances

    def _compute_distance_matrix_from_sequences_legacy(
        self, taxa: list[str], seqs: dict[str, str]
    ) -> np.ndarray:
        n = len(taxa)
        D = np.zeros((n, n))

        skip = {"-", "N", "X", "?", "n", "x"}

        for i in range(n):
            for j in range(i + 1, n):
                si = seqs[taxa[i]]
                sj = seqs[taxa[j]]
                matches = compared = 0
                for a, b in zip(si, sj):
                    if a in skip or b in skip:
                        continue
                    compared += 1
                    if a == b:
                        matches += 1
                if compared > 0:
                    p = 1 - matches / compared  # p-distance
                else:
                    p = 1.0

                if self.metric == "identity":
                    D[i, j] = D[j, i] = 1 - p
                elif self.metric == "jc":
                    # Jukes-Cantor correction
                    if p < 0.75:
                        D[i, j] = D[j, i] = -0.75 * np.log(1 - 4 / 3 * p)
                    else:
                        D[i, j] = D[j, i] = 10.0  # saturated
                else:
                    # p-distance (default)
                    D[i, j] = D[j, i] = p

        return D

    def _read_distance_matrix(self) -> tuple[list[str], np.ndarray]:
        """Read a pre-computed distance matrix from a CSV file."""
        path = Path(self.distance_matrix_path)
        if not path.exists():
            raise PhykitUserError(
                [
                    f"{self.distance_matrix_path} corresponds to no such file or directory.",
                    "Please check filename and pathing",
                ],
                code=2,
            )

        with open(path) as fh:
            header_line = fh.readline()
            if not header_line.strip() or '"' in header_line:
                return self._read_distance_matrix_csv(path)

            header = header_line.rstrip("\r\n").split(",")
            # First cell may be empty or a label
            first_cell = header[0].strip()
            if first_cell == "" or first_cell.lower() in (
                "taxon", "taxa", "label", "name",
            ):
                taxa = [h.strip() for h in header[1:]]
            else:
                taxa = [first_cell] + [h.strip() for h in header[1:]]

            n = len(taxa)
            D = np.zeros((n, n))
            for i, line in enumerate(fh):
                if '"' in line:
                    return self._read_distance_matrix_csv(path)
                stripped = line.strip()
                if not stripped:
                    parsed = np.fromstring("", sep=",", dtype=float)
                else:
                    field_count = stripped.count(",") + 1
                    if field_count > n:
                        values_text = stripped.split(",", 1)[1]
                        expected_count = field_count - 1
                    else:
                        values_text = stripped
                        expected_count = field_count
                    parsed = np.fromstring(values_text, sep=",", dtype=float)
                    if len(parsed) != expected_count:
                        return self._read_distance_matrix_csv(path)
                D[i, :len(parsed)] = parsed

        return taxa, D

    @staticmethod
    def _read_distance_matrix_csv(path) -> tuple[list[str], np.ndarray]:
        """Read distance matrices with Python's CSV parser fallback."""
        import csv

        with open(path) as fh:
            reader = csv.reader(fh)
            header = next(reader)
            # First cell may be empty or a label
            first_cell = header[0].strip()
            if first_cell == "" or first_cell.lower() in (
                "taxon", "taxa", "label", "name",
            ):
                taxa = [h.strip() for h in header[1:]]
            else:
                taxa = [first_cell] + [h.strip() for h in header[1:]]

            n = len(taxa)
            D = np.zeros((n, n))
            for i, row in enumerate(reader):
                # Row may start with a taxon label
                values = row[1:] if len(row) > n else row
                parsed = [float(val.strip()) for val in values]
                D[i, :len(parsed)] = parsed

        return taxa, D

    # ------------------------------------------------------------------
    # NJ circular ordering
    # ------------------------------------------------------------------

    @staticmethod
    def _terminal_names_direct(tree) -> list[str] | None:
        """Return terminal names from a standard Bio.Phylo tree, preserving order."""
        try:
            root = tree.root
            root.clades
        except AttributeError:
            return None

        names = []
        stack = [root]
        while stack:
            clade = stack.pop()
            children = clade.clades
            child_count = len(children)
            if child_count == 0:
                names.append(clade.name)
            elif child_count == 2:
                stack.append(children[1])
                stack.append(children[0])
            elif child_count == 1:
                stack.append(children[0])
            else:
                for child in reversed(children):
                    stack.append(child)

        return names

    @staticmethod
    def _nj_circular_ordering(D: np.ndarray, taxa: list[str]) -> list[str]:
        """Build a NJ tree and extract a circular ordering from leaf traversal."""
        from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

        n = len(taxa)
        # Convert to BioPython DistanceMatrix (lower triangle including diagonal)
        matrix = []
        for i, row in enumerate(D):
            matrix.append(row[:i + 1].tolist())

        dm = DistanceMatrix(taxa, matrix)
        constructor = DistanceTreeConstructor()
        nj_tree = constructor.nj(dm)

        # Extract leaf ordering from the tree
        ordering = NeighborNet._terminal_names_direct(nj_tree)
        if ordering is None:
            ordering = [tip.name for tip in nj_tree.get_terminals()]
        return ordering

    # ------------------------------------------------------------------
    # Circular splits
    # ------------------------------------------------------------------

    @staticmethod
    def _enumerate_circular_splits(ordering: list[str]) -> list[frozenset]:
        """Enumerate all non-trivial circular splits compatible with the ordering.

        A circular split divides the circular ordering into two contiguous arcs.
        We exclude trivial splits (single taxon or all-but-one).
        """
        n = len(ordering)
        all_taxa = frozenset(ordering[i] for i in range(n))
        seen = set()
        splits = []
        for start in range(n):
            side_set = {ordering[start], ordering[(start + 1) % n]}
            for size in range(2, n - 1):
                complement_size = n - size
                # Canonical form: smaller side (or lexicographic tiebreak)
                if size < complement_size:
                    canonical = frozenset(side_set)
                elif size > complement_size:
                    canonical = all_taxa - side_set
                else:
                    side = frozenset(side_set)
                    complement = all_taxa - side_set
                    # Disjoint equal-size sides differ at their smallest
                    # sorted element, so this preserves the sorted-list
                    # tiebreak without sorting both halves.
                    if min(side) < min(complement):
                        canonical = side
                    else:
                        canonical = complement
                if canonical not in seen:
                    seen.add(canonical)
                    splits.append(canonical)
                if size != n - 2:
                    side_set.add(ordering[(start + size) % n])
        return splits

    # ------------------------------------------------------------------
    # NNLS weight estimation
    # ------------------------------------------------------------------

    @staticmethod
    def _estimate_split_weights(
        D: np.ndarray,
        splits: list[frozenset],
        ordering: list[str],
        taxa: list[str],
    ) -> np.ndarray:
        """Estimate non-negative split weights via NNLS."""
        from scipy.optimize import nnls

        n = len(taxa)
        n_splits = len(splits)
        if n_splits == 0:
            return np.zeros(0, dtype=float)

        membership = np.zeros((n, n_splits), dtype=bool)
        taxon_to_idx = {name: i for i, name in enumerate(taxa)}
        for split_idx, split in enumerate(splits):
            for taxon in split:
                idx = taxon_to_idx.get(taxon)
                if idx is not None:
                    membership[idx, split_idx] = True

        pair_i, pair_j = np.triu_indices(n, k=1)
        A = np.logical_xor(membership[pair_i], membership[pair_j]).astype(
            np.float64,
            copy=False,
        )
        b = D[pair_i, pair_j]

        weights, _residual = nnls(A, b)
        return weights

    # ------------------------------------------------------------------
    # Visualization (reuses patterns from ConsensusNetwork)
    # ------------------------------------------------------------------

    @staticmethod
    def _canonical_split(taxa_side: frozenset, all_taxa: frozenset) -> frozenset:
        taxa_side_len = len(taxa_side)
        all_taxa_len = len(all_taxa)
        if taxa_side_len * 2 < all_taxa_len:
            return taxa_side
        complement = all_taxa - taxa_side
        if taxa_side_len * 2 > all_taxa_len:
            return complement
        if not taxa_side:
            return taxa_side
        if min(taxa_side) < min(complement):
            return taxa_side
        return complement

    @staticmethod
    def _circular_gap_positions(split, ordering):
        """Return the two circular boundary gap positions, or None."""
        n = len(ordering)
        if n == 0:
            return None
        gap_positions = []
        first_in_split = ordering[0] in split
        previous_in_split = first_in_split
        for i in range(1, n):
            current_in_split = ordering[i] in split
            if previous_in_split != current_in_split:
                gap_positions.append(i - 1)
                if len(gap_positions) > 2:
                    return None
            previous_in_split = current_in_split
        if previous_in_split != first_in_split:
            gap_positions.append(n - 1)
        if len(gap_positions) == 2:
            return (gap_positions[0], gap_positions[1])
        return None

    @staticmethod
    def _is_circular_split(split, ordering):
        """Check if a split has exactly 2 boundary gaps in the circular ordering."""
        return NeighborNet._circular_gap_positions(split, ordering) is not None

    @staticmethod
    def _compute_split_directions(ordering, circular_splits,
                                  gap_positions_by_split=None):
        """Compute 2D direction vectors for each circular split."""
        n = len(ordering)
        angles = {taxon: 2 * math.pi * i / n for i, taxon in enumerate(ordering)}
        directions = {}
        for split, weight in circular_splits:
            if gap_positions_by_split is None:
                gap_positions = NeighborNet._circular_gap_positions(split, ordering)
            else:
                gap_positions = gap_positions_by_split.get(split)
            if gap_positions is None:
                continue
            g1 = math.pi * (2 * gap_positions[0] + 1) / n
            g2 = math.pi * (2 * gap_positions[1] + 1) / n
            cx = math.cos(g2) - math.cos(g1)
            cy = math.sin(g2) - math.sin(g1)
            dx = -cy
            dy = cx
            length = math.sqrt(dx * dx + dy * dy)
            if length > 1e-10:
                dx /= length
                dy /= length
            else:
                dx, dy = 1.0, 0.0
            cx_split = sum(math.cos(angles[t]) for t in split) / len(split)
            cy_split = sum(math.sin(angles[t]) for t in split) / len(split)
            if dx * cx_split + dy * cy_split < 0:
                dx = -dx
                dy = -dy
            directions[split] = (dx, dy)
        return directions

    @staticmethod
    def _build_splits_graph(circular_splits, all_taxa):
        """Build the Buneman splits graph: valid sign vectors and edges."""
        splits_list = [s[0] for s in circular_splits]
        weights = [s[1] for s in circular_splits]
        n_splits = len(splits_list)
        if n_splits == 0:
            return {}, set(), [], [], []
        taxon_signs = {}
        for taxon in all_taxa:
            signs = tuple(1 if taxon in sp else -1 for sp in splits_list)
            taxon_signs[taxon] = signs
        forbidden = {}
        for i in range(n_splits):
            for j in range(i + 1, n_splits):
                si_pos = splits_list[i]
                si_neg = all_taxa - si_pos
                sj_pos = splits_list[j]
                sj_neg = all_taxa - sj_pos
                fb = set()
                if not (si_pos & sj_pos):
                    fb.add((1, 1))
                if not (si_pos & sj_neg):
                    fb.add((1, -1))
                if not (si_neg & sj_pos):
                    fb.add((-1, 1))
                if not (si_neg & sj_neg):
                    fb.add((-1, -1))
                if fb:
                    forbidden[(i, j)] = fb
        valid_nodes = set()
        constraints_by_split = [[] for _ in range(n_splits)]
        for (i, j), fb_pairs in forbidden.items():
            constraints_by_split[j].append((i, fb_pairs))

        signs = [-1] * n_splits

        def add_valid_signs(split_idx):
            if split_idx == n_splits:
                valid_nodes.add(tuple(signs))
                return
            for sign in (-1, 1):
                valid = True
                for prev_idx, fb_pairs in constraints_by_split[split_idx]:
                    if (signs[prev_idx], sign) in fb_pairs:
                        valid = False
                        break
                if valid:
                    signs[split_idx] = sign
                    add_valid_signs(split_idx + 1)

        add_valid_signs(0)
        for signs in taxon_signs.values():
            valid_nodes.add(signs)
        edges = []
        for node in valid_nodes:
            for split_idx in range(n_splits):
                flipped = (
                    node[:split_idx]
                    + (-node[split_idx],)
                    + node[split_idx + 1:]
                )
                if flipped in valid_nodes and node < flipped:
                    edges.append((node, flipped, split_idx))
        return taxon_signs, valid_nodes, edges, splits_list, weights

    def _draw_network(
        self,
        ordering: list[str],
        plot_splits: list[tuple[frozenset, float]],
        all_taxa: frozenset,
    ):
        """Draw a planar splits graph using matplotlib."""
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.collections import LineCollection

        n = len(ordering)
        angles = {taxon: 2 * math.pi * i / n for i, taxon in enumerate(ordering)}

        # Filter to circular splits only
        circular_splits = []
        gap_positions_by_split = {}
        for split, weight in plot_splits:
            gap_positions = self._circular_gap_positions(split, ordering)
            if gap_positions is not None:
                circular_splits.append((split, weight))
                gap_positions_by_split[split] = gap_positions

        # Compute direction vectors
        directions = self._compute_split_directions(
            ordering,
            circular_splits,
            gap_positions_by_split,
        )

        # Keep only splits that have computed directions
        circular_splits = [
            (split, weight)
            for split, weight in circular_splits
            if split in directions
        ]

        # Build splits graph
        taxon_signs, valid_nodes, edges, splits_list, weights = self._build_splits_graph(
            circular_splits, all_taxa
        )

        config = self.plot_config
        # Force square figure for network graphs
        if config.fig_width is None and config.fig_height is None:
            size = max(10, min(30, 8 + n * 0.05))
            config.fig_width = size
            config.fig_height = size
        config.resolve(n_rows=n, n_cols=None)
        fig, ax = plt.subplots(1, 1, figsize=(config.fig_width, config.fig_height))
        ax.set_aspect("equal")

        # Determine label fontsize
        if config.ylabel_fontsize is not None:
            label_fontsize = config.ylabel_fontsize
        elif n > 100:
            label_fontsize = 0
        elif n > 50:
            label_fontsize = max(3, 8 - (n - 50) * 0.1)
        else:
            label_fontsize = 10

        if not splits_list:
            point_x = []
            point_y = []
            for taxon in ordering:
                angle = angles[taxon]
                x = math.cos(angle)
                y = math.sin(angle)
                if label_fontsize > 0:
                    deg = math.degrees(angle)
                    ha = "left" if -90 < deg < 90 or deg > 270 else "right"
                    rotation = deg if -90 < deg < 90 or deg > 270 else deg + 180
                    ax.text(x, y, taxon, ha=ha, va="center",
                            fontsize=label_fontsize, rotation=rotation,
                            rotation_mode="anchor")
                else:
                    point_x.append(x)
                    point_y.append(y)
            if point_x:
                ax.scatter(point_x, point_y, color="black", s=4)
        else:
            node_positions = {}
            for node in valid_nodes:
                x, y = 0.0, 0.0
                for i, sp in enumerate(splits_list):
                    dx, dy = directions[sp]
                    x += node[i] * weights[i] * dx / 2
                    y += node[i] * weights[i] * dy / 2
                node_positions[node] = (x, y)

            extent = _position_extent(node_positions.values())
            pendant_len = max(0.15, extent * 0.3)

            internal_segments = []
            for s1, s2, split_idx in edges:
                x1, y1 = node_positions[s1]
                x2, y2 = node_positions[s2]
                internal_segments.append(((x1, y1), (x2, y2)))
            if internal_segments:
                ax.add_collection(
                    LineCollection(
                        internal_segments,
                        colors="black",
                        linewidths=1.5,
                        zorder=2,
                    ),
                    autolim=True,
                )

            boundary_taxa = set()
            for split, weight in circular_splits:
                for i in range(n):
                    curr = ordering[i]
                    nxt = ordering[(i + 1) % n]
                    if (curr in split) != (nxt in split):
                        boundary_taxa.add(curr)
                        boundary_taxa.add(nxt)

            pendant_segments = []
            pendant_widths = []
            pendant_colors = []
            for taxon in ordering:
                angle = angles[taxon]
                if taxon in taxon_signs and taxon_signs[taxon] in node_positions:
                    nx, ny = node_positions[taxon_signs[taxon]]
                else:
                    nx, ny = 0.0, 0.0

                if taxon in boundary_taxa or n <= 50:
                    plen = pendant_len
                else:
                    plen = pendant_len * 0.3

                tx = nx + plen * math.cos(angle)
                ty = ny + plen * math.sin(angle)
                is_boundary = taxon in boundary_taxa
                pendant_segments.append(((nx, ny), (tx, ty)))
                pendant_widths.append(1.5 if is_boundary else 0.8)
                pendant_colors.append((0.0, 0.0, 0.0, 1.0 if is_boundary else 0.3))

                if label_fontsize > 0 and (taxon in boundary_taxa or n <= 50):
                    lx = tx + 0.03 * math.cos(angle)
                    ly = ty + 0.03 * math.sin(angle)
                    deg = math.degrees(angle)
                    ha = "left" if -90 < deg < 90 or deg > 270 else "right"
                    rotation = deg if -90 < deg < 90 or deg > 270 else deg + 180
                    ax.text(lx, ly, taxon, ha=ha, va="center",
                            fontsize=label_fontsize, rotation=rotation,
                            rotation_mode="anchor", zorder=4)
            if pendant_segments:
                ax.add_collection(
                    LineCollection(
                        pendant_segments,
                        colors=pendant_colors,
                        linewidths=pendant_widths,
                        zorder=2,
                    ),
                    autolim=True,
                )
            ax.autoscale_view()

            labeled_splits = set()
            for s1, s2, split_idx in edges:
                if split_idx not in labeled_splits:
                    x1, y1 = node_positions[s1]
                    x2, y2 = node_positions[s2]
                    mx = (x1 + x2) / 2
                    my = (y1 + y2) / 2
                    w = weights[split_idx]
                    ax.text(
                        mx, my, f"{w:.4f}",
                        fontsize=7, ha="center", va="center",
                        color="gray", zorder=3,
                    )
                    labeled_splits.add(split_idx)

        ax.axis("off")
        if config.show_title:
            ax.set_title(
                config.title or "NeighborNet Splits Network",
                fontsize=config.title_fontsize or 14,
                pad=20,
            )

        plt.tight_layout()
        plt.savefig(self.output_path, dpi=config.dpi, bbox_inches="tight")
        plt.close()

    # ------------------------------------------------------------------
    # Output helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _format_split(split: frozenset) -> str:
        return "{" + ", ".join(sorted(split)) + "}"

    # ------------------------------------------------------------------
    # Main entry point
    # ------------------------------------------------------------------

    def run(self):
        # 1. Get distance matrix
        if self.alignment_path:
            taxa, D = self._compute_distances_from_alignment()
        else:
            taxa, D = self._read_distance_matrix()

        n = len(taxa)
        if n < 4:
            raise PhykitUserError(
                ["At least 4 taxa are required for NeighborNet analysis."], code=2
            )

        # 2. Build NJ tree to get circular ordering
        ordering = self._nj_circular_ordering(D, taxa)

        # 3. Enumerate circular splits
        splits = self._enumerate_circular_splits(ordering)

        # 4. Estimate split weights via NNLS
        weights = self._estimate_split_weights(D, splits, ordering, taxa)

        # 5. Filter to positive-weight splits
        positive = [
            (splits[i], weights[i])
            for i in range(len(splits))
            if weights[i] > 1e-8
        ]
        positive.sort(key=lambda x: -x[1])

        # 6. Cap for visualization
        if len(positive) > self.max_splits:
            print(
                f"Warning: {len(positive)} splits with positive weight; "
                f"using top {self.max_splits} for visualization.",
                file=sys.stderr,
            )
            plot_splits = positive[: self.max_splits]
        else:
            plot_splits = positive

        # 7. Draw network
        all_taxa = frozenset(taxa)
        self._draw_network(ordering, plot_splits, all_taxa)

        # 8. Output
        if self.json_output:
            splits_list = [
                {
                    "split": sorted(split),
                    "weight": round(float(weight), 6),
                }
                for split, weight in positive
            ]
            print_json(
                dict(
                    taxa=sorted(taxa),
                    taxa_count=n,
                    distance_metric=self.metric,
                    circular_ordering=ordering,
                    total_circular_splits=len(splits),
                    positive_weight_splits=len(positive),
                    splits=splits_list,
                    output=self.output_path,
                )
            )
        else:
            lines = [
                "NeighborNet Analysis",
                f"Taxa: {n}",
                f"Distance metric: {self.metric}",
                f"Circular splits with positive weight: {len(positive)}",
                "---",
            ]
            lines.extend(
                f"{self._format_split(split)}\t{weight:.6f}"
                for split, weight in positive
            )
            lines.append(f"Output: {self.output_path}")
            print("\n".join(lines))
