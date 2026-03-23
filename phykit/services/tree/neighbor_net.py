import math
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

from ...errors import PhykitUserError
from ...helpers.json_output import print_json
from ...helpers.plot_config import PlotConfig


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

    def process_args(self, args) -> Dict:
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

    def _compute_distances_from_alignment(self) -> Tuple[List[str], np.ndarray]:
        """Compute pairwise distance matrix from a FASTA alignment."""
        from Bio import SeqIO

        try:
            records = list(SeqIO.parse(self.alignment_path, "fasta"))
        except FileNotFoundError:
            raise PhykitUserError(
                [
                    f"{self.alignment_path} corresponds to no such file or directory.",
                    "Please check filename and pathing",
                ],
                code=2,
            )

        if len(records) == 0:
            raise PhykitUserError(
                ["No sequences found in the alignment file."], code=2
            )

        taxa = [r.id for r in records]
        seqs = {r.id: str(r.seq).upper() for r in records}
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

        return taxa, D

    def _read_distance_matrix(self) -> Tuple[List[str], np.ndarray]:
        """Read a pre-computed distance matrix from a CSV file."""
        import csv

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
            reader = csv.reader(fh)
            header = next(reader)
            # First cell may be empty or a label
            if header[0].strip() == "" or header[0].strip().lower() in ("taxon", "taxa", "label", "name"):
                taxa = [h.strip() for h in header[1:]]
            else:
                taxa = [h.strip() for h in header]

            n = len(taxa)
            D = np.zeros((n, n))
            for i, row in enumerate(reader):
                # Row may start with a taxon label
                values = row[1:] if len(row) > n else row
                for j, val in enumerate(values):
                    D[i, j] = float(val.strip())

        return taxa, D

    # ------------------------------------------------------------------
    # NJ circular ordering
    # ------------------------------------------------------------------

    @staticmethod
    def _nj_circular_ordering(D: np.ndarray, taxa: List[str]) -> List[str]:
        """Build a NJ tree and extract a circular ordering from leaf traversal."""
        from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

        n = len(taxa)
        # Convert to BioPython DistanceMatrix (lower triangle including diagonal)
        matrix = []
        for i in range(n):
            row = [D[i, j] for j in range(i + 1)]
            matrix.append(row)

        dm = DistanceMatrix(taxa, matrix)
        constructor = DistanceTreeConstructor()
        nj_tree = constructor.nj(dm)

        # Extract leaf ordering from the tree
        ordering = [tip.name for tip in nj_tree.get_terminals()]
        return ordering

    # ------------------------------------------------------------------
    # Circular splits
    # ------------------------------------------------------------------

    @staticmethod
    def _enumerate_circular_splits(ordering: List[str]) -> List[frozenset]:
        """Enumerate all non-trivial circular splits compatible with the ordering.

        A circular split divides the circular ordering into two contiguous arcs.
        We exclude trivial splits (single taxon or all-but-one).
        """
        n = len(ordering)
        seen = set()
        splits = []
        for start in range(n):
            for size in range(2, n - 1):
                indices = [(start + k) % n for k in range(size)]
                side = frozenset(ordering[idx] for idx in indices)
                complement = frozenset(ordering) - side
                # Canonical form: smaller side (or lexicographic tiebreak)
                if len(side) < len(complement):
                    canonical = side
                elif len(side) > len(complement):
                    canonical = complement
                else:
                    if sorted(side) < sorted(complement):
                        canonical = side
                    else:
                        canonical = complement
                if canonical not in seen:
                    seen.add(canonical)
                    splits.append(canonical)
        return splits

    # ------------------------------------------------------------------
    # NNLS weight estimation
    # ------------------------------------------------------------------

    @staticmethod
    def _estimate_split_weights(
        D: np.ndarray,
        splits: List[frozenset],
        ordering: List[str],
        taxa: List[str],
    ) -> np.ndarray:
        """Estimate non-negative split weights via NNLS."""
        from scipy.optimize import nnls

        n = len(taxa)
        name_to_idx = {name: i for i, name in enumerate(taxa)}

        pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
        n_pairs = len(pairs)
        n_splits = len(splits)

        A = np.zeros((n_pairs, n_splits))
        b = np.zeros(n_pairs)

        for pair_idx, (i, j) in enumerate(pairs):
            b[pair_idx] = D[i, j]
            ti = taxa[i]
            tj = taxa[j]
            for split_idx, split in enumerate(splits):
                # Does this split separate taxa i and j?
                if (ti in split) != (tj in split):
                    A[pair_idx, split_idx] = 1.0

        weights, _residual = nnls(A, b)
        return weights

    # ------------------------------------------------------------------
    # Visualization (reuses patterns from ConsensusNetwork)
    # ------------------------------------------------------------------

    @staticmethod
    def _canonical_split(taxa_side: frozenset, all_taxa: frozenset) -> frozenset:
        complement = all_taxa - taxa_side
        if len(taxa_side) < len(complement):
            return taxa_side
        if len(taxa_side) > len(complement):
            return complement
        if sorted(taxa_side) < sorted(complement):
            return taxa_side
        return complement

    @staticmethod
    def _is_circular_split(split, ordering):
        """Check if a split has exactly 2 boundary gaps in the circular ordering."""
        n = len(ordering)
        gaps = 0
        for i in range(n):
            curr = ordering[i]
            nxt = ordering[(i + 1) % n]
            if (curr in split) != (nxt in split):
                gaps += 1
        return gaps == 2

    @staticmethod
    def _compute_split_directions(ordering, circular_splits):
        """Compute 2D direction vectors for each circular split."""
        n = len(ordering)
        angles = {taxon: 2 * math.pi * i / n for i, taxon in enumerate(ordering)}
        directions = {}
        for split, weight in circular_splits:
            gap_positions = []
            for i in range(n):
                curr = ordering[i]
                nxt = ordering[(i + 1) % n]
                if (curr in split) != (nxt in split):
                    gap_positions.append(i)
            if len(gap_positions) != 2:
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
        for combo in range(2 ** n_splits):
            signs = tuple(1 if (combo >> j) & 1 else -1 for j in range(n_splits))
            valid = True
            for (i, j), fb_pairs in forbidden.items():
                if (signs[i], signs[j]) in fb_pairs:
                    valid = False
                    break
            if valid:
                valid_nodes.add(signs)
        for signs in taxon_signs.values():
            valid_nodes.add(signs)
        edges = []
        valid_list = list(valid_nodes)
        for i in range(len(valid_list)):
            for j in range(i + 1, len(valid_list)):
                s1 = valid_list[i]
                s2 = valid_list[j]
                diff = [k for k in range(n_splits) if s1[k] != s2[k]]
                if len(diff) == 1:
                    edges.append((s1, s2, diff[0]))
        return taxon_signs, valid_nodes, edges, splits_list, weights

    def _draw_network(
        self,
        ordering: List[str],
        plot_splits: List[Tuple[frozenset, float]],
        all_taxa: frozenset,
    ):
        """Draw a planar splits graph using matplotlib."""
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        n = len(ordering)
        angles = {taxon: 2 * math.pi * i / n for i, taxon in enumerate(ordering)}

        # Filter to circular splits only
        circular_splits = [
            (split, weight)
            for split, weight in plot_splits
            if self._is_circular_split(split, ordering)
        ]

        # Compute direction vectors
        directions = self._compute_split_directions(ordering, circular_splits)

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
                    ax.plot(x, y, "o", color="black", markersize=2)
        else:
            node_positions = {}
            for node in valid_nodes:
                x, y = 0.0, 0.0
                for i, sp in enumerate(splits_list):
                    dx, dy = directions[sp]
                    x += node[i] * weights[i] * dx / 2
                    y += node[i] * weights[i] * dy / 2
                node_positions[node] = (x, y)

            all_x = [p[0] for p in node_positions.values()]
            all_y = [p[1] for p in node_positions.values()]
            extent = max(
                max(all_x) - min(all_x) if all_x else 0,
                max(all_y) - min(all_y) if all_y else 0,
            )
            pendant_len = max(0.15, extent * 0.3)

            for s1, s2, split_idx in edges:
                x1, y1 = node_positions[s1]
                x2, y2 = node_positions[s2]
                ax.plot([x1, x2], [y1, y2], "-", color="black", linewidth=1.5, zorder=2)

            boundary_taxa = set()
            for split, weight in circular_splits:
                for i in range(n):
                    curr = ordering[i]
                    nxt = ordering[(i + 1) % n]
                    if (curr in split) != (nxt in split):
                        boundary_taxa.add(curr)
                        boundary_taxa.add(nxt)

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
                ax.plot([nx, tx], [ny, ty], "-", color="black",
                        linewidth=0.8 if taxon not in boundary_taxa else 1.5,
                        alpha=0.3 if taxon not in boundary_taxa else 1.0,
                        zorder=2)

                if label_fontsize > 0 and (taxon in boundary_taxa or n <= 50):
                    lx = tx + 0.03 * math.cos(angle)
                    ly = ty + 0.03 * math.sin(angle)
                    deg = math.degrees(angle)
                    ha = "left" if -90 < deg < 90 or deg > 270 else "right"
                    rotation = deg if -90 < deg < 90 or deg > 270 else deg + 180
                    ax.text(lx, ly, taxon, ha=ha, va="center",
                            fontsize=label_fontsize, rotation=rotation,
                            rotation_mode="anchor", zorder=4)

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
            print(f"NeighborNet Analysis")
            print(f"Taxa: {n}")
            print(f"Distance metric: {self.metric}")
            print(f"Circular splits with positive weight: {len(positive)}")
            print("---")
            for split, weight in positive:
                print(f"{self._format_split(split)}\t{weight:.6f}")
            print(f"Output: {self.output_path}")
