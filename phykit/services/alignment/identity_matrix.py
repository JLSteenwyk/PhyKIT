import sys
import json
from typing import Dict, List, Optional, Tuple

import numpy as np
from Bio import SeqIO, Phylo
from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list
from scipy.spatial.distance import squareform

from .base import Alignment
from ...helpers.json_output import print_json
from ...helpers.plot_config import PlotConfig


class IdentityMatrix(Alignment):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(
            alignment_file_path=parsed["alignment_path"],
        )
        self.output_path = parsed["output_path"]
        self.metric = parsed["metric"]
        self.tree_path = parsed["tree_path"]
        self.sort_method = parsed["sort_method"]
        self.partition_path = parsed["partition_path"]
        self.json_output = parsed["json_output"]
        self.plot_config = parsed["plot_config"]

    def process_args(self, args) -> Dict[str, object]:
        return dict(
            alignment_path=args.alignment,
            output_path=args.output,
            metric=getattr(args, "metric", "identity"),
            tree_path=getattr(args, "tree", None),
            sort_method=getattr(args, "sort", "cluster"),
            partition_path=getattr(args, "partition", None),
            json_output=getattr(args, "json", False),
            plot_config=PlotConfig.from_args(args),
        )

    def run(self) -> None:
        # 1. Parse alignment
        sequences, taxa_names, aln_length = self._parse_alignment()
        n_taxa = len(taxa_names)

        if n_taxa < 2:
            print("Error: alignment must contain at least 2 taxa.", file=sys.stderr)
            raise SystemExit(2)

        # 2. Compute pairwise identity matrix
        matrix = self._compute_identity_matrix(sequences, taxa_names)

        # 3. Apply metric transformation
        if self.metric == "p-distance":
            matrix = 1.0 - matrix

        # 4. Determine ordering
        order, ordered_labels = self._determine_order(
            matrix, taxa_names, n_taxa
        )

        # 5. Parse partitions if provided
        partitions = None
        if self.partition_path:
            partitions = self._parse_partitions(self.partition_path)

        # 6. Plot heatmap
        self._plot_heatmap(
            matrix, taxa_names, order, ordered_labels, n_taxa,
            partitions=partitions, sequences=sequences,
        )

        # 6. Compute summary stats for output
        # Use identity values (not p-distance) for min/max reporting
        if self.metric == "p-distance":
            identity_matrix = 1.0 - matrix
        else:
            identity_matrix = matrix

        # Extract upper triangle (excluding diagonal) for stats
        upper_indices = np.triu_indices(n_taxa, k=1)
        pairwise_values = identity_matrix[upper_indices]
        mean_val = float(np.mean(pairwise_values))
        min_val = float(np.min(pairwise_values))
        max_val = float(np.max(pairwise_values))

        # Find min and max pairs
        min_idx = np.argmin(pairwise_values)
        max_idx = np.argmax(pairwise_values)
        min_pair = [taxa_names[upper_indices[0][min_idx]],
                    taxa_names[upper_indices[1][min_idx]]]
        max_pair = [taxa_names[upper_indices[0][max_idx]],
                    taxa_names[upper_indices[1][max_idx]]]

        # 7. Output
        if self.json_output:
            reordered = matrix[np.ix_(order, order)]
            payload = {
                "n_taxa": n_taxa,
                "alignment_length": aln_length,
                "metric": self.metric,
                "mean_identity": round(mean_val, 4),
                "min_identity": {
                    "value": round(min_val, 4),
                    "taxa": min_pair,
                },
                "max_identity": {
                    "value": round(max_val, 4),
                    "taxa": max_pair,
                },
                "matrix": reordered.tolist(),
                "taxa_order": ordered_labels,
                "output_file": self.output_path,
            }
            print_json(payload)
            return

        try:
            print(f"Sequence Identity Matrix")
            print(f"Alignment: {self.alignment_file_path}")
            print(f"Taxa: {n_taxa}")
            print(f"Alignment length: {aln_length}")
            print(f"Metric: {self.metric}")
            print(f"Sort: {self.sort_method}")
            print(f"Mean pairwise identity: {round(mean_val, 4)}")
            print(f"Min: {round(min_val, 4)} ({min_pair[0]} vs {min_pair[1]})")
            print(f"Max: {round(max_val, 4)} ({max_pair[0]} vs {max_pair[1]})")
            print(f"Output: {self.output_path}")
        except BrokenPipeError:
            pass

    def _parse_partitions(self, path: str) -> List[Tuple[str, int, int]]:
        """Parse RAxML-style partition file.

        Format: AUTO, gene_name=start-end (1-indexed)
        Returns list of (name, start_0indexed, end_0indexed).
        """
        partitions = []
        with open(path) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                # Format: AUTO, name=start-end  OR  DNA, name=start-end
                if "=" in line and "-" in line:
                    parts = line.split(",", 1)
                    if len(parts) == 2:
                        name_range = parts[1].strip()
                    else:
                        name_range = parts[0].strip()
                    name, range_str = name_range.split("=", 1)
                    name = name.strip()
                    start_str, end_str = range_str.strip().split("-", 1)
                    start = int(start_str) - 1  # convert to 0-indexed
                    end = int(end_str)  # end is inclusive in RAxML, so keep as-is for slicing
                    partitions.append((name, start, end))
        return partitions

    def _compute_partition_identities(
        self,
        sequences: Dict[str, str],
        taxa_names: List[str],
        partitions: List[Tuple[str, int, int]],
    ) -> Tuple[List[str], np.ndarray]:
        """Compute per-partition pairwise identity.

        Returns (partition_names, identity_array) where identity_array
        is shape (n_pairs, n_partitions) with the upper-triangle pairs
        in row-major order.
        """
        n = len(taxa_names)
        n_parts = len(partitions)
        ambiguous = {"-", "?", "N", "X", "*", "n", "x"}

        # Build pair indices (upper triangle)
        pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
        n_pairs = len(pairs)

        part_names = [p[0] for p in partitions]
        result = np.zeros((n_pairs, n_parts), dtype=np.float64)

        for p_idx, (pname, start, end) in enumerate(partitions):
            for pair_idx, (i, j) in enumerate(pairs):
                seq_i = sequences[taxa_names[i]][start:end]
                seq_j = sequences[taxa_names[j]][start:end]
                matches = 0
                compared = 0
                for a, b in zip(seq_i, seq_j):
                    if a in ambiguous or b in ambiguous:
                        continue
                    compared += 1
                    if a == b:
                        matches += 1
                result[pair_idx, p_idx] = matches / compared if compared > 0 else 0.0

        return part_names, result

    def _parse_alignment(self) -> Tuple[Dict[str, str], List[str], int]:
        """Parse FASTA alignment and return sequences, taxa names, and alignment length."""
        alignment, _, is_protein = self.get_alignment_and_format()

        taxa_names = []
        sequences = {}
        for record in alignment:
            name = record.id
            taxa_names.append(name)
            sequences[name] = str(record.seq).upper()

        if not taxa_names:
            print("Error: no sequences found in alignment.", file=sys.stderr)
            raise SystemExit(2)

        # Validate all same length
        lengths = set(len(seq) for seq in sequences.values())
        if len(lengths) > 1:
            print("Error: sequences have different lengths. Is this an alignment?",
                  file=sys.stderr)
            raise SystemExit(2)

        aln_length = lengths.pop()
        return sequences, taxa_names, aln_length

    def _compute_identity_matrix(
        self,
        sequences: Dict[str, str],
        taxa_names: List[str],
    ) -> np.ndarray:
        """Compute pairwise identity matrix. Returns identity values (0-1)."""
        n = len(taxa_names)
        matrix = np.ones((n, n), dtype=np.float64)

        # Pre-convert sequences to arrays for speed
        seq_arrays = {}
        for name in taxa_names:
            seq_arrays[name] = np.array(list(sequences[name]), dtype="U1")

        ambiguous_chars = set(["-", "?", "N", "X", "*", "n", "x"])

        for i in range(n):
            seq_i = seq_arrays[taxa_names[i]]
            for j in range(i + 1, n):
                seq_j = seq_arrays[taxa_names[j]]

                # Create masks for valid (non-gap, non-ambiguous) positions
                valid_i = np.array(
                    [c not in ambiguous_chars for c in sequences[taxa_names[i]]],
                    dtype=bool,
                )
                valid_j = np.array(
                    [c not in ambiguous_chars for c in sequences[taxa_names[j]]],
                    dtype=bool,
                )
                valid = valid_i & valid_j

                compared = int(np.sum(valid))
                if compared > 0:
                    matches = int(np.sum((seq_i == seq_j) & valid))
                    identity = matches / compared
                else:
                    identity = 0.0

                matrix[i, j] = identity
                matrix[j, i] = identity

        return matrix

    def _determine_order(
        self,
        matrix: np.ndarray,
        taxa_names: List[str],
        n_taxa: int,
    ) -> Tuple[List[int], List[str]]:
        """Determine taxa ordering based on sort method."""
        if self.sort_method == "alpha":
            sorted_indices = sorted(range(n_taxa), key=lambda i: taxa_names[i])
            order = sorted_indices
        elif self.sort_method == "tree":
            if self.tree_path is None:
                print(
                    "Error: --sort tree requires --tree <file>.",
                    file=sys.stderr,
                )
                raise SystemExit(2)
            tree = Phylo.read(self.tree_path, "newick")
            tip_order = [tip.name for tip in tree.get_terminals()]
            # Map tree tip order to matrix indices
            name_to_idx = {name: idx for idx, name in enumerate(taxa_names)}
            order = []
            for tip_name in tip_order:
                if tip_name in name_to_idx:
                    order.append(name_to_idx[tip_name])
            # Add any taxa not in the tree at the end
            in_tree = set(order)
            for idx in range(n_taxa):
                if idx not in in_tree:
                    order.append(idx)
        elif self.sort_method == "cluster":
            if n_taxa >= 3:
                # For clustering, always use distance = 1 - identity
                if self.metric == "p-distance":
                    # matrix already holds p-distance; identity = 1 - matrix
                    distance_matrix = np.copy(matrix)
                else:
                    distance_matrix = 1.0 - matrix
                distance_matrix = np.clip(distance_matrix, 0.0, 1.0)
                np.fill_diagonal(distance_matrix, 0.0)
                condensed = squareform(distance_matrix, checks=False)
                Z = linkage(condensed, method="average")
                order = list(leaves_list(Z))
            else:
                order = list(range(n_taxa))
        else:
            order = list(range(n_taxa))

        ordered_labels = [taxa_names[i] for i in order]
        return order, ordered_labels

    def _plot_heatmap(
        self,
        matrix: np.ndarray,
        taxa_names: List[str],
        order: List[int],
        ordered_labels: List[str],
        n_taxa: int,
        partitions: Optional[List[Tuple[str, int, int]]] = None,
        sequences: Optional[Dict[str, str]] = None,
    ) -> None:
        """Plot the identity/p-distance heatmap with optional dendrograms."""
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        except ImportError:
            print(
                "matplotlib is required for identity_matrix. "
                "Install matplotlib and retry.",
                file=sys.stderr,
            )
            raise SystemExit(2)

        config = self.plot_config
        config.resolve(n_rows=n_taxa, n_cols=n_taxa)

        # Reorder matrix
        reordered = matrix[np.ix_(order, order)]

        if self.sort_method == "cluster" and n_taxa >= 3:
            # Layout with dendrograms
            fig_w = config.fig_width or max(8, min(20, n_taxa * 0.4))
            fig_h = config.fig_height or fig_w

            fig = plt.figure(figsize=(fig_w, fig_h))

            # Axes positions: [left, bottom, width, height]
            ax_dendro_top = fig.add_axes([0.10, 0.80, 0.65, 0.15])
            ax_dendro_left = fig.add_axes([0.01, 0.10, 0.08, 0.68])
            ax_heat = fig.add_axes([0.10, 0.10, 0.65, 0.68])
            ax_cbar = fig.add_axes([0.78, 0.10, 0.02, 0.68])

            # Compute linkage for dendrograms
            if self.metric == "p-distance":
                dist_for_linkage = np.copy(matrix)
            else:
                dist_for_linkage = 1.0 - matrix
            dist_for_linkage = np.clip(dist_for_linkage, 0.0, 1.0)
            np.fill_diagonal(dist_for_linkage, 0.0)
            condensed = squareform(dist_for_linkage, checks=False)
            Z = linkage(condensed, method="average")

            # Top dendrogram
            dendrogram(
                Z, ax=ax_dendro_top,
                leaf_rotation=90, leaf_font_size=0,
                color_threshold=0, above_threshold_color="#333333",
                no_labels=True,
            )
            ax_dendro_top.axis("off")

            # Left dendrogram (rotated)
            dendrogram(
                Z, ax=ax_dendro_left,
                orientation="left", leaf_font_size=0,
                color_threshold=0, above_threshold_color="#333333",
                no_labels=True,
            )
            ax_dendro_left.axis("off")

        else:
            # Simple layout without dendrograms
            fig_w = config.fig_width or max(6, min(20, n_taxa * 0.4))
            fig_h = config.fig_height or fig_w

            fig = plt.figure(figsize=(fig_w, fig_h))
            ax_heat = fig.add_axes([0.08, 0.10, 0.72, 0.80])
            ax_cbar = fig.add_axes([0.83, 0.10, 0.02, 0.80])

        # Draw heatmap
        cmap = "viridis" if self.metric == "identity" else "viridis_r"
        im = ax_heat.imshow(
            reordered, cmap=cmap, aspect="auto",
            interpolation="nearest", vmin=0, vmax=1,
        )

        # Colorbar
        cbar_label = "Pairwise identity" if self.metric == "identity" else "p-distance"
        cbar = fig.colorbar(im, cax=ax_cbar)
        cbar.set_label(cbar_label)

        # Labels
        label_fontsize = config.ylabel_fontsize
        if label_fontsize and label_fontsize > 0:
            ax_heat.set_xticks(np.arange(n_taxa))
            ax_heat.set_yticks(np.arange(n_taxa))
            ax_heat.set_xticklabels(
                ordered_labels, rotation=90,
                fontsize=config.xlabel_fontsize or label_fontsize,
            )
            ax_heat.set_yticklabels(
                ordered_labels,
                fontsize=label_fontsize,
            )
        else:
            ax_heat.set_xticks([])
            ax_heat.set_yticks([])

        # Axis labels
        sort_label = f"Taxa ({self.sort_method})"
        ax_heat.set_xlabel(sort_label)
        ax_heat.set_ylabel(sort_label)

        if config.axis_fontsize:
            ax_heat.xaxis.label.set_fontsize(config.axis_fontsize)
            ax_heat.yaxis.label.set_fontsize(config.axis_fontsize)

        # Partition annotation panel
        if partitions and sequences:
            part_names, part_identities = self._compute_partition_identities(
                sequences, taxa_names, partitions
            )
            n_parts = len(part_names)

            # Build a taxa x taxa matrix per partition using the same order
            # part_identities is (n_pairs, n_partitions) for upper-triangle pairs
            # We need to map pairs back to the reordered matrix indices
            n = n_taxa
            pair_to_idx = {}
            idx = 0
            for i in range(n):
                for j in range(i + 1, n):
                    pair_to_idx[(i, j)] = idx
                    idx += 1

            # Build reordered per-partition strip: rows = reordered taxa pairs
            # For the strip, show average identity per partition per taxon
            # as a narrow heatmap: rows = taxa (same order), cols = partitions
            # Value = mean identity of that taxon against all others for that partition
            strip = np.zeros((n_taxa, n_parts), dtype=np.float64)
            for p_idx in range(n_parts):
                for row_pos, orig_i in enumerate(order):
                    vals = []
                    for orig_j in range(n):
                        if orig_i == orig_j:
                            continue
                        key = (min(orig_i, orig_j), max(orig_i, orig_j))
                        if key in pair_to_idx:
                            vals.append(part_identities[pair_to_idx[key], p_idx])
                    strip[row_pos, p_idx] = np.mean(vals) if vals else 0.0

            # Add partition strip panel to the right of the heatmap
            heat_right = ax_heat.get_position().x1
            strip_width = min(0.15, 0.01 * n_parts + 0.03)
            ax_strip = fig.add_axes([
                heat_right + 0.02, ax_heat.get_position().y0,
                strip_width, ax_heat.get_position().height,
            ])
            strip_im = ax_strip.imshow(
                strip, cmap="viridis", aspect="auto",
                interpolation="nearest", vmin=0, vmax=1,
            )
            ax_strip.set_xticks(range(n_parts))
            ax_strip.set_xticklabels(
                part_names, rotation=90,
                fontsize=max(3, min(6, 7 - n_parts * 0.05)),
            )
            ax_strip.set_yticks([])
            ax_strip.set_xlabel("Per-gene\nidentity", fontsize=7)

        # Title
        if config.show_title:
            title_text = config.title or (
                "Pairwise Identity Matrix" if self.metric == "identity"
                else "p-distance Matrix"
            )
            fig.suptitle(title_text, fontsize=config.title_fontsize)

        fig.savefig(self.output_path, dpi=config.dpi, bbox_inches="tight")
        plt.close(fig)
