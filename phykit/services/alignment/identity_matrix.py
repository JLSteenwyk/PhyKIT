from __future__ import annotations

import math
import sys

from .base import Alignment


class _LazyNumpy:
    def __getattr__(self, name):
        import numpy as _np

        return getattr(_np, name)


np = _LazyNumpy()

_LINKAGE = None
_DENDROGRAM = None
_LEAVES_LIST = None
_SQUAREFORM = None
_IDENTITY_INVALID_BYTES = b"-?NX*nx"
_IDENTITY_INVALID_CHARS = frozenset("-?NX*nx")
_NO_INVALID_DIRECT_MIN_LENGTH = 2048
_NO_INVALID_DIRECT_MIN_TAXA = 100
_NO_INVALID_DIRECT_SHORT_ALIGNMENT_MIN_TAXA = 150
_NO_INVALID_DIRECT_MAX_TAXA = 700


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


def linkage(*args, **kwargs):
    global _LINKAGE
    if _LINKAGE is None:
        from scipy.cluster.hierarchy import linkage as _linkage

        _LINKAGE = _linkage

    return _LINKAGE(*args, **kwargs)


def dendrogram(*args, **kwargs):
    global _DENDROGRAM
    if _DENDROGRAM is None:
        from scipy.cluster.hierarchy import dendrogram as _dendrogram

        _DENDROGRAM = _dendrogram

    return _DENDROGRAM(*args, **kwargs)


def leaves_list(*args, **kwargs):
    global _LEAVES_LIST
    if _LEAVES_LIST is None:
        from scipy.cluster.hierarchy import leaves_list as _leaves_list

        _LEAVES_LIST = _leaves_list

    return _LEAVES_LIST(*args, **kwargs)


def squareform(*args, **kwargs):
    global _SQUAREFORM
    if _SQUAREFORM is None:
        from scipy.spatial.distance import squareform as _squareform

        _SQUAREFORM = _squareform

    return _SQUAREFORM(*args, **kwargs)


def _condensed_index_to_pair(index: int, n_items: int) -> tuple[int, int]:
    row = int(
        n_items
        - 2
        - math.floor(
            math.sqrt(-8 * index + 4 * n_items * (n_items - 1) - 7) / 2.0
            - 0.5
        )
    )
    col = int(index + row + 1 - n_items * row + row * (row + 1) // 2)
    return row, col


def _identical_sequence_identity_value(sequence: str) -> float:
    try:
        return (
            1.0
            if sequence.encode("ascii").translate(None, _IDENTITY_INVALID_BYTES)
            else 0.0
        )
    except UnicodeEncodeError:
        return (
            1.0
            if any(char not in _IDENTITY_INVALID_CHARS for char in sequence)
            else 0.0
        )


def _all_sequences_identical(sequences: dict[str, str], taxa_names: list[str]) -> bool:
    iterator = iter(taxa_names)
    first_name = next(iterator, None)
    if first_name is None:
        return False
    first_sequence = sequences[first_name]
    for name in iterator:
        if sequences[name] != first_sequence:
            return False
    return True


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

    def process_args(self, args) -> dict[str, object]:
        from ...helpers.plot_config import PlotConfig

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
        order, ordered_labels, cluster_linkage = self._determine_order(
            matrix, taxa_names, n_taxa, return_linkage=True,
        )

        # 5. Parse partitions if provided
        partitions = None
        if self.partition_path:
            partitions = self._parse_partitions(self.partition_path)

        # 6. Plot heatmap
        self._plot_heatmap(
            matrix, taxa_names, order, ordered_labels, n_taxa,
            partitions=partitions, sequences=sequences,
            cluster_linkage=cluster_linkage,
        )

        # 6. Compute summary stats for output
        # Use identity values (not p-distance) for min/max reporting
        if self.metric == "p-distance":
            identity_matrix = 1.0 - matrix
        else:
            identity_matrix = matrix

        mean_val, min_val, min_pair, max_val, max_pair = (
            self._summarize_identity_matrix(identity_matrix, taxa_names)
        )

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
            self._print_text_output(
                n_taxa,
                aln_length,
                mean_val,
                min_val,
                min_pair,
                max_val,
                max_pair,
            )
        except BrokenPipeError:
            pass

    def _print_text_output(
        self,
        n_taxa: int,
        aln_length: int,
        mean_val: float,
        min_val: float,
        min_pair,
        max_val: float,
        max_pair,
    ) -> None:
        print(
            (
                "Sequence Identity Matrix\n"
                "Alignment: %s\n"
                "Taxa: %s\n"
                "Alignment length: %s\n"
                "Metric: %s\n"
                "Sort: %s\n"
                "Mean pairwise identity: %s\n"
                "Min: %s (%s vs %s)\n"
                "Max: %s (%s vs %s)\n"
                "Output: %s"
            )
            % (
                self.alignment_file_path,
                n_taxa,
                aln_length,
                self.metric,
                self.sort_method,
                round(mean_val, 4),
                round(min_val, 4),
                min_pair[0],
                min_pair[1],
                round(max_val, 4),
                max_pair[0],
                max_pair[1],
                self.output_path,
            )
        )

    @staticmethod
    def _summarize_identity_matrix(identity_matrix, taxa_names):
        pairwise_values = squareform(identity_matrix, checks=False)
        mean_val = float(pairwise_values.mean())
        min_index = int(pairwise_values.argmin())
        max_index = int(pairwise_values.argmax())
        min_val = float(pairwise_values[min_index])
        max_val = float(pairwise_values[max_index])

        min_row, min_col = _condensed_index_to_pair(min_index, len(taxa_names))
        max_row, max_col = _condensed_index_to_pair(max_index, len(taxa_names))
        min_pair = [taxa_names[min_row], taxa_names[min_col]]
        max_pair = [taxa_names[max_row], taxa_names[max_col]]
        return mean_val, min_val, min_pair, max_val, max_pair

    def _parse_partitions(self, path: str) -> list[tuple[str, int, int]]:
        """Parse RAxML-style partition file.

        Format: AUTO, gene_name=start-end (1-indexed)
        Returns list of (name, start_0indexed, end_0indexed).
        """
        partitions = []
        append = partitions.append
        with open(path) as f:
            for line in f:
                line = line.strip()
                if not line or line[0] == "#":
                    continue
                # Format: AUTO, name=start-end  OR  DNA, name=start-end
                if "=" in line and "-" in line:
                    name_range, sep, rest = line.partition(",")
                    if sep:
                        name_range = rest.strip()
                    else:
                        name_range = name_range.strip()
                    name, _sep, range_str = name_range.partition("=")
                    name = name.strip()
                    start_str, _sep, end_str = range_str.strip().partition("-")
                    start = int(start_str) - 1  # convert to 0-indexed
                    end = int(end_str)  # end is inclusive in RAxML, so keep as-is for slicing
                    append((name, start, end))
        return partitions

    def _compute_partition_identities(
        self,
        sequences: dict[str, str],
        taxa_names: list[str],
        partitions: list[tuple[str, int, int]],
    ) -> tuple[list[str], np.ndarray]:
        """Compute per-partition pairwise identity.

        Returns (partition_names, identity_array) where identity_array
        is shape (n_pairs, n_partitions) with the upper-triangle pairs
        in row-major order.
        """
        n = len(taxa_names)
        n_parts = len(partitions)
        if _all_sequences_identical(sequences, taxa_names):
            first_sequence = sequences[taxa_names[0]]
            pair_count = n * (n - 1) // 2
            part_values = np.fromiter(
                (
                    _identical_sequence_identity_value(first_sequence[start:end])
                    for _name, start, end in partitions
                ),
                dtype=np.float64,
                count=n_parts,
            )
            return [partition[0] for partition in partitions], np.tile(
                part_values,
                (pair_count, 1),
            )

        try:
            seq_matrix = self._sequence_byte_matrix(sequences, taxa_names)
        except UnicodeEncodeError:
            return self._compute_partition_identities_pairwise(
                sequences, taxa_names, partitions
            )

        invalid = np.frombuffer(b"-?NX*nx", dtype=np.uint8)
        valid_masks = ~np.isin(seq_matrix, invalid)
        valid_symbols = np.unique(seq_matrix[valid_masks])

        n_pairs = n * (n - 1) // 2

        part_names = [p[0] for p in partitions]
        result = np.zeros((n_pairs, n_parts), dtype=np.float64)
        valid_float = valid_masks.astype(np.float64, copy=False)

        for p_idx, (pname, start, end) in enumerate(partitions):
            part_seq_matrix = seq_matrix[:, start:end]
            part_valid_masks = valid_masks[:, start:end]
            # Float masks route dense products through optimized BLAS while
            # preserving exact integer counts for practical alignment lengths.
            valid_i = valid_float[:, start:end]
            compared = valid_i @ valid_i.T

            matches = np.zeros((n, n), dtype=np.float64)
            for symbol in valid_symbols:
                symbol_mask = (part_seq_matrix == symbol) & part_valid_masks
                symbol_i = symbol_mask.astype(np.float64, copy=False)
                matches += symbol_i @ symbol_i.T

            pair_compared = squareform(compared, checks=False)
            np.divide(
                squareform(matches, checks=False),
                pair_compared,
                out=result[:, p_idx],
                where=pair_compared > 0,
            )

        return part_names, result

    @staticmethod
    def _sequence_byte_matrix(
        sequences: dict[str, str], taxa_names: list[str],
    ) -> np.ndarray:
        seq_len = len(sequences[taxa_names[0]]) if taxa_names else 0
        return np.frombuffer(
            "".join(sequences[name] for name in taxa_names).encode("ascii"),
            dtype=np.uint8,
        ).reshape(len(taxa_names), seq_len)

    def _compute_partition_identities_pairwise(
        self,
        sequences: dict[str, str],
        taxa_names: list[str],
        partitions: list[tuple[str, int, int]],
    ) -> tuple[list[str], np.ndarray]:
        n = len(taxa_names)
        n_parts = len(partitions)
        seq_matrix = np.array([list(sequences[name]) for name in taxa_names], dtype="U1")
        valid_masks = ~np.isin(
            seq_matrix,
            np.array(["-", "?", "N", "X", "*", "n", "x"], dtype="U1"),
        )

        pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
        part_names = [p[0] for p in partitions]
        result = np.zeros((len(pairs), n_parts), dtype=np.float64)

        for p_idx, (pname, start, end) in enumerate(partitions):
            part_seq_matrix = seq_matrix[:, start:end]
            part_valid_masks = valid_masks[:, start:end]
            for pair_idx, (i, j) in enumerate(pairs):
                valid = part_valid_masks[i] & part_valid_masks[j]
                compared = int(np.count_nonzero(valid))
                matches = int(
                    np.count_nonzero((part_seq_matrix[i] == part_seq_matrix[j]) & valid)
                )
                result[pair_idx, p_idx] = matches / compared if compared > 0 else 0.0

        return part_names, result

    def _parse_alignment(self) -> tuple[dict[str, str], list[str], int]:
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

        seq_iter = iter(sequences.values())
        aln_length = len(next(seq_iter))
        for seq in seq_iter:
            if len(seq) != aln_length:
                print(
                    "Error: sequences have different lengths. Is this an alignment?",
                    file=sys.stderr,
                )
                raise SystemExit(2)
        return sequences, taxa_names, aln_length

    def _compute_identity_matrix(
        self,
        sequences: dict[str, str],
        taxa_names: list[str],
    ) -> np.ndarray:
        """Compute pairwise identity matrix. Returns identity values (0-1)."""
        n = len(taxa_names)
        if _all_sequences_identical(sequences, taxa_names):
            identity = _identical_sequence_identity_value(sequences[taxa_names[0]])
            matrix = np.full((n, n), identity, dtype=np.float64)
            np.fill_diagonal(matrix, 1.0)
            return matrix

        try:
            seq_len = len(sequences[taxa_names[0]]) if taxa_names else 0
            joined_bytes = "".join(sequences[name] for name in taxa_names).encode(
                "ascii"
            )
            seq_matrix = np.frombuffer(
                joined_bytes,
                dtype=np.uint8,
            ).reshape(n, seq_len)
        except UnicodeEncodeError:
            return self._compute_identity_matrix_pairwise(sequences, taxa_names)

        direct_clean_ascii = (
            seq_len > 0
            and n <= _NO_INVALID_DIRECT_MAX_TAXA
            and (
                n >= _NO_INVALID_DIRECT_SHORT_ALIGNMENT_MIN_TAXA
                or (
                    n >= _NO_INVALID_DIRECT_MIN_TAXA
                    and seq_len >= _NO_INVALID_DIRECT_MIN_LENGTH
                )
            )
            and not any(code in joined_bytes for code in _IDENTITY_INVALID_BYTES)
        )
        if direct_clean_ascii:
            target_bytes = 16 * 1024 * 1024
            block_size = max(
                1,
                min(64, target_bytes // max(1, n * max(1, seq_len))),
            )
            matrix = np.empty((n, n), dtype=np.float64)
            denominator = float(seq_len)
            for start in range(0, n, block_size):
                stop = min(n, start + block_size)
                match_counts = (
                    seq_matrix[start:stop, None, :] == seq_matrix[None, :, :]
                ).sum(axis=2, dtype=np.int32)
                matrix[start:stop] = match_counts / denominator
            np.fill_diagonal(matrix, 1.0)
            return matrix

        invalid = np.frombuffer(_IDENTITY_INVALID_BYTES, dtype=np.uint8)
        valid_masks = ~np.isin(seq_matrix, invalid)
        # Float masks route dense products through optimized BLAS while
        # preserving exact integer counts for practical alignment lengths.
        valid_i = valid_masks.astype(np.float64, copy=False)
        compared = valid_i @ valid_i.T

        matches = np.zeros((n, n), dtype=np.float64)
        for symbol in np.unique(seq_matrix[valid_masks]):
            symbol_mask = (seq_matrix == symbol) & valid_masks
            symbol_i = symbol_mask.astype(np.float64, copy=False)
            matches += symbol_i @ symbol_i.T

        matrix = np.zeros((n, n), dtype=np.float64)
        np.divide(matches, compared, out=matrix, where=compared > 0)
        np.fill_diagonal(matrix, 1.0)
        return matrix

    def _compute_identity_matrix_pairwise(
        self,
        sequences: dict[str, str],
        taxa_names: list[str],
    ) -> np.ndarray:
        n = len(taxa_names)
        matrix = np.ones((n, n), dtype=np.float64)

        seq_matrix = np.array([list(sequences[name]) for name in taxa_names], dtype="U1")
        valid_masks = ~np.isin(
            seq_matrix,
            np.array(["-", "?", "N", "X", "*", "n", "x"], dtype="U1"),
        )

        for i in range(n):
            seq_i = seq_matrix[i]
            valid_i = valid_masks[i]
            for j in range(i + 1, n):
                valid = valid_i & valid_masks[j]

                compared = int(np.count_nonzero(valid))
                if compared > 0:
                    matches = int(
                        np.count_nonzero((seq_i == seq_matrix[j]) & valid)
                    )
                    identity = matches / compared
                else:
                    identity = 0.0

                matrix[i, j] = identity
                matrix[j, i] = identity

        return matrix

    @staticmethod
    def _partition_identity_strip(
        part_identities: np.ndarray,
        order: list[int],
        n_taxa: int,
        n_parts: int,
    ) -> np.ndarray:
        if n_taxa <= 1:
            return np.zeros((n_taxa, n_parts), dtype=np.float64)

        sums = np.zeros((n_taxa, n_parts), dtype=np.float64)
        cursor = 0
        for idx in range(n_taxa - 1):
            n_remaining = n_taxa - idx - 1
            block = part_identities[cursor:cursor + n_remaining]
            sums[idx] += block.sum(axis=0)
            sums[idx + 1:] += block
            cursor += n_remaining
        return (sums / (n_taxa - 1))[order]

    def _determine_order(
        self,
        matrix: np.ndarray,
        taxa_names: list[str],
        n_taxa: int,
        return_linkage: bool = False,
    ) -> tuple[list[int], list[str]]:
        """Determine taxa ordering based on sort method."""
        cluster_linkage = None
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
            from Bio import Phylo
            from ..tree.base import Tree

            tree = Phylo.read(self.tree_path, "newick")
            tip_order = Tree.calculate_terminal_names_fast(tree)
            if tip_order is None:
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
                cluster_linkage = linkage(condensed, method="average")
                order = list(leaves_list(cluster_linkage))
            else:
                order = list(range(n_taxa))
        else:
            order = list(range(n_taxa))

        ordered_labels = [taxa_names[i] for i in order]
        if return_linkage:
            return order, ordered_labels, cluster_linkage
        return order, ordered_labels

    def _plot_heatmap(
        self,
        matrix: np.ndarray,
        taxa_names: list[str],
        order: list[int],
        ordered_labels: list[str],
        n_taxa: int,
        partitions: list[tuple[str, int, int]] | None = None,
        sequences: dict[str, str] | None = None,
        cluster_linkage: np.ndarray | None = None,
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

            if cluster_linkage is None:
                # Compute linkage for dendrograms when called directly.
                if self.metric == "p-distance":
                    dist_for_linkage = np.copy(matrix)
                else:
                    dist_for_linkage = 1.0 - matrix
                dist_for_linkage = np.clip(dist_for_linkage, 0.0, 1.0)
                np.fill_diagonal(dist_for_linkage, 0.0)
                condensed = squareform(dist_for_linkage, checks=False)
                cluster_linkage = linkage(condensed, method="average")

            # Top dendrogram
            dendrogram(
                cluster_linkage, ax=ax_dendro_top,
                leaf_rotation=90, leaf_font_size=0,
                color_threshold=0, above_threshold_color="#333333",
                no_labels=True,
            )
            ax_dendro_top.axis("off")

            # Left dendrogram (rotated)
            dendrogram(
                cluster_linkage, ax=ax_dendro_left,
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

            # Build reordered per-partition strip: rows = reordered taxa pairs
            # For the strip, show average identity per partition per taxon
            # as a narrow heatmap: rows = taxa (same order), cols = partitions
            # Value = mean identity of that taxon against all others for that partition
            strip = self._partition_identity_strip(
                part_identities, order, n_taxa, n_parts
            )

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
