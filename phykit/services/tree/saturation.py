from __future__ import annotations

import itertools
import os

from .base import Tree


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


def get_alignment_and_format_helper(*args, **kwargs):
    from ...helpers.files import get_alignment_and_format

    return get_alignment_and_format(*args, **kwargs)


class _LazyNumpy:
    def __getattr__(self, name):
        import numpy as _np

        return getattr(_np, name)


np = _LazyNumpy()


def _all_sequences_identical(sequences) -> bool:
    if len(sequences) < 2:
        return True
    first_sequence = sequences[0]
    for idx in range(1, len(sequences)):
        if sequences[idx] != first_sequence:
            return False
    return True


class _LazyMultiprocessing:
    def cpu_count(self):
        import multiprocessing as _mp

        return _mp.cpu_count()

    def Pool(self, *args, **kwargs):
        import multiprocessing as _mp

        return _mp.Pool(*args, **kwargs)


mp = _LazyMultiprocessing()


class Saturation(Tree):
    MP_MIN_COMBOS = 2000
    MAX_MP_WORKERS = 8
    MAX_MATRIX_DISTANCE_TAXA = 1000
    _BYTE_GAP_LOOKUP_CACHE = {}

    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(
            tree_file_path=parsed["tree_file_path"],
            alignment_file_path=parsed["alignment_file_path"],
            exclude_gaps=parsed["exclude_gaps"],
            verbose=parsed["verbose"],
        )
        self.json_output = parsed["json_output"]
        self.plot = parsed["plot"]
        self.plot_output = parsed["plot_output"]
        self.plot_config = parsed["plot_config"]

    def _should_use_multiprocessing(self, n_combos: int) -> bool:
        if os.environ.get("PHYKIT_DISABLE_MP", "0") == "1":
            return False
        if os.environ.get("PHYKIT_FORCE_MP", "0") == "1":
            return True
        return n_combos >= self.MP_MIN_COMBOS

    def run(self) -> None:
        alignment, _, is_protein = get_alignment_and_format_helper(
            self.alignment_file_path
        )

        tree = self.read_tree_file_unmodified()

        tips = self.get_tip_names_from_tree(tree)
        combos = list(itertools.combinations(tips, 2))

        (
            patristic_distances,
            uncorrected_distances,
        ) = self.loop_through_combos_and_calculate_pds_and_pis(
            combos, alignment, tree, self.exclude_gaps, is_protein
        )

        # calculate slope while fitting the y-intercept to zero.
        # This follows Jeffroy et al. (fig 2):
        # https://www.cell.com/trends/genetics/fulltext/S0168-9525(06)00051-5
        x = np.asarray(patristic_distances, dtype=float)
        y = np.asarray(uncorrected_distances, dtype=float)
        finite_mask = np.isfinite(x) & np.isfinite(y)
        x = x[finite_mask]
        y = y[finite_mask]

        denom = float(np.dot(x, x))
        slope = float(np.dot(x, y) / denom) if denom != 0.0 else 0.0

        if self.plot:
            self._plot_saturation_scatter(x, y, slope)

        self.print_res(
            self.verbose, combos, uncorrected_distances, patristic_distances, slope
        )

    def process_args(self, args) -> dict[str, str]:
        from ...helpers.plot_config import PlotConfig

        return dict(
            tree_file_path=args.tree,
            alignment_file_path=args.alignment,
            exclude_gaps=args.exclude_gaps,
            verbose=args.verbose,
            json_output=getattr(args, "json", False),
            plot=getattr(args, "plot", False),
            plot_output=getattr(args, "plot_output", "saturation_plot.png"),
            plot_config=PlotConfig.from_args(args),
        )

    def _plot_saturation_scatter(
        self,
        patristic_distances: np.ndarray,
        uncorrected_distances: np.ndarray,
        slope: float,
    ) -> None:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        except ImportError:
            print("matplotlib is required for --plot in saturation. Install matplotlib and retry.")
            raise SystemExit(2)

        config = self.plot_config
        config.resolve(n_rows=len(patristic_distances), n_cols=None)
        default_colors = ["#2b8cbe", "#000000"]
        colors = config.merge_colors(default_colors)

        fig, ax = plt.subplots(figsize=(config.fig_width, config.fig_height))
        ax.scatter(
            patristic_distances,
            uncorrected_distances,
            s=14,
            alpha=0.6,
            color=colors[0],
            edgecolors="none",
        )

        if patristic_distances.size > 0:
            x_line = np.linspace(0.0, float(np.max(patristic_distances)), 200)
            y_line = slope * x_line
            ax.plot(
                x_line,
                y_line,
                color=colors[1],
                linestyle="--",
                linewidth=2.0,
                label=f"Fit through origin (slope={slope:.4f})",
            )
            legend_loc = config.legend_position or "best"
            if legend_loc != "none":
                ax.legend(loc=legend_loc, frameon=False)

        if config.show_title:
            ax.set_title(config.title or "Saturation: Patristic vs Uncorrected Distance", fontsize=config.title_fontsize)
        ax.set_xlabel("Patristic distance")
        ax.set_ylabel("Uncorrected distance")
        if config.axis_fontsize:
            ax.xaxis.label.set_fontsize(config.axis_fontsize)
            ax.yaxis.label.set_fontsize(config.axis_fontsize)
        fig.tight_layout()
        fig.savefig(self.plot_output, dpi=config.dpi, bbox_inches="tight")
        plt.close(fig)

    def _process_combo_batch(self, tree, seq_arrays, gap_mask, exclude_gaps, combo_batch):
        """Process a batch of combinations in parallel."""
        results = []
        for combo in combo_batch:
            # Calculate patristic distance
            pd = tree.distance(combo[0], combo[1])

            # Calculate uncorrected distance using numpy operations
            seq1_arr = seq_arrays[combo[0]]
            seq2_arr = seq_arrays[combo[1]]

            if exclude_gaps:
                # Use pre-computed gap masks
                gap_mask1 = gap_mask[combo[0]]
                gap_mask2 = gap_mask[combo[1]]
                valid_positions = ~(gap_mask1 | gap_mask2)
                adjusted_len = int(np.count_nonzero(valid_positions))

                if adjusted_len:
                    identities = int(
                        np.count_nonzero((seq1_arr == seq2_arr) & valid_positions)
                    )
                    ud = 1 - (identities / adjusted_len)
                else:
                    ud = float('nan')
            else:
                identities = int(np.count_nonzero(seq1_arr == seq2_arr))
                ud = 1 - (identities / len(seq1_arr))

            results.append((pd, ud))
        return results

    @staticmethod
    def _sequence_to_array(sequence):
        sequence = str(sequence).upper()
        try:
            return np.frombuffer(sequence.encode("ascii"), dtype=np.uint8)
        except UnicodeEncodeError:
            return np.array(list(sequence), dtype="U1")

    @staticmethod
    def _gap_values_for_array(seq_arr, gap_chars):
        gap_chars = {char.upper() for char in gap_chars}
        if seq_arr.dtype == np.uint8:
            return np.fromiter((ord(char) for char in gap_chars), dtype=np.uint8)
        if seq_arr.dtype.kind == "S":
            return [char.encode("ascii") for char in gap_chars]
        return list(gap_chars)

    @classmethod
    def _gap_mask_for_array(cls, seq_arr, gap_chars):
        gap_chars = frozenset(char.upper() for char in gap_chars)
        if seq_arr.dtype == np.uint8:
            lookup = cls._BYTE_GAP_LOOKUP_CACHE.get(gap_chars)
            if lookup is None:
                lookup = np.zeros(256, dtype=bool)
                lookup[
                    np.fromiter((ord(char) for char in gap_chars), dtype=np.uint8)
                ] = True
                cls._BYTE_GAP_LOOKUP_CACHE[gap_chars] = lookup
            return lookup[seq_arr]
        return np.isin(seq_arr, cls._gap_values_for_array(seq_arr, gap_chars))

    @staticmethod
    def _valid_site_count_for_identical_sequence(sequence: str, gap_chars) -> int:
        gap_chars = tuple({char.upper() for char in gap_chars})
        try:
            sequence_bytes = sequence.encode("ascii")
            gap_bytes = bytes(ord(char) for char in gap_chars)
            return len(sequence_bytes.translate(None, gap_bytes))
        except UnicodeEncodeError:
            return len(sequence) - sum(sequence.count(char) for char in gap_chars)

    @classmethod
    def _constant_uncorrected_distance_for_identical_sequences(
        cls,
        alignment,
        combo_tips: list[str],
        exclude_gaps: bool,
        gap_chars,
    ) -> float | None:
        if not combo_tips:
            return None

        sequences_by_tip = {
            record.name: str(record.seq).upper()
            for record in alignment
        }
        try:
            sequences = [sequences_by_tip[tip] for tip in combo_tips]
        except KeyError:
            return None

        if not _all_sequences_identical(sequences):
            return None

        first_sequence = sequences[0]
        if exclude_gaps:
            valid_sites = cls._valid_site_count_for_identical_sequence(
                first_sequence,
                gap_chars,
            )
            return float("nan") if valid_sites == 0 else 0.0

        return float("nan") if len(first_sequence) == 0 else 0.0

    @classmethod
    def _calculate_uncorrected_distances_matrix(
        cls,
        combo_tips: list[str],
        combos: list[tuple[str, str]],
        seq_arrays: dict[str, np.ndarray],
        gap_mask: dict[str, np.ndarray],
        exclude_gaps: bool,
    ) -> list[float] | None:
        """Compute all requested uncorrected distances from sequence masks."""
        n_tips = len(combo_tips)
        if n_tips == 0:
            return []
        if n_tips > cls.MAX_MATRIX_DISTANCE_TAXA:
            return None

        try:
            arrays = [seq_arrays[tip] for tip in combo_tips]
        except KeyError:
            return None
        if any(seq_arr.dtype != np.uint8 for seq_arr in arrays):
            return None

        seq_len = len(arrays[0])
        if any(len(seq_arr) != seq_len for seq_arr in arrays):
            return None

        seq_matrix = np.vstack(arrays)
        if not exclude_gaps:
            return cls._calculate_uncorrected_distances_no_gap_matrix(
                combo_tips,
                combos,
                seq_matrix,
            )

        try:
            if not any(gap_mask[tip].any() for tip in combo_tips):
                return cls._calculate_uncorrected_distances_no_gap_matrix(
                    combo_tips,
                    combos,
                    seq_matrix,
                )
        except KeyError:
            return None

        valid_matrix = ~np.vstack([gap_mask[tip] for tip in combo_tips])

        valid_float = valid_matrix.astype(np.float64)
        adjusted_lengths = valid_float @ valid_float.T
        identity_counts = np.zeros(adjusted_lengths.shape, dtype=np.float64)
        for symbol in np.unique(seq_matrix[valid_matrix]):
            symbol_mask = ((seq_matrix == symbol) & valid_matrix).astype(np.float64)
            identity_counts += symbol_mask @ symbol_mask.T

        standard_distances = cls._standard_upper_triangle_gappy_distances(
            combo_tips,
            combos,
            identity_counts,
            adjusted_lengths,
        )
        if standard_distances is not None:
            return standard_distances

        tip_indices = {tip: idx for idx, tip in enumerate(combo_tips)}
        uncorrected_distances = []
        for tip_a, tip_b in combos:
            idx_a = tip_indices[tip_a]
            idx_b = tip_indices[tip_b]
            adjusted_len = adjusted_lengths[idx_a, idx_b]
            if adjusted_len == 0.0:
                uncorrected_distances.append(float("nan"))
            else:
                uncorrected_distances.append(
                    1.0 - (identity_counts[idx_a, idx_b] / adjusted_len)
                )
        return uncorrected_distances

    @classmethod
    def _calculate_uncorrected_distances_no_gap_matrix(
        cls,
        combo_tips: list[str],
        combos: list[tuple[str, str]],
        seq_matrix,
    ) -> list[float]:
        n_tips, seq_len = seq_matrix.shape
        if seq_len == 0:
            tip_indices = {tip: idx for idx, tip in enumerate(combo_tips)}
            return [float("nan") for _ in combos]

        target_bytes = 16 * 1024 * 1024
        block_size = max(
            1,
            min(64, target_bytes // max(1, n_tips * max(1, seq_len))),
        )
        distances = np.empty((n_tips, n_tips), dtype=np.float64)
        denominator = float(seq_len)
        for start in range(0, n_tips, block_size):
            stop = min(n_tips, start + block_size)
            identity_counts = (
                seq_matrix[start:stop, None, :] == seq_matrix[None, :, :]
            ).sum(axis=2, dtype=np.int32)
            distances[start:stop] = 1.0 - (identity_counts / denominator)

        standard_distances = cls._standard_upper_triangle_values(
            combo_tips,
            combos,
            distances,
        )
        if standard_distances is not None:
            return standard_distances

        tip_indices = {tip: idx for idx, tip in enumerate(combo_tips)}
        return [
            float(distances[tip_indices[tip_a], tip_indices[tip_b]])
            for tip_a, tip_b in combos
        ]

    @staticmethod
    def _combos_are_standard_upper_triangle(
        combo_tips: list[str],
        combos: list[tuple[str, str]],
    ) -> bool:
        n_tips = len(combo_tips)
        if len(combos) != n_tips * (n_tips - 1) // 2:
            return False

        cursor = 0
        for idx, tip_a in enumerate(combo_tips[:-1]):
            for tip_b in combo_tips[idx + 1:]:
                if combos[cursor] != (tip_a, tip_b):
                    return False
                cursor += 1
        return True

    @classmethod
    def _standard_upper_triangle_values(
        cls,
        combo_tips: list[str],
        combos: list[tuple[str, str]],
        values,
    ) -> list[float] | None:
        if not cls._combos_are_standard_upper_triangle(combo_tips, combos):
            return None

        distances = []
        extend = distances.extend
        for idx in range(len(combo_tips) - 1):
            extend(values[idx, idx + 1:].tolist())
        return distances

    @classmethod
    def _standard_upper_triangle_gappy_distances(
        cls,
        combo_tips: list[str],
        combos: list[tuple[str, str]],
        identity_counts,
        adjusted_lengths,
    ) -> list[float] | None:
        if not cls._combos_are_standard_upper_triangle(combo_tips, combos):
            return None

        if (adjusted_lengths != 0.0).all():
            distances = []
            extend = distances.extend
            for idx in range(len(combo_tips) - 1):
                extend(
                    (
                        1.0
                        - (
                            identity_counts[idx, idx + 1:]
                            / adjusted_lengths[idx, idx + 1:]
                        )
                    ).tolist()
                )
            return distances

        valid = adjusted_lengths != 0.0
        values = np.empty(adjusted_lengths.shape, dtype=np.float64)
        values[valid] = 1.0 - (identity_counts[valid] / adjusted_lengths[valid])
        values[~valid] = np.nan

        distances = []
        extend = distances.extend
        for idx in range(len(combo_tips) - 1):
            extend(values[idx, idx + 1:].tolist())
        return distances

    def loop_through_combos_and_calculate_pds_and_pis(
        self,
        combos: list[tuple[str, str]],
        alignment: Align.MultipleSeqAlignment,
        tree: Newick.Tree,
        exclude_gaps: bool,
        is_protein: bool = False,
    ) -> tuple[
        list[float],
        list[float],
    ]:
        """
        loop through all taxon combinations and determine
        their patristic distance and pairwise identity
        """
        gap_chars = self.get_gap_chars(is_protein)
        combo_tips = []
        seen_tips = set()
        for tip_a, tip_b in combos:
            if tip_a not in seen_tips:
                combo_tips.append(tip_a)
                seen_tips.add(tip_a)
            if tip_b not in seen_tips:
                combo_tips.append(tip_b)
                seen_tips.add(tip_b)
        fast_pair_distances = None
        direct_pair_distances = None
        try:
            fast_result = self.calculate_pairwise_tip_distances_fast(tree, combo_tips)
        except (AttributeError, TypeError, KeyError):
            fast_result = None
        if fast_result is not None:
            fast_combos, fast_distances = fast_result
            if fast_combos == combos:
                direct_pair_distances = fast_distances
            else:
                fast_pair_distances = {
                    frozenset((tip_a, tip_b)): distance
                    for (tip_a, tip_b), distance in zip(fast_combos, fast_distances)
                }

        if direct_pair_distances is not None or fast_pair_distances is not None:
            constant_uncorrected_distance = (
                self._constant_uncorrected_distance_for_identical_sequences(
                    alignment,
                    combo_tips,
                    exclude_gaps,
                    gap_chars,
                )
            )
            if constant_uncorrected_distance is not None:
                if direct_pair_distances is not None:
                    patristic_distances = list(direct_pair_distances)
                else:
                    patristic_distances = [
                        fast_pair_distances[frozenset(combo)]
                        for combo in combos
                    ]
                return (
                    patristic_distances,
                    [constant_uncorrected_distance] * len(combos),
                )

        # Convert sequences to numpy arrays for vectorized operations
        seq_arrays = {}
        gap_mask = {}
        for record in alignment:
            seq_arr = self._sequence_to_array(record.seq)
            seq_arrays[record.name] = seq_arr
            if exclude_gaps:
                gap_mask[record.name] = self._gap_mask_for_array(
                    seq_arr, gap_chars
                )

        if direct_pair_distances is not None:
            matrix_distances = self._calculate_uncorrected_distances_matrix(
                combo_tips,
                combos,
                seq_arrays,
                gap_mask,
                exclude_gaps,
            )
            if matrix_distances is not None:
                return list(direct_pair_distances), matrix_distances

        # For small/medium workloads, multiprocessing overhead dominates.
        if (
            direct_pair_distances is not None
            or fast_pair_distances is not None
            or not self._should_use_multiprocessing(len(combos))
        ):
            patristic_distances = []
            uncorrected_distances = []
            for combo_idx, combo in enumerate(combos):
                if fast_pair_distances is None and direct_pair_distances is None:
                    pd = tree.distance(combo[0], combo[1])
                elif direct_pair_distances is not None:
                    pd = direct_pair_distances[combo_idx]
                else:
                    pd = fast_pair_distances[frozenset(combo)]
                patristic_distances.append(pd)

                seq1_arr = seq_arrays[combo[0]]
                seq2_arr = seq_arrays[combo[1]]

                if exclude_gaps:
                    gap_mask1 = gap_mask[combo[0]]
                    gap_mask2 = gap_mask[combo[1]]
                    valid_positions = ~(gap_mask1 | gap_mask2)
                    adjusted_len = int(np.count_nonzero(valid_positions))

                    if adjusted_len:
                        identities = int(
                            np.count_nonzero((seq1_arr == seq2_arr) & valid_positions)
                        )
                        ud = 1 - (identities / adjusted_len)
                    else:
                        ud = float('nan')
                else:
                    identities = int(np.count_nonzero(seq1_arr == seq2_arr))
                    ud = 1 - (identities / len(seq1_arr))

                uncorrected_distances.append(ud)
        else:
            # Use multiprocessing for larger datasets
            from functools import partial

            num_workers = min(mp.cpu_count(), self.MAX_MP_WORKERS)
            chunk_size = max(1, len(combos) // (num_workers * 4))
            combo_chunks = [combos[i:i + chunk_size] for i in range(0, len(combos), chunk_size)]

            # Create partial function
            process_func = partial(
                self._process_combo_batch,
                tree,
                seq_arrays,
                gap_mask,
                exclude_gaps
            )

            # Process in parallel
            with mp.Pool(processes=num_workers) as pool:
                chunk_results = pool.map(process_func, combo_chunks)

            # Flatten results
            patristic_distances = []
            uncorrected_distances = []
            for chunk_result in chunk_results:
                for pd, ud in chunk_result:
                    patristic_distances.append(pd)
                    uncorrected_distances.append(ud)

        return patristic_distances, uncorrected_distances

    def print_res(
        self,
        verbose: bool,
        combos: list[str],
        uncorrected_distances: list[float],
        patristic_distances: list[float],
        slope: float,
    ) -> None:
        """
        print results to stdout
        """
        try:
            if self.json_output:
                payload = dict(verbose=verbose, exclude_gaps=self.exclude_gaps)
                if verbose:
                    rows = [
                        dict(
                            taxon_a=cbo[0],
                            taxon_b=cbo[1],
                            uncorrected_distance=round(dist, 4),
                            patristic_distance=round(pd, 4),
                        )
                        for cbo, dist, pd in zip(
                            combos, uncorrected_distances, patristic_distances
                        )
                    ]
                    payload["rows"] = rows
                    payload["pairs"] = rows
                else:
                    payload["summary"] = dict(
                        slope=round(slope, 4),
                        one_minus_slope_abs=abs(round(1 - slope, 4)),
                    )
                if self.plot:
                    payload["plot_output"] = self.plot_output
                print_json(payload)
                return

            if verbose:
                lines = [
                    f"{cbo[0]}\t{cbo[1]}\t{round(dist,4)}\t{round(pd, 4)}"
                    for cbo, dist, pd in zip(
                        combos, uncorrected_distances, patristic_distances
                    )
                ]
                if lines:
                    print("\n".join(lines))
            else:
                print(f"{round(slope, 4)}\t{abs(round(1-slope, 4))}")
            if self.plot:
                print(f"Saved saturation plot: {self.plot_output}")
        except BrokenPipeError:
            pass
