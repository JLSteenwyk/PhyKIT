from __future__ import annotations

from .base import Alignment, _all_sequences_identical


class _LazyNumpy:
    def __getattr__(self, name):
        import numpy as _np

        return getattr(_np, name)


np = _LazyNumpy()
_SITE_ENTROPY_DIRECT_SUM_MAX_SIZE = 100_000


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


def _row_l2_norms(matrix):
    return np.sqrt(np.einsum("ij,ij->i", matrix, matrix))


def _row_sums(matrix):
    return matrix.sum(axis=1)


def _column_dot(left, right):
    return np.einsum("ij,ij->j", left, right)


def _average_site_entropy(site_entropies, aln_len: int) -> float:
    if site_entropies.size <= _SITE_ENTROPY_DIRECT_SUM_MAX_SIZE:
        return float(site_entropies.sum() / aln_len)
    return float(np.sum(site_entropies) / aln_len)


def _rcvt_from_composition_matrix(comp_matrix, valid_lengths, n_taxa: int):
    if comp_matrix.shape[1] == 0:
        return np.zeros(comp_matrix.shape[0], dtype=np.float64)

    average_counts = comp_matrix.sum(axis=0) / n_taxa
    deviations = np.abs(comp_matrix - average_counts)
    seq_sums = deviations.sum(axis=1)
    denom = n_taxa * valid_lengths
    return np.divide(
        seq_sums,
        denom,
        out=np.zeros_like(seq_sums, dtype=np.float64),
        where=denom > 0,
    )


class AlignmentOutlierTaxa(Alignment):
    _INVALID_LOOKUP_CACHE = {}
    _PAIRWISE_MATRIX_MAX_CELLS = 4_000_000
    _BLOCKED_SYMBOL_MASK_MAX_CELLS = 20_000_000

    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(alignment_file_path=parsed["alignment_file_path"])
        self.gap_z = parsed["gap_z"]
        self.composition_z = parsed["composition_z"]
        self.distance_z = parsed["distance_z"]
        self.rcvt_z = parsed["rcvt_z"]
        self.occupancy_z = parsed["occupancy_z"]
        self.entropy_z = parsed["entropy_z"]
        self.json_output = parsed["json_output"]

    def process_args(self, args) -> dict[str, str]:
        return dict(
            alignment_file_path=args.alignment,
            gap_z=args.gap_z,
            composition_z=args.composition_z,
            distance_z=args.distance_z,
            rcvt_z=args.rcvt_z,
            occupancy_z=args.occupancy_z,
            entropy_z=args.entropy_z,
            json_output=getattr(args, "json", False),
        )

    def _high_outlier_threshold(self, values: np.ndarray, z: float) -> float:
        finite_values = values[np.isfinite(values)]
        if finite_values.size < 2:
            return float("inf")

        median = float(np.median(finite_values))
        mad = float(np.median(np.abs(finite_values - median)))
        if mad > 0.0:
            robust_sigma = 1.4826 * mad
            return median + (z * robust_sigma)

        # Fallback when variation is very low: flag values above second-highest.
        unique_vals = np.unique(finite_values)
        if unique_vals.size > 1:
            return float(unique_vals[-2])
        return float("inf")

    def _low_outlier_threshold(self, values: np.ndarray, z: float) -> float:
        finite_values = values[np.isfinite(values)]
        if finite_values.size < 2:
            return float("-inf")

        median = float(np.median(finite_values))
        mad = float(np.median(np.abs(finite_values - median)))
        if mad > 0.0:
            robust_sigma = 1.4826 * mad
            return median - (z * robust_sigma)

        # Fallback when variation is very low: flag values below second-lowest.
        unique_vals = np.unique(finite_values)
        if unique_vals.size > 1:
            return float(unique_vals[1])
        return float("-inf")

    @classmethod
    def _invalid_lookup_for_chars(cls, invalid_chars):
        invalid_chars = frozenset(char.upper() for char in invalid_chars)
        lookup = cls._INVALID_LOOKUP_CACHE.get(invalid_chars)
        if lookup is None:
            lookup = np.zeros(256, dtype=bool)
            lookup[
                np.fromiter((ord(char) for char in invalid_chars), dtype=np.uint8)
            ] = True
            cls._INVALID_LOOKUP_CACHE[invalid_chars] = lookup
        return lookup

    @staticmethod
    def _constant_composition_result(
        alignment,
        valid_length: float,
        aln_len: int,
    ) -> dict[str, object]:
        n_taxa = len(alignment)
        if aln_len > 0:
            gap_rate = 1.0 - (valid_length / float(aln_len))
            occupancy = valid_length / float(aln_len)
        else:
            gap_rate = 1.0
            occupancy = 0.0
        branch_proxy = 0.0 if n_taxa > 1 else None
        rows = [
            {
                "taxon": record.id,
                "gap_rate": round(float(gap_rate), 4),
                "occupancy": round(float(occupancy), 4),
                "composition_distance": 0.0,
                "long_branch_proxy": branch_proxy,
                "rcvt": 0.0,
                "entropy_burden": 0.0,
                "flagged": False,
                "reasons": [],
            }
            for record in alignment
        ]
        return dict(
            features=[
                "gap_rate",
                "occupancy",
                "composition_distance",
                "long_branch_proxy",
                "rcvt",
                "entropy_burden",
            ],
            thresholds=dict(
                gap_rate=None,
                occupancy=None,
                composition_distance=None,
                long_branch_proxy=None,
                rcvt=None,
                entropy_burden=None,
            ),
            rows=rows,
            outliers=[],
        )

    @staticmethod
    def _valid_length_for_identical_sequence(
        sequence: str,
        invalid_chars: list[str],
    ) -> int:
        try:
            sequence_bytes = sequence.encode("ascii")
            invalid_bytes = bytes(ord(char) for char in invalid_chars)
            return len(sequence_bytes.translate(None, invalid_bytes))
        except UnicodeEncodeError:
            return len(sequence) - sum(sequence.count(char) for char in invalid_chars)

    @staticmethod
    def _symbol_counts_by_row(alignment_array, symbols):
        if alignment_array.dtype == np.uint8 and symbols.size >= 16:
            n_rows, aln_len = alignment_array.shape
            max_code = int(symbols.max()) + 1
            if alignment_array.size <= 2_000_000 or (
                n_rows >= 15_000 and aln_len <= 256
            ):
                encoded = alignment_array.astype(np.int64)
                encoded += (np.arange(n_rows, dtype=np.int64) * max_code)[:, None]
                counts = np.bincount(
                    encoded.ravel(),
                    minlength=n_rows * max_code,
                ).reshape(n_rows, max_code)
                return counts[:, symbols].astype(np.float64, copy=False)

            counts = np.empty((n_rows, symbols.size), dtype=np.float64)
            for row_idx, row in enumerate(alignment_array):
                counts[row_idx] = np.bincount(row, minlength=max_code)[symbols]
            return counts

        return np.array(
            [
                np.sum(alignment_array == symbol, axis=1)
                for symbol in symbols
            ],
            dtype=np.float64,
        ).T

    @staticmethod
    def _symbol_counts_by_site(alignment_array, symbols):
        if alignment_array.dtype == np.uint8 and alignment_array.size <= 32_000_000:
            n_sites = alignment_array.shape[1]
            max_code = int(symbols.max()) + 1
            encoded = alignment_array.astype(np.int64)
            encoded += np.arange(n_sites, dtype=np.int64) * max_code
            counts = np.bincount(
                encoded.ravel(),
                minlength=n_sites * max_code,
            ).reshape(n_sites, max_code)
            return counts[:, symbols].T.astype(np.float64, copy=False)

        return np.array(
            [
                np.sum(alignment_array == symbol, axis=0)
                for symbol in symbols
            ],
            dtype=np.float64,
        )

    @staticmethod
    def _all_valid_long_branch_proxy(alignment_array, symbols, site_counts):
        n_taxa, aln_len = alignment_array.shape
        long_branch_proxy = np.full(n_taxa, np.nan, dtype=np.float64)
        if n_taxa <= 1 or aln_len == 0:
            return long_branch_proxy

        symbol_lookup = np.empty(256, dtype=np.int16)
        symbol_lookup.fill(-1)
        symbol_lookup[symbols] = np.arange(symbols.size, dtype=np.int16)
        symbol_indices = symbol_lookup[alignment_array]
        same_symbol_counts = site_counts[
            symbol_indices,
            np.arange(aln_len, dtype=np.intp),
        ]
        if aln_len >= 1000:
            other_matches = same_symbol_counts.sum(axis=1) - aln_len
        else:
            other_matches = np.sum(same_symbol_counts, axis=1) - aln_len
        denom = float(aln_len * (n_taxa - 1))
        long_branch_proxy[:] = 1.0 - (other_matches / denom)
        return long_branch_proxy

    @classmethod
    def _blocked_long_branch_proxy(cls, alignment_array, valid_mask, symbols):
        n_taxa, aln_len = alignment_array.shape
        long_branch_proxy = np.full(n_taxa, np.nan, dtype=np.float64)
        if n_taxa <= 1 or aln_len == 0:
            return long_branch_proxy

        symbol_mask_cells = int(symbols.size) * int(alignment_array.size)
        if symbol_mask_cells > cls._BLOCKED_SYMBOL_MASK_MAX_CELLS:
            return None

        valid_float = valid_mask.astype(np.float64)
        symbol_masks = [
            ((alignment_array == symbol) & valid_mask).astype(np.float64)
            for symbol in symbols
        ]
        block_size = max(1, min(n_taxa, cls._PAIRWISE_MATRIX_MAX_CELLS // n_taxa))
        for start in range(0, n_taxa, block_size):
            stop = min(start + block_size, n_taxa)
            overlap_counts = valid_float[start:stop] @ valid_float.T
            mismatch_counts = overlap_counts.copy()
            for symbol_mask in symbol_masks:
                mismatch_counts -= symbol_mask[start:stop] @ symbol_mask.T

            comparable = overlap_counts > 0
            rows = np.arange(stop - start, dtype=np.intp)
            comparable[rows, start + rows] = False
            comparable_counts = np.count_nonzero(comparable, axis=1)
            mismatch_counts[~comparable] = 0.0
            np.divide(
                mismatch_counts,
                overlap_counts,
                out=mismatch_counts,
                where=comparable,
            )
            long_branch_proxy[start:stop] = np.divide(
                np.sum(mismatch_counts, axis=1),
                comparable_counts,
                out=np.full(stop - start, np.nan, dtype=np.float64),
                where=comparable_counts > 0,
            )
        return long_branch_proxy

    def calculate_outliers(self, alignment, is_protein: bool) -> dict[str, object]:
        taxa = [record.id for record in alignment]
        n_taxa = len(taxa)
        if n_taxa == 0:
            return dict(
                features=[
                    "gap_rate",
                    "occupancy",
                    "composition_distance",
                    "long_branch_proxy",
                    "rcvt",
                    "entropy_burden",
                ],
                thresholds=dict(
                    gap_rate=float("inf"),
                    occupancy=float("-inf"),
                    composition_distance=float("inf"),
                    long_branch_proxy=float("inf"),
                    rcvt=float("inf"),
                    entropy_burden=float("inf"),
                ),
                rows=[],
                outliers=[],
            )

        invalid_chars = (
            ["-", "?", "*", "X"]
            if is_protein
            else ["-", "?", "*", "X", "N"]
        )
        aln_len = alignment.get_alignment_length()
        if n_taxa == 1:
            valid_length = self._valid_length_for_identical_sequence(
                str(alignment[0].seq).upper(),
                invalid_chars,
            )
            return self._constant_composition_result(
                alignment,
                float(valid_length),
                aln_len,
            )

        sequences = [str(record.seq).upper() for record in alignment]
        if _all_sequences_identical(sequences):
            first_sequence = sequences[0]
            valid_length = self._valid_length_for_identical_sequence(
                first_sequence,
                invalid_chars,
            )
            if valid_length > 0:
                return self._constant_composition_result(
                    alignment,
                    float(valid_length),
                    aln_len,
                )

        all_valid_ascii = False
        try:
            joined_bytes = "".join(sequences).encode("ascii")
            alignment_array = np.frombuffer(
                joined_bytes,
                dtype=np.uint8,
            ).reshape(n_taxa, aln_len)
            invalid_bytes = bytes(ord(char) for char in invalid_chars)
            all_valid_ascii = not any(code in joined_bytes for code in invalid_bytes)
            if all_valid_ascii:
                valid_lengths = np.full(n_taxa, aln_len, dtype=np.float64)
            else:
                invalid_lookup = self._invalid_lookup_for_chars(invalid_chars)
                valid_mask = ~invalid_lookup[alignment_array]
        except UnicodeEncodeError:
            alignment_array = np.array([list(seq) for seq in sequences], dtype="U1")
            invalid_chars_array = np.array(invalid_chars, dtype="U1")
            valid_mask = ~np.isin(alignment_array, invalid_chars_array)
        if not all_valid_ascii:
            valid_lengths = np.count_nonzero(valid_mask, axis=1).astype(np.float64)

        if aln_len > 0:
            gap_rates = 1.0 - (valid_lengths / float(aln_len))
            occupancies = valid_lengths / float(aln_len)
        else:
            gap_rates = np.ones(n_taxa, dtype=np.float64)
            occupancies = np.zeros(n_taxa, dtype=np.float64)

        if all_valid_ascii:
            valid_chars = alignment_array.reshape(-1)
        else:
            valid_chars = alignment_array[valid_mask]
        if valid_chars.size > 0:
            symbols = np.unique(valid_chars)
            if symbols.size == 1 and (valid_lengths == valid_lengths[0]).all():
                return self._constant_composition_result(
                    alignment,
                    float(valid_lengths[0]),
                    aln_len,
                )
            symbol_counts = self._symbol_counts_by_row(alignment_array, symbols)
            comp_matrix = np.divide(
                symbol_counts,
                valid_lengths[:, None],
                out=np.zeros_like(symbol_counts, dtype=np.float64),
                where=valid_lengths[:, None] > 0,
            )
            center = np.median(comp_matrix, axis=0)
            composition_distances = _row_l2_norms(comp_matrix - center)
        else:
            composition_distances = np.zeros(n_taxa, dtype=np.float64)
            comp_matrix = np.zeros((n_taxa, 0), dtype=np.float64)

        long_branch_proxy = np.full(n_taxa, np.nan, dtype=np.float64)
        pairwise_cells = n_taxa * n_taxa
        site_counts = None
        if valid_chars.size > 0 and all_valid_ascii:
            site_counts = self._symbol_counts_by_site(alignment_array, symbols)
            long_branch_proxy = self._all_valid_long_branch_proxy(
                alignment_array,
                symbols,
                site_counts,
            )
        elif (
            valid_chars.size > 0
            and pairwise_cells <= self._PAIRWISE_MATRIX_MAX_CELLS
        ):
            valid_float = valid_mask.astype(np.float64)
            overlap_counts = valid_float @ valid_float.T
            match_counts = np.zeros_like(overlap_counts, dtype=np.float64)
            for symbol in symbols:
                symbol_mask = ((alignment_array == symbol) & valid_mask).astype(
                    np.float64
                )
                match_counts += symbol_mask @ symbol_mask.T

            comparable = overlap_counts > 0
            np.fill_diagonal(comparable, False)
            comparable_counts = np.count_nonzero(comparable, axis=1)
            np.subtract(overlap_counts, match_counts, out=match_counts)
            match_counts[~comparable] = 0.0
            np.divide(
                match_counts,
                overlap_counts,
                out=match_counts,
                where=comparable,
            )
            long_branch_proxy = np.divide(
                _row_sums(match_counts),
                comparable_counts,
                out=np.full(n_taxa, np.nan, dtype=np.float64),
                where=comparable_counts > 0,
            )
        else:
            blocked_proxy = None
            if valid_chars.size > 0 and not all_valid_ascii:
                blocked_proxy = self._blocked_long_branch_proxy(
                    alignment_array,
                    valid_mask,
                    symbols,
                )
            if blocked_proxy is not None:
                long_branch_proxy = blocked_proxy
            else:
                for i in range(n_taxa):
                    if all_valid_ascii:
                        comparable = np.ones(n_taxa, dtype=bool)
                        comparable[i] = False
                        if comparable.any():
                            mismatches = np.sum(
                                alignment_array != alignment_array[i],
                                axis=1,
                            )
                            long_branch_proxy[i] = float(
                                np.mean(mismatches[comparable] / float(aln_len))
                            )
                    else:
                        overlap = valid_mask & valid_mask[i]
                        overlap_counts = np.sum(overlap, axis=1)
                        mismatches = np.sum(
                            (alignment_array != alignment_array[i]) & overlap,
                            axis=1,
                        )
                        comparable = overlap_counts > 0
                        comparable[i] = False
                        if comparable.any():
                            distances = (
                                mismatches[comparable] / overlap_counts[comparable]
                            )
                            long_branch_proxy[i] = float(np.mean(distances))

        gap_threshold = self._high_outlier_threshold(gap_rates, self.gap_z)
        composition_threshold = self._high_outlier_threshold(
            composition_distances, self.composition_z
        )
        distance_threshold = self._high_outlier_threshold(
            long_branch_proxy, self.distance_z
        )
        rcvt_values = _rcvt_from_composition_matrix(
            comp_matrix,
            valid_lengths,
            n_taxa,
        )
        rcvt_threshold = self._high_outlier_threshold(rcvt_values, self.rcvt_z)
        occupancy_threshold = self._low_outlier_threshold(occupancies, self.occupancy_z)

        # Entropy burden: average site entropy over valid positions for each taxon.
        if aln_len > 0:
            site_entropies = np.zeros(aln_len, dtype=np.float64)
            if valid_chars.size > 0:
                if site_counts is None:
                    site_counts = self._symbol_counts_by_site(
                        alignment_array,
                        symbols,
                    )
                site_totals = np.sum(site_counts, axis=0)
                site_probs = np.divide(
                    site_counts,
                    site_totals,
                    out=np.zeros_like(site_counts, dtype=np.float64),
                    where=site_totals > 0,
                )
                log_probs = np.zeros_like(site_probs, dtype=np.float64)
                positive_probs = site_probs > 0
                log_probs[positive_probs] = np.log2(site_probs[positive_probs])
                site_entropies = -_column_dot(site_probs, log_probs)

            if all_valid_ascii:
                entropy_burden = np.full(
                    n_taxa,
                    _average_site_entropy(site_entropies, aln_len),
                    dtype=np.float64,
                )
            else:
                entropy_sums = np.dot(valid_mask.astype(np.float64), site_entropies)
                entropy_burden = np.divide(
                    entropy_sums,
                    valid_lengths,
                    out=np.zeros(n_taxa, dtype=np.float64),
                    where=valid_lengths > 0,
                )
        else:
            entropy_burden = np.zeros(n_taxa, dtype=np.float64)
        entropy_threshold = self._high_outlier_threshold(entropy_burden, self.entropy_z)

        rows = []
        outliers = []
        for (
            taxon,
            gap_rate,
            occupancy,
            composition_distance,
            branch_proxy,
            rcvt,
            entropy,
        ) in zip(
            taxa,
            gap_rates,
            occupancies,
            composition_distances,
            long_branch_proxy,
            rcvt_values,
            entropy_burden,
        ):
            reasons: list[dict[str, float]] = []

            if np.isfinite(gap_rate) and gap_rate > gap_threshold:
                reasons.append(
                    {
                        "feature": "gap_rate",
                        "value": float(gap_rate),
                        "threshold": float(gap_threshold),
                        "explanation": (
                            "High fraction of gap/ambiguous symbols compared to other taxa."
                        ),
                    }
                )
            if np.isfinite(occupancy) and occupancy < occupancy_threshold:
                reasons.append(
                    {
                        "feature": "occupancy",
                        "value": float(occupancy),
                        "threshold": float(occupancy_threshold),
                        "explanation": "Low fraction of valid symbols compared to other taxa.",
                    }
                )
            if (
                np.isfinite(composition_distance)
                and composition_distance > composition_threshold
            ):
                reasons.append(
                    {
                        "feature": "composition_distance",
                        "value": float(composition_distance),
                        "threshold": float(composition_threshold),
                        "explanation": (
                            "Unusual sequence composition profile relative to other taxa."
                        ),
                    }
                )
            branch_proxy_is_finite = np.isfinite(branch_proxy)
            if branch_proxy_is_finite and branch_proxy > distance_threshold:
                reasons.append(
                    {
                        "feature": "long_branch_proxy",
                        "value": float(branch_proxy),
                        "threshold": float(distance_threshold),
                        "explanation": "High mean pairwise sequence distance to other taxa.",
                    }
                )
            if np.isfinite(rcvt) and rcvt > rcvt_threshold:
                reasons.append(
                    {
                        "feature": "rcvt",
                        "value": float(rcvt),
                        "threshold": float(rcvt_threshold),
                        "explanation": "High relative composition variability for this taxon.",
                    }
                )
            if np.isfinite(entropy) and entropy > entropy_threshold:
                reasons.append(
                    {
                        "feature": "entropy_burden",
                        "value": float(entropy),
                        "threshold": float(entropy_threshold),
                        "explanation": (
                            "High average site entropy across this taxon's valid positions."
                        ),
                    }
                )

            row = {
                "taxon": taxon,
                "gap_rate": round(float(gap_rate), 4),
                "occupancy": round(float(occupancy), 4),
                "composition_distance": round(float(composition_distance), 4),
                "long_branch_proxy": (
                    None
                    if not branch_proxy_is_finite
                    else round(float(branch_proxy), 4)
                ),
                "rcvt": round(float(rcvt), 4),
                "entropy_burden": round(float(entropy), 4),
                "flagged": bool(reasons),
                "reasons": [
                    {
                        "feature": r["feature"],
                        "value": round(float(r["value"]), 4),
                        "threshold": round(float(r["threshold"]), 4),
                        "direction": (
                            "low" if r["feature"] == "occupancy" else "high"
                        ),
                        "explanation": r["explanation"],
                    }
                    for r in reasons
                ],
            }
            rows.append(row)
            if reasons:
                outliers.append(row)

        return dict(
            features=[
                "gap_rate",
                "occupancy",
                "composition_distance",
                "long_branch_proxy",
                "rcvt",
                "entropy_burden",
            ],
            thresholds=dict(
                gap_rate=round(float(gap_threshold), 4)
                if np.isfinite(gap_threshold)
                else None,
                occupancy=round(float(occupancy_threshold), 4)
                if np.isfinite(occupancy_threshold)
                else None,
                composition_distance=round(float(composition_threshold), 4)
                if np.isfinite(composition_threshold)
                else None,
                long_branch_proxy=round(float(distance_threshold), 4)
                if np.isfinite(distance_threshold)
                else None,
                rcvt=round(float(rcvt_threshold), 4)
                if np.isfinite(rcvt_threshold)
                else None,
                entropy_burden=round(float(entropy_threshold), 4)
                if np.isfinite(entropy_threshold)
                else None,
            ),
            rows=rows,
            outliers=outliers,
        )

    def run(self) -> None:
        alignment, _, is_protein = self.get_alignment_and_format()
        result = self.calculate_outliers(alignment, is_protein)

        if self.json_output:
            payload = dict(
                features_evaluated=result["features"],
                thresholds=result["thresholds"],
                outliers=result["outliers"],
                rows=result["rows"],
                taxa=result["rows"],
            )
            print_json(payload)
            return

        lines = [
            "features_evaluated\t"
            "gap_rate,occupancy,composition_distance,long_branch_proxy,rcvt,entropy_burden"
        ]
        lines.append(
            "thresholds\t"
            f"gap_rate>{result['thresholds']['gap_rate']};"
            f"occupancy<{result['thresholds']['occupancy']};"
            f"composition_distance>{result['thresholds']['composition_distance']};"
            f"long_branch_proxy>{result['thresholds']['long_branch_proxy']};"
            f"rcvt>{result['thresholds']['rcvt']};"
            f"entropy_burden>{result['thresholds']['entropy_burden']}"
        )

        if not result["outliers"]:
            lines.append("No outlier taxa detected.")
            print("\n".join(lines))
            return

        for row in result["outliers"]:
            reason_string = ";".join(
                f"{reason['feature']}={reason['value']}"
                f"{'<' if reason.get('direction') == 'low' else '>'}{reason['threshold']}"
                for reason in row["reasons"]
            )
            explanation_string = " | ".join(
                reason["explanation"] for reason in row["reasons"]
            )
            lines.append(f"{row['taxon']}\t{reason_string}\t{explanation_string}")

        print("\n".join(lines))
