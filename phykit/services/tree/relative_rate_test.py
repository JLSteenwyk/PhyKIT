import itertools
import math
import os

from .base import Tree
from ...errors import PhykitUserError


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


def get_alignment_and_format_helper(*args, **kwargs):
    from ...helpers.files import get_alignment_and_format

    return get_alignment_and_format(*args, **kwargs)


def Path(*args, **kwargs):
    from pathlib import Path as _Path

    return _Path(*args, **kwargs)


def _chi2_sf(x, df):
    """Chi-squared survival function (lazy scipy import)."""
    if df == 1:
        return math.erfc(math.sqrt(x / 2.0))

    from scipy.stats import chi2
    return float(chi2.sf(x, df))


_NO_SINGLETON_OUTGROUP = object()
_TAJIMA_VECTOR_MIN_LENGTH = 256
_TAJIMA_SKIP_BYTES = b"-?NXnx"
_TAJIMA_SKIP_SCAN_BYTES = 4096
_CORRECTION_VECTOR_MIN_LENGTH = 2048
_path_isabs = os.path.isabs


class _LazyNumpy:
    def __getattr__(self, name):
        import numpy as _np

        return getattr(_np, name)


np = _LazyNumpy()


class RelativeRateTest(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(
            tree_file_path=parsed["tree_file_path"],
            alignment_file_path=parsed.get("alignment_file_path"),
            verbose=parsed["verbose"],
        )
        self.alignment_list = parsed.get("alignment_list")
        self.json_output = parsed["json_output"]
        self.plot_output = parsed.get("plot_output")
        self.plot_config = parsed["plot_config"]

    def process_args(self, args) -> dict:
        from ...helpers.plot_config import PlotConfig

        return dict(
            tree_file_path=args.tree,
            alignment_file_path=getattr(args, "alignment", None),
            alignment_list=getattr(args, "alignment_list", None),
            verbose=getattr(args, "verbose", False),
            json_output=getattr(args, "json", False),
            plot_output=getattr(args, "plot_output", None),
            plot_config=PlotConfig.from_args(args),
        )

    # ------------------------------------------------------------------
    # Core Tajima's test
    # ------------------------------------------------------------------

    @staticmethod
    def _tajima_test(seq_x: str, seq_y: str, seq_out: str) -> dict:
        """Tajima's relative rate test (Equation 4).

        Counts unique substitutions on each ingroup lineage relative
        to the outgroup.  Returns m1, m2, chi-squared, and p-value.
        """
        n_sites = min(len(seq_x), len(seq_y), len(seq_out))
        if n_sites < _TAJIMA_VECTOR_MIN_LENGTH:
            return RelativeRateTest._tajima_test_scalar(seq_x, seq_y, seq_out)

        try:
            x_bytes = seq_x[:n_sites].encode("ascii")
            y_bytes = seq_y[:n_sites].encode("ascii")
            out_bytes = seq_out[:n_sites].encode("ascii")
            x = np.frombuffer(x_bytes, dtype=np.uint8)
            y = np.frombuffer(y_bytes, dtype=np.uint8)
            out = np.frombuffer(out_bytes, dtype=np.uint8)
        except UnicodeEncodeError:
            return RelativeRateTest._tajima_test_scalar(seq_x, seq_y, seq_out)

        has_skip_code = any(
            code in x_bytes[:_TAJIMA_SKIP_SCAN_BYTES]
            or code in y_bytes[:_TAJIMA_SKIP_SCAN_BYTES]
            or code in out_bytes[:_TAJIMA_SKIP_SCAN_BYTES]
            or code in x_bytes[-_TAJIMA_SKIP_SCAN_BYTES:]
            or code in y_bytes[-_TAJIMA_SKIP_SCAN_BYTES:]
            or code in out_bytes[-_TAJIMA_SKIP_SCAN_BYTES:]
            for code in _TAJIMA_SKIP_BYTES
        )
        if not has_skip_code:
            has_skip_code = any(
                code in x_bytes or code in y_bytes or code in out_bytes
                for code in _TAJIMA_SKIP_BYTES
            )

        x_upper = RelativeRateTest._ascii_upper_array(x)
        y_upper = RelativeRateTest._ascii_upper_array(y)
        out_upper = RelativeRateTest._ascii_upper_array(out)

        x_diff = x_upper != out_upper
        y_diff = y_upper != out_upper
        if has_skip_code:
            valid = RelativeRateTest._tajima_valid_mask(
                x_upper,
                y_upper,
                out_upper,
            )
            m1 = int(np.count_nonzero(valid & x_diff & ~y_diff))
            m2 = int(np.count_nonzero(valid & y_diff & ~x_diff))
        else:
            m1 = int(np.count_nonzero(x_diff & ~y_diff))
            m2 = int(np.count_nonzero(y_diff & ~x_diff))

        return RelativeRateTest._tajima_result(m1, m2)

    @staticmethod
    def _tajima_valid_mask(x_upper, y_upper, out_upper):
        return (
            (x_upper != 45) & (x_upper != 63)
            & (x_upper != 78) & (x_upper != 88)
            & (y_upper != 45) & (y_upper != 63)
            & (y_upper != 78) & (y_upper != 88)
            & (out_upper != 45) & (out_upper != 63)
            & (out_upper != 78) & (out_upper != 88)
        )

    @staticmethod
    def _ascii_upper_array(values):
        upper = values.copy()
        lowercase = (upper >= 97) & (upper <= 122)
        upper[lowercase] -= 32
        return upper

    @staticmethod
    def _tajima_test_scalar(seq_x: str, seq_y: str, seq_out: str) -> dict:
        m1 = 0  # unique changes on lineage x
        m2 = 0  # unique changes on lineage y
        skip_chars = {"-", "?", "N", "X", "n", "x"}

        for cx, cy, co in zip(seq_x, seq_y, seq_out):
            if cx in skip_chars or cy in skip_chars or co in skip_chars:
                continue
            x_diff = cx.upper() != co.upper()
            y_diff = cy.upper() != co.upper()
            if x_diff and not y_diff:
                m1 += 1
            elif y_diff and not x_diff:
                m2 += 1
            # Both differ or neither differs -> uninformative

        return RelativeRateTest._tajima_result(m1, m2)

    @staticmethod
    def _tajima_result(m1: int, m2: int) -> dict:
        if m1 + m2 == 0:
            return {"m1": m1, "m2": m2, "chi2": 0.0, "p_value": 1.0}

        chi2 = (m1 - m2) ** 2 / (m1 + m2)
        p_value = _chi2_sf(chi2, df=1)

        return {"m1": m1, "m2": m2, "chi2": chi2, "p_value": p_value}

    # ------------------------------------------------------------------
    # Outgroup identification
    # ------------------------------------------------------------------

    @staticmethod
    def _single_terminal_name_if_singleton(clade):
        stack = [clade]
        terminal_count = 0
        terminal_name = _NO_SINGLETON_OUTGROUP

        while stack:
            node = stack.pop()
            children = node.clades
            if children:
                stack.extend(children)
                continue

            terminal_count += 1
            if terminal_count == 1:
                terminal_name = node.name
            else:
                return _NO_SINGLETON_OUTGROUP

        return terminal_name if terminal_count == 1 else _NO_SINGLETON_OUTGROUP

    @staticmethod
    def _generic_single_terminal_name_if_singleton(clade):
        terminal_count = 0
        terminal_name = _NO_SINGLETON_OUTGROUP

        for terminal in clade.get_terminals():
            terminal_count += 1
            if terminal_count == 1:
                terminal_name = terminal.name
            else:
                return _NO_SINGLETON_OUTGROUP

        return terminal_name if terminal_count == 1 else _NO_SINGLETON_OUTGROUP

    @staticmethod
    def _identify_outgroup(tree) -> str:
        """Identify outgroup as the earliest-diverging taxon in a rooted tree.

        The outgroup is the single taxon on the smaller side of the root split.
        """
        root = tree.root
        children = root.clades
        if len(children) < 2:
            raise PhykitUserError(
                ["Tree does not appear to be rooted (root has < 2 children)."]
            )

        try:
            for child in children:
                terminal_name = (
                    RelativeRateTest._single_terminal_name_if_singleton(child)
                )
                if terminal_name is not _NO_SINGLETON_OUTGROUP:
                    return terminal_name
        except AttributeError:
            pass
        else:
            raise PhykitUserError(
                [
                    "Cannot determine a single outgroup taxon from the rooted tree. "
                    "The smallest clade at the root contains multiple taxa. "
                    "Please ensure the tree is rooted with a single outgroup taxon."
                ]
            )

        # Generic fallback for nonstandard tree-like objects.
        for child in children:
            terminal_name = (
                RelativeRateTest._generic_single_terminal_name_if_singleton(child)
            )
            if terminal_name is not _NO_SINGLETON_OUTGROUP:
                return terminal_name

        raise PhykitUserError(
            [
                "Cannot determine a single outgroup taxon from the rooted tree. "
                "The smallest clade at the root contains multiple taxa. "
                "Please ensure the tree is rooted with a single outgroup taxon."
            ]
        )

    # ------------------------------------------------------------------
    # Multiple testing correction
    # ------------------------------------------------------------------

    @staticmethod
    def _bonferroni(p_values: list[float]) -> list[float]:
        n = len(p_values)
        if n == 0:
            return []
        if n < _CORRECTION_VECTOR_MIN_LENGTH:
            return [min(p_value * n, 1.0) for p_value in p_values]
        return np.minimum(np.asarray(p_values, dtype=float) * n, 1.0).tolist()

    @staticmethod
    def _fdr(p_values: list[float]) -> list[float]:
        """Benjamini-Hochberg FDR correction."""
        n = len(p_values)
        if n == 0:
            return []
        if n < _CORRECTION_VECTOR_MIN_LENGTH:
            indexed = sorted(enumerate(p_values), key=lambda item: item[1])
            corrected = [0.0] * n
            previous = 1.0
            for rank_minus_1 in range(n - 1, -1, -1):
                original_idx, p_value = indexed[rank_minus_1]
                rank = rank_minus_1 + 1
                adjusted = min(p_value * n / rank, previous)
                adjusted = min(adjusted, 1.0)
                corrected[original_idx] = adjusted
                previous = adjusted
            return corrected
        p_arr = np.asarray(p_values, dtype=float)
        order = np.argsort(p_arr)
        sorted_p = p_arr[order]
        ranks = np.arange(1, n + 1, dtype=float)
        adjusted = np.minimum.accumulate((sorted_p * n / ranks)[::-1])[::-1]
        adjusted = np.minimum(adjusted, 1.0)
        corrected = np.empty(n, dtype=float)
        corrected[order] = adjusted
        return corrected.tolist()

    # ------------------------------------------------------------------
    # Single alignment analysis
    # ------------------------------------------------------------------

    def _run_single(self, alignment_path: str, tree, outgroup: str) -> list[dict]:
        """Run Tajima's test for all ingroup pairs in one alignment."""
        alignment, _, _ = get_alignment_and_format_helper(alignment_path)

        # Build sequence dict
        seq_dict = {}
        for record in alignment:
            seq_dict[record.id] = str(record.seq)

        if outgroup not in seq_dict:
            raise PhykitUserError(
                [
                    f"Outgroup taxon '{outgroup}' from tree not found in alignment. "
                    f"Alignment taxa: {sorted(seq_dict.keys())}"
                ]
            )

        seq_out = seq_dict[outgroup]
        ingroup = sorted(t for t in seq_dict if t != outgroup)

        if len(ingroup) < 2:
            raise PhykitUserError(
                ["At least 2 ingroup taxa are required (3 total including outgroup)."]
            )

        results = self._run_pairwise_tests(seq_dict, ingroup, seq_out)

        self._add_multiple_testing_corrections(results)
        return results

    def _run_pairwise_tests(
        self, seq_dict: dict[str, str], ingroup: list[str], seq_out: str
    ) -> list[dict]:
        lengths = {len(seq_dict[taxon]) for taxon in ingroup}
        lengths.add(len(seq_out))
        if len(lengths) == 1:
            try:
                return self._run_pairwise_tests_vectorized(seq_dict, ingroup, seq_out)
            except UnicodeEncodeError:
                pass
        return self._run_pairwise_tests_legacy(seq_dict, ingroup, seq_out)

    def _run_pairwise_tests_legacy(
        self, seq_dict: dict[str, str], ingroup: list[str], seq_out: str
    ) -> list[dict]:
        results = []
        for t1, t2 in itertools.combinations(ingroup, 2):
            test_result = self._tajima_test(seq_dict[t1], seq_dict[t2], seq_out)
            test_result["taxon1"] = t1
            test_result["taxon2"] = t2
            results.append(test_result)
        return results

    def _run_pairwise_tests_vectorized(
        self, seq_dict: dict[str, str], ingroup: list[str], seq_out: str
    ) -> list[dict]:
        out = np.frombuffer(seq_out.upper().encode("ascii"), dtype=np.uint8)
        ingroup_matrix = self._ingroup_ascii_matrix(seq_dict, ingroup, out.size)

        skip = np.frombuffer(b"-?NX", dtype=np.uint8)
        out_valid = ~np.isin(out, skip)
        ingroup_valid = ~np.isin(ingroup_matrix, skip)
        valid = ingroup_valid & out_valid
        differs = valid & (ingroup_matrix != out)
        matches_out = valid & ~differs

        differs_i = differs.astype(np.int32, copy=False)
        matches_i = matches_out.astype(np.int32, copy=False)
        m1_matrix = differs_i @ matches_i.T

        pair_count = len(ingroup) * (len(ingroup) - 1) // 2
        if pair_count <= 1024:
            return self._pairwise_results_small(m1_matrix, ingroup)
        return self._pairwise_results_large(m1_matrix, ingroup)

    @staticmethod
    def _ingroup_ascii_matrix(
        seq_dict: dict[str, str],
        ingroup: list[str],
        seq_len: int,
    ):
        return np.frombuffer(
            "".join(seq_dict[taxon].upper() for taxon in ingroup).encode("ascii"),
            dtype=np.uint8,
        ).reshape(len(ingroup), seq_len)

    @staticmethod
    def _pairwise_results_small(m1_matrix, ingroup: list[str]) -> list[dict]:
        results = []
        for i, taxon1 in enumerate(ingroup[:-1]):
            for j in range(i + 1, len(ingroup)):
                m1 = int(m1_matrix[i, j])
                m2 = int(m1_matrix[j, i])
                result = RelativeRateTest._tajima_result(m1, m2)
                result["taxon1"] = taxon1
                result["taxon2"] = ingroup[j]
                results.append(result)
        return results

    @staticmethod
    def _pairwise_results_large(m1_matrix, ingroup: list[str]) -> list[dict]:
        from scipy.special import erfc

        results = []
        append = results.append
        n_ingroup = len(ingroup)
        for i, taxon1 in enumerate(ingroup[:-1]):
            j_start = i + 1
            m1_values = m1_matrix[i, j_start:]
            m2_values = m1_matrix[j_start:, i]
            totals = m1_values + m2_values

            chi2_values = np.zeros_like(totals, dtype=np.float64)
            informative = totals > 0
            if informative.any():
                diff_values = m1_values - m2_values
                chi2_values[informative] = (
                    diff_values[informative] * diff_values[informative]
                ) / totals[informative]

            p_values = np.ones_like(chi2_values, dtype=np.float64)
            if informative.any():
                p_values[informative] = erfc(
                    np.sqrt(chi2_values[informative] * 0.5)
                )

            for taxon2, m1, m2, chi2, p_value in zip(
                ingroup[j_start:n_ingroup],
                m1_values,
                m2_values,
                chi2_values,
                p_values,
            ):
                append(
                    {
                        "m1": int(m1),
                        "m2": int(m2),
                        "chi2": float(chi2),
                        "p_value": float(p_value),
                        "taxon1": taxon1,
                        "taxon2": taxon2,
                    }
                )
        return results

    def _add_multiple_testing_corrections(self, results: list[dict]) -> None:
        p_values = [r["p_value"] for r in results]
        p_bonf = self._bonferroni(p_values)
        p_fdr = self._fdr(p_values)
        for r, bonferroni_value, fdr_value in zip(results, p_bonf, p_fdr):
            r["p_bonf"] = bonferroni_value
            r["p_fdr"] = fdr_value

    # ------------------------------------------------------------------
    # Main entry point
    # ------------------------------------------------------------------

    def run(self):
        tree = self.read_tree_file_unmodified()

        outgroup = self._identify_outgroup(tree)

        if self.alignment_list:
            # Batch mode
            list_path = Path(self.alignment_list)
            if not list_path.exists():
                raise PhykitUserError([f"{self.alignment_list} not found."])
            with list_path.open() as handle:
                aln_paths = [
                    stripped
                    for line in handle
                    if (stripped := line.strip())
                    and stripped[0] != "#"
                ]
            if not aln_paths:
                raise PhykitUserError(["No alignment paths found in list file."])

            # Run each alignment
            all_results = {}  # (t1, t2) -> list of per-gene results
            n_alignments = 0
            parent_str = str(list_path.parent)
            parent_prefix = "" if parent_str == "." else parent_str + os.sep
            for aln_path in aln_paths:
                resolved = (
                    aln_path if _path_isabs(aln_path) else parent_prefix + aln_path
                )
                try:
                    results = self._run_single(resolved, tree, outgroup)
                    n_alignments += 1
                    for r in results:
                        key = (r["taxon1"], r["taxon2"])
                        all_results.setdefault(key, []).append(r)
                except (SystemExit, ValueError, FileNotFoundError, KeyError):
                    continue  # skip problematic alignments

            self._output_batch(outgroup, n_alignments, all_results)

        elif self.alignment_file_path:
            # Single mode
            results = self._run_single(self.alignment_file_path, tree, outgroup)
            self._output_single(outgroup, results)
            if self.plot_output:
                self._plot_heatmap(results, self.plot_output)
        else:
            raise PhykitUserError(
                ["Provide -a (single alignment) or -l (alignment list)."]
            )

    # ------------------------------------------------------------------
    # Plotting
    # ------------------------------------------------------------------

    @staticmethod
    def _build_heatmap_matrices(results: list[dict]):
        taxa = sorted(
            set(r["taxon1"] for r in results) | set(r["taxon2"] for r in results)
        )
        n = len(taxa)
        idx = {t: i for i, t in enumerate(taxa)}

        matrix = np.full((n, n), np.nan)
        p_fdr_matrix = np.full((n, n), np.nan)
        for r in results:
            i, j = idx[r["taxon1"]], idx[r["taxon2"]]
            p_fdr = r["p_fdr"]
            val = -math.log10(p_fdr) if p_fdr > 0 else 10.0
            matrix[i, j] = val
            matrix[j, i] = val
            p_fdr_matrix[i, j] = p_fdr
            p_fdr_matrix[j, i] = p_fdr

        return taxa, matrix, p_fdr_matrix

    @staticmethod
    def _significant_heatmap_cells(p_fdr_matrix, alpha=0.05):
        if p_fdr_matrix.size == 0 or p_fdr_matrix.shape[0] <= 1:
            return (
                np.empty(0, dtype=np.intp),
                np.empty(0, dtype=np.intp),
            )
        if np.nanmin(p_fdr_matrix) >= alpha:
            return (
                np.empty(0, dtype=np.intp),
                np.empty(0, dtype=np.intp),
            )
        flat_indices = np.flatnonzero(p_fdr_matrix < alpha)
        if flat_indices.size == 0:
            return (
                np.empty(0, dtype=np.intp),
                np.empty(0, dtype=np.intp),
            )
        return np.unravel_index(flat_indices, p_fdr_matrix.shape)

    def _plot_heatmap(self, results: list[dict], output_path: str) -> None:
        """Plot a symmetric heatmap of -log10(p_fdr) for all pairwise tests."""
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            import numpy as np
        except ImportError:
            print("matplotlib and numpy are required for --plot-output. "
                  "Install them and retry.")
            raise SystemExit(2)

        config = self.plot_config
        config.resolve(n_rows=len(results), n_cols=None)

        taxa, matrix, p_fdr_matrix = self._build_heatmap_matrices(results)
        n = len(taxa)

        fig, ax = plt.subplots(figsize=(config.fig_width, config.fig_height))
        cmap = plt.cm.YlOrRd.copy()
        cmap.set_bad(color="white")

        # Fix scale so non-significant data always looks light:
        # vmin=0, vmax at least -log10(0.05) ≈ 1.3
        sig_threshold = -math.log10(0.05)
        data_max = float(np.nanmax(matrix)) if not np.isnan(matrix).all() else 0.0
        vmax = max(sig_threshold, data_max) * 1.1

        im = ax.imshow(
            matrix, cmap=cmap, aspect="equal", interpolation="nearest",
            vmin=0.0, vmax=vmax,
        )

        ax.set_xticks(range(n))
        ax.set_yticks(range(n))
        ax.set_xticklabels(taxa, rotation=45, ha="right", fontsize=9)
        ax.set_yticklabels(taxa, fontsize=9)

        sig_rows, sig_cols = self._significant_heatmap_cells(p_fdr_matrix)
        if sig_rows.size:
            ax.scatter(
                sig_cols,
                sig_rows,
                marker="$*$",
                s=170,
                c="black",
                linewidths=0,
                zorder=3,
            )

        cbar = fig.colorbar(im, ax=ax, shrink=0.8)
        cbar.set_label("$-\\log_{10}$(FDR-corrected p-value)", fontsize=10)

        cbar.ax.axhline(y=sig_threshold, color="black", linewidth=1.5, linestyle="--")

        if config.show_title:
            ax.set_title(config.title or "Relative Rate Test", fontsize=config.title_fontsize)
        if config.axis_fontsize:
            ax.xaxis.label.set_fontsize(config.axis_fontsize)
            ax.yaxis.label.set_fontsize(config.axis_fontsize)

        fig.tight_layout()
        fig.savefig(output_path, dpi=config.dpi, bbox_inches="tight")
        plt.close(fig)

    # ------------------------------------------------------------------
    # Output
    # ------------------------------------------------------------------

    def _output_single(self, outgroup: str, results: list[dict]):
        if self.json_output:
            taxa = set()
            output_results = []
            append_result = output_results.append
            add_taxon = taxa.add
            for r in results:
                taxon1 = r["taxon1"]
                taxon2 = r["taxon2"]
                add_taxon(taxon1)
                add_taxon(taxon2)
                append_result({
                    "taxon1": taxon1, "taxon2": taxon2,
                    "m1": r["m1"], "m2": r["m2"],
                    "chi2": round(r["chi2"], 4),
                    "p_value": r["p_value"],
                    "p_bonferroni": r["p_bonf"],
                    "p_fdr": r["p_fdr"],
                })
            print_json({
                "outgroup": outgroup,
                "n_ingroup_taxa": len(taxa),
                "n_tests": len(results),
                "results": output_results,
            })
        else:
            try:
                taxa = set()
                body_lines = []
                append_body = body_lines.append
                add_taxon = taxa.add
                for r in results:
                    taxon1 = r["taxon1"]
                    taxon2 = r["taxon2"]
                    add_taxon(taxon1)
                    add_taxon(taxon2)
                    sig = "*" if r["p_fdr"] < 0.05 else ""
                    append_body(
                        f"{taxon1}\t{taxon2}\t{r['m1']}\t{r['m2']}\t"
                        f"{r['chi2']:.4f}\t{r['p_value']:.4e}\t{r['p_bonf']:.4e}\t"
                        f"{r['p_fdr']:.4e}\t{sig}"
                    )
                lines = [
                    f"Outgroup: {outgroup}",
                    f"Number of ingroup taxa: {len(taxa)}",
                    f"Number of pairwise tests: {len(results)}",
                    "---",
                    "taxon1\ttaxon2\tm1\tm2\tchi2\tp_value\tp_bonf\tp_fdr\tsignificant",
                ]
                lines.extend(body_lines)
                print("\n".join(lines))
            except BrokenPipeError:
                pass

    def _output_batch(self, outgroup: str, n_alignments: int, all_results: dict):
        if self.json_output:
            pairs = []
            for (t1, t2), gene_results in sorted(all_results.items()):
                n_reject, n_total, pct_reject, median_chi2 = (
                    self._summarize_batch_gene_results(gene_results)
                )
                pairs.append(dict(
                    taxon1=t1, taxon2=t2,
                    n_reject=n_reject, n_total=n_total,
                    pct_reject=round(pct_reject, 1),
                    median_chi2=round(median_chi2, 4),
                ))
            print_json(dict(
                outgroup=outgroup,
                n_alignments=n_alignments,
                pairs=pairs,
            ))
        else:
            try:
                lines = [
                    f"Number of alignments: {n_alignments}",
                    f"Outgroup: {outgroup}",
                    "---",
                    "taxon1\ttaxon2\tn_reject\tn_total\tpct_reject\tmedian_chi2",
                ]
                append = lines.append
                for (t1, t2), gene_results in sorted(all_results.items()):
                    n_reject, n_total, pct, median_chi2 = (
                        self._summarize_batch_gene_results(gene_results)
                    )
                    append(
                        f"{t1}\t{t2}\t{n_reject}\t{n_total}\t"
                        f"{pct:.1f}\t{median_chi2:.4f}"
                    )
                print("\n".join(lines))
            except BrokenPipeError:
                pass

    @staticmethod
    def _summarize_batch_gene_results(gene_results: list[dict]):
        n_reject = 0
        chi2_values = []
        append_chi2 = chi2_values.append
        for result in gene_results:
            if result["p_value"] < 0.05:
                n_reject += 1
            append_chi2(result["chi2"])
        chi2_values.sort()
        n_total = len(gene_results)
        median_chi2 = chi2_values[n_total // 2] if chi2_values else 0.0
        pct_reject = 100 * n_reject / n_total
        return n_reject, n_total, pct_reject, median_chi2
