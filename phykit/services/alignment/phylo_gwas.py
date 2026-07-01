"""Phylogenetic GWAS following the Pease et al. (2016) approach.

Performs per-site association tests between alignment columns and a
phenotype, applies Benjamini-Hochberg FDR correction, optionally
classifies significant associations as monophyletic or polyphyletic
using a phylogenetic tree, and produces a Manhattan plot.
"""

from __future__ import annotations

import csv
import heapq
import math
import sys
from bisect import bisect_right
from collections import Counter
from operator import itemgetter

from ._fasta import read_fasta_first_token_upper
from .base import Alignment


class _LazyNumpy:
    def __getattr__(self, name):
        import numpy as _np

        return getattr(_np, name)


np = _LazyNumpy()
_P_VALUE_GETTER = itemgetter("p_value")
_PARTITION_ROW_PATTERN = None


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


def _log_comb(n: int, k: int) -> float:
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def _get_partition_row_pattern():
    global _PARTITION_ROW_PATTERN
    if _PARTITION_ROW_PATTERN is None:
        import re

        _PARTITION_ROW_PATTERN = re.compile(
            r"[^,]*,\s*(\S+)\s*=\s*(\d+)\s*-\s*(\d+)"
        )
    return _PARTITION_ROW_PATTERN


def _log_comb_cached(n: int, k: int, cache: dict) -> float:
    key = ("log_comb", n, k)
    cached = cache.get(key)
    if cached is not None:
        return cached

    value = _log_comb(n, k)
    cache[key] = value
    return value


def _fisher_exact_2x2_two_sided(
    table: np.ndarray, cache: Optional[dict] = None
) -> float:
    """Return SciPy-compatible two-sided Fisher exact p-value for a 2x2 table."""
    return _fisher_exact_2x2_counts(
        int(table[0, 0]),
        int(table[0, 1]),
        int(table[1, 0]),
        int(table[1, 1]),
        cache,
    )


def _fisher_exact_2x2_counts(
    a: int,
    b: int,
    c: int,
    d: int,
    cache: Optional[dict] = None,
) -> float:
    """Return SciPy-compatible two-sided Fisher exact p-value for 2x2 counts."""
    if cache is not None:
        key = (a, b, c, d)
        cached = cache.get(key)
        if cached is not None:
            return cached

    row0 = a + b
    row1 = c + d
    col0 = a + c
    total = row0 + row1
    if total == 0:
        if cache is not None:
            cache[key] = 1.0
        return 1.0

    lower = max(0, col0 - row1)
    upper = min(row0, col0)
    if cache is None:
        log_comb = _log_comb
    else:
        log_comb = lambda n, k: _log_comb_cached(n, k, cache)

    denominator = log_comb(total, col0)

    def log_probability(x: int) -> float:
        return log_comb(row0, x) + log_comb(row1, col0 - x) - denominator

    cutoff_log_p = log_probability(a) + math.log1p(1e-12)
    log_p = log_probability(lower)
    log_total = None
    x = lower

    while True:
        if log_p <= cutoff_log_p:
            if log_total is None:
                log_total = log_p
            elif log_total >= log_p:
                log_total += math.log1p(math.exp(log_p - log_total))
            else:
                log_total = log_p + math.log1p(math.exp(log_total - log_p))

        if x == upper:
            break

        log_p += (
            math.log(row0 - x)
            + math.log(col0 - x)
            - math.log(x + 1)
            - math.log(row1 - col0 + x + 1)
        )
        x += 1

    if log_total is None:
        if cache is not None:
            cache[key] = 0.0
        return 0.0

    p_value = math.exp(log_total)
    p_value = min(1.0, p_value)
    if cache is not None:
        cache[key] = p_value
    return p_value


class PhyloGwas(Alignment):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__()
        self.alignment_path = parsed["alignment_path"]
        self.phenotype_path = parsed["phenotype_path"]
        self.output_path = parsed["output_path"]
        self.tree_path = parsed["tree_path"]
        self.partition_path = parsed["partition_path"]
        self.alpha = parsed["alpha"]
        self.exclude_monophyletic = parsed["exclude_monophyletic"]
        self.csv_output = parsed["csv_output"]
        self.dot_size = parsed["dot_size"]
        self.json_output = parsed["json_output"]
        self.plot_config = parsed["plot_config"]

    def process_args(self, args) -> Dict[str, object]:
        from ...helpers.plot_config import PlotConfig

        return dict(
            alignment_path=args.alignment,
            phenotype_path=args.phenotype,
            output_path=args.output,
            tree_path=getattr(args, "tree", None),
            partition_path=getattr(args, "partition", None),
            alpha=getattr(args, "alpha", 0.05),
            exclude_monophyletic=getattr(args, "exclude_monophyletic", False),
            csv_output=getattr(args, "csv", None),
            dot_size=getattr(args, "dot_size", 1.0),
            json_output=getattr(args, "json", False),
            plot_config=PlotConfig.from_args(args),
        )

    # ------------------------------------------------------------------
    # Input parsing
    # ------------------------------------------------------------------

    @staticmethod
    def _read_fasta(path: str) -> Dict[str, str]:
        """Read FASTA alignment into {taxon: sequence} dict."""
        return read_fasta_first_token_upper(path)

    @staticmethod
    def _read_phenotype(path: str) -> Dict[str, str]:
        """Read two-column TSV phenotype file into {taxon: value} dict."""
        pheno = {}
        with open(path, "rb") as fh:
            for line in fh:
                line = line.strip()
                if not line or line[0] == 35:
                    continue
                parts = line.split(b"\t", 2)
                if len(parts) < 2:
                    continue
                pheno[parts[0].decode()] = parts[1].decode()
        return pheno

    @staticmethod
    def _detect_phenotype_type(values: List[str]) -> str:
        """Return 'continuous' if all values castable to float, else 'categorical'."""
        for v in values:
            try:
                float(v)
            except ValueError:
                return "categorical"
        return "continuous"

    @staticmethod
    def _parse_partition_file(path: str) -> List[Tuple[str, int, int]]:
        """Parse RAxML-style partition file.

        Expected format per line (0- or 1-based):
            DNA, gene_name = start-end
        Returns list of (name, start_0based, end_0based_inclusive).
        """
        partitions = []
        append = partitions.append
        match = _get_partition_row_pattern().match
        with open(path) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                # Match: TYPE, name = start-end
                m = match(line)
                if m:
                    name = m.group(1)
                    start = int(m.group(2)) - 1  # convert to 0-based
                    end = int(m.group(3)) - 1
                    append((name, start, end))
        return partitions

    # ------------------------------------------------------------------
    # Per-site association tests
    # ------------------------------------------------------------------

    @staticmethod
    def _is_ambiguous(char: str) -> bool:
        return char in {"-", "N", "X", "?", "*"}

    @staticmethod
    def _build_ascii_alignment_matrix(
        sequences: List[str], aln_len: int
    ) -> Optional[np.ndarray]:
        """Build a byte matrix for fast column extraction, or None for non-ASCII."""
        try:
            data = "".join(sequences).encode("ascii")
        except UnicodeEncodeError:
            return None
        return np.frombuffer(data, dtype=np.uint8).reshape(len(sequences), aln_len)

    @staticmethod
    def _ascii_ambiguity_lookup() -> np.ndarray:
        lookup = np.zeros(256, dtype=bool)
        for char in b"-NX?*":
            lookup[char] = True
        return lookup

    @staticmethod
    def _valid_ascii_columns(
        alignment_matrix: np.ndarray, ambiguous_lookup: np.ndarray
    ) -> np.ndarray:
        return ~ambiguous_lookup[alignment_matrix].any(axis=0)

    @staticmethod
    def _biallelic_valid_ascii_columns(
        alignment_matrix: np.ndarray, ambiguous_lookup: np.ndarray
    ) -> np.ndarray:
        lower = alignment_matrix.min(axis=0)
        upper = alignment_matrix.max(axis=0)
        lower_counts = np.count_nonzero(alignment_matrix == lower, axis=0)
        upper_counts = np.count_nonzero(alignment_matrix == upper, axis=0)
        return (
            ~ambiguous_lookup[alignment_matrix].any(axis=0)
            & (lower != upper)
            & ((lower_counts + upper_counts) == alignment_matrix.shape[0])
        )

    @staticmethod
    def _should_prefilter_biallelic_ascii_columns(
        alignment_matrix: np.ndarray, valid_ascii_columns: np.ndarray
    ) -> bool:
        valid_indices = np.flatnonzero(valid_ascii_columns)
        if valid_indices.size == 0:
            return False

        if valid_indices.size > 512:
            sample_positions = np.linspace(
                0,
                valid_indices.size - 1,
                num=512,
                dtype=np.intp,
            )
            valid_indices = valid_indices[sample_positions]

        for col_idx in valid_indices:
            if PhyloGwas._major_minor_bytes(alignment_matrix[:, col_idx]) is None:
                return True
        return False

    @staticmethod
    def _extract_column_alleles(
        col_idx: int,
        sequences: List[str],
        alignment_matrix: Optional[np.ndarray],
        ambiguous_lookup: Optional[np.ndarray],
    ) -> Optional[List[str]]:
        if alignment_matrix is not None and ambiguous_lookup is not None:
            column = alignment_matrix[:, col_idx]
            if ambiguous_lookup[column].any():
                return None
            return list(column.tobytes().decode("ascii"))

        alleles = []
        for sequence in sequences:
            char = sequence[col_idx]
            if PhyloGwas._is_ambiguous(char):
                return None
            alleles.append(char)
        return alleles

    @staticmethod
    def _major_minor_bytes(allele_array: np.ndarray):
        if allele_array.size == 0:
            return None

        lower = int(allele_array.min())
        upper = int(allele_array.max())
        if lower == upper:
            return None

        lower_count = int(np.count_nonzero(allele_array == lower))
        upper_count = int(np.count_nonzero(allele_array == upper))
        if lower_count + upper_count != len(allele_array):
            return None

        if lower_count < upper_count:
            return upper, lower, lower_count
        return lower, upper, upper_count

    @staticmethod
    def _major_minor_two_group_counts(allele_array: np.ndarray, group_codes):
        """Return major/minor byte alleles and minor counts for 2 group codes."""
        if allele_array.size == 0:
            return None

        lower = int(allele_array.min())
        upper = int(allele_array.max())
        if lower == upper:
            return None

        lower_mask = allele_array == lower
        lower_count = int(np.count_nonzero(lower_mask))
        upper_mask = allele_array == upper
        upper_count = int(np.count_nonzero(upper_mask))
        if lower_count + upper_count != len(allele_array):
            return None

        if lower_count < upper_count:
            allele_0_byte = upper
            allele_1_byte = lower
            minor_count = lower_count
            derived_mask = lower_mask
        else:
            allele_0_byte = lower
            allele_1_byte = upper
            minor_count = upper_count
            derived_mask = upper_mask

        derived_1 = int(np.count_nonzero(group_codes & derived_mask))
        return (
            allele_0_byte,
            allele_1_byte,
            (minor_count - derived_1, derived_1),
        )

    def _test_site_categorical(
        self,
        alleles: List[str],
        groups: List[str],
        unique_groups: List[str],
        group_idx: Optional[Dict[str, int]] = None,
        group_codes: Optional[np.ndarray] = None,
        group_row_counts: Optional[np.ndarray] = None,
        group_code_offsets: Optional[np.ndarray] = None,
        chdtrc_func=None,
        fisher_cache: Optional[dict] = None,
    ) -> Optional[Tuple[float, str, str, Dict[str, float]]]:
        """Fisher's exact or chi2 test for a categorical phenotype at one site.

        Returns (p_value, allele_0, allele_1, group_freqs) or None if site
        should be skipped.
        """
        n_groups = len(unique_groups)
        allele_array = None
        used_two_group_byte_counts = False
        if (
            n_groups == 2
            and isinstance(alleles, np.ndarray)
            and alleles.dtype == np.uint8
            and isinstance(group_codes, np.ndarray)
        ):
            byte_summary = self._major_minor_two_group_counts(
                alleles, group_codes
            )
            if byte_summary is None:
                return None
            allele_0_byte, allele_1_byte, derived_counts = byte_summary
            allele_0 = chr(int(allele_0_byte))
            allele_1 = chr(int(allele_1_byte))
            used_two_group_byte_counts = True
        elif isinstance(alleles, np.ndarray) and alleles.dtype == np.uint8:
            allele_array = alleles
            byte_summary = self._major_minor_bytes(allele_array)
            if byte_summary is None:
                return None
            allele_0_byte, allele_1_byte, _ = byte_summary
            allele_0 = chr(int(allele_0_byte))
            allele_1 = chr(int(allele_1_byte))
        else:
            if self._starts_with_non_ascii_allele(alleles):
                allele_pair = self._binary_alleles_from_unicode_site(alleles)
                if allele_pair is None:
                    return None
                allele_0, allele_1 = allele_pair
            else:
                try:
                    allele_bytes = "".join(alleles).encode("ascii")
                    allele_array = np.frombuffer(allele_bytes, dtype=np.uint8)
                    byte_summary = self._major_minor_bytes(allele_array)
                    if byte_summary is None:
                        return None
                    allele_0_byte, allele_1_byte, _ = byte_summary
                    allele_0 = chr(int(allele_0_byte))
                    allele_1 = chr(int(allele_1_byte))
                except UnicodeEncodeError:
                    allele_pair = self._binary_alleles_from_unicode_site(alleles)
                    if allele_pair is None:
                        return None
                    allele_0, allele_1 = allele_pair

        if group_codes is None:
            if n_groups == 2 and group_idx is None:
                group_codes = np.fromiter(
                    (0 if group == unique_groups[0] else 1 for group in groups),
                    dtype=np.intp,
                    count=len(groups),
                )
            else:
                if group_idx is None:
                    group_idx = {g: i for i, g in enumerate(unique_groups)}
                group_codes = np.fromiter(
                    (group_idx[group] for group in groups),
                    dtype=np.intp,
                    count=len(groups),
                )
        if group_row_counts is None:
            group_row_counts = np.bincount(group_codes, minlength=n_groups)

        if (
            not used_two_group_byte_counts
            and allele_array is not None
            and group_code_offsets is not None
        ):
            group_allele_counts = np.bincount(
                group_code_offsets + allele_array,
                minlength=n_groups * 256,
            )
            derived_counts = group_allele_counts[int(allele_1_byte)::256]
        elif (
            not used_two_group_byte_counts
            and allele_array is not None
        ):
            derived = allele_array != ord(allele_0)
            derived_counts = np.bincount(group_codes[derived], minlength=n_groups)
        elif not used_two_group_byte_counts:
            derived = np.fromiter(
                (allele != allele_0 for allele in alleles),
                dtype=bool,
                count=len(alleles),
            )
            derived_counts = np.bincount(group_codes[derived], minlength=n_groups)
        if n_groups == 2:
            derived_0 = int(derived_counts[0])
            derived_1 = int(derived_counts[1])
            row_count_0 = int(group_row_counts[0])
            row_count_1 = int(group_row_counts[1])
            ancestral_0 = row_count_0 - derived_0
            ancestral_1 = row_count_1 - derived_1
            p_value = _fisher_exact_2x2_counts(
                ancestral_0,
                derived_0,
                ancestral_1,
                derived_1,
                fisher_cache,
            )
            group_freqs = {
                unique_groups[0]: (float(derived_0) / row_count_0 if row_count_0 else 0.0),
                unique_groups[1]: (float(derived_1) / row_count_1 if row_count_1 else 0.0),
            }
        else:
            table = np.empty((n_groups, 2), dtype=int)
            table[:, 1] = derived_counts
            table[:, 0] = group_row_counts - derived_counts
            if chdtrc_func is None:
                from scipy.special import chdtrc as chdtrc_func

            # Check for zero-sum rows/columns
            col_sums = table.sum(axis=0)
            if (group_row_counts == 0).any() or (col_sums == 0).any():
                return None
            degrees_of_freedom = n_groups - 1
            if degrees_of_freedom == 0:
                p_value = 1.0
            else:
                total = group_row_counts.sum()
                expected = group_row_counts[:, None] * col_sums[None, :] / total
                chi2_stat = float(((table - expected) ** 2 / expected).sum())
                p_value = float(chdtrc_func(degrees_of_freedom, chi2_stat))

            group_freqs = {
                group: (
                    float(derived_counts[i]) / group_row_counts[i]
                    if group_row_counts[i]
                    else 0.0
                )
                for i, group in enumerate(unique_groups)
            }

        return p_value, allele_0, allele_1, group_freqs

    @staticmethod
    def _starts_with_non_ascii_allele(alleles):
        try:
            first_allele = alleles[0]
        except (IndexError, TypeError):
            return False

        try:
            first_allele.encode("ascii")
        except UnicodeEncodeError:
            return True
        return False

    @staticmethod
    def _binary_alleles_from_unicode_site(alleles):
        allele_counts = {}
        for allele in alleles:
            if allele in allele_counts:
                allele_counts[allele] += 1
                continue
            if len(allele_counts) == 2:
                return None
            allele_counts[allele] = 1

        if len(allele_counts) != 2:
            return None

        allele_0, allele_1 = sorted(allele_counts)
        if allele_counts[allele_0] < allele_counts[allele_1]:
            allele_0, allele_1 = allele_1, allele_0
        return allele_0, allele_1

    def _test_site_continuous(
        self,
        alleles: List[str],
        phenotype_values: List[float],
        phenotype_centered: Optional[np.ndarray] = None,
        phenotype_ss: Optional[float] = None,
        stdtr_func=None,
    ) -> Optional[Tuple[float, str, str, float]]:
        """Point-biserial correlation for a continuous phenotype at one site.

        Returns (p_value, allele_0, allele_1, correlation_r) or None.
        """
        allele_array = None
        allele_0_byte = None
        derived_count = None
        if isinstance(alleles, np.ndarray) and alleles.dtype == np.uint8:
            allele_array = alleles
        else:
            if self._starts_with_non_ascii_allele(alleles):
                allele_pair = self._binary_alleles_from_unicode_site(alleles)
                if allele_pair is None:
                    return None
                allele_0, allele_1 = allele_pair
            else:
                try:
                    allele_bytes = "".join(alleles).encode("ascii")
                    allele_array = np.frombuffer(allele_bytes, dtype=np.uint8)
                except UnicodeEncodeError:
                    allele_pair = self._binary_alleles_from_unicode_site(alleles)
                    if allele_pair is None:
                        return None
                    allele_0, allele_1 = allele_pair

        if allele_array is not None:
            byte_summary = self._major_minor_bytes(allele_array)
            if byte_summary is None:
                return None
            allele_0_byte, allele_1_byte, derived_count = byte_summary
            allele_0 = chr(int(allele_0_byte))
            allele_1 = chr(int(allele_1_byte))

        if phenotype_centered is None or phenotype_ss is None:
            values = np.array(phenotype_values, dtype=float)
            phenotype_centered = values - values.mean()
            phenotype_ss = float(np.dot(phenotype_centered, phenotype_centered))

        # Need variance in both arrays
        if allele_array is not None:
            n_observations = len(allele_array)
            binary_ss = (
                derived_count * (n_observations - derived_count) / n_observations
            )
            phenotype_dot = float(
                np.dot(phenotype_centered, allele_array != allele_0_byte)
            )
        else:
            binary = np.array([0.0 if a == allele_0 else 1.0 for a in alleles])
            n_observations = len(binary)
            binary_centered = binary - binary.mean()
            binary_ss = float(np.dot(binary_centered, binary_centered))
            phenotype_dot = float(np.dot(binary_centered, phenotype_centered))
        if binary_ss == 0.0 or phenotype_ss == 0.0:
            return None

        r = phenotype_dot / math.sqrt(binary_ss * phenotype_ss)
        if math.isnan(r):
            return np.nan, allele_0, allele_1, np.nan
        r = max(-1.0, min(1.0, r))
        df = n_observations - 2
        if df <= 0:
            p_value = 1.0
        elif abs(r) == 1.0:
            p_value = 0.0
        else:
            if stdtr_func is None:
                from scipy.special import stdtr as stdtr_func

            t_stat = abs(r) * math.sqrt(df / ((1.0 - r) * (1.0 + r)))
            p_value = float(2.0 * stdtr_func(df, -t_stat))
        return p_value, allele_0, allele_1, float(r)

    # ------------------------------------------------------------------
    # Multiple testing correction
    # ------------------------------------------------------------------

    @staticmethod
    def _benjamini_hochberg(p_values: List[float]) -> np.ndarray:
        """Return BH-adjusted p-values."""
        n = len(p_values)
        if n == 0:
            return np.array([])
        p_arr = np.array(p_values, dtype=float)
        sorted_indices = np.argsort(p_arr)
        adjusted = p_arr[sorted_indices]
        adjusted *= n
        adjusted /= np.arange(1, n + 1, dtype=float)
        adjusted_reversed = adjusted[::-1]
        np.minimum.accumulate(adjusted_reversed, out=adjusted_reversed)
        np.minimum(adjusted, 1.0, out=adjusted)

        result = np.empty(n, dtype=float)
        result[sorted_indices] = adjusted
        return result

    # ------------------------------------------------------------------
    # Phylogenetic pattern classification
    # ------------------------------------------------------------------

    @staticmethod
    def _minor_allele_taxa(
        shared_taxa: list[str],
        seqs: dict[str, str],
        position: int,
        allele: str,
    ) -> list[str]:
        return [taxon for taxon in shared_taxa if seqs[taxon][position] == allele]

    @staticmethod
    def _build_phylo_pattern_index(tree) -> set:
        """Return descendant taxon sets for all clades in one postorder pass."""
        direct_result = PhyloGwas._build_phylo_pattern_index_direct(tree)
        if direct_result is not None:
            return direct_result

        clade_tips = {}
        monophyletic_sets = set()
        for clade in tree.find_clades(order="postorder"):
            if clade.is_terminal():
                tips = frozenset([clade.name])
            else:
                tips = frozenset().union(
                    *(clade_tips[id(child)] for child in clade.clades)
                )
            clade_tips[id(clade)] = tips
            monophyletic_sets.add(tips)
        return monophyletic_sets

    @staticmethod
    def _build_phylo_pattern_index_direct(tree):
        try:
            root = tree.root
            root.clades
        except AttributeError:
            return None

        preorder = []
        stack = [root]
        append = preorder.append
        pop = stack.pop
        extend = stack.extend
        try:
            while stack:
                clade = pop()
                append(clade)
                children = clade.clades
                if children:
                    extend(children)
        except AttributeError:
            return None

        clade_tips = {}
        monophyletic_sets = set()
        try:
            for clade in reversed(preorder):
                children = clade.clades
                if children:
                    child_count = len(children)
                    if child_count == 1:
                        tips = clade_tips[id(children[0])]
                    elif child_count == 2:
                        tips = (
                            clade_tips[id(children[0])]
                            | clade_tips[id(children[1])]
                        )
                    else:
                        tips = frozenset().union(
                            *(clade_tips[id(child)] for child in children)
                        )
                else:
                    tips = frozenset([clade.name])
                clade_tips[id(clade)] = tips
                monophyletic_sets.add(tips)
        except (AttributeError, KeyError, TypeError):
            return None

        return monophyletic_sets

    @staticmethod
    def _classify_phylo_pattern(
        tree,
        taxa_with_derived_allele: List[str],
        monophyletic_sets: Optional[set] = None,
    ) -> str:
        """Check if taxa with the derived allele form a monophyletic group."""
        if len(taxa_with_derived_allele) <= 1:
            return "monophyletic"

        derived_set = set(taxa_with_derived_allele)
        if monophyletic_sets is not None:
            return (
                "monophyletic"
                if frozenset(derived_set) in monophyletic_sets
                else "polyphyletic"
            )

        try:
            mrca = tree.common_ancestor(taxa_with_derived_allele)
        except Exception:
            return "polyphyletic"

        mrca_tips = {t.name for t in mrca.get_terminals()}

        if derived_set == mrca_tips:
            return "monophyletic"
        else:
            return "polyphyletic"

    @staticmethod
    def _prune_tree_to_taxa(tree, taxa: set) -> None:
        """Prune tips whose names are not in ``taxa`` using one terminal pass."""
        for tip in PhyloGwas._terminal_clades(tree):
            if tip.name not in taxa:
                tree.prune(tip)

    @staticmethod
    def _terminal_clades(tree):
        try:
            root = tree.root
            root.clades
        except AttributeError:
            return list(tree.get_terminals())

        terminals = []
        stack = [root]
        pop = stack.pop
        append = stack.append
        append_terminal = terminals.append
        try:
            while stack:
                clade = pop()
                children = clade.clades
                if children:
                    child_count = len(children)
                    if child_count == 2:
                        append(children[1])
                        append(children[0])
                    else:
                        for index in range(child_count - 1, -1, -1):
                            append(children[index])
                else:
                    append_terminal(clade)
        except AttributeError:
            return list(tree.get_terminals())
        return terminals

    # ------------------------------------------------------------------
    # Partition gene lookup
    # ------------------------------------------------------------------

    @staticmethod
    def _position_to_gene(
        position: int, partitions: List[Tuple[str, int, int]]
    ) -> Optional[str]:
        """Return gene name for a given 0-based position, or None."""
        for name, start, end in partitions:
            if start <= position <= end:
                return name
        return None

    @staticmethod
    def _build_partition_lookup(
        partitions: List[Tuple[str, int, int]]
    ) -> Optional[Tuple[List[int], List[int], List[str]]]:
        """Build a binary-searchable lookup for sorted, non-overlapping ranges."""
        starts = []
        ends = []
        names = []
        previous_end = -1
        for name, start, end in partitions:
            if start <= previous_end or end < start:
                return None
            starts.append(start)
            ends.append(end)
            names.append(name)
            previous_end = end
        return starts, ends, names

    @staticmethod
    def _position_to_gene_from_lookup(
        position: int, lookup: Tuple[List[int], List[int], List[str]]
    ) -> Optional[str]:
        """Return gene name using a lookup from _build_partition_lookup."""
        starts, ends, names = lookup
        index = bisect_right(starts, position) - 1
        if index >= 0 and position <= ends[index]:
            return names[index]
        return None

    @staticmethod
    def _sorted_shared_taxa(seqs: dict, phenotypes: dict) -> list[str]:
        return sorted(set(seqs).intersection(phenotypes))

    # ------------------------------------------------------------------
    # Main run
    # ------------------------------------------------------------------

    def run(self) -> None:
        # 1. Parse inputs
        seqs = self._read_fasta(self.alignment_path)
        phenotypes = self._read_phenotype(self.phenotype_path)

        shared_taxa = self._sorted_shared_taxa(seqs, phenotypes)
        if len(shared_taxa) < 4:
            print(
                f"Error: need at least 4 shared taxa between alignment and phenotype "
                f"file, found {len(shared_taxa)}."
            )
            sys.exit(2)

        # Validate alignment lengths
        aln_len = len(next(iter(seqs.values())))
        for taxon, seq in seqs.items():
            if len(seq) != aln_len:
                print(f"Error: sequences have unequal lengths.")
                sys.exit(2)
        shared_seqs = [seqs[t] for t in shared_taxa]
        alignment_matrix = self._build_ascii_alignment_matrix(shared_seqs, aln_len)
        ambiguous_lookup = (
            self._ascii_ambiguity_lookup() if alignment_matrix is not None else None
        )
        valid_ascii_columns = (
            self._valid_ascii_columns(alignment_matrix, ambiguous_lookup)
            if alignment_matrix is not None and ambiguous_lookup is not None
            else None
        )
        if (
            valid_ascii_columns is not None
            and self._should_prefilter_biallelic_ascii_columns(
                alignment_matrix, valid_ascii_columns
            )
        ):
            valid_ascii_columns = self._biallelic_valid_ascii_columns(
                alignment_matrix, ambiguous_lookup
            )

        # Detect phenotype type
        pheno_values = [phenotypes[t] for t in shared_taxa]
        pheno_type = self._detect_phenotype_type(pheno_values)
        is_continuous = pheno_type == "continuous"

        if is_continuous:
            from scipy.special import stdtr

            pheno_float = [float(phenotypes[t]) for t in shared_taxa]
            pheno_array = np.array(pheno_float, dtype=float)
            pheno_centered = pheno_array - pheno_array.mean()
            pheno_ss = float(np.dot(pheno_centered, pheno_centered))
        else:
            from scipy.special import chdtrc

            unique_groups = sorted(set(pheno_values))
            group_counts = Counter(pheno_values)
            groups = [phenotypes[t] for t in shared_taxa]
            group_idx = {g: i for i, g in enumerate(unique_groups)}
            group_codes = np.fromiter(
                (group_idx[group] for group in groups),
                dtype=np.intp,
                count=len(groups),
            )
            group_row_counts = np.bincount(
                group_codes, minlength=len(unique_groups)
            )
            group_code_offsets = group_codes * 256
            fisher_cache = {} if len(unique_groups) == 2 else None

        # Parse tree if provided
        tree = None
        if self.tree_path:
            from Bio import Phylo

            tree = Phylo.read(self.tree_path, "newick")
            # Prune tree to shared taxa
            self._prune_tree_to_taxa(tree, set(shared_taxa))
        monophyletic_sets = (
            self._build_phylo_pattern_index(tree) if tree is not None else None
        )

        # Parse partitions if provided
        partitions = []
        if self.partition_path:
            partitions = self._parse_partition_file(self.partition_path)
        partition_lookup = (
            self._build_partition_lookup(partitions) if partitions else None
        )

        # 2. Per-site association tests
        site_results = []
        raw_p_values = []

        column_indices = (
            np.flatnonzero(valid_ascii_columns)
            if valid_ascii_columns is not None
            else range(aln_len)
        )
        for col_idx in column_indices:
            # Extract alleles for shared taxa
            if alignment_matrix is not None and ambiguous_lookup is not None:
                alleles = alignment_matrix[:, col_idx]
            else:
                alleles = self._extract_column_alleles(
                    col_idx, shared_seqs, alignment_matrix, ambiguous_lookup
                )
            if alleles is None:
                continue

            if is_continuous:
                result = self._test_site_continuous(
                    alleles, pheno_float, pheno_centered, pheno_ss, stdtr
                )
                if result is None:
                    continue
                p_value, allele_0, allele_1, corr_r = result
                site_results.append(
                    {
                        "position": col_idx + 1,  # 1-based
                        "allele_0": allele_0,
                        "allele_1": allele_1,
                        "p_value": p_value,
                        "correlation_r": corr_r,
                    }
                )
            else:
                result = self._test_site_categorical(
                    alleles,
                    groups,
                    unique_groups,
                    group_idx,
                    group_codes,
                    group_row_counts,
                    group_code_offsets,
                    chdtrc,
                    fisher_cache,
                )
                if result is None:
                    continue
                p_value, allele_0, allele_1, group_freqs = result
                site_results.append(
                    {
                        "position": col_idx + 1,  # 1-based
                        "allele_0": allele_0,
                        "allele_1": allele_1,
                        "p_value": p_value,
                        "group_freqs": group_freqs,
                    }
                )

            raw_p_values.append(p_value)

        # 3. Multiple testing correction
        if raw_p_values:
            adjusted = self._benjamini_hochberg(raw_p_values)
            for r, fdr_p_value in zip(site_results, adjusted):
                fdr_p_value = float(fdr_p_value)
                r["fdr_p_value"] = fdr_p_value
                r["fdr_significant"] = fdr_p_value < self.alpha
        else:
            for r in site_results:
                r["fdr_p_value"] = 1.0
                r["fdr_significant"] = False

        # 4. Gene annotation and phylo pattern
        for r in site_results:
            position = r["position"] - 1
            if partition_lookup is not None:
                r["gene"] = self._position_to_gene_from_lookup(
                    position, partition_lookup
                )
            else:
                r["gene"] = self._position_to_gene(position, partitions)

            if tree is not None and r["fdr_significant"]:
                taxa_minor = self._minor_allele_taxa(
                    shared_taxa,
                    seqs,
                    position,
                    r["allele_1"],
                )
                r["phylo_pattern"] = self._classify_phylo_pattern(
                    tree, taxa_minor, monophyletic_sets
                )
            elif tree is not None:
                r["phylo_pattern"] = None
            else:
                r["phylo_pattern"] = None

        # 5. Exclude monophyletic if requested
        if self.exclude_monophyletic:
            site_results = [
                r
                for r in site_results
                if not (r["fdr_significant"] and r.get("phylo_pattern") == "monophyletic")
            ]

        # Compute summary stats
        significant, n_polyphyletic, n_monophyletic = (
            self._summarize_significant_results(site_results)
        )

        # 6. Manhattan plot
        self._create_manhattan_plot(
            site_results, partitions, tree is not None
        )

        # 7. CSV output
        if self.csv_output:
            self._write_csv(site_results, is_continuous)

        # 8. Text / JSON output
        if self.json_output:
            payload = self._build_json_payload(
                len(shared_taxa),
                aln_len,
                pheno_type,
                site_results,
                len(significant),
                group_counts if not is_continuous else None,
                n_polyphyletic,
                n_monophyletic,
                tree is not None,
            )
            print_json(payload)
            return

        # Text output
        self._print_text_output(
            len(shared_taxa),
            aln_len,
            pheno_type,
            is_continuous,
            unique_groups,
            group_counts,
            site_results,
            significant,
            n_polyphyletic,
            n_monophyletic,
            tree is not None,
        )

    def _build_json_payload(
        self,
        n_taxa: int,
        aln_len: int,
        pheno_type: str,
        site_results: List[dict],
        n_significant: int,
        group_counts,
        n_polyphyletic: int,
        n_monophyletic: int,
        has_tree: bool,
    ) -> dict:
        payload = dict(
            n_taxa=n_taxa,
            alignment_length=aln_len,
            phenotype_type=pheno_type,
            sites_tested=len(site_results),
            significant_sites=n_significant,
            alpha=self.alpha,
            results=site_results,
        )
        if group_counts is not None:
            payload["groups"] = dict(group_counts)
        if has_tree:
            payload["polyphyletic"] = n_polyphyletic
            payload["monophyletic"] = n_monophyletic
        return payload

    @staticmethod
    def _summarize_significant_results(
        site_results: List[dict],
    ) -> Tuple[List[dict], int, int]:
        significant = []
        append = significant.append
        n_polyphyletic = 0
        n_monophyletic = 0
        for result in site_results:
            if result["fdr_significant"]:
                append(result)
                pattern = result.get("phylo_pattern")
                if pattern == "polyphyletic":
                    n_polyphyletic += 1
                elif pattern == "monophyletic":
                    n_monophyletic += 1
        return significant, n_polyphyletic, n_monophyletic

    @staticmethod
    def _top_significant_sites(significant: List[dict]) -> List[dict]:
        if len(significant) <= 10:
            return sorted(significant, key=_P_VALUE_GETTER)
        return heapq.nsmallest(10, significant, key=_P_VALUE_GETTER)

    @staticmethod
    def _prepare_manhattan_series(
        results: List[dict],
        color_nonsig: str,
        color_polyphyletic: str,
        color_monophyletic: str,
    ):
        n_results = len(results)
        positions = np.empty(n_results, dtype=np.int64)
        p_values = np.empty(n_results, dtype=np.float64)
        point_colors = [color_nonsig] * n_results
        max_sig_p = None

        for idx, result in enumerate(results):
            p_value = result["p_value"]
            positions[idx] = result["position"]
            p_values[idx] = p_value if p_value > 1e-300 else 1e-300

            if result["fdr_significant"]:
                if result.get("phylo_pattern") == "monophyletic":
                    point_colors[idx] = color_monophyletic
                else:
                    point_colors[idx] = color_polyphyletic
                if max_sig_p is None or p_value > max_sig_p:
                    max_sig_p = p_value

        neg_log_p = -np.log10(p_values)
        threshold = -np.log10(max_sig_p) if max_sig_p is not None else None
        return positions, neg_log_p, point_colors, threshold

    def _print_text_output(
        self,
        n_taxa: int,
        aln_len: int,
        pheno_type: str,
        is_continuous: bool,
        unique_groups: List[str],
        group_counts,
        site_results: List[dict],
        significant: List[dict],
        n_polyphyletic: int,
        n_monophyletic: int,
        has_tree: bool,
    ) -> None:
        lines = []
        lines.append("Phylogenetic GWAS")
        lines.append(f"Alignment: {self.alignment_path}")
        lines.append(f"Taxa: {n_taxa}")
        lines.append(f"Alignment length: {aln_len}")
        lines.append(f"Phenotype type: {pheno_type}")
        if not is_continuous:
            group_str = ", ".join(
                f"{g} ({group_counts[g]})" for g in unique_groups
            )
            lines.append(f"Groups: {group_str}")
        lines.append(f"Biallelic sites tested: {len(site_results)}")
        lines.append(f"Significant sites (FDR < {self.alpha}): {len(significant)}")
        if has_tree:
            lines.append(f"  Polyphyletic: {n_polyphyletic}")
            lines.append(f"  Monophyletic: {n_monophyletic}")

        if significant:
            lines.append("")
            top_sig = self._top_significant_sites(significant)
            if is_continuous:
                lines.append("Top significant sites:")
                lines.append(
                    f"  {'Position':<10}{'Gene':<12}{'Allele':<10}"
                    f"{'r':<10}{'p-value':<12}{'FDR_p':<12}"
                    f"{'Pattern':<14}"
                )
                for r in top_sig:
                    gene = r["gene"] or "."
                    pattern = r["phylo_pattern"] or "."
                    lines.append(
                        f"  {r['position']:<10}{gene:<12}"
                        f"{r['allele_0']}>{r['allele_1']:<8}"
                        f"{r.get('correlation_r', 0):<10.4f}"
                        f"{r['p_value']:<12.4g}"
                        f"{r['fdr_p_value']:<12.4g}"
                        f"{pattern:<14}"
                    )
            else:
                lines.append("Top significant sites:")
                freq_header = "/".join(unique_groups)
                lines.append(
                    f"  {'Position':<10}{'Gene':<12}{'Allele':<10}"
                    f"{freq_header:<20}{'p-value':<12}{'FDR_p':<12}"
                    f"{'Pattern':<14}"
                )
                for r in top_sig:
                    gene = r["gene"] or "."
                    pattern = r["phylo_pattern"] or "."
                    freqs = "/".join(
                        f"{r['group_freqs'].get(g, 0):.2f}" for g in unique_groups
                    )
                    lines.append(
                        f"  {r['position']:<10}{gene:<12}"
                        f"{r['allele_0']}>{r['allele_1']:<8}"
                        f"{freqs:<20}"
                        f"{r['p_value']:<12.4g}"
                        f"{r['fdr_p_value']:<12.4g}"
                        f"{pattern:<14}"
                    )

        try:
            print("\n".join(lines))
        except BrokenPipeError:
            pass

    # ------------------------------------------------------------------
    # Manhattan plot
    # ------------------------------------------------------------------

    def _create_manhattan_plot(
        self,
        results: List[dict],
        partitions: List[Tuple[str, int, int]],
        has_tree: bool,
    ) -> None:
        try:
            import matplotlib

            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            from matplotlib.collections import PatchCollection
            from matplotlib.patches import Patch, Rectangle
        except ImportError:
            print(
                "matplotlib is required for phylo_gwas. Install matplotlib and retry."
            )
            sys.exit(2)

        if not results:
            # Create an empty plot
            fig, ax = plt.subplots(figsize=(10, 4))
            ax.set_xlabel("Alignment position")
            ax.set_ylabel("-log10(p-value)")
            ax.text(
                0.5, 0.5, "No biallelic sites tested",
                transform=ax.transAxes, ha="center", va="center",
            )
            fig.savefig(self.output_path, dpi=self.plot_config.dpi, bbox_inches="tight")
            plt.close(fig)
            return

        config = self.plot_config
        config.resolve(n_rows=None, n_cols=None)

        default_colors = ["#377eb8", "#e41a1c", "#999999"]
        colors = config.merge_colors(default_colors)
        color_nonsig = colors[0]
        color_polyphyletic = colors[1]
        color_monophyletic = colors[2]

        positions, neg_log_p, point_colors, threshold = (
            self._prepare_manhattan_series(
                results,
                color_nonsig,
                color_polyphyletic,
                color_monophyletic,
            )
        )

        fig, ax = plt.subplots(figsize=(config.fig_width, config.fig_height))
        ax.scatter(
            positions, neg_log_p, c=point_colors,
            s=8 * self.dot_size, alpha=0.6, edgecolors="none",
        )

        # Significance threshold line
        if threshold is not None:
            ax.axhline(
                y=threshold,
                color="red",
                linestyle="--",
                lw=0.8,
                label=f"FDR {self.alpha}",
            )

        # Partition annotations
        if partitions:
            span_patches = [
                Rectangle(
                    (start + 1, 0.0),
                    (end + 1) - (start + 1),
                    1.0,
                )
                for i, (_, start, end) in enumerate(partitions)
                if i % 2 == 0
            ]
            if span_patches:
                ax.add_collection(
                    PatchCollection(
                        span_patches,
                        facecolor="gray",
                        edgecolor="none",
                        alpha=0.05,
                        transform=ax.get_xaxis_transform(),
                        zorder=0,
                    )
                )

            for i, (name, start, end) in enumerate(partitions):
                mid = (start + end) / 2 + 1
                if len(partitions) <= 50:
                    ax.text(
                        mid,
                        ax.get_ylim()[1] * 0.98,
                        name,
                        ha="center",
                        fontsize=4,
                        rotation=90,
                    )

        ax.set_xlabel("Alignment position")
        ax.set_ylabel(r"-log$_{10}$(p-value)")

        # Legend
        legend_elements = [
            Patch(facecolor=color_nonsig, label="Not significant"),
        ]
        if has_tree:
            legend_elements.append(
                Patch(facecolor=color_polyphyletic, label="Significant (polyphyletic)")
            )
            legend_elements.append(
                Patch(facecolor=color_monophyletic, label="Significant (monophyletic)")
            )
        else:
            legend_elements.append(
                Patch(facecolor=color_polyphyletic, label="Significant")
            )
        legend_loc = config.legend_position or "upper right"
        if legend_loc != "none":
            ax.legend(handles=legend_elements, loc=legend_loc, fontsize=7)

        if config.show_title:
            ax.set_title(
                config.title or "Phylogenetic GWAS", fontsize=config.title_fontsize
            )
        if config.axis_fontsize:
            ax.xaxis.label.set_fontsize(config.axis_fontsize)
            ax.yaxis.label.set_fontsize(config.axis_fontsize)

        fig.tight_layout()
        fig.savefig(self.output_path, dpi=config.dpi, bbox_inches="tight")
        plt.close(fig)

    # ------------------------------------------------------------------
    # CSV output
    # ------------------------------------------------------------------

    def _write_csv(self, results: List[dict], is_continuous: bool) -> None:
        fieldnames = [
            "position",
            "gene",
            "allele_0",
            "allele_1",
        ]
        if is_continuous:
            fieldnames.append("correlation_r")
        fieldnames.extend(
            [
                "p_value",
                "fdr_p_value",
                "significant",
                "phylo_pattern",
            ]
        )

        with open(self.csv_output, "w", newline="") as fh:
            writer = csv.writer(fh)
            writer.writerow(fieldnames)
            writerow = writer.writerow
            if is_continuous:
                for r in results:
                    try:
                        row = (
                            r["position"],
                            r.get("gene") or "",
                            r["allele_0"],
                            r["allele_1"],
                            r["correlation_r"],
                            r["p_value"],
                            r["fdr_p_value"],
                            "yes" if r.get("fdr_significant") else "no",
                            r.get("phylo_pattern") or "",
                        )
                    except KeyError:
                        row = (
                            r.get("position", ""),
                            r.get("gene") or "",
                            r.get("allele_0", ""),
                            r.get("allele_1", ""),
                            r.get("correlation_r", ""),
                            r.get("p_value", ""),
                            r.get("fdr_p_value", ""),
                            "yes" if r.get("fdr_significant") else "no",
                            r.get("phylo_pattern") or "",
                        )
                    writerow(row)
            else:
                for r in results:
                    try:
                        row = (
                            r["position"],
                            r.get("gene") or "",
                            r["allele_0"],
                            r["allele_1"],
                            r["p_value"],
                            r["fdr_p_value"],
                            "yes" if r.get("fdr_significant") else "no",
                            r.get("phylo_pattern") or "",
                        )
                    except KeyError:
                        row = (
                            r.get("position", ""),
                            r.get("gene") or "",
                            r.get("allele_0", ""),
                            r.get("allele_1", ""),
                            r.get("p_value", ""),
                            r.get("fdr_p_value", ""),
                            "yes" if r.get("fdr_significant") else "no",
                            r.get("phylo_pattern") or "",
                        )
                    writerow(row)
