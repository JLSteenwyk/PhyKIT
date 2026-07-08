from __future__ import annotations

import math
import os
from io import StringIO

from .base import Tree
from ...errors import PhykitUserError

_path_isabs = os.path.isabs


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


def Path(*args, **kwargs):
    from pathlib import Path as _Path

    return _Path(*args, **kwargs)


def _mannwhitneyu(x, y, *, alternative):
    from scipy.stats import mannwhitneyu

    return mannwhitneyu(x, y, alternative=alternative)


class _LazyNumpy:
    _module = None

    def __getattr__(self, name):
        module = self._module
        if module is None:
            import numpy as _np

            module = _np
            self._module = module

        value = getattr(module, name)
        setattr(self, name, value)
        return value


class _LazyPhylo:
    def read(self, *args, **kwargs):
        from Bio import Phylo as _Phylo

        self.read = _Phylo.read
        return self.read(*args, **kwargs)


np = _LazyNumpy()
Phylo = _LazyPhylo()
_FDR_VECTOR_MIN_LENGTH = 32
_MANN_WHITNEY_EXACT_MAX_MIN_N = 8
_MANN_WHITNEY_EXACT_COUNTS = {}
_MANN_WHITNEY_EXACT_CUMULATIVE = {}


def _median(values: list[float]) -> float:
    ordered = sorted(values)
    mid = len(ordered) // 2
    if len(ordered) % 2:
        return float(ordered[mid])
    return float((ordered[mid - 1] + ordered[mid]) / 2.0)


def _sample_std(values: list[float], mean: float) -> float | None:
    n = len(values)
    if n < 2:
        return None
    sum_squared_deviations = 0.0
    for value in values:
        delta = value - mean
        sum_squared_deviations += delta * delta
    variance = sum_squared_deviations / (n - 1)
    return math.sqrt(variance)


def _summarize_lengths(
    values: list[float],
) -> tuple[float | None, float | None, float | None]:
    if len(values) == 0:
        return None, None, None
    mean = sum(values) / len(values)
    return float(mean), _median(values), _sample_std(values, mean)


def _summarize_mean_median(
    values: list[float],
) -> tuple[float | None, float | None]:
    if len(values) == 0:
        return None, None
    return float(sum(values) / len(values)), _median(values)


def _mann_whitney_exact_counts(n_x: int, n_y: int) -> list[int]:
    key = (n_x, n_y)
    counts = _MANN_WHITNEY_EXACT_COUNTS.get(key)
    if counts is not None:
        return counts

    dp = [[None] * (n_y + 1) for _ in range(n_x + 1)]
    dp[0][0] = [1]
    for i in range(n_x + 1):
        for j in range(n_y + 1):
            if i == 0 and j == 0:
                continue
            max_u = i * j
            current = [0] * (max_u + 1)
            if i:
                previous = dp[i - 1][j]
                shift = j
                for u_value, count in enumerate(previous):
                    current[u_value + shift] += count
            if j:
                previous = dp[i][j - 1]
                for u_value, count in enumerate(previous):
                    current[u_value] += count
            dp[i][j] = current

    counts = dp[n_x][n_y]
    _MANN_WHITNEY_EXACT_COUNTS[key] = counts
    return counts


def _mann_whitney_exact_cumulative(n_x: int, n_y: int) -> tuple[list[int], int]:
    key = (n_x, n_y)
    cached = _MANN_WHITNEY_EXACT_CUMULATIVE.get(key)
    if cached is not None:
        return cached

    cumulative = []
    running_total = 0
    for count in _mann_whitney_exact_counts(n_x, n_y):
        running_total += count
        cumulative.append(running_total)

    result = (cumulative, running_total)
    _MANN_WHITNEY_EXACT_CUMULATIVE[key] = result
    return result


def _mannwhitneyu_no_ties(
    x: list[float], y: list[float]
) -> tuple[float, float] | None:
    n_x = len(x)
    n_y = len(y)
    if n_x == 0 or n_y == 0:
        return None

    combined = [(float(value), 0) for value in x]
    combined.extend((float(value), 1) for value in y)
    if any(not math.isfinite(value) for value, _ in combined):
        return None
    if len({value for value, _ in combined}) != len(combined):
        return None

    rank_sum_x = 0.0
    for rank, (_, group) in enumerate(sorted(combined), start=1):
        if group == 0:
            rank_sum_x += rank
    u_statistic = rank_sum_x - (n_x * (n_x + 1) / 2.0)

    if min(n_x, n_y) <= _MANN_WHITNEY_EXACT_MAX_MIN_N:
        u_index = int(u_statistic)
        cumulative, total = _mann_whitney_exact_cumulative(n_x, n_y)
        lower = cumulative[u_index]
        upper = total - (cumulative[u_index - 1] if u_index else 0)
        p_value = min(1.0, 2.0 * min(lower, upper) / total)
        return float(u_statistic), float(p_value)

    mean = n_x * n_y / 2.0
    variance = n_x * n_y * (n_x + n_y + 1) / 12.0
    z_value = (abs(u_statistic - mean) - 0.5) / math.sqrt(variance)
    p_value = math.erfc(z_value / math.sqrt(2.0))
    return float(u_statistic), float(p_value)


class EvoTempoMap(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"])
        self.gene_trees_path = parsed["gene_trees_path"]
        self.verbose = parsed["verbose"]
        self.json_output = parsed["json_output"]
        self.plot_output = parsed["plot_output"]
        self.plot_config = parsed["plot_config"]

    def process_args(self, args) -> dict:
        from ...helpers.plot_config import PlotConfig

        return dict(
            tree_file_path=args.tree,
            gene_trees_path=args.gene_trees,
            verbose=getattr(args, "verbose", False),
            json_output=getattr(args, "json", False),
            plot_output=getattr(args, "plot_output", None),
            plot_config=PlotConfig.from_args(args),
        )

    def run(self) -> None:
        species_tree = self.read_tree_file_unmodified()
        gene_trees = self._parse_and_validate_gene_trees()
        classification = self._classify_gene_trees(species_tree, gene_trees)

        # Test each branch and collect results
        branch_results = []
        for branch_key in sorted(classification.keys()):
            data = classification[branch_key]
            test_result = self._test_branch(
                data["concordant_lengths"],
                data["discordant_lengths"],
            )
            entry = dict(
                split=data["split"],
                n_concordant=data["n_concordant"],
                n_discordant=data["n_discordant"],
                _concordant_lengths=data["concordant_lengths"],
                _discordant_lengths=data["discordant_lengths"],
            )
            entry.update(test_result)
            branch_results.append(entry)

        # FDR correction across testable p-values
        testable_indices = []
        testable_pvals = []
        for i, entry in enumerate(branch_results):
            if entry["mann_whitney_p"] is not None:
                testable_indices.append(i)
                testable_pvals.append(entry["mann_whitney_p"])

        fdr_corrected = self._fdr(testable_pvals)
        for idx, fdr_p in zip(testable_indices, fdr_corrected):
            branch_results[idx]["fdr_p"] = fdr_p

        # Set fdr_p to None for untestable branches
        for entry in branch_results:
            if "fdr_p" not in entry:
                entry["fdr_p"] = None

        # Global treeness
        global_stats = self._compute_global_treeness(species_tree, gene_trees)
        global_stats["n_gene_trees"] = len(gene_trees)
        global_stats["n_branches_tested"] = len(testable_indices)
        global_stats["n_significant_fdr05"] = sum(
            1 for entry in branch_results
            if entry["fdr_p"] is not None and entry["fdr_p"] < 0.05
        )

        # Output
        if self.json_output:
            self._output_json(branch_results, global_stats)
        else:
            self._output_text(branch_results, global_stats)

        if self.plot_output:
            self._plot(branch_results, self.plot_output)

    def _output_text(self, branch_results, global_stats) -> None:
        try:
            lines = []
            append = lines.append
            # Header
            header = (
                f"{'branch':<30}"
                f"{'n_conc':>8}"
                f"{'n_disc':>8}"
                f"{'med_conc':>12}"
                f"{'med_disc':>12}"
                f"{'U_pval':>12}"
                f"{'perm_pval':>12}"
                f"{'fdr_p':>12}"
            )
            append(header)
            append("-" * len(header))

            row_format = (
                "{:<30}{:>8}{:>8}{:>12}{:>12}{:>12}{:>12}{:>12}"
            ).format
            verbose_lines = [] if self.verbose else None
            verbose_append = (
                verbose_lines.append if verbose_lines is not None else None
            )
            for entry in branch_results:
                branch_label = ",".join(entry["split"])
                concordant_median = entry["concordant_median"]
                discordant_median = entry["discordant_median"]
                mann_whitney_p = entry["mann_whitney_p"]
                permutation_p = entry["permutation_p"]
                fdr_value = entry["fdr_p"]
                med_conc = (
                    f"{concordant_median:.6f}"
                    if concordant_median is not None
                    else "NA"
                )
                med_disc = (
                    f"{discordant_median:.6f}"
                    if discordant_median is not None
                    else "NA"
                )
                u_pval = (
                    f"{mann_whitney_p:.6f}"
                    if mann_whitney_p is not None
                    else "NA"
                )
                perm_pval = (
                    f"{permutation_p:.6f}"
                    if permutation_p is not None
                    else "NA"
                )
                fdr_p = (
                    f"{fdr_value:.6f}"
                    if fdr_value is not None
                    else "NA"
                )
                append(
                    row_format(
                        branch_label,
                        entry["n_concordant"],
                        entry["n_discordant"],
                        med_conc,
                        med_disc,
                        u_pval,
                        perm_pval,
                        fdr_p,
                    )
                )
                if verbose_append is not None:
                    verbose_append(f"Branch: {branch_label}")
                    verbose_append(
                        f"  Concordant lengths: "
                        f"{entry['_concordant_lengths']}"
                    )
                    verbose_append(
                        f"  Discordant lengths: "
                        f"{entry['_discordant_lengths']}"
                    )

            append("---")

            # Global treeness summary
            conc_t = global_stats["treeness_concordant"]
            disc_t = global_stats["treeness_discordant"]
            conc_med = (
                f"{conc_t['median']:.6f}" if conc_t["median"] is not None else "NA"
            )
            disc_med = (
                f"{disc_t['median']:.6f}" if disc_t["median"] is not None else "NA"
            )
            append(
                f"Global treeness: concordant={conc_med} (n={conc_t['n']}), "
                f"discordant={disc_med} (n={disc_t['n']})"
            )
            append(
                f"Branches tested: {global_stats['n_branches_tested']}, "
                f"significant (FDR<0.05): {global_stats['n_significant_fdr05']}"
            )

            # Verbose: print per-branch raw lengths
            if verbose_lines is not None:
                append("")
                lines.extend(verbose_lines)
            print("\n".join(lines))
        except BrokenPipeError:
            pass

    def _output_json(self, branch_results, global_stats) -> None:
        result = dict(
            branches=[
                {k: v for k, v in entry.items() if not k.startswith("_")}
                for entry in branch_results
            ],
            global_=global_stats,
        )
        print_json(result)

    def _plot(self, branch_results, output_path):
        """Grouped box/strip plot of concordant vs discordant branch lengths."""
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np

        labels = []
        conc_data = []
        disc_data = []
        sig_flags = []

        for entry in branch_results:
            conc = entry.get("_concordant_lengths", [])
            disc = entry.get("_discordant_lengths", [])
            if not conc and not disc:
                continue
            split_str = "{" + ",".join(entry["split"]) + "}"
            labels.append(split_str)
            conc_data.append(conc if conc else [0])  # need at least 1 value for boxplot
            disc_data.append(disc if disc else [0])
            sig_flags.append(
                entry.get("fdr_p") is not None and entry["fdr_p"] < 0.05
            )

        if not labels:
            return

        n = len(labels)
        config = self.plot_config
        config.resolve(n_rows=n, n_cols=None)
        fig, ax = plt.subplots(figsize=(config.fig_width, config.fig_height))

        positions = np.arange(n)
        width = 0.35

        # Box plots for concordant
        bp_conc = ax.boxplot(
            conc_data, positions=positions - width / 2, widths=width * 0.8,
            patch_artist=True, showfliers=False,
            boxprops=dict(facecolor="#4C72B0", alpha=0.7),
            medianprops=dict(color="black"),
        )
        # Box plots for discordant
        bp_disc = ax.boxplot(
            disc_data, positions=positions + width / 2, widths=width * 0.8,
            patch_artist=True, showfliers=False,
            boxprops=dict(facecolor="#DD8452", alpha=0.7),
            medianprops=dict(color="black"),
        )

        # Strip (jitter) points
        rng = np.random.default_rng(42)
        conc_x = []
        conc_y = []
        disc_x = []
        disc_y = []
        sig_positions = []
        for i in range(n):
            if conc_data[i] and conc_data[i] != [0]:
                jitter = rng.uniform(-0.05, 0.05, len(conc_data[i]))
                conc_x.extend(positions[i] - width / 2 + jitter)
                conc_y.extend(conc_data[i])
            if disc_data[i] and disc_data[i] != [0]:
                jitter = rng.uniform(-0.05, 0.05, len(disc_data[i]))
                disc_x.extend(positions[i] + width / 2 + jitter)
                disc_y.extend(disc_data[i])
            if sig_flags[i]:
                max_val = max(
                    max(conc_data[i]) if conc_data[i] else 0,
                    max(disc_data[i]) if disc_data[i] else 0,
                )
                sig_positions.append((positions[i], max_val * 1.05))

        if conc_x:
            ax.scatter(
                conc_x, conc_y, color="#4C72B0", alpha=0.6, s=20, zorder=3,
            )
        if disc_x:
            ax.scatter(
                disc_x, disc_y, color="#DD8452", alpha=0.6, s=20, zorder=3,
            )
        if sig_positions:
            sig_x, sig_y = zip(*sig_positions)
            ax.scatter(
                sig_x,
                sig_y,
                marker="$*$",
                s=170,
                c="black",
                linewidths=0,
                zorder=4,
            )

        ax.set_xticks(positions)
        ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=9)
        ax.set_ylabel("Branch length (subs/site)", fontsize=12)
        ax.set_xlabel("Species tree branch", fontsize=12)
        if config.show_title:
            ax.set_title(config.title or "Evolutionary Tempo Map", fontsize=config.title_fontsize)
        ax.legend(
            [bp_conc["boxes"][0], bp_disc["boxes"][0]],
            ["Concordant", "Discordant"],
            loc="upper right",
        )
        if config.axis_fontsize:
            ax.xaxis.label.set_fontsize(config.axis_fontsize)
            ax.yaxis.label.set_fontsize(config.axis_fontsize)

        fig.savefig(output_path, dpi=config.dpi, bbox_inches="tight")
        plt.close(fig)

    # ------------------------------------------------------------------
    # Gene tree parsing
    # ------------------------------------------------------------------

    def _parse_gene_trees(self, path: str) -> list:
        source = Path(path)
        try:
            with source.open() as handle:
                cleaned = [
                    stripped
                    for line in handle
                    if (stripped := line.strip())
                    and stripped[0] != "#"
                ]
        except FileNotFoundError:
            raise PhykitUserError(
                [
                    f"{path} corresponds to no such file or directory.",
                    "Please check filename and pathing",
                ],
                code=2,
            )

        trees = []
        parent_str = str(source.parent)
        parent_prefix = "" if parent_str == "." else parent_str + os.sep
        for line in cleaned:
            if line.startswith("("):
                trees.append(Phylo.read(StringIO(line), "newick"))
            else:
                tree_path = line if _path_isabs(line) else parent_prefix + line
                trees.append(Phylo.read(tree_path, "newick"))
        return trees

    def _parse_and_validate_gene_trees(self) -> list:
        gene_trees = self._parse_gene_trees(self.gene_trees_path)

        if len(gene_trees) < 2:
            raise PhykitUserError(
                [
                    "At least 2 gene trees are required for evolutionary",
                    "tempo mapping. Please provide more gene trees.",
                ],
                code=2,
            )

        for i, gt in enumerate(gene_trees):
            clades = self._iter_preorder_clades(gt)
            if clades is None:
                clades = gt.find_clades()
            root = gt.root
            for clade in clades:
                if clade is root:
                    continue
                if clade.branch_length is None:
                    raise PhykitUserError(
                        [
                            f"Gene tree {i + 1} has branches without branch lengths.",
                            "All gene tree branches must have branch lengths.",
                        ],
                        code=2,
                    )

        return gene_trees

    # ------------------------------------------------------------------
    # Bipartition extraction and concordance classification
    # ------------------------------------------------------------------

    @staticmethod
    def _canonical_split(taxa_side, all_taxa):
        """Normalize a bipartition to canonical form.

        Returns the smaller side as a frozenset; ties are broken
        lexicographically.
        """
        taxa_side_len = len(taxa_side)
        all_taxa_len = len(all_taxa)
        if taxa_side_len * 2 < all_taxa_len:
            return frozenset(taxa_side)
        complement = all_taxa - taxa_side
        if taxa_side_len * 2 > all_taxa_len:
            return frozenset(complement)
        if not taxa_side:
            return frozenset(taxa_side)
        if min(taxa_side) <= min(complement):
            return frozenset(taxa_side)
        return frozenset(complement)

    @staticmethod
    def _build_parent_map(tree) -> dict:
        """Build a dict mapping child id -> parent clade."""
        parent_map = {}
        try:
            root = tree.root
            root.clades
        except AttributeError:
            root = None

        if root is not None:
            stack = [root]
            pop = stack.pop
            extend = stack.extend
            try:
                while stack:
                    clade = pop()
                    children = clade.clades
                    for child in children:
                        parent_map[id(child)] = clade
                    if children:
                        extend(children)
            except AttributeError:
                parent_map = {}
            else:
                return parent_map

        for clade in tree.find_clades(order="preorder"):
            for child in clade.clades:
                parent_map[id(child)] = clade
        return parent_map

    @staticmethod
    def _collect_clade_taxa(tree) -> dict[int, frozenset]:
        clade_taxa: dict[int, frozenset] = {}
        postorder_clades = EvoTempoMap._iter_postorder_clades(tree)
        if postorder_clades is None:
            postorder_clades = tree.find_clades(order="postorder")
        empty = frozenset()
        for clade in postorder_clades:
            children = clade.clades
            if not children:
                clade_taxa[id(clade)] = frozenset((clade.name,))
            else:
                child_count = len(children)
                if child_count == 2:
                    clade_taxa[id(clade)] = (
                        clade_taxa.get(id(children[0]), empty)
                        | clade_taxa.get(id(children[1]), empty)
                    )
                else:
                    taxa = set()
                    for child in children:
                        taxa.update(clade_taxa.get(id(child), empty))
                    clade_taxa[id(clade)] = frozenset(taxa)
        return clade_taxa

    @staticmethod
    def _iter_preorder_clades(tree):
        try:
            root = tree.root
            root.clades
        except AttributeError:
            return None

        clades = []
        stack = [root]
        try:
            pop = stack.pop
            append = stack.append
            append_clade = clades.append
            while stack:
                clade = pop()
                append_clade(clade)
                children = clade.clades
                if children:
                    child_count = len(children)
                    if child_count == 2:
                        append(children[1])
                        append(children[0])
                    else:
                        for index in range(child_count - 1, -1, -1):
                            append(children[index])
        except AttributeError:
            return None
        return clades

    @staticmethod
    def _iter_postorder_clades(tree):
        try:
            root = tree.root
            root.clades
        except AttributeError:
            return None

        clades = []
        stack = [root]
        try:
            while stack:
                clade = stack.pop()
                clades.append(clade)
                children = clade.clades
                if children:
                    stack.extend(children)
        except AttributeError:
            return None
        clades.reverse()
        return clades

    def _split_set_from_clade_taxa(
        self, tree, all_taxa_fs, clade_taxa, clades=None
    ):
        if clades is None:
            clades = self._iter_preorder_clades(tree)
            if clades is None:
                clades = tree.get_nonterminals()

        splits = set()
        for clade in clades:
            if not clade.clades:
                continue
            tips = clade_taxa.get(id(clade), frozenset())
            if len(tips) <= 1 or tips == all_taxa_fs:
                continue
            splits.add(self._canonical_split(tips, all_taxa_fs))
        return splits

    def _get_four_groups(self, tree, node, parent_map, all_taxa_fs, clade_taxa=None):
        """Identify the four subtree groups around an internal branch.

        For the branch connecting *node* to its parent:
          C1 = tips of node's first child
          C2 = tips of node's second child (extra children merged for polytomies)
          S  = tips of node's sibling under parent
          D  = remaining tips (everything above parent)

        Returns (C1, C2, S, D) as frozensets, or None if decomposition
        is not possible (e.g., node is root, leaf, or has <2 children).
        """
        if node.is_terminal() or len(node.clades) < 2:
            return None

        if clade_taxa is None:
            clade_taxa = self._collect_clade_taxa(tree)

        empty = frozenset()
        children = node.clades
        C1 = clade_taxa.get(id(children[0]), empty)
        C2 = clade_taxa.get(id(children[1]), empty)
        # If node has >2 children (polytomy), merge extras into C2
        if len(children) > 2:
            C2 = C2.union(
                *(
                    clade_taxa.get(id(extra_child), empty)
                    for extra_child in children[2:]
                )
            )

        parent = parent_map.get(id(node))
        if parent is None:
            # node is root — no branch above it
            return None

        # Get siblings of node under parent
        siblings = [c for c in parent.clades if id(c) != id(node)]
        if not siblings:
            return None

        S = clade_taxa.get(id(siblings[0]), empty)
        # D = everything else (other siblings + above parent)
        D = all_taxa_fs - C1 - C2 - S

        return C1, C2, S, D

    def _extract_bipartitions_with_lengths(self, tree, all_taxa_fs, clade_taxa=None):
        """Extract non-trivial bipartitions from a tree with branch lengths.

        Returns a dict mapping canonical_split (frozenset) -> branch_length.
        """
        if clade_taxa is None:
            clade_taxa = {}
            bp_to_length = {}
            postorder_clades = self._iter_postorder_clades(tree)
            if postorder_clades is None:
                postorder_clades = tree.find_clades(order="postorder")
            empty = frozenset()
            for clade in postorder_clades:
                children = clade.clades
                if not children:
                    clade_taxa[id(clade)] = frozenset((clade.name,))
                    continue

                child_count = len(children)
                if child_count == 2:
                    taxa = (
                        clade_taxa.get(id(children[0]), empty)
                        | clade_taxa.get(id(children[1]), empty)
                    )
                else:
                    taxa_set = set()
                    for child in children:
                        taxa_set.update(clade_taxa.get(id(child), empty))
                    taxa = frozenset(taxa_set)
                clade_taxa[id(clade)] = taxa

                if len(taxa) <= 1 or taxa == all_taxa_fs:
                    continue
                bp = self._canonical_split(taxa, all_taxa_fs)
                bl = clade.branch_length if clade.branch_length else 0.0
                bp_to_length[bp] = bl
            return bp_to_length

        bp_to_length = {}
        clades = self._iter_preorder_clades(tree)
        if clades is None:
            clades = tree.get_nonterminals()
        for clade in clades:
            if clade.is_terminal():
                continue
            tips = clade_taxa.get(id(clade), frozenset())
            if len(tips) <= 1 or tips == all_taxa_fs:
                continue
            bp = self._canonical_split(tips, all_taxa_fs)
            bl = clade.branch_length if clade.branch_length else 0.0
            bp_to_length[bp] = bl
        return bp_to_length

    def _classify_gene_trees(self, species_tree, gene_trees) -> dict:
        """Classify each gene tree as concordant or discordant at each
        species tree branch, and extract the homologous branch lengths.

        Returns a dict keyed by branch label (comma-joined sorted taxa) with:
          split: list of taxon names in the smaller partition side
          n_concordant: int
          n_discordant: int
          concordant_lengths: list of floats
          discordant_lengths: list of floats
        """
        species_tips = set(self.get_tip_names_from_tree(species_tree))
        gene_tip_sets = [
            set(self.get_tip_names_from_tree(gt)) for gt in gene_trees
        ]
        all_taxa = sorted(species_tips & set().union(*gene_tip_sets))
        all_taxa_fs = frozenset(all_taxa)
        parent_map = self._build_parent_map(species_tree)
        species_clade_taxa = self._collect_clade_taxa(species_tree)

        # Extract bipartitions + branch lengths from all gene trees
        gene_tree_bp_lengths = []
        for gt, gt_taxa in zip(gene_trees, gene_tip_sets):
            # Prune to shared taxa if needed
            if gt_taxa != all_taxa_fs:
                taxa_to_remove = gt_taxa - all_taxa_fs
                terminals = gt.get_terminals()
                for tip in terminals:
                    if tip.name in taxa_to_remove:
                        gt.prune(tip)
            gene_tree_bp_lengths.append(
                self._extract_bipartitions_with_lengths(gt, all_taxa_fs)
            )

        result = {}
        species_preorder = self._iter_preorder_clades(species_tree)
        if species_preorder is None:
            species_preorder = species_tree.find_clades(order="preorder")
        for clade in species_preorder:
            if clade.is_terminal():
                continue

            groups = self._get_four_groups(
                species_tree, clade, parent_map, all_taxa_fs, species_clade_taxa
            )
            if groups is None:
                continue

            C1, C2, S, D = groups

            concordant_bp = self._canonical_split(C1 | C2, all_taxa_fs)
            nni_alt1_bp = self._canonical_split(S | C2, all_taxa_fs)
            nni_alt2_bp = self._canonical_split(C1 | S, all_taxa_fs)

            concordant_lengths = []
            discordant_lengths = []

            for bp_lengths in gene_tree_bp_lengths:
                if concordant_bp in bp_lengths:
                    concordant_lengths.append(bp_lengths[concordant_bp])
                elif nni_alt1_bp in bp_lengths:
                    discordant_lengths.append(bp_lengths[nni_alt1_bp])
                elif nni_alt2_bp in bp_lengths:
                    discordant_lengths.append(bp_lengths[nni_alt2_bp])

            # Label by the smaller side of the species tree split
            node_tips = species_clade_taxa.get(id(clade), frozenset())
            split_label = (
                sorted(node_tips)
                if len(node_tips) <= len(all_taxa_fs) - len(node_tips)
                else sorted(all_taxa_fs - node_tips)
            )

            branch_key = ",".join(split_label)
            result[branch_key] = dict(
                split=split_label,
                n_concordant=len(concordant_lengths),
                n_discordant=len(discordant_lengths),
                concordant_lengths=concordant_lengths,
                discordant_lengths=discordant_lengths,
            )

        return result

    # ------------------------------------------------------------------
    # Statistical testing
    # ------------------------------------------------------------------

    def _test_branch(
        self,
        concordant_lengths: list[float],
        discordant_lengths: list[float],
        n_permutations: int = 1000,
    ) -> dict | None:
        """Compare branch length distributions between concordant and
        discordant gene trees using Mann-Whitney U and a permutation test.

        Returns a dict with summary statistics and p-values, or a dict with
        None values where tests cannot be computed.
        """
        result = dict(
            concordant_mean=None,
            concordant_median=None,
            concordant_std=None,
            discordant_mean=None,
            discordant_median=None,
            discordant_std=None,
            mann_whitney_U=None,
            mann_whitney_p=None,
            permutation_p=None,
        )

        (
            result["concordant_mean"],
            result["concordant_median"],
            result["concordant_std"],
        ) = _summarize_lengths(concordant_lengths)
        (
            result["discordant_mean"],
            result["discordant_median"],
            result["discordant_std"],
        ) = _summarize_lengths(discordant_lengths)

        # Need at least 2 observations in each group for tests
        if len(concordant_lengths) < 2 or len(discordant_lengths) < 2:
            return result

        # Mann-Whitney U test (two-sided)
        exact_result = _mannwhitneyu_no_ties(
            concordant_lengths, discordant_lengths
        )
        conc = np.array(concordant_lengths, dtype=float)
        disc = np.array(discordant_lengths, dtype=float)
        if exact_result is None:
            U, mw_p = _mannwhitneyu(conc, disc, alternative="two-sided")
        else:
            U, mw_p = exact_result
        result["mann_whitney_U"] = float(U)
        result["mann_whitney_p"] = float(mw_p)

        # Permutation test on difference in medians
        median = np.median
        observed_diff = abs(float(median(conc)) - float(median(disc)))
        combined = np.concatenate([conc, disc])
        n_conc = len(conc)
        rng = np.random.default_rng(42)
        shuffle = rng.shuffle

        count = 0
        for _ in range(n_permutations):
            shuffle(combined)
            perm_conc = combined[:n_conc]
            perm_disc = combined[n_conc:]
            perm_diff = abs(float(median(perm_conc)) - float(median(perm_disc)))
            if perm_diff >= observed_diff:
                count += 1

        result["permutation_p"] = count / n_permutations
        return result

    @staticmethod
    def _fdr(p_values: list[float]) -> list[float]:
        """Benjamini-Hochberg FDR correction."""
        n = len(p_values)
        if n == 0:
            return []
        if n < _FDR_VECTOR_MIN_LENGTH:
            indexed = sorted(enumerate(p_values), key=lambda item: item[1])
            corrected = [0.0] * n
            previous = 1.0
            for rank_index in range(n - 1, -1, -1):
                original_index, p_value = indexed[rank_index]
                rank = rank_index + 1
                adjusted = min(p_value * n / rank, previous)
                adjusted = min(adjusted, 1.0)
                corrected[original_index] = adjusted
                previous = adjusted
            return corrected
        p_arr = np.asarray(p_values, dtype=float)
        order = np.argsort(p_arr)
        ranks = np.arange(1, n + 1, dtype=float)
        adjusted = p_arr[order].copy()
        adjusted *= n
        adjusted /= ranks
        adjusted_reversed = adjusted[::-1]
        np.minimum.accumulate(adjusted_reversed, out=adjusted_reversed)
        np.minimum(adjusted, 1.0, out=adjusted)
        corrected = np.empty(n, dtype=float)
        corrected[order] = adjusted
        return corrected.tolist()

    # ------------------------------------------------------------------
    # Global treeness comparison
    # ------------------------------------------------------------------

    @staticmethod
    def _compute_treeness(tree) -> float:
        """Compute treeness = sum(internal branch lengths) / sum(all branch lengths).

        Produces identical results to Tree.calculate_treeness().
        """
        inter_len, total_len = Tree.calculate_internal_and_total_branch_length_fast(
            tree
        )
        if total_len == 0:
            return 0.0
        return inter_len / total_len

    def _compute_global_treeness(self, species_tree, gene_trees) -> dict:
        """Classify gene trees as globally concordant or discordant,
        compute treeness for each group, and test for a difference.

        A gene tree is globally concordant if ALL of its non-trivial
        bipartitions match the species tree bipartitions exactly.
        """
        species_clade_taxa = self._collect_clade_taxa(species_tree)
        all_taxa_fs = species_clade_taxa.get(id(species_tree.root), frozenset())

        # Species tree bipartitions
        sp_splits = self._split_set_from_clade_taxa(
            species_tree,
            all_taxa_fs,
            species_clade_taxa,
        )

        concordant_treeness = []
        discordant_treeness = []

        for gt in gene_trees:
            gt_clade_taxa = self._collect_clade_taxa(gt)
            gt_splits = self._split_set_from_clade_taxa(
                gt,
                all_taxa_fs,
                gt_clade_taxa,
            )

            treeness = self._compute_treeness(gt)

            if gt_splits == sp_splits:
                concordant_treeness.append(treeness)
            else:
                discordant_treeness.append(treeness)

        concordant_mean, concordant_median = _summarize_mean_median(
            concordant_treeness
        )
        discordant_mean, discordant_median = _summarize_mean_median(
            discordant_treeness
        )

        result = dict(
            treeness_concordant=dict(
                mean=concordant_mean,
                median=concordant_median,
                n=len(concordant_treeness),
            ),
            treeness_discordant=dict(
                mean=discordant_mean,
                median=discordant_median,
                n=len(discordant_treeness),
            ),
        )

        if len(concordant_treeness) >= 2 and len(discordant_treeness) >= 2:
            exact_result = _mannwhitneyu_no_ties(
                concordant_treeness, discordant_treeness
            )
            if exact_result is None:
                _, p = _mannwhitneyu(
                    concordant_treeness, discordant_treeness,
                    alternative="two-sided",
                )
            else:
                _, p = exact_result
            result["treeness_U_p"] = float(p)
        else:
            result["treeness_U_p"] = None

        return result
