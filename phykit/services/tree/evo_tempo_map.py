from io import StringIO
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
from Bio import Phylo
from scipy.stats import mannwhitneyu

from .base import Tree
from ...helpers.json_output import print_json
from ...helpers.plot_config import PlotConfig
from ...errors import PhykitUserError


class EvoTempoMap(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"])
        self.gene_trees_path = parsed["gene_trees_path"]
        self.verbose = parsed["verbose"]
        self.json_output = parsed["json_output"]
        self.plot_output = parsed["plot_output"]
        self.plot_config = parsed["plot_config"]

    def process_args(self, args) -> Dict:
        return dict(
            tree_file_path=args.tree,
            gene_trees_path=args.gene_trees,
            verbose=getattr(args, "verbose", False),
            json_output=getattr(args, "json", False),
            plot_output=getattr(args, "plot_output", None),
            plot_config=PlotConfig.from_args(args),
        )

    def run(self) -> None:
        species_tree = self.read_tree_file()
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
            print(header)
            print("-" * len(header))

            for entry in branch_results:
                branch_label = ",".join(entry["split"])
                med_conc = (
                    f"{entry['concordant_median']:.6f}"
                    if entry["concordant_median"] is not None
                    else "NA"
                )
                med_disc = (
                    f"{entry['discordant_median']:.6f}"
                    if entry["discordant_median"] is not None
                    else "NA"
                )
                u_pval = (
                    f"{entry['mann_whitney_p']:.6f}"
                    if entry["mann_whitney_p"] is not None
                    else "NA"
                )
                perm_pval = (
                    f"{entry['permutation_p']:.6f}"
                    if entry["permutation_p"] is not None
                    else "NA"
                )
                fdr_p = (
                    f"{entry['fdr_p']:.6f}"
                    if entry["fdr_p"] is not None
                    else "NA"
                )
                print(
                    f"{branch_label:<30}"
                    f"{entry['n_concordant']:>8}"
                    f"{entry['n_discordant']:>8}"
                    f"{med_conc:>12}"
                    f"{med_disc:>12}"
                    f"{u_pval:>12}"
                    f"{perm_pval:>12}"
                    f"{fdr_p:>12}"
                )

            print("---")

            # Global treeness summary
            conc_t = global_stats["treeness_concordant"]
            disc_t = global_stats["treeness_discordant"]
            conc_med = (
                f"{conc_t['median']:.6f}" if conc_t["median"] is not None else "NA"
            )
            disc_med = (
                f"{disc_t['median']:.6f}" if disc_t["median"] is not None else "NA"
            )
            print(
                f"Global treeness: concordant={conc_med} (n={conc_t['n']}), "
                f"discordant={disc_med} (n={disc_t['n']})"
            )
            print(
                f"Branches tested: {global_stats['n_branches_tested']}, "
                f"significant (FDR<0.05): {global_stats['n_significant_fdr05']}"
            )

            # Verbose: print per-branch raw lengths
            if self.verbose:
                print()
                for entry in branch_results:
                    branch_label = ",".join(entry["split"])
                    print(f"Branch: {branch_label}")
                    print(
                        f"  Concordant lengths: "
                        f"{entry['_concordant_lengths']}"
                    )
                    print(
                        f"  Discordant lengths: "
                        f"{entry['_discordant_lengths']}"
                    )
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
        for i in range(n):
            if conc_data[i] and conc_data[i] != [0]:
                jitter = rng.uniform(-0.05, 0.05, len(conc_data[i]))
                ax.scatter(
                    positions[i] - width / 2 + jitter,
                    conc_data[i], color="#4C72B0", alpha=0.6, s=20, zorder=3,
                )
            if disc_data[i] and disc_data[i] != [0]:
                jitter = rng.uniform(-0.05, 0.05, len(disc_data[i]))
                ax.scatter(
                    positions[i] + width / 2 + jitter,
                    disc_data[i], color="#DD8452", alpha=0.6, s=20, zorder=3,
                )
            if sig_flags[i]:
                max_val = max(
                    max(conc_data[i]) if conc_data[i] else 0,
                    max(disc_data[i]) if disc_data[i] else 0,
                )
                ax.text(positions[i], max_val * 1.05, "*", ha="center", fontsize=14, fontweight="bold")

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

        fig.tight_layout()
        fig.savefig(output_path, dpi=config.dpi, bbox_inches="tight")
        plt.close(fig)

    # ------------------------------------------------------------------
    # Gene tree parsing
    # ------------------------------------------------------------------

    def _parse_gene_trees(self, path: str) -> list:
        try:
            lines = Path(path).read_text().splitlines()
        except FileNotFoundError:
            raise PhykitUserError(
                [
                    f"{path} corresponds to no such file or directory.",
                    "Please check filename and pathing",
                ],
                code=2,
            )

        cleaned = [l.strip() for l in lines if l.strip() and not l.strip().startswith("#")]
        trees = []
        for line in cleaned:
            if line.startswith("("):
                trees.append(Phylo.read(StringIO(line), "newick"))
            else:
                tree_path = Path(path).parent / line
                trees.append(Phylo.read(str(tree_path), "newick"))
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
            for clade in gt.find_clades():
                if clade == gt.root:
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
        complement = all_taxa - taxa_side
        if len(taxa_side) < len(complement):
            return frozenset(taxa_side)
        elif len(taxa_side) > len(complement):
            return frozenset(complement)
        else:
            return min(frozenset(taxa_side), frozenset(complement),
                       key=lambda s: sorted(s))

    @staticmethod
    def _build_parent_map(tree) -> Dict:
        """Build a dict mapping child id -> parent clade."""
        parent_map = {}
        for clade in tree.find_clades(order="preorder"):
            for child in clade.clades:
                parent_map[id(child)] = clade
        return parent_map

    def _get_four_groups(self, tree, node, parent_map, all_taxa_fs):
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

        C1 = frozenset(t.name for t in node.clades[0].get_terminals())
        C2 = frozenset(t.name for t in node.clades[1].get_terminals())
        # If node has >2 children (polytomy), merge extras into C2
        for extra_child in node.clades[2:]:
            C2 = C2 | frozenset(t.name for t in extra_child.get_terminals())

        parent = parent_map.get(id(node))
        if parent is None:
            # node is root — no branch above it
            return None

        # Get siblings of node under parent
        siblings = [c for c in parent.clades if id(c) != id(node)]
        if not siblings:
            return None

        S = frozenset(t.name for t in siblings[0].get_terminals())
        # D = everything else (other siblings + above parent)
        D = all_taxa_fs - C1 - C2 - S

        return C1, C2, S, D

    def _extract_bipartitions_with_lengths(self, tree, all_taxa_fs):
        """Extract non-trivial bipartitions from a tree with branch lengths.

        Returns a dict mapping canonical_split (frozenset) -> branch_length.
        """
        bp_to_length = {}
        for clade in tree.get_nonterminals():
            tips = frozenset(t.name for t in clade.get_terminals())
            if len(tips) <= 1 or tips == all_taxa_fs:
                continue
            bp = self._canonical_split(tips, all_taxa_fs)
            bl = clade.branch_length if clade.branch_length else 0.0
            bp_to_length[bp] = bl
        return bp_to_length

    def _classify_gene_trees(self, species_tree, gene_trees) -> Dict:
        """Classify each gene tree as concordant or discordant at each
        species tree branch, and extract the homologous branch lengths.

        Returns a dict keyed by branch label (comma-joined sorted taxa) with:
          split: list of taxon names in the smaller partition side
          n_concordant: int
          n_discordant: int
          concordant_lengths: list of floats
          discordant_lengths: list of floats
        """
        all_taxa = sorted(
            set(t.name for t in species_tree.get_terminals())
            & set().union(*(
                set(t.name for t in gt.get_terminals()) for gt in gene_trees
            ))
        )
        all_taxa_fs = frozenset(all_taxa)
        parent_map = self._build_parent_map(species_tree)

        # Extract bipartitions + branch lengths from all gene trees
        gene_tree_bp_lengths = []
        for gt in gene_trees:
            # Prune to shared taxa if needed
            gt_taxa = set(t.name for t in gt.get_terminals())
            if gt_taxa != set(all_taxa):
                taxa_to_remove = gt_taxa - set(all_taxa)
                for taxon in taxa_to_remove:
                    gt.prune(taxon)
            gene_tree_bp_lengths.append(
                self._extract_bipartitions_with_lengths(gt, all_taxa_fs)
            )

        result = {}
        for clade in species_tree.find_clades(order="preorder"):
            if clade.is_terminal():
                continue

            groups = self._get_four_groups(
                species_tree, clade, parent_map, all_taxa_fs
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
            node_tips = frozenset(
                t.name for t in clade.get_terminals()
            )
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
        concordant_lengths: List[float],
        discordant_lengths: List[float],
        n_permutations: int = 1000,
    ) -> Optional[Dict]:
        """Compare branch length distributions between concordant and
        discordant gene trees using Mann-Whitney U and a permutation test.

        Returns a dict with summary statistics and p-values, or a dict with
        None values where tests cannot be computed.
        """
        conc = np.array(concordant_lengths, dtype=float)
        disc = np.array(discordant_lengths, dtype=float)

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

        # Summary stats for concordant
        if len(conc) > 0:
            result["concordant_mean"] = float(np.mean(conc))
            result["concordant_median"] = float(np.median(conc))
            result["concordant_std"] = float(np.std(conc, ddof=1)) if len(conc) > 1 else None

        # Summary stats for discordant
        if len(disc) > 0:
            result["discordant_mean"] = float(np.mean(disc))
            result["discordant_median"] = float(np.median(disc))
            result["discordant_std"] = float(np.std(disc, ddof=1)) if len(disc) > 1 else None

        # Need at least 2 observations in each group for tests
        if len(conc) < 2 or len(disc) < 2:
            return result

        # Mann-Whitney U test (two-sided)
        U, mw_p = mannwhitneyu(conc, disc, alternative="two-sided")
        result["mann_whitney_U"] = float(U)
        result["mann_whitney_p"] = float(mw_p)

        # Permutation test on difference in medians
        observed_diff = abs(float(np.median(conc)) - float(np.median(disc)))
        combined = np.concatenate([conc, disc])
        n_conc = len(conc)
        rng = np.random.default_rng(42)

        count = 0
        for _ in range(n_permutations):
            rng.shuffle(combined)
            perm_conc = combined[:n_conc]
            perm_disc = combined[n_conc:]
            perm_diff = abs(float(np.median(perm_conc)) - float(np.median(perm_disc)))
            if perm_diff >= observed_diff:
                count += 1

        result["permutation_p"] = count / n_permutations
        return result

    @staticmethod
    def _fdr(p_values: List[float]) -> List[float]:
        """Benjamini-Hochberg FDR correction."""
        n = len(p_values)
        if n == 0:
            return []
        indexed = sorted(enumerate(p_values), key=lambda x: x[1])
        corrected = [0.0] * n
        prev = 1.0
        for rank_minus_1 in range(n - 1, -1, -1):
            orig_idx, p = indexed[rank_minus_1]
            rank = rank_minus_1 + 1
            adjusted = min(p * n / rank, prev)
            adjusted = min(adjusted, 1.0)
            corrected[orig_idx] = adjusted
            prev = adjusted
        return corrected

    # ------------------------------------------------------------------
    # Global treeness comparison
    # ------------------------------------------------------------------

    @staticmethod
    def _compute_treeness(tree) -> float:
        """Compute treeness = sum(internal branch lengths) / sum(all branch lengths).

        Produces identical results to Tree.calculate_treeness().
        """
        inter_len = 0.0
        for internal in tree.get_nonterminals():
            if internal.branch_length is not None:
                inter_len += internal.branch_length
        total_len = tree.total_branch_length()
        if total_len == 0:
            return 0.0
        return inter_len / total_len

    def _compute_global_treeness(self, species_tree, gene_trees) -> Dict:
        """Classify gene trees as globally concordant or discordant,
        compute treeness for each group, and test for a difference.

        A gene tree is globally concordant if ALL of its non-trivial
        bipartitions match the species tree bipartitions exactly.
        """
        all_taxa = sorted(t.name for t in species_tree.get_terminals())
        all_taxa_fs = frozenset(all_taxa)

        # Species tree bipartitions
        sp_splits = set()
        for clade in species_tree.get_nonterminals():
            tips = frozenset(t.name for t in clade.get_terminals())
            if len(tips) <= 1 or tips == all_taxa_fs:
                continue
            sp_splits.add(self._canonical_split(tips, all_taxa_fs))

        concordant_treeness = []
        discordant_treeness = []

        for gt in gene_trees:
            gt_splits = set()
            for clade in gt.get_nonterminals():
                tips = frozenset(t.name for t in clade.get_terminals())
                if len(tips) <= 1 or tips == all_taxa_fs:
                    continue
                gt_splits.add(self._canonical_split(tips, all_taxa_fs))

            treeness = self._compute_treeness(gt)

            if gt_splits == sp_splits:
                concordant_treeness.append(treeness)
            else:
                discordant_treeness.append(treeness)

        result = dict(
            treeness_concordant=dict(
                mean=float(np.mean(concordant_treeness)) if concordant_treeness else None,
                median=float(np.median(concordant_treeness)) if concordant_treeness else None,
                n=len(concordant_treeness),
            ),
            treeness_discordant=dict(
                mean=float(np.mean(discordant_treeness)) if discordant_treeness else None,
                median=float(np.median(discordant_treeness)) if discordant_treeness else None,
                n=len(discordant_treeness),
            ),
        )

        if len(concordant_treeness) >= 2 and len(discordant_treeness) >= 2:
            _, p = mannwhitneyu(
                concordant_treeness, discordant_treeness,
                alternative="two-sided",
            )
            result["treeness_U_p"] = float(p)
        else:
            result["treeness_U_p"] = None

        return result
