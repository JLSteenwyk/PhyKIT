from io import StringIO
from pathlib import Path
from typing import Dict, List, Optional

from Bio import Phylo
from scipy.stats import binomtest

from .base import Tree
from ...helpers.json_output import print_json
from ...errors import PhykitUserError


class DiscordanceAsymmetry(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"])
        self.gene_trees_path = parsed["gene_trees_path"]
        self.verbose = parsed["verbose"]
        self.json_output = parsed["json_output"]
        self.plot_output = parsed["plot_output"]

    def process_args(self, args) -> Dict:
        return dict(
            tree_file_path=args.tree,
            gene_trees_path=args.gene_trees,
            verbose=getattr(args, "verbose", False),
            json_output=getattr(args, "json", False),
            plot_output=getattr(args, "plot_output", None),
        )

    def run(self) -> None:
        species_tree = self.read_tree_file()
        gene_trees = self._parse_gene_trees(self.gene_trees_path)

        topology_counts = self._count_topologies(species_tree, gene_trees)

        # Test each branch and collect results
        branch_results = []
        for branch_key in sorted(topology_counts.keys()):
            data = topology_counts[branch_key]
            test_result = self._test_asymmetry(data["n_alt1"], data["n_alt2"])
            entry = dict(
                split=data["split"],
                n_concordant=data["n_concordant"],
                n_alt1=data["n_alt1"],
                n_alt2=data["n_alt2"],
            )
            entry.update(test_result)
            branch_results.append(entry)

        # FDR correction across testable p-values
        testable_indices = []
        testable_pvals = []
        for i, entry in enumerate(branch_results):
            if entry["p_value"] is not None:
                testable_indices.append(i)
                testable_pvals.append(entry["p_value"])

        fdr_corrected = self._fdr(testable_pvals)
        for idx, fdr_p in zip(testable_indices, fdr_corrected):
            branch_results[idx]["fdr_p"] = fdr_p

        # Set fdr_p to None for untestable branches
        for entry in branch_results:
            if "fdr_p" not in entry:
                entry["fdr_p"] = None

        # Summary
        summary = dict(
            n_gene_trees=len(gene_trees),
            n_branches_tested=len(testable_indices),
            n_significant_fdr05=sum(
                1 for entry in branch_results
                if entry["fdr_p"] is not None and entry["fdr_p"] < 0.05
            ),
        )

        # Output
        if self.json_output:
            self._output_json(branch_results, summary)
        else:
            self._output_text(branch_results, summary)

        if self.plot_output:
            self._plot(species_tree, branch_results, self.plot_output)

    # ------------------------------------------------------------------
    # Output methods
    # ------------------------------------------------------------------

    def _output_text(self, branch_results, summary) -> None:
        try:
            header = (
                f"{'branch':<30}"
                f"{'n_conc':>8}"
                f"{'n_alt1':>8}"
                f"{'n_alt2':>8}"
                f"{'asym_ratio':>12}"
                f"{'binom_p':>12}"
                f"{'fdr_p':>12}"
                f"{'gene_flow':>12}"
            )
            print(header)
            print("-" * len(header))

            for entry in branch_results:
                branch_label = ",".join(entry["split"])
                asym = (
                    f"{entry['asymmetry_ratio']:.3f}"
                    if entry["asymmetry_ratio"] is not None
                    else "NA"
                )
                binom_p = (
                    f"{entry['p_value']:.4f}"
                    if entry["p_value"] is not None
                    else "NA"
                )
                fdr_p = (
                    f"{entry['fdr_p']:.4f}"
                    if entry["fdr_p"] is not None
                    else "NA"
                )
                gene_flow = "-"
                if (entry["fdr_p"] is not None and entry["fdr_p"] < 0.05
                        and entry["favored_alt"] is not None):
                    gene_flow = entry["favored_alt"]
                print(
                    f"{branch_label:<30}"
                    f"{entry['n_concordant']:>8}"
                    f"{entry['n_alt1']:>8}"
                    f"{entry['n_alt2']:>8}"
                    f"{asym:>12}"
                    f"{binom_p:>12}"
                    f"{fdr_p:>12}"
                    f"{gene_flow:>12}"
                )

            print("---")
            print(
                f"Summary: {summary['n_branches_tested']} branches tested, "
                f"{summary['n_significant_fdr05']} significant (FDR<0.05)"
            )
        except BrokenPipeError:
            pass

    def _output_json(self, branch_results, summary) -> None:
        result = dict(
            branches=branch_results,
            summary=summary,
        )
        print_json(result)

    def _plot(self, species_tree, branch_results, output_path) -> None:
        """Phylogram colored by asymmetry ratio at each branch."""
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.colors import Normalize
        import numpy as np

        # Build lookup from split label -> branch result
        branch_lookup = {}
        for entry in branch_results:
            key = ",".join(entry["split"])
            branch_lookup[key] = entry

        parent_map = self._build_parent_map(species_tree)
        tips = list(species_tree.get_terminals())
        all_taxa_fs = frozenset(t.name for t in tips)

        # Compute node positions
        node_x = {}
        node_y = {}

        for i, tip in enumerate(tips):
            node_y[id(tip)] = i

        root = species_tree.root
        for clade in species_tree.find_clades(order="preorder"):
            if clade == root:
                node_x[id(clade)] = 0.0
            else:
                if id(clade) in parent_map:
                    parent = parent_map[id(clade)]
                    t = clade.branch_length if clade.branch_length else 0.0
                    node_x[id(clade)] = node_x.get(id(parent), 0.0) + t

        for clade in species_tree.find_clades(order="postorder"):
            if not clade.is_terminal() and id(clade) not in node_y:
                child_ys = [
                    node_y[id(c)] for c in clade.clades if id(c) in node_y
                ]
                if child_ys:
                    node_y[id(clade)] = np.mean(child_ys)
                else:
                    node_y[id(clade)] = 0.0

        # Map internal nodes to their branch result
        node_to_result = {}
        for clade in species_tree.find_clades(order="preorder"):
            if clade.is_terminal():
                continue
            node_tips = frozenset(t.name for t in clade.get_terminals())
            split_label = (
                sorted(node_tips)
                if len(node_tips) <= len(all_taxa_fs) - len(node_tips)
                else sorted(all_taxa_fs - node_tips)
            )
            key = ",".join(split_label)
            if key in branch_lookup:
                node_to_result[id(clade)] = branch_lookup[key]

        # Color setup: diverging from blue (0.5 = symmetric) to red (1.0 = asymmetric)
        cmap = plt.cm.RdYlBu_r
        norm = Normalize(vmin=0.5, vmax=1.0)

        fig, ax = plt.subplots(figsize=(10, max(4, len(tips) * 0.4)))

        # Draw branches
        for clade in species_tree.find_clades(order="preorder"):
            if clade == root:
                continue
            if id(clade) not in parent_map:
                continue
            parent = parent_map[id(clade)]
            if id(parent) not in node_x or id(clade) not in node_x:
                continue

            x0 = node_x[id(parent)]
            x1 = node_x[id(clade)]
            y0 = node_y.get(id(parent), 0)
            y1 = node_y.get(id(clade), 0)

            # Color the horizontal branch by asymmetry ratio if this is an internal node
            color = "gray"
            lw = 2
            if id(clade) in node_to_result:
                entry = node_to_result[id(clade)]
                if entry["asymmetry_ratio"] is not None:
                    color = cmap(norm(entry["asymmetry_ratio"]))
                    lw = 3

            ax.plot([x0, x1], [y1, y1], color=color, lw=lw)
            ax.plot([x0, x0], [y0, y1], color="gray", lw=1.5)

        # Annotate internal nodes
        for clade in species_tree.find_clades(order="preorder"):
            if clade.is_terminal():
                continue
            if id(clade) not in node_to_result:
                continue
            entry = node_to_result[id(clade)]
            x = node_x.get(id(clade), 0)
            y = node_y.get(id(clade), 0)

            # Show gCF value
            total = entry["n_concordant"] + entry["n_alt1"] + entry["n_alt2"]
            gcf = entry["n_concordant"] / total if total > 0 else 1.0
            ax.annotate(
                f"gCF={gcf:.2f}",
                (x, y),
                textcoords="offset points",
                xytext=(5, 5),
                fontsize=7,
            )

            # Mark significant branches (FDR < 0.05)
            if (entry["fdr_p"] is not None and entry["fdr_p"] < 0.05
                    and entry["favored_alt"] is not None):
                ax.scatter(x, y, s=100, c="red", marker="*", zorder=5)

        # Tip labels
        max_x = max(node_x.values()) if node_x else 0
        offset = max_x * 0.02
        for tip in tips:
            ax.text(
                node_x[id(tip)] + offset, node_y[id(tip)],
                tip.name, va="center", fontsize=9,
            )

        # Colorbar
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax, pad=0.15)
        cbar.set_label("Asymmetry ratio")

        ax.set_xlabel("Branch length (subs/site)")
        ax.set_yticks([])
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)
        ax.set_title("Discordance Asymmetry")
        fig.tight_layout()
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
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

    # ------------------------------------------------------------------
    # Bipartition extraction and topology counting
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

    def _count_topologies(self, species_tree, gene_trees) -> Dict:
        """Count concordant and two NNI-alternative topologies for each
        internal branch of the species tree across gene trees.

        Returns a dict keyed by branch label (comma-joined sorted taxa
        in the smaller partition side) with:
          split: list of sorted taxon names
          n_concordant: int
          n_alt1: int
          n_alt2: int
        """
        all_taxa = sorted(
            set(t.name for t in species_tree.get_terminals())
            & set().union(*(
                set(t.name for t in gt.get_terminals()) for gt in gene_trees
            ))
        )
        all_taxa_fs = frozenset(all_taxa)
        parent_map = self._build_parent_map(species_tree)

        # Extract bipartitions from all gene trees (topology only, no lengths)
        gene_tree_splits = []
        for gt in gene_trees:
            gt_taxa = set(t.name for t in gt.get_terminals())
            if gt_taxa != set(all_taxa):
                taxa_to_remove = gt_taxa - set(all_taxa)
                for taxon in taxa_to_remove:
                    gt.prune(taxon)
            splits = set()
            for clade in gt.get_nonterminals():
                tips = frozenset(t.name for t in clade.get_terminals())
                if len(tips) <= 1 or tips == all_taxa_fs:
                    continue
                splits.add(self._canonical_split(tips, all_taxa_fs))
            gene_tree_splits.append(splits)

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

            n_concordant = sum(1 for splits in gene_tree_splits if concordant_bp in splits)
            n_alt1 = sum(1 for splits in gene_tree_splits if nni_alt1_bp in splits)
            n_alt2 = sum(1 for splits in gene_tree_splits if nni_alt2_bp in splits)

            node_tips = frozenset(t.name for t in clade.get_terminals())
            split_label = (
                sorted(node_tips)
                if len(node_tips) <= len(all_taxa_fs) - len(node_tips)
                else sorted(all_taxa_fs - node_tips)
            )
            branch_key = ",".join(split_label)
            result[branch_key] = dict(
                split=split_label,
                n_concordant=n_concordant,
                n_alt1=n_alt1,
                n_alt2=n_alt2,
            )
        return result

    # ------------------------------------------------------------------
    # Statistical testing
    # ------------------------------------------------------------------

    def _test_asymmetry(self, n_alt1: int, n_alt2: int) -> Dict:
        """Run a two-sided binomial test on the two NNI alternatives.

        Returns a dict with asymmetry_ratio, p_value, and favored_alt.
        """
        total = n_alt1 + n_alt2
        if total == 0:
            return dict(
                asymmetry_ratio=None,
                p_value=None,
                favored_alt=None,
            )

        asymmetry_ratio = max(n_alt1, n_alt2) / total
        result = binomtest(n_alt1, total, p=0.5, alternative='two-sided')
        p_value = result.pvalue

        if n_alt1 > n_alt2:
            favored_alt = "alt1"
        elif n_alt2 > n_alt1:
            favored_alt = "alt2"
        else:
            favored_alt = None

        return dict(
            asymmetry_ratio=asymmetry_ratio,
            p_value=p_value,
            favored_alt=favored_alt,
        )

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
