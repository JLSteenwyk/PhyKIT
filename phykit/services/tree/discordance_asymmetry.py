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
        pass

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
