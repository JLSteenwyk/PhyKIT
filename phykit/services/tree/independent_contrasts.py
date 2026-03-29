"""
Phylogenetically independent contrasts (PIC).

Computes Felsenstein's (1985) phylogenetically independent contrasts
for a continuous trait on a phylogeny. Each internal node yields one
contrast (standardized difference), producing n-1 contrasts for n tips.

Cross-validated against R's ape::pic().
"""
import copy
import sys
from typing import Dict, List, Tuple

import numpy as np

from .base import Tree
from ...helpers.json_output import print_json
from ...errors import PhykitUserError


class IndependentContrasts(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"])
        self.trait_data_path = parsed["trait_data_path"]
        self.json_output = parsed["json_output"]

    def run(self) -> None:
        tree = self.read_tree_file()
        tree = copy.deepcopy(tree)
        self.validate_tree(tree, min_tips=3, assign_default_branch_length=1e-8, context="independent contrasts")

        tree_tips = [t.name for t in tree.get_terminals()]
        tip_traits = self._parse_trait_data(self.trait_data_path, tree_tips)

        # Prune tree to shared taxa
        shared = set(tip_traits.keys())
        tips_to_prune = [t for t in tree_tips if t not in shared]
        if tips_to_prune:
            tree = self.prune_tree_using_taxa_list(tree, tips_to_prune)

        # Resolve multifurcations (PIC requires fully dichotomous tree)
        self._resolve_polytomies(tree)

        contrasts, node_labels = self._compute_pic(tree, tip_traits)

        if self.json_output:
            self._print_json(contrasts, node_labels, tip_traits)
        else:
            self._print_text(contrasts, node_labels)

    def process_args(self, args) -> Dict:
        return dict(
            tree_file_path=args.tree,
            trait_data_path=args.trait_data,
            json_output=getattr(args, "json", False),
        )

    def _parse_trait_data(
        self, path: str, tree_tips: List[str]
    ) -> Dict[str, float]:
        """Parse two-column trait file (taxon<tab>value)."""
        try:
            with open(path) as f:
                lines = f.readlines()
        except FileNotFoundError:
            raise PhykitUserError(
                [f"{path} not found. Check filename and path."], code=2
            )

        traits = {}
        for line in lines:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            parts = stripped.split("\t")
            if len(parts) != 2:
                continue
            try:
                traits[parts[0]] = float(parts[1])
            except ValueError:
                continue

        shared = set(tree_tips) & set(traits.keys())
        if len(shared) < 3:
            raise PhykitUserError(
                [
                    f"Only {len(shared)} shared taxa between tree and trait file.",
                    "At least 3 shared taxa are required.",
                ],
                code=2,
            )

        return {t: traits[t] for t in shared}

    def _resolve_polytomies(self, tree) -> None:
        """Resolve multifurcations by adding zero-length branches."""
        for clade in tree.find_clades(order="postorder"):
            while len(clade.clades) > 2:
                # Take the last two children and merge them
                from Bio.Phylo import Newick
                child1 = clade.clades.pop()
                child2 = clade.clades.pop()
                new_internal = Newick.Clade(branch_length=0.0)
                new_internal.clades = [child1, child2]
                clade.clades.append(new_internal)

    def _compute_pic(
        self, tree, tip_traits: Dict[str, float]
    ) -> Tuple[List[float], List[List[str]]]:
        """Compute phylogenetically independent contrasts.

        Felsenstein (1985) algorithm:
        - Postorder traversal through tree
        - At each internal node with children L, R:
          contrast = (x_L - x_R) / sqrt(v_L + v_R)
          x_node = (x_L/v_L + x_R/v_R) / (1/v_L + 1/v_R)
          v_node += v_L * v_R / (v_L + v_R)

        Returns (contrasts, node_tip_labels) where node_tip_labels[i]
        is the list of tip names descending from the node that produced
        contrast[i].
        """
        # Store trait values and effective branch lengths per node
        node_val = {}
        node_bl = {}  # effective branch length (adjusted during traversal)

        # Initialize tips
        for tip in tree.get_terminals():
            if tip.name in tip_traits:
                node_val[id(tip)] = tip_traits[tip.name]
                node_bl[id(tip)] = tip.branch_length if tip.branch_length else 1e-8

        contrasts = []
        node_labels = []

        for clade in tree.find_clades(order="postorder"):
            if clade.is_terminal():
                continue
            if len(clade.clades) != 2:
                continue

            left, right = clade.clades
            lid, rid = id(left), id(right)

            if lid not in node_val or rid not in node_val:
                continue

            x_l = node_val[lid]
            x_r = node_val[rid]
            v_l = node_bl[lid]
            v_r = node_bl[rid]

            # Standardized contrast
            contrast = (x_l - x_r) / np.sqrt(v_l + v_r)
            contrasts.append(float(contrast))

            # Tip labels for this node
            tips_here = sorted(t.name for t in clade.get_terminals())
            node_labels.append(tips_here)

            # Weighted average trait value for this node
            node_val[id(clade)] = (x_l / v_l + x_r / v_r) / (1.0 / v_l + 1.0 / v_r)

            # Adjusted branch length
            parent_bl = clade.branch_length if clade.branch_length else 0.0
            node_bl[id(clade)] = parent_bl + (v_l * v_r) / (v_l + v_r)

        return contrasts, node_labels

    def _print_text(self, contrasts, node_labels):
        print(f"Number of contrasts: {len(contrasts)}")
        print()
        print(f"{'Node':<6}{'Contrast':>12}  Tips")
        print("-" * 60)
        for i, (c, tips) in enumerate(zip(contrasts, node_labels), 1):
            tips_str = ", ".join(tips[:3])
            if len(tips) > 3:
                tips_str += f", ... ({len(tips)} total)"
            print(f"{i:<6}{c:>12.6f}  {tips_str}")
        print()
        print(f"Mean absolute contrast: {np.mean(np.abs(contrasts)):.6f}")
        print(f"Variance of contrasts:  {np.var(contrasts, ddof=1):.6f}")

    def _print_json(self, contrasts, node_labels, tip_traits):
        nodes = []
        for i, (c, tips) in enumerate(zip(contrasts, node_labels)):
            nodes.append({
                "node": i + 1,
                "contrast": round(c, 6),
                "tips": tips,
            })
        payload = {
            "n_taxa": len(tip_traits),
            "n_contrasts": len(contrasts),
            "contrasts": nodes,
            "mean_absolute_contrast": round(float(np.mean(np.abs(contrasts))), 6),
            "variance_of_contrasts": round(float(np.var(contrasts, ddof=1)), 6),
        }
        print_json(payload)
