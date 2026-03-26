import copy
import json
import os
from io import StringIO
from typing import Dict, List, Tuple

from Bio import Phylo
from Bio.Phylo.Newick import Clade

from .base import Tree
from ...errors import PhykitUserError
from ...helpers.json_output import print_json


class Spr(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"])
        self.subtree_taxa = parsed["subtree_taxa"]
        self.output_path = parsed["output_path"]
        self.json_output = parsed["json_output"]

    def process_args(self, args) -> Dict[str, object]:
        taxa_str = args.subtree
        # Support both a file (one taxon per line) and comma-separated taxa
        if os.path.isfile(taxa_str):
            with open(taxa_str) as f:
                taxa = [line.strip() for line in f if line.strip()]
        else:
            taxa = [t.strip() for t in taxa_str.split(",")]
        return dict(
            tree_file_path=args.tree,
            subtree_taxa=taxa,
            output_path=getattr(args, "output", None),
            json_output=getattr(args, "json", False),
        )

    def run(self) -> None:
        tree = self.read_tree_file()

        # Validate taxa exist
        tree_tips = {t.name for t in tree.get_terminals()}
        for taxon in self.subtree_taxa:
            if taxon not in tree_tips:
                raise PhykitUserError(
                    [
                        f"Taxon '{taxon}' not found in tree. "
                        f"Available: {', '.join(sorted(tree_tips))}"
                    ],
                    code=2,
                )

        # Find MRCA
        if len(self.subtree_taxa) == 1:
            subtree = None
            for tip in tree.get_terminals():
                if tip.name == self.subtree_taxa[0]:
                    subtree = tip
                    break
        else:
            subtree = tree.common_ancestor(self.subtree_taxa)

        if subtree is None or subtree == tree.root:
            raise PhykitUserError(
                ["Cannot prune the root of the tree."], code=2
            )

        # Generate SPR trees
        spr_trees = self._generate_spr_trees(tree, subtree)

        # Output
        if self.json_output:
            trees_payload = []
            for desc, spr_tree in spr_trees:
                s = StringIO()
                Phylo.write(spr_tree, s, "newick")
                trees_payload.append(
                    {"description": desc, "newick": s.getvalue().strip()}
                )
            print_json(
                {
                    "subtree_taxa": self.subtree_taxa,
                    "n_spr_trees": len(spr_trees),
                    "trees": trees_payload,
                },
                sort_keys=False,
            )
        elif self.output_path:
            with open(self.output_path, "w") as f:
                for desc, spr_tree in spr_trees:
                    Phylo.write(spr_tree, f, "newick")
            print(f"Generated {len(spr_trees)} SPR trees")
            print(f"Subtree: ({', '.join(self.subtree_taxa)})")
            print(f"Regraft positions: {len(spr_trees)}")
            print(f"Output: {self.output_path}")
        else:
            for desc, spr_tree in spr_trees:
                s = StringIO()
                Phylo.write(spr_tree, s, "newick")
                try:
                    print(s.getvalue().strip())
                except BrokenPipeError:
                    return

    @staticmethod
    def _generate_spr_trees(
        tree, subtree_clade
    ) -> List[Tuple[str, object]]:
        """Generate all SPR rearrangements for a given subtree.

        Returns list of (regraft_description, tree) tuples.
        """
        results = []

        # Step 1: Build parent map
        parent_map = {}
        for clade in tree.find_clades(order="preorder"):
            for child in clade.clades:
                parent_map[id(child)] = clade

        # Step 2: Identify the subtree and its parent
        subtree_parent = parent_map.get(id(subtree_clade))
        if subtree_parent is None:
            return []  # Can't prune the root

        # Step 3: Get all node ids within the subtree
        subtree_node_ids = set()
        for node in subtree_clade.find_clades(order="preorder"):
            subtree_node_ids.add(id(node))

        # Collect regraft targets: all branches (parent->child) not in the subtree
        regraft_targets = []
        for clade in tree.find_clades(order="preorder"):
            if clade == tree.root:
                continue
            if id(clade) in subtree_node_ids:
                continue
            if id(clade) == id(subtree_clade):
                continue
            parent = parent_map.get(id(clade))
            if parent is not None and id(parent) not in subtree_node_ids:
                regraft_targets.append(clade)

        # Step 4: For each regraft target, create a new tree
        for target in regraft_targets:
            new_tree = copy.deepcopy(tree)

            # Build new parent map
            new_parent_map = {}
            for c in new_tree.find_clades(order="preorder"):
                for child in c.clades:
                    new_parent_map[id(child)] = c

            # Find subtree in new tree by matching taxa
            subtree_taxa = frozenset(
                t.name for t in subtree_clade.get_terminals()
            )
            new_subtree = None

            if subtree_clade.is_terminal():
                for c in new_tree.get_terminals():
                    if c.name in subtree_taxa:
                        new_subtree = c
                        break
            else:
                for c in new_tree.find_clades(order="preorder"):
                    if not c.is_terminal():
                        tips = frozenset(t.name for t in c.get_terminals())
                        if tips == subtree_taxa:
                            new_subtree = c
                            break

            if new_subtree is None:
                continue

            # Find target in new tree by matching taxa
            if target.is_terminal():
                target_taxa = frozenset([target.name])
            else:
                target_taxa = frozenset(
                    t.name for t in target.get_terminals()
                )
            new_target = None
            for c in new_tree.find_clades(order="preorder"):
                if c.is_terminal():
                    if frozenset([c.name]) == target_taxa:
                        new_target = c
                        break
                else:
                    tips = frozenset(t.name for t in c.get_terminals())
                    if tips == target_taxa:
                        new_target = c
                        break

            if new_target is None:
                continue

            # Rebuild parent map for new tree after potential changes
            new_parent_map = {}
            for c in new_tree.find_clades(order="preorder"):
                for child in c.clades:
                    new_parent_map[id(child)] = c

            # PRUNE: remove subtree from its parent
            subtree_parent_new = new_parent_map.get(id(new_subtree))
            if subtree_parent_new is None:
                continue

            subtree_bl = new_subtree.branch_length or 0.0
            subtree_parent_new.clades.remove(new_subtree)

            # If parent now has only 1 child, collapse it
            if len(subtree_parent_new.clades) == 1:
                only_child = subtree_parent_new.clades[0]
                grandparent = new_parent_map.get(id(subtree_parent_new))
                if grandparent is not None:
                    only_child.branch_length = (
                        only_child.branch_length or 0.0
                    ) + (subtree_parent_new.branch_length or 0.0)
                    idx = grandparent.clades.index(subtree_parent_new)
                    grandparent.clades[idx] = only_child
                else:
                    # Parent is root - make only_child the new root
                    new_tree.root = only_child
                    only_child.branch_length = None

            # Rebuild parent map after pruning
            new_parent_map = {}
            for c in new_tree.find_clades(order="preorder"):
                for child in c.clades:
                    new_parent_map[id(child)] = c

            # REGRAFT: insert subtree onto the target branch
            target_parent = new_parent_map.get(id(new_target))
            if target_parent is None:
                continue

            target_bl = new_target.branch_length or 0.0

            # Create new internal node
            new_node = Clade(branch_length=target_bl / 2)
            new_target.branch_length = target_bl / 2
            new_subtree.branch_length = subtree_bl

            # new_node gets the target and the subtree as children
            new_node.clades = [new_target, new_subtree]

            # Replace target with new_node in target's parent
            idx = target_parent.clades.index(new_target)
            target_parent.clades[idx] = new_node

            # Generate description
            if new_target.is_terminal():
                desc = (
                    f"regraft onto branch leading to {new_target.name}"
                )
            else:
                tip_names = sorted(
                    t.name for t in new_target.get_terminals()
                )
                if len(tip_names) <= 3:
                    desc = (
                        f"regraft onto branch leading to "
                        f"({', '.join(tip_names)})"
                    )
                else:
                    desc = (
                        f"regraft onto branch leading to "
                        f"({tip_names[0]}, ..., {tip_names[-1]}; "
                        f"{len(tip_names)} taxa)"
                    )

            results.append((desc, new_tree))

        return results
