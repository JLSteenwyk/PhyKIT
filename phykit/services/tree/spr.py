from __future__ import annotations

import os
from io import StringIO

from .base import Tree
from ...errors import PhykitUserError


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


class _LazyPickle:
    def dumps(self, *args, **kwargs):
        import pickle as _pickle

        return _pickle.dumps(*args, **kwargs)

    def loads(self, *args, **kwargs):
        import pickle as _pickle

        return _pickle.loads(*args, **kwargs)

    def __getattr__(self, name):
        import pickle as _pickle

        return getattr(_pickle, name)


pickle = _LazyPickle()


class _LazyPhylo:
    def write(self, *args, **kwargs):
        from Bio import Phylo as _Phylo

        self.write = _Phylo.write
        return self.write(*args, **kwargs)


class _LazyClade:
    def __call__(self, *args, **kwargs):
        from Bio.Phylo.Newick import Clade as _Clade

        return _Clade(*args, **kwargs)


Phylo = _LazyPhylo()
Clade = _LazyClade()


class Spr(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"])
        self.subtree_taxa = parsed["subtree_taxa"]
        self.output_path = parsed["output_path"]
        self.json_output = parsed["json_output"]

    def process_args(self, args) -> dict[str, object]:
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
        terminal_by_name = self._terminal_by_name_fast(tree)
        if terminal_by_name is not None:
            tree_tips = set(terminal_by_name)
        else:
            terminals = tree.get_terminals()
            terminal_by_name = {tip.name: tip for tip in terminals}
            tree_tips = set(terminal_by_name)
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
            subtree = terminal_by_name.get(self.subtree_taxa[0])
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
            self._print_output_summary(len(spr_trees))
        else:
            self._print_spr_trees(spr_trees)

    def _print_output_summary(self, n_spr_trees: int) -> None:
        print(
            "\n".join(
                [
                    f"Generated {n_spr_trees} SPR trees",
                    f"Subtree: ({', '.join(self.subtree_taxa)})",
                    f"Regraft positions: {n_spr_trees}",
                    f"Output: {self.output_path}",
                ]
            )
        )

    @staticmethod
    def _print_spr_trees(spr_trees) -> None:
        lines = []
        for desc, spr_tree in spr_trees:
            s = StringIO()
            Phylo.write(spr_tree, s, "newick")
            lines.append(s.getvalue().strip())
        if not lines:
            return
        try:
            print("\n".join(lines))
        except BrokenPipeError:
            return

    @staticmethod
    def _collect_clade_taxa(tree) -> dict[int, frozenset]:
        clade_taxa: dict[int, frozenset] = {}
        for clade in Spr._iter_postorder(tree.root):
            children = clade.clades
            if not children:
                clade_taxa[id(clade)] = frozenset((clade.name,))
            else:
                child_count = len(children)
                if child_count == 2:
                    clade_taxa[id(clade)] = (
                        clade_taxa[id(children[0])]
                        | clade_taxa[id(children[1])]
                    )
                else:
                    taxa = set()
                    for child in children:
                        taxa.update(clade_taxa[id(child)])
                    clade_taxa[id(clade)] = frozenset(taxa)
        return clade_taxa

    @staticmethod
    def _iter_preorder(root):
        stack = [root]
        pop = stack.pop
        append = stack.append
        while stack:
            clade = pop()
            yield clade
            children = clade.clades
            if children:
                child_count = len(children)
                if child_count == 2:
                    append(children[1])
                    append(children[0])
                else:
                    for index in range(child_count - 1, -1, -1):
                        append(children[index])

    @staticmethod
    def _iter_postorder(root):
        clades = []
        stack = [root]
        pop = stack.pop
        extend = stack.extend
        append_clade = clades.append
        while stack:
            clade = pop()
            append_clade(clade)
            children = clade.clades
            if children:
                extend(children)
        yield from reversed(clades)

    @staticmethod
    def _clades_and_parent_map(root):
        clades = []
        parent_map = {}
        stack = [root]
        pop = stack.pop
        append = stack.append
        append_clade = clades.append
        while stack:
            clade = pop()
            append_clade(clade)
            children = clade.clades
            if children:
                if len(children) == 2:
                    left, right = children
                    parent_map[id(left)] = clade
                    parent_map[id(right)] = clade
                    append(right)
                    append(left)
                else:
                    for index in range(len(children) - 1, -1, -1):
                        child = children[index]
                        parent_map[id(child)] = clade
                        append(child)
        return clades, parent_map

    @staticmethod
    def _describe_regraft_target(target, target_taxa) -> str:
        if target.is_terminal():
            return f"regraft onto branch leading to {target.name}"

        tip_names = sorted(target_taxa)
        if len(tip_names) <= 3:
            return (
                f"regraft onto branch leading to "
                f"({', '.join(tip_names)})"
            )
        return (
            f"regraft onto branch leading to "
            f"({tip_names[0]}, ..., {tip_names[-1]}; "
            f"{len(tip_names)} taxa)"
        )

    @staticmethod
    def _generate_spr_trees(
        tree, subtree_clade
    ) -> list[tuple[str, object]]:
        """Generate all SPR rearrangements for a given subtree.

        Returns list of (regraft_description, tree) tuples.
        """
        results = []

        # Step 1: Build parent map
        original_clades, parent_map = Spr._clades_and_parent_map(tree.root)

        # Step 2: Identify the subtree and its parent
        subtree_parent = parent_map.get(id(subtree_clade))
        if subtree_parent is None:
            return []  # Can't prune the root

        # Step 3: Get all node ids within the subtree
        subtree_node_ids = set()
        for node in Spr._iter_preorder(subtree_clade):
            subtree_node_ids.add(id(node))

        clade_index = {id(clade): idx for idx, clade in enumerate(original_clades)}
        subtree_index = clade_index[id(subtree_clade)]
        clade_taxa = Spr._collect_clade_taxa(tree)

        # Collect regraft targets: all branches (parent->child) not in the subtree
        regraft_targets = []
        target_descriptions = {}
        for clade in original_clades:
            if clade == tree.root:
                continue
            if id(clade) in subtree_node_ids:
                continue
            if id(clade) == id(subtree_clade):
                continue
            parent = parent_map.get(id(clade))
            if parent is not None and id(parent) not in subtree_node_ids:
                regraft_targets.append(clade)
                target_descriptions[id(clade)] = Spr._describe_regraft_target(
                    clade,
                    clade_taxa.get(id(clade), frozenset()),
                )

        # Step 4: For each regraft target, create a new tree
        for target in regraft_targets:
            new_tree = pickle.loads(pickle.dumps(tree, protocol=pickle.HIGHEST_PROTOCOL))
            copy_clades, new_parent_map = Spr._clades_and_parent_map(new_tree.root)
            new_subtree = copy_clades[subtree_index]
            new_target = copy_clades[clade_index[id(target)]]

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
            _, new_parent_map = Spr._clades_and_parent_map(new_tree.root)

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

            desc = target_descriptions[id(target)]

            results.append((desc, new_tree))

        return results
