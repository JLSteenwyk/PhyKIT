from __future__ import annotations

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


class _LazyPhylo:
    def write(self, *args, **kwargs):
        from Bio import Phylo as _Phylo

        self.write = _Phylo.write
        return self.write(*args, **kwargs)


Phylo = _LazyPhylo()
pickle = _LazyPickle()


BranchSpec = tuple[str, list[str]]


def _split_branch_fields(value: str) -> list[str]:
    return [
        part.strip()
        for part in value.replace("\t", ",").split(",")
        if part.strip()
    ]


class NearestNeighborInterchange(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(
            tree_file_path=parsed["tree_file_path"],
            output_file_path=parsed["output_file_path"],
        )
        self.json_output = parsed["json_output"]
        self.branch = parsed["branch"]
        self.branches_file = parsed["branches_file"]
        self.no_input_tree = parsed["no_input_tree"]

    def run(self) -> None:
        tree = self.read_tree_file()

        branch_specs = self._resolve_branch_specs()

        if branch_specs:
            output_trees, branch_report = self._generate_targeted_nnis(
                tree, branch_specs
            )
        else:
            output_trees = self.get_neighbors(tree)
            branch_report = None

        trees_to_write = output_trees if self.no_input_tree else [tree, *output_trees]
        Phylo.write(trees_to_write, self.output_file_path, "newick")

        if self.json_output:
            payload = dict(
                input_tree=self.tree_file_path,
                total_trees=len(trees_to_write),
                nni_neighbors=len(output_trees),
                output_file=self.output_file_path,
            )
            if branch_report is not None:
                payload["branches"] = branch_report
            print_json(payload)

    def process_args(self, args) -> dict[str, str]:
        tree_file_path = args.tree
        output_file_path = \
            f"{args.output}" if args.output else f"{tree_file_path}.nnis"

        branch = None
        raw_branch = getattr(args, "branch", None)
        if raw_branch:
            parts = _split_branch_fields(raw_branch)
            if len(parts) != 2:
                raise PhykitUserError(
                    [
                        "--branch must specify exactly two taxa separated by a "
                        "comma (e.g., --branch A,B)"
                    ],
                    code=2,
                )
            branch = tuple(parts)

        return dict(
            tree_file_path=tree_file_path,
            output_file_path=output_file_path,
            branch=branch,
            branches_file=getattr(args, "branches", None),
            no_input_tree=getattr(args, "no_input_tree", False),
            json_output=getattr(args, "json", False),
        )

    def _resolve_branch_specs(self) -> list[BranchSpec]:
        specs: list[BranchSpec] = []
        if self.branch:
            t1, t2 = self.branch
            specs.append((f"{t1}|{t2}", [t1, t2]))
        if self.branches_file:
            with open(self.branches_file) as fh:
                for i, raw in enumerate(fh, start=1):
                    line = raw.strip()
                    if not line or line.startswith("#"):
                        continue
                    parts = _split_branch_fields(line)
                    if len(parts) < 2:
                        raise PhykitUserError(
                            [
                                f"--branches file line {i}: need at least two "
                                "taxa per line (comma- or tab-separated)"
                            ],
                            code=2,
                        )
                    if len(parts) >= 3:
                        label, t1, t2 = parts[0], parts[1], parts[2]
                    else:
                        t1, t2 = parts[0], parts[1]
                        label = f"{t1}|{t2}"
                    specs.append((label, [t1, t2]))
        return specs

    def _generate_targeted_nnis(
        self,
        tree: Newick.Tree,
        branch_specs: list[BranchSpec],
    ) -> tuple[list[Newick.Tree], list[dict]]:
        tip_names = self.calculate_terminal_names_fast(tree)
        if tip_names is None:
            tip_names = [term.name for term in tree.get_terminals()]
        tip_names = set(tip_names)
        output_trees: list[Newick.Tree] = []
        report: list[dict] = []

        for label, taxa in branch_specs:
            missing = [t for t in taxa if t not in tip_names]
            if missing:
                raise PhykitUserError(
                    [
                        f"branch '{label}': taxa not found in tree: "
                        + ", ".join(missing)
                    ],
                    code=2,
                )

            working = self._fast_tree_copy(tree)
            parents = self._build_parent_map(working)
            target = working.common_ancestor(taxa)

            if target is working.root:
                raise PhykitUserError(
                    [
                        f"branch '{label}': MRCA of {taxa[0]} and {taxa[1]} is "
                        "the root; no internal branch above to rearrange"
                    ],
                    code=2,
                )
            if target.is_terminal():
                raise PhykitUserError(
                    [
                        f"branch '{label}': MRCA of {taxa[0]} and {taxa[1]} is "
                        "a terminal; NNI requires an internal branch"
                    ],
                    code=2,
                )

            nnis = self._nnis_around_branch(working, target, parents)
            output_trees.extend(nnis)
            report.append(
                dict(label=label, taxa=list(taxa), n_nnis=len(nnis))
            )

        return output_trees, report

    def _fast_tree_copy(self, tree: Newick.Tree) -> Newick.Tree:
        """Fast tree copying using pickle instead of deep copy."""
        return pickle.loads(pickle.dumps(tree, protocol=pickle.HIGHEST_PROTOCOL))

    def _build_parent_map(self, tree: Newick.Tree) -> dict:
        parents = self._build_parent_map_direct(tree)
        if parents is not None:
            return parents

        parents = {}
        for clade in tree.find_clades():
            for child in clade.clades:
                parents[child] = clade
        return parents

    @staticmethod
    def _build_parent_map_direct(tree: Newick.Tree):
        try:
            root = tree.root
            root.clades
        except AttributeError:
            return None

        parents = {}
        stack = [root]
        try:
            pop = stack.pop
            extend = stack.extend
            while stack:
                clade = pop()
                children = clade.clades
                for child in children:
                    parents[child] = clade
                if children:
                    extend(children)
        except AttributeError:
            return None
        return parents

    @staticmethod
    def _nonterminals_level_order_direct(tree: Newick.Tree):
        try:
            root = tree.root
            root.clades
        except AttributeError:
            return None

        nonterminals = []
        queue = [root]
        index = 0
        append = nonterminals.append
        extend = queue.extend
        try:
            while index < len(queue):
                clade = queue[index]
                index += 1
                children = clade.clades
                if children:
                    append(clade)
                    extend(children)
        except AttributeError:
            return None
        return nonterminals

    def _nnis_around_branch(
        self,
        tree: Newick.Tree,
        clade,
        parents: dict,
    ) -> list[Newick.Tree]:
        """Generate the 2 NNI rearrangements around the branch leading to clade.

        Mutates the tree in place but restores the original topology before
        returning so the caller can request additional rearrangements on the
        same tree.
        """
        if clade.is_terminal():
            return []

        root_children = tree.root.clades
        is_root_child = clade in root_children
        if is_root_child:
            other = (
                root_children[1] if clade is root_children[0] else root_children[0]
            )
            if other.is_terminal():
                # Implicit root edge can only be rearranged if both root
                # children are internal — otherwise there is no subtree on
                # the sister side to swap.
                return []
            return self._nnis_root_edge(tree)

        parent = parents.get(clade)
        if parent is None:
            return []
        return self._nnis_internal_edge(tree, clade, parent)

    def _nnis_root_edge(self, tree: Newick.Tree) -> list[Newick.Tree]:
        neighbors: list[Newick.Tree] = []
        left = tree.root.clades[0]
        right = tree.root.clades[1]
        left_right = left.clades[1]
        right_left = right.clades[0]
        right_right = right.clades[1]

        # neighbor 1: swap left.clades[1] with right.clades[1]
        del left.clades[1]
        del right.clades[1]
        left.clades.append(right_right)
        right.clades.append(left_right)
        neighbors.append(self._fast_tree_copy(tree))

        # neighbor 2: swap left.clades[1] with right.clades[0]
        del left.clades[1]
        del right.clades[0]
        left.clades.append(right_left)
        right.clades.append(right_right)
        neighbors.append(self._fast_tree_copy(tree))

        # restore original
        del left.clades[1]
        del right.clades[0]
        left.clades.append(left_right)
        right.clades.insert(0, right_left)
        return neighbors

    def _nnis_internal_edge(
        self,
        tree: Newick.Tree,
        clade,
        parent,
    ) -> list[Newick.Tree]:
        neighbors: list[Newick.Tree] = []
        left = clade.clades[0]
        right = clade.clades[1]

        if clade is parent.clades[0]:
            sister = parent.clades[1]
            # neighbor 1: swap clade.clades[1] with sister
            del parent.clades[1]
            del clade.clades[1]
            parent.clades.append(right)
            clade.clades.append(sister)
            neighbors.append(self._fast_tree_copy(tree))
            # neighbor 2: swap clade.clades[0] with sister (current parent.clades[1] is `right`)
            del parent.clades[1]
            del clade.clades[0]
            parent.clades.append(left)
            clade.clades.append(right)
            neighbors.append(self._fast_tree_copy(tree))
            # restore
            del parent.clades[1]
            del clade.clades[0]
            parent.clades.append(sister)
            clade.clades.insert(0, left)
        else:
            sister = parent.clades[0]
            # neighbor 1
            del parent.clades[0]
            del clade.clades[1]
            parent.clades.insert(0, right)
            clade.clades.append(sister)
            neighbors.append(self._fast_tree_copy(tree))
            # neighbor 2
            del parent.clades[0]
            del clade.clades[0]
            parent.clades.insert(0, left)
            clade.clades.append(right)
            neighbors.append(self._fast_tree_copy(tree))
            # restore
            del parent.clades[0]
            del clade.clades[0]
            parent.clades.insert(0, sister)
            clade.clades.insert(0, left)

        return neighbors

    def get_neighbors(
        self,
        tree: Newick.Tree
    ) -> list[Newick.Tree]:
        ### This code is from BioPython (so is this comment)
        # Get all neighbor trees of the given tree (PRIVATE).
        # Currently only for binary rooted trees.
        ###
        parents = self._build_parent_map(tree)
        neighbors: list[Newick.Tree] = []
        root_childs = []
        nonterminals = self._nonterminals_level_order_direct(tree)
        if nonterminals is None:
            nonterminals = tree.get_nonterminals(order="level")
        for clade in nonterminals:
            if clade == tree.root:
                left = clade.clades[0]
                right = clade.clades[1]
                root_childs.append(left)
                root_childs.append(right)
                if not left.is_terminal() and not right.is_terminal():
                    neighbors.extend(self._nnis_root_edge(tree))
            elif clade in root_childs:
                continue
            else:
                parent = parents[clade]
                neighbors.extend(self._nnis_internal_edge(tree, clade, parent))

        return neighbors
