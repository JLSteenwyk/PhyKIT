import re
from typing import Dict, List, Optional, Tuple
import pickle

from Bio import Phylo
from Bio.Phylo import Newick

from .base import Tree
from ...errors import PhykitUserError
from ...helpers.json_output import print_json


BranchSpec = Tuple[str, List[str]]


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

    def process_args(self, args) -> Dict[str, str]:
        tree_file_path = args.tree
        output_file_path = \
            f"{args.output}" if args.output else f"{tree_file_path}.nnis"

        branch = None
        raw_branch = getattr(args, "branch", None)
        if raw_branch:
            parts = [p.strip() for p in re.split(r"[,\t]", raw_branch) if p.strip()]
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

    def _resolve_branch_specs(self) -> List[BranchSpec]:
        specs: List[BranchSpec] = []
        if self.branch:
            t1, t2 = self.branch
            specs.append((f"{t1}|{t2}", [t1, t2]))
        if self.branches_file:
            with open(self.branches_file) as fh:
                for i, raw in enumerate(fh, start=1):
                    line = raw.strip()
                    if not line or line.startswith("#"):
                        continue
                    parts = [p.strip() for p in re.split(r"[,\t]", line) if p.strip()]
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
        branch_specs: List[BranchSpec],
    ) -> Tuple[List[Newick.Tree], List[Dict]]:
        tip_names = {term.name for term in tree.get_terminals()}
        output_trees: List[Newick.Tree] = []
        report: List[Dict] = []

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

    def _build_parent_map(self, tree: Newick.Tree) -> Dict:
        parents = {}
        for clade in tree.find_clades():
            if clade != tree.root:
                node_path = tree.get_path(clade)
                if len(node_path) == 1:
                    parents[clade] = tree.root
                else:
                    parents[clade] = node_path[-2]
        return parents

    def _nnis_around_branch(
        self,
        tree: Newick.Tree,
        clade,
        parents: Dict,
    ) -> List[Newick.Tree]:
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

    def _nnis_root_edge(self, tree: Newick.Tree) -> List[Newick.Tree]:
        neighbors: List[Newick.Tree] = []
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
    ) -> List[Newick.Tree]:
        neighbors: List[Newick.Tree] = []
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
    ) -> List[Newick.Tree]:
        ### This code is from BioPython (so is this comment)
        # Get all neighbor trees of the given tree (PRIVATE).
        # Currently only for binary rooted trees.
        ###
        parents = self._build_parent_map(tree)
        neighbors: List[Newick.Tree] = []
        root_childs = []
        for clade in tree.get_nonterminals(order="level"):
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
