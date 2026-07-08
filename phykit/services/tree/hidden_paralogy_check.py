from __future__ import annotations

import sys

from .base import Tree


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


class _LazyPhylo:
    _module = None

    def _load(self):
        module = self._module
        if module is None:
            from Bio import Phylo as _Phylo

            module = _Phylo
            self._module = module
        return module

    def read(self, *args, **kwargs):
        read = self._load().read
        self.read = read
        return read(*args, **kwargs)


Phylo = _LazyPhylo()


class _LazyMultiprocessing:
    _module = None

    def _load(self):
        module = self._module
        if module is None:
            import multiprocessing as module

            self._module = module
        return module

    def cpu_count(self):
        return self._load().cpu_count()

    def Pool(self, *args, **kwargs):
        return self._load().Pool(*args, **kwargs)


mp = _LazyMultiprocessing()

class HiddenParalogyCheck(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"], clade=parsed["clade"])
        self.json_output = parsed["json_output"]

    @staticmethod
    def _process_clade_batch(
        clade_batch, tree_file_path, master_tree_tips, exact_clades=None
    ):
        """Process a batch of clades in parallel."""
        batch_results = []
        if exact_clades is None:
            try:
                batch_tree = Phylo.read(tree_file_path, "newick")
                exact_clades = HiddenParalogyCheck._build_exact_clade_index(batch_tree)
            except (ValueError, AttributeError, FileNotFoundError):
                exact_clades = None

        for clade in clade_batch:
            clade_of_interest = master_tree_tips.intersection(clade)

            if len(clade_of_interest) <= 1:
                batch_results.append(["insufficient_taxon_representation"])
                continue

            if exact_clades is not None and frozenset(clade_of_interest) in exact_clades:
                batch_results.append(["monophyletic", []])
                continue

            # Read a fresh copy for each non-exact clade to preserve rooting behavior.
            tree = Phylo.read(tree_file_path, "newick")
            diff_tips = master_tree_tips - clade_of_interest

            # Root and find common ancestor
            try:
                tree.root_with_outgroup(list(diff_tips))
                subtree = tree.common_ancestor(clade_of_interest)

                common_ancestor_tips = (
                    HiddenParalogyCheck._terminal_names_direct(subtree)
                )
                if common_ancestor_tips is None:
                    common_ancestor_tips = set(
                        tip.name for tip in subtree.get_terminals()
                    )

                diff_tips_between_clade_and_curr_tree = \
                    clade_of_interest.symmetric_difference(common_ancestor_tips)

                batch_results.append([
                    "monophyletic" if not diff_tips_between_clade_and_curr_tree else "not_monophyletic",
                    list(diff_tips_between_clade_and_curr_tree),
                ])
            except (ValueError, AttributeError):
                # Handle edge cases where rooting fails
                batch_results.append(["processing_error"])

        return batch_results

    def run(self) -> None:
        # Read the master tree once to get all tip names
        master_tree = self.read_tree_file_unmodified()
        master_tree_tips = frozenset(self.get_tip_names_from_tree(master_tree))
        exact_clades = self._build_exact_clade_index(master_tree)

        # Read clades
        clades = self.read_clades_file(self.clade)

        # For small datasets, process sequentially
        if len(clades) < 10:
            res_arr = []
            for clade in clades:
                clade_of_interest = master_tree_tips.intersection(clade)

                if len(clade_of_interest) <= 1:
                    res_arr.append(["insufficient_taxon_representation"])
                    continue

                if exact_clades is not None and frozenset(clade_of_interest) in exact_clades:
                    res_arr.append(["monophyletic", []])
                    continue

                # Read a fresh tree for each non-exact clade to preserve rooting behavior.
                tree = Phylo.read(self.tree_file_path, "newick")
                diff_tips = master_tree_tips - clade_of_interest
                tree.root_with_outgroup(list(diff_tips))

                subtree = tree.common_ancestor(clade_of_interest)
                common_ancestor_tips = self._terminal_names_direct(subtree)
                if common_ancestor_tips is None:
                    common_ancestor_tips = set(self.get_tip_names_from_tree(subtree))

                diff_tips_between_clade_and_curr_tree = \
                    clade_of_interest.symmetric_difference(common_ancestor_tips)

                res_arr.append([
                    "monophyletic" if not diff_tips_between_clade_and_curr_tree else "not_monophyletic",
                    list(diff_tips_between_clade_and_curr_tree),
                ])
        else:
            # Use multiprocessing for larger datasets
            from functools import partial

            num_workers = min(mp.cpu_count(), 8)
            batch_size = max(1, len(clades) // num_workers)

            # Create clade batches
            clade_batches = [clades[i:i + batch_size]
                           for i in range(0, len(clades), batch_size)]

            # Process batches in parallel
            process_func = partial(
                self._process_clade_batch,
                tree_file_path=self.tree_file_path,
                master_tree_tips=master_tree_tips,
                exact_clades=exact_clades,
            )

            with mp.Pool(processes=num_workers) as pool:
                batch_results = pool.map(process_func, clade_batches)

            # Flatten results
            res_arr = []
            for batch_result in batch_results:
                res_arr.extend(batch_result)

        self.print_results(res_arr)

    @staticmethod
    def _terminal_names_direct(clade):
        try:
            clade.clades
        except AttributeError:
            return None

        names = set()
        stack = [clade]
        try:
            add_name = names.add
            pop_stack = stack.pop
            extend_stack = stack.extend
            while stack:
                node = pop_stack()
                children = node.clades
                if children:
                    extend_stack(children)
                else:
                    add_name(node.name)
        except (AttributeError, TypeError):
            return None
        return names

    @staticmethod
    def _build_exact_clade_index(tree):
        try:
            clade_tips = {}
            exact_clades = set()
            postorder = HiddenParalogyCheck._postorder_clades_direct(tree)
            if postorder is None:
                postorder = tree.find_clades(order="postorder")

            for clade in postorder:
                if clade.is_terminal():
                    tips = frozenset([clade.name])
                else:
                    children = clade.clades
                    if len(children) == 2:
                        tips = (
                            clade_tips[id(children[0])]
                            | clade_tips[id(children[1])]
                        )
                    else:
                        tips = frozenset().union(
                            *(clade_tips[id(child)] for child in children)
                        )
                clade_tips[id(clade)] = tips
                exact_clades.add(tips)
            return exact_clades
        except (AttributeError, KeyError, TypeError):
            return None

    @staticmethod
    def _postorder_clades_direct(tree):
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
                stack.extend(clade.clades)
        except AttributeError:
            return None
        clades.reverse()
        return clades

    def process_args(self, args) -> dict[str, str]:
        return dict(
            tree_file_path=args.tree,
            clade=args.clade,
            json_output=getattr(args, "json", False),
        )

    def read_clades_file(self, clades: str) -> list[list[str]]:
        try:
            with open(clades) as file:
                return [line.split() for line in file]
        except FileNotFoundError:
            print("Clade file not found. Please check the path.")
            sys.exit(2)

    def print_results(self, res_arr: list[list[list | str]]) -> None:
        if self.json_output:
            rows = []
            append_row = rows.append
            for idx, res in enumerate(res_arr, start=1):
                unexpected = []
                if len(res) > 1 and isinstance(res[1], list):
                    unexpected = sorted(res[1])
                append_row(
                    {
                        "clade_index": idx,
                        "status": res[0],
                        "unexpected_taxa": unexpected,
                    }
                )
            print_json(dict(rows=rows, clades=rows))
            return

        lines = [res[0] for res in res_arr]
        if lines:
            print("\n".join(lines))
