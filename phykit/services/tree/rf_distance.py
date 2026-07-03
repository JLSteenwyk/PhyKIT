from __future__ import annotations

from .base import Tree


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


def ProcessPoolExecutor(*args, **kwargs):
    from concurrent.futures import ProcessPoolExecutor as _ProcessPoolExecutor

    return _ProcessPoolExecutor(*args, **kwargs)


pickle = _LazyPickle()


class RobinsonFouldsDistance(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(
            tree_file_path=parsed["tree_file_path"],
            tree1_file_path=parsed["tree1_file_path"],
        )
        self.json_output = parsed["json_output"]

    def run(self):
        if self.tree_file_path == self.tree1_file_path:
            tree_zero = self.read_tree_file_unmodified()
            tip_count = self._terminal_count_direct(tree_zero)
            if tip_count is None:
                tip_count = tree_zero.count_terminals()
            self._output_result(0, 0 / (2 * (tip_count - 3)))
            return

        tree_zero = self.read_tree_file()
        tree_one = self.read_tree1_file()

        # get shared tree tip names - use sets for efficiency
        tree_zero_tips = set(self.get_tip_names_from_tree(tree_zero))
        tree_one_tips = set(self.get_tip_names_from_tree(tree_one))
        shared_tree_tips = tree_zero_tips & tree_one_tips

        # prune to common set - already have sets
        tree_zero_tips_to_prune = list(tree_zero_tips - shared_tree_tips)
        tree_one_tips_to_prune = list(tree_one_tips - shared_tree_tips)

        if tree_zero_tips_to_prune:
            tree_zero = self.prune_tree_using_taxa_list(tree_zero, tree_zero_tips_to_prune)
        if tree_one_tips_to_prune:
            tree_one = self.prune_tree_using_taxa_list(tree_one, tree_one_tips_to_prune)

        # Get first terminal for rooting
        tip_for_rooting = self._first_terminal_name(tree_zero)
        tree_zero.root_with_outgroup(tip_for_rooting)
        tree_one.root_with_outgroup(tip_for_rooting)

        plain_rf, normalized_rf = self.calculate_robinson_foulds_distance(
            tree_zero, tree_one
        )

        self._output_result(plain_rf, normalized_rf)

    def _output_result(self, plain_rf, normalized_rf):
        if self.json_output:
            print_json(
                dict(
                    plain_rf=plain_rf,
                    normalized_rf=round(normalized_rf, 4),
                )
            )
            return

        print(f"{plain_rf}\t{round(normalized_rf, 4)}")

    def process_args(self, args) -> dict[str, str]:
        return dict(
            tree_file_path=args.tree_zero,
            tree1_file_path=args.tree_one,
            json_output=getattr(args, "json", False),
        )

    def calculate_robinson_foulds_distance(self, tree_zero, tree_one):
        if tree_zero is tree_one:
            tip_count = self._terminal_count_direct(tree_zero)
            if tip_count is None:
                tip_count = tree_zero.count_terminals()
            return 0, 0 / (2 * (tip_count - 3))

        try:
            id_result_zero = self._get_all_bipartition_id_sets_direct(tree_zero)
            if id_result_zero is None:
                bipartitions_zero = self.get_all_bipartitions(tree_zero)
                bipartitions_one = self.get_all_bipartitions(tree_one)
                plain_rf = len(bipartitions_zero ^ bipartitions_one)
            else:
                bipartitions_zero, tip_index = id_result_zero
                id_result_one = self._get_all_bipartition_id_sets_direct(
                    tree_one,
                    tip_index,
                )
                if id_result_one is None:
                    bipartitions_zero = self.get_all_bipartitions(tree_zero)
                    bipartitions_one = self.get_all_bipartitions(tree_one)
                    plain_rf = len(bipartitions_zero ^ bipartitions_one)
                else:
                    bipartitions_one, _ = id_result_one
                    plain_rf = len(bipartitions_zero ^ bipartitions_one)
        except (AttributeError, TypeError):
            plain_rf = 0
            plain_rf = self.compare_trees_optimized(plain_rf, tree_zero, tree_one)
            plain_rf = self.compare_trees_optimized(plain_rf, tree_one, tree_zero)

        tip_count = self._terminal_count_direct(tree_zero)
        if tip_count is None:
            tip_count = tree_zero.count_terminals()
        normalized_rf = plain_rf / (2 * (tip_count - 3))

        return plain_rf, normalized_rf

    @staticmethod
    def _terminal_count_direct(tree: Newick.Tree):
        try:
            root = tree.root
            root.clades
        except AttributeError:
            return None

        count = 0
        stack = [root]
        try:
            pop = stack.pop
            extend = stack.extend
            while stack:
                clade = pop()
                children = clade.clades
                if children:
                    extend(children)
                else:
                    count += 1
        except AttributeError:
            return None
        return count

    def get_all_bipartitions(self, tree: Newick.Tree) -> set[frozenset]:
        """Return rooted non-root internal descendant tip sets.

        This intentionally matches the historical RF calculation, which compares
        descendant taxa below each non-root internal clade rather than canonical
        unrooted splits.
        """
        bipartitions = self._get_all_bipartitions_direct(tree)
        if bipartitions is not None:
            return bipartitions

        clade_tips: dict[int, frozenset] = {}
        bipartitions: set[frozenset] = set()

        for clade in tree.find_clades(order="postorder"):
            if clade.is_terminal():
                tip_set = frozenset([clade.name])
            else:
                child_sets = [clade_tips[id(child)] for child in clade.clades]
                tip_set = frozenset().union(*child_sets) if child_sets else frozenset()
                if clade is not tree.root:
                    bipartitions.add(tip_set)
            clade_tips[id(clade)] = tip_set

        return bipartitions

    @staticmethod
    def _get_all_bipartitions_direct(tree: Newick.Tree):
        try:
            root = tree.root
            root.clades
        except AttributeError:
            return None

        clade_tips: dict[int, frozenset] = {}
        bipartitions: set[frozenset] = set()
        preorder = []
        stack = [root]
        append_preorder = preorder.append
        append_stack = stack.append
        pop_stack = stack.pop

        try:
            while stack:
                clade = pop_stack()
                append_preorder(clade)
                children = clade.clades
                child_count = len(children)
                if child_count == 2:
                    append_stack(children[0])
                    append_stack(children[1])
                elif child_count:
                    for child in children:
                        append_stack(child)

            for clade in reversed(preorder):
                children = clade.clades
                if children:
                    child_count = len(children)
                    if child_count == 2:
                        tip_set = (
                            clade_tips[id(children[0])]
                            | clade_tips[id(children[1])]
                        )
                    elif child_count == 1:
                        tip_set = clade_tips[id(children[0])]
                    else:
                        tip_set = frozenset().union(
                            *(clade_tips[id(child)] for child in children)
                        )
                    if clade is not root:
                        bipartitions.add(tip_set)
                else:
                    tip_set = frozenset([clade.name])
                clade_tips[id(clade)] = tip_set
        except (AttributeError, TypeError):
            return None

        return bipartitions

    @staticmethod
    def _get_all_bipartition_id_sets_direct(tree: Newick.Tree, tip_index=None):
        try:
            root = tree.root
            root.clades
        except AttributeError:
            return None

        if tip_index is None:
            tip_index = {}
        clade_tips: dict[int, frozenset] = {}
        bipartitions: set[frozenset] = set()
        preorder = []
        stack = [root]

        try:
            while stack:
                clade = stack.pop()
                preorder.append(clade)
                stack.extend(clade.clades)

            for clade in reversed(preorder):
                children = clade.clades
                if children:
                    child_count = len(children)
                    if child_count == 2:
                        tip_set = (
                            clade_tips[id(children[0])]
                            | clade_tips[id(children[1])]
                        )
                    elif child_count == 1:
                        tip_set = clade_tips[id(children[0])]
                    else:
                        tip_set = frozenset().union(
                            *(clade_tips[id(child)] for child in children)
                        )
                    if clade is not root:
                        bipartitions.add(tip_set)
                else:
                    tip_id = tip_index.get(clade.name)
                    if tip_id is None:
                        tip_id = len(tip_index)
                        tip_index[clade.name] = tip_id
                    tip_set = frozenset([tip_id])
                clade_tips[id(clade)] = tip_set
        except (AttributeError, TypeError):
            return None

        return bipartitions, tip_index

    @staticmethod
    def _first_terminal_name(tree: Newick.Tree) -> str:
        try:
            root = tree.root
            root.clades
        except AttributeError:
            return tree.get_terminals()[0].name

        clade = root
        while True:
            children = getattr(clade, "clades", None)
            if not isinstance(children, list):
                return tree.get_terminals()[0].name
            if not children:
                return clade.name
            clade = children[0]

    def compare_trees_optimized(
        self,
        plain_rf: int,
        tree_zero: Newick.Tree,
        tree_one: Newick.Tree
    ) -> int:
        # Cache tip names for clades to avoid recomputation
        tip_names_cache = {}

        def get_cached_tips(clade):
            clade_id = id(clade)
            if clade_id not in tip_names_cache:
                tip_names_cache[clade_id] = frozenset(self.get_tip_names_from_tree(clade))
            return tip_names_cache[clade_id]

        # loop through tree_zero and find similar clade in tree_one
        for clade_zero in tree_zero.get_nonterminals()[1:]:
            # Get tip names from tree_zero clade
            tip_names_zero = get_cached_tips(clade_zero)
            # get common ancestor of tree_zero tip names in tree_one
            clade_one = tree_one.common_ancestor(list(tip_names_zero))
            # Get tip names from tree_one clade
            tip_names_one = get_cached_tips(clade_one)
            # compare the list of tip names
            if tip_names_zero != tip_names_one:
                plain_rf += 1

        return plain_rf

    @staticmethod
    def _calculate_rf_batch(tree_pairs_pickle):
        """Calculate RF distance for a batch of tree pairs in parallel."""
        tree_pairs = pickle.loads(tree_pairs_pickle)
        results = []

        for tree_zero, tree_one in tree_pairs:
            rf_calc = RobinsonFouldsDistance.__new__(RobinsonFouldsDistance)
            rf_calc.__dict__.update({'tree_format': 'newick'})

            plain_rf, normalized_rf = rf_calc.calculate_robinson_foulds_distance(
                tree_zero,
                tree_one,
            )
            results.append((plain_rf, normalized_rf))

        return results

    def calculate_multiple_rf_distances(
        self,
        tree_pairs: list[tuple],
    ) -> list[tuple[int, float]]:
        """Calculate RF distances for multiple tree pairs in parallel."""
        if len(tree_pairs) < 5:
            # Sequential for small datasets
            results = []
            for tree_zero, tree_one in tree_pairs:
                plain_rf, normalized_rf = self.calculate_robinson_foulds_distance(tree_zero, tree_one)
                results.append((plain_rf, normalized_rf))
            return results

        # Parallel processing for larger datasets
        batch_size = max(2, len(tree_pairs) // 4)
        batches = [tree_pairs[i:i + batch_size] for i in range(0, len(tree_pairs), batch_size)]

        with ProcessPoolExecutor(max_workers=min(4, len(batches))) as executor:
            futures = []
            for batch in batches:
                batch_pickle = pickle.dumps(batch)
                futures.append(executor.submit(self._calculate_rf_batch, batch_pickle))

            all_results = []
            for future in futures:
                all_results.extend(future.result())

        return all_results
