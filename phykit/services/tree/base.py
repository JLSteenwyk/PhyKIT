from __future__ import annotations

from functools import lru_cache
import os

from ..base import BaseService
from ...errors import PhykitUserError


class _LazyPickle:
    _module = None

    def _load(self):
        module = self._module
        if module is None:
            import pickle as module

            self._module = module
        return module

    def dumps(self, *args, **kwargs):
        dumps = self._load().dumps
        self.dumps = dumps
        return dumps(*args, **kwargs)

    def loads(self, *args, **kwargs):
        loads = self._load().loads
        self.loads = loads
        return loads(*args, **kwargs)

    def __getattr__(self, name):
        attr = getattr(self._load(), name)
        setattr(self, name, attr)
        return attr


pickle = _LazyPickle()


class Tree(BaseService):
    _PAIRWISE_LCA_DEPTH_THRESHOLD = 64
    _PAIRWISE_PREALLOC_MIN_PAIRS = 100_000
    _PAIRWISE_PREALLOC_WITH_COMBOS_MIN_PAIRS = 500_000
    _ORDERED_MAPPING_PRUNE_MIN_SIZE = 200_000
    _ORDERED_NAMES_PRUNE_MIN_SIZE = 50_000

    def __init__(
        self,
        *args,
        tree_file_path=None,
        idmap=None,
        alignment_file_path=None,
        tree1_file_path=None,
        outgroup_taxa_file_path=None,
        output_file_path=None,
        factor=None,
        remove=None,
        verbose=None,
        reference=None,
        list_of_taxa=None,
        trees=None,
        groups=None,
        support=None,
        tip_1=None,
        tip_2=None,
        clade=None,
        keep=None,
        exclude_gaps=None,
    ):
        self.tree_file_path = tree_file_path
        self.tree1_file_path = tree1_file_path
        self.alignment_file_path = alignment_file_path
        self.output_file_path = output_file_path
        self.outgroup_taxa_file_path = outgroup_taxa_file_path
        self.tree_format = "newick"
        self.verbose = verbose
        self.factor = factor
        self.remove = remove
        self.idmap = idmap
        self.reference = reference
        self.list_of_taxa = list_of_taxa
        self.trees = trees
        self.groups = groups
        self.support = support
        self.tip_1 = tip_1
        self.tip_2 = tip_2
        self.clade = clade
        self.keep = keep
        self.exclude_gaps = exclude_gaps

    @staticmethod
    @lru_cache(maxsize=32)
    def _cached_tree_read(file_path: str, tree_format: str, file_hash: str):
        """Cached tree reading with file hash for cache invalidation."""
        from Bio import Phylo

        return Phylo.read(file_path, tree_format)

    @staticmethod
    @lru_cache(maxsize=32)
    def _cached_simple_newick_summary(file_path: str, file_hash: str):
        return Tree._scan_simple_newick_summary(file_path)

    @staticmethod
    @lru_cache(maxsize=32)
    def _cached_simple_newick_tip_names(file_path: str, file_hash: str):
        return Tree._scan_simple_newick_tip_names(file_path)

    @staticmethod
    @lru_cache(maxsize=32)
    def _cached_simple_newick_terminal_distance_stats(
        file_path: str,
        file_hash: str,
    ):
        return Tree._scan_simple_newick_terminal_distance_stats(file_path)

    @staticmethod
    def _get_file_hash(file_path: str) -> str:
        """Get a hash based on file path, size, and modification time."""
        try:
            stat = os.stat(file_path)
            return f"{file_path}_{stat.st_size}_{stat.st_mtime_ns}"
        except OSError:
            return ""

    def read_tree_file(self):
        return self._read_tree_with_error(self.tree_file_path, "tree_file_path")

    def read_tree_file_unmodified(self):
        return self._read_tree_with_error(
            self.tree_file_path,
            "tree_file_path",
            copy_tree=False,
        )

    def read_tree1_file(self):
        return self._read_tree_with_error(self.tree1_file_path, "tree1_file_path")

    def read_tree1_file_unmodified(self):
        return self._read_tree_with_error(
            self.tree1_file_path,
            "tree1_file_path",
            copy_tree=False,
        )

    def read_reference_tree_file(self):
        return self._read_tree_with_error(self.reference, "reference")

    def read_reference_tree_file_unmodified(self):
        return self._read_tree_with_error(
            self.reference,
            "reference",
            copy_tree=False,
        )

    @staticmethod
    def _fast_copy(tree):
        """Copy a tree using pickle instead of deepcopy.

        Avoids RecursionError on deeply nested trees (e.g., ladder-like
        topologies with hundreds of cascading bifurcations) where
        copy.deepcopy exceeds Python's default recursion limit.
        Falls back to deepcopy for objects that can't be pickled (e.g.,
        mocks in unit tests).
        """
        try:
            return pickle.loads(pickle.dumps(tree, protocol=pickle.HIGHEST_PROTOCOL))
        except (pickle.PicklingError, TypeError, AttributeError):
            import copy

            return copy.deepcopy(tree)

    def _read_tree_with_error(
        self,
        tree_path: str,
        attr_name: str,
        copy_tree: bool = True,
    ):
        try:
            file_hash = self._get_file_hash(tree_path)
            tree = self._cached_tree_read(tree_path, self.tree_format, file_hash)
            if not copy_tree:
                return tree
            # Return a copy to prevent modifications to the cached tree.
            # Uses pickle instead of deepcopy to avoid RecursionError on
            # deeply nested trees.
            return self._fast_copy(tree)
        except FileNotFoundError:
            path = getattr(self, attr_name)
            raise PhykitUserError(
                [
                    f"{path} corresponds to no such file or directory.",
                    "Please check filename and pathing",
                ],
                code=2,
            )

    @staticmethod
    def _scan_simple_newick_summary(file_path: str):
        with open(file_path) as handle:
            text = handle.read()

        if not text or any(char in text for char in "'\"[]"):
            return None

        tip_names = []
        total_len = 0.0
        internal_len = 0.0
        prev_sig = None
        i = 0
        text_len = len(text)
        delimiters = "(),:;"
        whitespace = " \t\r\n"

        while i < text_len:
            char = text[i]
            if char in whitespace:
                i += 1
                continue

            if char == "(":
                prev_sig = "("
                i += 1
                continue
            if char == ",":
                prev_sig = ","
                i += 1
                continue
            if char == ")":
                prev_sig = ")"
                i += 1
                continue
            if char == ";":
                return tuple(tip_names), total_len, internal_len

            if char == ":":
                i += 1
                start = i
                while i < text_len and text[i] not in ",);":
                    if text[i] in whitespace:
                        break
                    i += 1
                number = text[start:i]
                if not number:
                    return None
                try:
                    branch_length = float(number)
                except ValueError:
                    return None
                total_len += branch_length
                if prev_sig in (")", "internal_label"):
                    internal_len += branch_length
                prev_sig = "branch_length"
                continue

            start = i
            while i < text_len and text[i] not in delimiters and text[i] not in whitespace:
                i += 1
            token = text[start:i]
            if not token:
                return None
            if prev_sig in ("(", ",", None):
                tip_names.append(token)
                prev_sig = "terminal_label"
            elif prev_sig == ")":
                prev_sig = "internal_label"
            else:
                return None

        return None

    @staticmethod
    def _scan_simple_newick_tip_names(file_path: str):
        with open(file_path) as handle:
            text = handle.read()

        if not text or any(char in text for char in "'\"[]"):
            return None

        tip_names = []
        prev_sig = None
        i = 0
        text_len = len(text)
        delimiters = "(),:;"
        whitespace = " \t\r\n"

        while i < text_len:
            char = text[i]
            if char in whitespace:
                i += 1
                continue

            if char == "(":
                prev_sig = "("
                i += 1
                continue
            if char == ",":
                prev_sig = ","
                i += 1
                continue
            if char == ")":
                prev_sig = ")"
                i += 1
                continue
            if char == ";":
                return tuple(tip_names)

            if char == ":":
                i += 1
                start = i
                while i < text_len and text[i] not in ",);" and text[i] not in whitespace:
                    i += 1
                if start == i:
                    return None
                prev_sig = "branch_length"
                continue

            start = i
            while i < text_len and text[i] not in delimiters and text[i] not in whitespace:
                i += 1
            token = text[start:i]
            if not token:
                return None
            if prev_sig in ("(", ",", None):
                tip_names.append(token)
                prev_sig = "terminal_label"
            elif prev_sig == ")":
                prev_sig = "internal_label"
            else:
                return None

        return None

    def _get_simple_newick_summary(self, tree_path: str, attr_name: str):
        try:
            file_hash = self._get_file_hash(tree_path)
            return self._cached_simple_newick_summary(tree_path, file_hash)
        except FileNotFoundError:
            path = getattr(self, attr_name)
            raise PhykitUserError(
                [
                    f"{path} corresponds to no such file or directory.",
                    "Please check filename and pathing",
                ],
                code=2,
            )

    def _get_simple_newick_tip_names(self, tree_path: str, attr_name: str):
        try:
            file_hash = self._get_file_hash(tree_path)
            return self._cached_simple_newick_tip_names(tree_path, file_hash)
        except FileNotFoundError:
            path = getattr(self, attr_name)
            raise PhykitUserError(
                [
                    f"{path} corresponds to no such file or directory.",
                    "Please check filename and pathing",
                ],
                code=2,
            )

    @staticmethod
    def _scan_simple_newick_terminal_distance_stats(file_path: str):
        with open(file_path) as handle:
            text = handle.read()

        if not text or any(char in text for char in "'\"[]"):
            return None

        text_len = len(text)
        whitespace = " \t\r\n"
        delimiters = "(),:;"

        def skip_ws(index):
            while index < text_len and text[index] in whitespace:
                index += 1
            return index

        def parse_label(index):
            index = skip_ws(index)
            start = index
            while (
                index < text_len
                and text[index] not in delimiters
                and text[index] not in whitespace
            ):
                index += 1
            return index, index != start

        def parse_length(index):
            index = skip_ws(index)
            if index >= text_len or text[index] != ":":
                return index, 0.0

            index += 1
            index = skip_ws(index)
            start = index
            while index < text_len and text[index] not in ",);" and text[index] not in whitespace:
                index += 1
            number = text[start:index]
            if not number:
                return None
            try:
                length = float(number)
            except ValueError:
                return None
            return skip_ws(index), length

        def parse_node(index):
            index = skip_ws(index)
            if index >= text_len:
                return None

            if text[index] == "(":
                index += 1
                distances = []
                while True:
                    parsed = parse_node(index)
                    if parsed is None:
                        return None
                    child_distances, _child_length, index = parsed
                    distances.extend(child_distances)
                    index = skip_ws(index)
                    if index >= text_len:
                        return None
                    char = text[index]
                    if char == ",":
                        index += 1
                        continue
                    if char == ")":
                        index += 1
                        break
                    return None

                index, _has_label = parse_label(index)
                parsed_length = parse_length(index)
                if parsed_length is None:
                    return None
                index, branch_length = parsed_length
                if branch_length:
                    distances = [distance + branch_length for distance in distances]
                return distances, branch_length, index

            index, has_label = parse_label(index)
            if not has_label:
                return None
            parsed_length = parse_length(index)
            if parsed_length is None:
                return None
            index, branch_length = parsed_length
            return [branch_length], branch_length, index

        try:
            parsed = parse_node(0)
        except RecursionError:
            return None
        if parsed is None:
            return None

        distances, root_branch_length, index = parsed
        index = skip_ws(index)
        if index >= text_len or text[index] != ";":
            return None
        index = skip_ws(index + 1)
        if index != text_len:
            return None

        if root_branch_length:
            distances = [
                distance - root_branch_length
                for distance in distances
            ]
        count = len(distances)
        if count < 2:
            return None

        total = 0.0
        total_sq = 0.0
        for distance in distances:
            total += distance
            total_sq += distance * distance
        return count, total, total_sq

    def _get_simple_newick_terminal_distance_stats(
        self,
        tree_path: str,
        attr_name: str,
    ):
        try:
            file_hash = self._get_file_hash(tree_path)
            return self._cached_simple_newick_terminal_distance_stats(
                tree_path,
                file_hash,
            )
        except FileNotFoundError:
            path = getattr(self, attr_name)
            raise PhykitUserError(
                [
                    f"{path} corresponds to no such file or directory.",
                    "Please check filename and pathing",
                ],
                code=2,
            )

    def write_tree_file(self, tree, output_file_path):
        from Bio import Phylo

        return Phylo.write(tree, output_file_path, self.tree_format)

    def get_tip_names_from_tree(self, tree) -> list:
        """
        get tip names from a tree
        """
        names = self.calculate_terminal_names_fast(tree)
        if names is not None:
            return names
        return [tip.name for tip in tree.get_terminals()]

    @staticmethod
    def _tips_to_prune_for_ordered_mapping(
        tree_tips,
        values,
        min_ordered_size=None,
    ):
        if min_ordered_size is None:
            min_ordered_size = Tree._ORDERED_MAPPING_PRUNE_MIN_SIZE

        n_values = len(values)
        if n_values >= min_ordered_size and len(tree_tips) >= n_values:
            for tip, key in zip(tree_tips, values):
                if tip != key:
                    break
            else:
                return tree_tips[n_values:]

        return [tip for tip in tree_tips if tip not in values]

    @staticmethod
    def _tips_to_prune_for_ordered_names(
        tree_tips,
        ordered_names,
        min_ordered_size=None,
    ):
        if min_ordered_size is None:
            min_ordered_size = Tree._ORDERED_NAMES_PRUNE_MIN_SIZE

        n_names = len(ordered_names)
        tree_tip_count = len(tree_tips)
        if n_names >= min_ordered_size and tree_tip_count >= n_names:
            if (
                n_names == 0
                or (
                    tree_tips[0] == ordered_names[0]
                    and tree_tips[n_names - 1] == ordered_names[-1]
                )
            ):
                if tree_tip_count == n_names:
                    if tree_tips == ordered_names:
                        return []
                else:
                    for index, name in enumerate(ordered_names):
                        if tree_tips[index] != name:
                            break
                    else:
                        return tree_tips[n_names:]

        ordered_name_set = set(ordered_names)
        return [tip for tip in tree_tips if tip not in ordered_name_set]

    def get_first_tip_name_from_tree(self, tree) -> str:
        name = self.calculate_first_terminal_name_fast(tree)
        if name is not None:
            return name
        return tree.get_terminals()[0].name

    @staticmethod
    def calculate_first_terminal_name_fast(tree):
        root = getattr(tree, "root", tree)
        try:
            root.clades
        except AttributeError:
            return None

        clade = root
        children = root.clades
        while True:
            if not isinstance(children, list):
                return None
            if not children:
                return clade.name
            clade = children[0]
            try:
                children = clade.clades
            except AttributeError:
                return None

    @staticmethod
    def calculate_terminal_names_fast(tree):
        root = getattr(tree, "root", tree)
        try:
            root.clades
        except AttributeError:
            return None

        names = []
        stack = [root]
        pop = stack.pop
        append = stack.append
        names_append = names.append
        while stack:
            clade = pop()
            children = clade.clades
            if children:
                child_count = len(children)
                if child_count == 2:
                    append(children[1])
                    append(children[0])
                else:
                    for child in reversed(children):
                        append(child)
            else:
                names_append(clade.name)
        return names

    @staticmethod
    def calculate_terminal_clades_fast(tree):
        root = getattr(tree, "root", tree)
        try:
            root.clades
        except AttributeError:
            return None

        terminals = []
        stack = [root]
        try:
            pop = stack.pop
            append = stack.append
            terminals_append = terminals.append
            while stack:
                clade = pop()
                children = clade.clades
                if children:
                    child_count = len(children)
                    if child_count == 2:
                        append(children[1])
                        append(children[0])
                    else:
                        for child in reversed(children):
                            append(child)
                else:
                    terminals_append(clade)
        except AttributeError:
            return None
        return terminals

    @staticmethod
    def calculate_terminal_count_fast(tree):
        root = getattr(tree, "root", tree)
        try:
            root.clades
        except AttributeError:
            return None

        count = 0
        stack = [root]
        pop = stack.pop
        extend = stack.extend
        try:
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

    @staticmethod
    def calculate_total_branch_length_fast(tree):
        try:
            root = tree.root
        except AttributeError:
            return tree.total_branch_length()

        total_len = 0.0
        stack = [root]
        pop = stack.pop
        extend = stack.extend
        while stack:
            clade = pop()
            branch_length = clade.branch_length
            if branch_length:
                total_len += branch_length
            children = clade.clades
            if children:
                extend(children)
        return total_len

    @staticmethod
    def calculate_total_branch_length_and_terminal_count_fast(tree):
        try:
            root = tree.root
        except AttributeError:
            return tree.total_branch_length(), tree.count_terminals()

        total_len = 0.0
        terminal_count = 0
        stack = [root]
        pop = stack.pop
        extend = stack.extend
        while stack:
            clade = pop()
            branch_length = clade.branch_length
            if branch_length:
                total_len += branch_length

            children = clade.clades
            if children:
                extend(children)
            else:
                terminal_count += 1

        return total_len, terminal_count

    @staticmethod
    def calculate_internal_and_total_branch_length_fast(tree):
        try:
            root = tree.root
        except AttributeError:
            internal_len = 0.0
            for internal in tree.get_nonterminals():
                if internal.branch_length is not None:
                    internal_len += internal.branch_length
            return internal_len, tree.total_branch_length()

        internal_len = 0.0
        total_len = 0.0
        stack = [root]
        pop = stack.pop
        extend = stack.extend
        while stack:
            clade = pop()
            branch_length = clade.branch_length
            children = clade.clades

            if branch_length:
                total_len += branch_length
                if children:
                    internal_len += branch_length

            if children:
                extend(children)

        return internal_len, total_len

    @staticmethod
    def calculate_terminal_root_distances_fast(tree):
        import numpy as np

        try:
            root = tree.root
        except AttributeError:
            return None

        distances = []
        stack = [(root, 0.0)]
        pop = stack.pop
        append = stack.append
        distances_append = distances.append
        while stack:
            clade, distance = pop()
            children = clade.clades
            if children:
                child_count = len(children)
                if child_count == 2:
                    child = children[1]
                    append((child, distance + (child.branch_length or 0.0)))
                    child = children[0]
                    append((child, distance + (child.branch_length or 0.0)))
                else:
                    for idx in range(child_count - 1, -1, -1):
                        child = children[idx]
                        append((child, distance + (child.branch_length or 0.0)))
            else:
                distances_append(distance)

        return np.asarray(distances, dtype=float)

    def calculate_pairwise_tip_distances_fast(
        self,
        tree,
        tips: list[str],
        include_combos: bool = True,
    ) -> tuple[list[tuple[str, str]] | None, list[float]] | None:
        """Calculate all pairwise tip distances from cached root depths.

        Returns None for non-Biopython test doubles so callers can preserve
        their existing fallback behavior.
        """
        try:
            root = tree.root
            root.clades
        except AttributeError:
            return None

        parent_indices = []
        levels = []
        depths = []
        terminal_by_name = {}
        stack = [(root, -1, 0, 0.0)]
        try:
            pop = stack.pop
            append = stack.append
            while stack:
                clade, parent_idx, level, depth = pop()
                clade_idx = len(parent_indices)
                parent_indices.append(parent_idx)
                levels.append(level)
                depths.append(depth)
                children = clade.clades
                if children:
                    next_level = level + 1
                    child_count = len(children)
                    if child_count == 2:
                        child = children[1]
                        append(
                            (
                                child,
                                clade_idx,
                                next_level,
                                depth + (child.branch_length or 0.0),
                            )
                        )
                        child = children[0]
                        append(
                            (
                                child,
                                clade_idx,
                                next_level,
                                depth + (child.branch_length or 0.0),
                            )
                        )
                    else:
                        for idx in range(child_count - 1, -1, -1):
                            child = children[idx]
                            append(
                                (
                                    child,
                                    clade_idx,
                                    next_level,
                                    depth + (child.branch_length or 0.0),
                                )
                            )
                else:
                    terminal_by_name[clade.name] = clade_idx
        except (AttributeError, TypeError):
            return None

        try:
            tip_indices = [terminal_by_name[tip] for tip in tips]
        except (KeyError, TypeError):
            return None

        if not tip_indices:
            return [], []

        max_tip_level = max(levels[tip_idx] for tip_idx in tip_indices)
        if max_tip_level > self._PAIRWISE_LCA_DEPTH_THRESHOLD:
            return self._pairwise_tip_distances_from_lca_index(
                tips,
                tip_indices,
                parent_indices,
                levels,
                depths,
                max_tip_level,
                include_combos=include_combos,
            )

        return self._pairwise_tip_distances_from_paths(
            tips,
            tip_indices,
            parent_indices,
            depths,
            include_combos=include_combos,
        )

    @staticmethod
    def _pairwise_tip_distances_from_paths(
        tips: list[str],
        tip_indices,
        parent_indices,
        depths,
        include_combos: bool = True,
    ) -> tuple[list[tuple[str, str]] | None, list[float]]:
        tip_paths = []
        for tip_idx in tip_indices:
            path = [tip_idx]
            current = tip_idx
            while parent_indices[current] != -1:
                current = parent_indices[current]
                path.append(current)
            path.reverse()
            tip_paths.append(tuple(path))

        pair_count = len(tips) * (len(tips) - 1) // 2
        prealloc_min_pairs = (
            Tree._PAIRWISE_PREALLOC_WITH_COMBOS_MIN_PAIRS
            if include_combos
            else Tree._PAIRWISE_PREALLOC_MIN_PAIRS
        )
        if pair_count >= prealloc_min_pairs:
            combos = [None] * pair_count if include_combos else None
            distances = [0.0] * pair_count
            out_idx = 0
            for i in range(len(tips) - 1):
                tip_a = tips[i]
                path_a = tip_paths[i]
                depth_a = depths[tip_indices[i]]
                for j in range(i + 1, len(tips)):
                    mrca = 0
                    for clade_a, clade_b in zip(path_a, tip_paths[j]):
                        if clade_a != clade_b:
                            break
                        mrca = clade_a
                    if combos is not None:
                        combos[out_idx] = (tip_a, tips[j])
                    distances[out_idx] = (
                        depth_a + depths[tip_indices[j]] - 2 * depths[mrca]
                    )
                    out_idx += 1
            return combos, distances

        combos: list[tuple[str, str]] | None = [] if include_combos else None
        distances: list[float] = []
        for i in range(len(tips) - 1):
            tip_a = tips[i]
            path_a = tip_paths[i]
            depth_a = depths[tip_indices[i]]
            for j in range(i + 1, len(tips)):
                mrca = 0
                for clade_a, clade_b in zip(path_a, tip_paths[j]):
                    if clade_a != clade_b:
                        break
                    mrca = clade_a
                if combos is not None:
                    combos.append((tip_a, tips[j]))
                distances.append(depth_a + depths[tip_indices[j]] - 2 * depths[mrca])

        return combos, distances

    @staticmethod
    def _pairwise_tip_distances_from_lca_index(
        tips: list[str],
        tip_indices,
        parent_indices,
        levels,
        depths,
        max_tip_level,
        include_combos: bool = True,
    ) -> tuple[list[tuple[str, str]] | None, list[float]]:
        jump_count = max(1, max_tip_level.bit_length())
        ancestors = [parent_indices]
        for _ in range(1, jump_count):
            previous = ancestors[-1]
            ancestors.append([
                previous[parent_idx] if parent_idx != -1 else -1
                for parent_idx in previous
            ])

        def lca_index(node_a, node_b):
            if levels[node_a] < levels[node_b]:
                node_a, node_b = node_b, node_a

            level_diff = levels[node_a] - levels[node_b]
            bit = 0
            while level_diff:
                if level_diff & 1:
                    node_a = ancestors[bit][node_a]
                level_diff >>= 1
                bit += 1

            if node_a == node_b:
                return node_a

            for bit in range(jump_count - 1, -1, -1):
                ancestor_a = ancestors[bit][node_a]
                ancestor_b = ancestors[bit][node_b]
                if ancestor_a != ancestor_b:
                    node_a = ancestor_a
                    node_b = ancestor_b

            return parent_indices[node_a]

        pair_count = len(tips) * (len(tips) - 1) // 2
        prealloc_min_pairs = (
            Tree._PAIRWISE_PREALLOC_WITH_COMBOS_MIN_PAIRS
            if include_combos
            else Tree._PAIRWISE_PREALLOC_MIN_PAIRS
        )
        if pair_count >= prealloc_min_pairs:
            combos = [None] * pair_count if include_combos else None
            distances = [0.0] * pair_count
            out_idx = 0
            for i in range(len(tips) - 1):
                tip_a = tips[i]
                tip_a_idx = tip_indices[i]
                depth_a = depths[tip_a_idx]
                for j in range(i + 1, len(tips)):
                    tip_b_idx = tip_indices[j]
                    mrca = lca_index(tip_a_idx, tip_b_idx)
                    if combos is not None:
                        combos[out_idx] = (tip_a, tips[j])
                    distances[out_idx] = depth_a + depths[tip_b_idx] - 2 * depths[mrca]
                    out_idx += 1
            return combos, distances

        combos: list[tuple[str, str]] | None = [] if include_combos else None
        distances: list[float] = []
        for i in range(len(tips) - 1):
            tip_a = tips[i]
            tip_a_idx = tip_indices[i]
            depth_a = depths[tip_a_idx]
            for j in range(i + 1, len(tips)):
                tip_b_idx = tip_indices[j]
                mrca = lca_index(tip_a_idx, tip_b_idx)
                if combos is not None:
                    combos.append((tip_a, tips[j]))
                distances.append(depth_a + depths[tip_b_idx] - 2 * depths[mrca])

        return combos, distances

    def validate_tree(
        self,
        tree,
        min_tips: int = 3,
        require_branch_lengths: bool = False,
        assign_default_branch_length: float = None,
        context: str = "",
    ) -> None:
        """Validate a tree for common issues.

        Parameters
        ----------
        tree : Bio.Phylo tree
        min_tips : minimum number of tips required (default 3)
        require_branch_lengths : if True, raise error when branches lack lengths
        assign_default_branch_length : if set, assign this value to missing lengths
            instead of raising an error (e.g. 1e-8)
        context : description for error messages (e.g. "phylogenetic regression")
        """
        ctx = f" for {context}" if context else ""
        if self._validate_standard_tree(
            tree,
            min_tips,
            require_branch_lengths,
            assign_default_branch_length,
            ctx,
        ):
            return

        if not self._has_minimum_terminals(tree, min_tips):
            raise PhykitUserError(
                [f"Tree must have at least {min_tips} tips{ctx}."],
                code=2,
            )

        if assign_default_branch_length is not None:
            for clade in tree.find_clades():
                if clade.branch_length is None and clade != tree.root:
                    clade.branch_length = assign_default_branch_length
        elif require_branch_lengths:
            for clade in tree.find_clades():
                if clade.branch_length is None and clade != tree.root:
                    raise PhykitUserError(
                        [f"All branches in the tree must have lengths{ctx}."],
                        code=2,
                    )

    @staticmethod
    def _has_minimum_terminals(tree, min_tips: int) -> bool:
        tip_count = 0
        for _terminal in tree.get_terminals():
            tip_count += 1
            if tip_count >= min_tips:
                return True
        return False

    @staticmethod
    def _validate_standard_tree(
        tree,
        min_tips: int,
        require_branch_lengths: bool,
        assign_default_branch_length: float,
        ctx: str,
    ) -> bool:
        try:
            root = tree.root
            root.clades
        except AttributeError:
            return False

        if assign_default_branch_length is None and not require_branch_lengths:
            if min_tips <= 0:
                return True

            tip_count = 0
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
                        tip_count += 1
                        if tip_count >= min_tips:
                            return True
            except AttributeError:
                return False

            raise PhykitUserError(
                [f"Tree must have at least {min_tips} tips{ctx}."],
                code=2,
            )

        tip_count = 0
        missing_branch_length_clades = []
        stack = [root]
        try:
            pop = stack.pop
            extend = stack.extend
            missing_append = missing_branch_length_clades.append
            while stack:
                clade = pop()
                children = clade.clades
                if children:
                    extend(children)
                else:
                    tip_count += 1

                if clade is not root and clade.branch_length is None:
                    missing_append(clade)
        except AttributeError:
            return False

        if tip_count < min_tips:
            raise PhykitUserError(
                [f"Tree must have at least {min_tips} tips{ctx}."],
                code=2,
            )

        if assign_default_branch_length is not None:
            for clade in missing_branch_length_clades:
                clade.branch_length = assign_default_branch_length
        elif require_branch_lengths and missing_branch_length_clades:
            raise PhykitUserError(
                [f"All branches in the tree must have lengths{ctx}."],
                code=2,
            )

        return True

    def shared_tips(self, a, b):
        """
        Determines what tips are shared between two trees
        -------------------------------------------------
        argv: a
            list of tips from one tree
        argv: b
            list of tips from a second tree
        """

        shared = set(a).intersection(b)
        if shared:
            return list(shared)
        raise PhykitUserError(["no common tips"], code=2)

    def prune_tree_using_taxa_list(self, tree, taxa_to_prune: list):
        """
        prune taxa from tree
        """
        if not taxa_to_prune:
            return tree

        target_result = self._terminal_targets_and_count_fast(tree, taxa_to_prune)
        if target_result is not None:
            targets, terminal_count = target_result
            if targets is not None and all(hasattr(target, "clades") for target in targets):
                if (
                    len(targets) < terminal_count
                    and len(set(taxa_to_prune)) == len(taxa_to_prune)
                ):
                    target_ids = {id(target) for target in targets}
                    if self._prune_terminal_objects_batch_standard_tree(
                        tree,
                        target_ids,
                    ):
                        return tree
                for target in targets:
                    tree.prune(target)
                return tree

        terminal_by_name = None
        if target_result is None:
            terminal_by_name = self._terminal_by_name_fast(tree)

        if terminal_by_name is None:
            try:
                terminals = tree.get_terminals()
            except AttributeError:
                terminals = None

            if terminals is not None:
                terminal_by_name = {tip.name: tip for tip in terminals}

        if terminal_by_name is not None:
            try:
                targets = [terminal_by_name[taxon] for taxon in taxa_to_prune]
            except KeyError:
                targets = None
            else:
                if all(hasattr(target, "clades") for target in targets):
                    if (
                        len(targets) < len(terminal_by_name)
                        and len(set(taxa_to_prune)) == len(taxa_to_prune)
                    ):
                        target_ids = {id(target) for target in targets}
                        if self._prune_terminal_objects_batch_standard_tree(
                            tree,
                            target_ids,
                        ):
                            return tree
                    for target in targets:
                        tree.prune(target)
                    return tree

        for taxon in taxa_to_prune:
            tree.prune(taxon)

        return tree

    @staticmethod
    def _terminal_targets_and_count_fast(tree, taxa_to_prune):
        try:
            selected_names = set(taxa_to_prune)
            root = tree.root
            root.clades
        except (AttributeError, TypeError):
            return None

        terminal_by_name = {}
        terminal_count = 0
        stack = [root]
        try:
            pop = stack.pop
            extend = stack.extend
            while stack:
                node = pop()
                children = node.clades
                if children:
                    extend(children)
                else:
                    terminal_count += 1
                    try:
                        if node.name in selected_names:
                            terminal_by_name[node.name] = node
                    except TypeError:
                        pass
            targets = [terminal_by_name[taxon] for taxon in taxa_to_prune]
        except (AttributeError, KeyError):
            targets = None
        return targets, terminal_count

    @staticmethod
    def _terminal_by_name_fast(tree):
        try:
            root = tree.root
            root.clades
        except AttributeError:
            return None

        terminal_by_name = {}
        stack = [root]
        try:
            pop = stack.pop
            extend = stack.extend
            while stack:
                node = pop()
                children = node.clades
                if children:
                    extend(children)
                else:
                    terminal_by_name[node.name] = node
        except AttributeError:
            return None
        return terminal_by_name

    @staticmethod
    def _prune_terminal_objects_batch_standard_tree(tree, target_ids):
        try:
            root = tree.root
            root.clades
        except AttributeError:
            return False

        kept_by_id = {}
        clades = []
        stack = [root]
        try:
            pop = stack.pop
            extend = stack.extend
            append_clade = clades.append
            while stack:
                clade = pop()
                append_clade(clade)
                children = clade.clades
                if children:
                    extend(children)

            clades.reverse()
            for clade in clades:
                children = clade.clades
                if not children:
                    kept_by_id[id(clade)] = (
                        None if id(clade) in target_ids else clade
                    )
                    continue

                kept_children = []
                for child in children:
                    kept_child = kept_by_id[id(child)]
                    if kept_child is not None:
                        kept_children.append(kept_child)

                clade.clades = kept_children
                if not kept_children:
                    kept_by_id[id(clade)] = None
                elif len(kept_children) == 1:
                    child = kept_children[0]
                    if child.branch_length is not None:
                        child.branch_length += clade.branch_length or 0.0
                    kept_by_id[id(clade)] = child
                else:
                    kept_by_id[id(clade)] = clade
        except AttributeError:
            return False

        new_root = kept_by_id[id(root)]
        if new_root is None:
            root.clades = []
        elif new_root is not root:
            tree.root = new_root
        return True

    def calculate_treeness(self, tree=None, print_value=False):
        if not tree:
            tree = self.read_tree_file()

        inter_len, total_len = self.calculate_internal_and_total_branch_length_fast(
            tree
        )

        try:
            treeness = float(inter_len / total_len)
            try:
                if print_value:
                    print(f"{treeness}")
                return treeness
            except BrokenPipeError:
                pass
        except ZeroDivisionError:
            try:
                print("Invalid tree. Tree should contain branch lengths")
                return None
            except BrokenPipeError:
                pass

    @staticmethod
    def get_gap_chars(is_protein: bool = False) -> list[str]:
        if is_protein:
            return ["-", "?", "*", "X", "x"]
        else:
            return ["-", "?", "*", "X", "x", "N", "n"]
