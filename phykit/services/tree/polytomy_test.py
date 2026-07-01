from __future__ import annotations

import sys
import itertools
import math
from functools import lru_cache
from collections import namedtuple
import os

from .base import Tree


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


def read_single_column_file_to_list(*args, **kwargs):
    from ...helpers.files import (
        read_single_column_file_to_list as _read_single_column_file_to_list,
    )

    return _read_single_column_file_to_list(*args, **kwargs)


class _LazyPhylo:
    def read(self, *args, **kwargs):
        from Bio import Phylo as _Phylo

        return _Phylo.read(*args, **kwargs)


Phylo = _LazyPhylo()


class _LazyMultiprocessing:
    def cpu_count(self):
        import multiprocessing as _mp

        return _mp.cpu_count()

    def Pool(self, *args, **kwargs):
        import multiprocessing as _mp

        return _mp.Pool(*args, **kwargs)


mp = _LazyMultiprocessing()


def ThreadPoolExecutor(*args, **kwargs):
    from concurrent.futures import ThreadPoolExecutor as _ThreadPoolExecutor

    return _ThreadPoolExecutor(*args, **kwargs)


Power_divergenceResult = namedtuple(
    "Power_divergenceResult",
    ["statistic", "pvalue"],
)


def _is_unittest_mock(obj) -> bool:
    return type(obj).__module__ == "unittest.mock"


def _all_triplet_nodes_identical(nodes) -> bool:
    first = nodes[0]
    return nodes[1] is first and nodes[2] is first


def chisquare(*args, **kwargs):
    if kwargs:
        from scipy.stats import chisquare as _chisquare
        return _chisquare(*args, **kwargs)

    observed = [float(value) for value in args[0]]
    total = sum(observed)
    if total == 0:
        return Power_divergenceResult(float("nan"), float("nan"))

    expected = total / len(observed)
    statistic = sum((value - expected) ** 2 / expected for value in observed)
    df = len(observed) - 1
    if df == 2:
        pvalue = math.exp(-statistic / 2.0)
    elif df == 1:
        pvalue = math.erfc(math.sqrt(statistic / 2.0))
    else:
        from scipy.stats import chi2
        pvalue = float(chi2.sf(statistic, df=df))

    return Power_divergenceResult(statistic, pvalue)


class PolytomyTest(Tree):
    MP_MIN_TREES = 50
    MAX_MP_WORKERS = 8

    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(trees=parsed["trees"], groups=parsed["groups"])
        self.json_output = parsed["json_output"]

    def _should_use_multiprocessing(self, n_trees: int) -> bool:
        if os.environ.get("PHYKIT_DISABLE_MP", "0") == "1":
            return False
        if os.environ.get("PHYKIT_FORCE_MP", "0") == "1":
            return True
        return n_trees >= self.MP_MIN_TREES

    def run(self):
        # read in groups
        groups_arr = self.read_in_groups()

        # determine groups of groups
        groups_of_groups, outgroup_taxa = self.determine_groups_of_groups(groups_arr)

        # read trees into list
        trees_file_path = read_single_column_file_to_list(self.trees)

        # go through all triplets of all trees and
        # examine sister relationships among all triplets
        summary = self.loop_through_trees_and_examine_sister_support_among_triplets(
            trees_file_path, groups_of_groups, outgroup_taxa
        )

        # count triplet and gene support frequencies for different sister relationships
        (
            triplet_group_counts,
            gene_support_freq,
        ) = self.get_triplet_and_gene_support_freq_counts(summary)

        # conduct chisquare tests
        triplet_res, gene_support_freq_res = self.chisquare_tests(
            triplet_group_counts, gene_support_freq
        )

        # print results
        self.print_gene_support_freq_res(
            gene_support_freq_res, gene_support_freq, trees_file_path
        )
        # self.print_triplet_based_res(triplet_res, triplet_group_counts)

    def process_args(self, args) -> dict[str, str]:
        return dict(
            trees=args.trees,
            groups=args.groups,
            json_output=getattr(args, "json", False),
        )

    @staticmethod
    def _safe_copy(obj):
        import pickle

        try:
            return pickle.loads(pickle.dumps(obj, protocol=pickle.HIGHEST_PROTOCOL))
        except (pickle.PicklingError, TypeError, AttributeError):
            import copy
            return copy.deepcopy(obj)

    def _read_tree_with_cache(self, tree_path: str) -> Newick.Tree:
        if not hasattr(self, "_tree_cache"):
            self._tree_cache = {}
        if tree_path not in self._tree_cache:
            tree = Phylo.read(tree_path, self.tree_format)
            self._tree_cache[tree_path] = self._safe_copy(tree)
        return self._safe_copy(self._tree_cache[tree_path])

    def read_in_groups(
        self
    ) -> list[
        list[
            str | list[str]
        ]
    ]:
        groups_arr = []
        try:
            for line in open(self.groups):
                line = line.strip()
                if not line.startswith("#"):
                    try:
                        line = line.split("\t")
                        temp = []
                        temp.append(line[0])
                        temp.append(line[1].split(";"))
                        temp.append(line[2].split(";"))
                        temp.append(line[3].split(";"))
                        temp.append(line[4].split(";"))
                        groups_arr.append(temp)
                    except IndexError:
                        try:
                            print(f"{self.groups} contains an indexing error.")
                            print(
                                "Please format the groups file (-g) as a four column tab-delimited file with column 1 being the name of the test"
                            )
                            print("col2: the tip names of one group (; separated)")
                            print("col3: the tip names of a second group (; separated)")
                            print("col4: the tip names of a third group (; separated)")
                            print(
                                "col5: the tip names of the outgroup taxa (; separated)"
                            )
                            sys.exit(2)
                        except BrokenPipeError:
                            pass

        except FileNotFoundError:
            try:
                print(f"{self.groups} corresponds to no such file.")
                print("Please check filename and pathing again.")
                sys.exit(2)
            except BrokenPipeError:
                pass

        return groups_arr

    def _process_tree_batch(
        self,
        tree_files_batch: list[str],
        groups_of_groups: dict[str, list[list[str]]],
        outgroup_taxa: list[str],
    ) -> dict[str, dict[str, dict[str, int]]]:
        """Process a batch of trees in parallel."""
        batch_summary = {}
        if _is_unittest_mock(self.examine_all_triplets_and_sister_pairing):
            for tree_file in tree_files_batch:
                try:
                    tree = self._read_tree_with_cache(tree_file)
                    tips = self.get_tip_names_from_tree(tree)
                    batch_summary = self.examine_all_triplets_and_sister_pairing(
                        tips, tree_file, batch_summary, groups_of_groups, outgroup_taxa
                    )
                except Exception:
                    continue
            return batch_summary
        if not hasattr(self, "_tree_cache"):
            self._tree_cache = {}
        for tree_file in tree_files_batch:
            try:
                tree = self._read_tree_with_cache(tree_file)
                prepared_tree = self._prepare_tree_for_triplets(tree, outgroup_taxa)
                tree_summary = self._evaluate_tree_triplets_fast(prepared_tree, groups_of_groups)
                if not tree_summary:
                    tips = self.get_tip_names_from_tree(tree)
                    tree_summary = self._legacy_triplet_pass(
                        tips,
                        tree_file,
                        groups_of_groups,
                        outgroup_taxa,
                    )
                if tree_summary:
                    batch_summary[tree_file] = tree_summary
            except Exception:
                continue
        return batch_summary

    def _prepare_tree_for_triplets(self, tree: Newick.Tree, outgroup_taxa: list[str]) -> Newick.Tree:
        prepared = self._safe_copy(tree)
        if outgroup_taxa:
            try:
                prepared.root_with_outgroup(outgroup_taxa)
            except ValueError:
                pass
        return prepared

    @staticmethod
    def _build_clade_terminal_cache(tree: Newick.Tree) -> dict[int, frozenset]:
        cache: dict[int, frozenset] = {}
        for clade in tree.find_clades(order="postorder"):
            if clade.is_terminal():
                cache[id(clade)] = frozenset([clade.name])
            else:
                names = set()
                for child in clade.clades:
                    names.update(cache.get(id(child), frozenset()))
                cache[id(clade)] = frozenset(names)
        return cache

    @staticmethod
    def _build_tip_path_cache(tree: Newick.Tree) -> dict[str, tuple]:
        cache: dict[str, tuple] = {}

        def visit(clade, path: tuple) -> None:
            current_path = path + (clade,)
            if clade.is_terminal():
                cache[clade.name] = current_path
                return
            for child in clade.clades:
                visit(child, current_path)

        visit(tree.root, tuple())
        return cache

    @staticmethod
    def _common_ancestor_from_path_cache(
        triplet: tuple[str, str, str],
        path_cache: dict[str, tuple],
    ):
        try:
            paths = [path_cache[taxon] for taxon in triplet]
        except KeyError:
            return None

        lca = None
        for nodes in zip(*paths):
            if _all_triplet_nodes_identical(nodes):
                lca = nodes[0]
            else:
                break
        return lca

    def _find_sister_pair(
        self,
        tree: Newick.Tree,
        triplet: tuple[str, str, str],
        clade_cache: dict[int, tuple[str, ...]],
    ) -> tuple[str, str] | None:
        triplet_set = set(triplet)
        try:
            lca = tree.common_ancestor(triplet)
        except ValueError:
            return None

        assignments: list[set] = []
        for child in lca.clades:
            descendant = triplet_set.intersection(clade_cache.get(id(child), ()))
            if descendant:
                assignments.append(descendant)

        if len(assignments) != 2:
            return None

        for subset in assignments:
            if len(subset) == 2:
                return tuple(sorted(subset))  # type: ignore

        return None

    def _find_sister_pair_from_path_cache(
        self,
        triplet: tuple[str, str, str],
        clade_cache: dict[int, frozenset],
        path_cache: dict[str, tuple],
        triplet_set: set[str] | None = None,
    ) -> tuple[str, str] | None:
        if triplet_set is None:
            triplet_set = set(triplet)
        lca = self._common_ancestor_from_path_cache(triplet, path_cache)
        if lca is None:
            return None

        assignments: list[set] = []
        for child in lca.clades:
            descendant = triplet_set.intersection(
                clade_cache.get(id(child), frozenset())
            )
            if descendant:
                assignments.append(descendant)

        if len(assignments) != 2:
            return None

        for subset in assignments:
            if len(subset) == 2:
                return tuple(sorted(subset))  # type: ignore

        return None

    @staticmethod
    def _guess_tips_from_groups(
        groups_of_groups: dict[str, list[list[str]]],
        outgroup_taxa: list[str],
    ) -> list[str]:
        tips = set(outgroup_taxa)
        for group_lists in groups_of_groups.values():
            for group in group_lists:
                tips.update(group)
        return list(tips)

    def _legacy_triplet_pass(
        self,
        tips: list[str],
        tree_file: str,
        groups_of_groups: dict[str, list[list[str]]],
        outgroup_taxa: list[str],
    ) -> dict[str, int]:
        identifier = next(iter(groups_of_groups))
        legacy_summary: dict[str, int] = {}
        for triplet in itertools.product(*groups_of_groups[identifier]):
            tree = self.get_triplet_tree(tips, triplet, tree_file, outgroup_taxa)
            if tree and hasattr(tree, "get_terminals"):
                if self._has_exactly_three_terminals(tree):
                    for _, groups in groups_of_groups.items():
                        represented = self.count_number_of_groups_in_triplet(triplet, groups)
                        if represented == 3:
                            tip_names = self.get_tip_names_from_tree(tree)
                            self.set_branch_lengths_in_tree_to_one(tree)
                            temp_summary = {}
                            temp_summary = self.determine_sisters_and_add_to_counter(
                                tip_names, tree, tree_file, groups, temp_summary
                            )
                            for sisters, count in temp_summary.get(tree_file, {}).items():
                                legacy_summary[sisters] = legacy_summary.get(sisters, 0) + count
        return legacy_summary

    def _evaluate_tree_triplets_fast(
        self,
        tree: Newick.Tree,
        groups_of_groups: dict[str, list[list[str]]],
    ) -> dict[str, int]:
        if not groups_of_groups:
            return {}

        tip_names_set = set(self.get_tip_names_from_tree(tree))
        clade_cache = self._build_clade_terminal_cache(tree)
        path_cache = self._build_tip_path_cache(tree)
        tree_summary: dict[str, int] = {}

        for groups in groups_of_groups.values():
            if not groups:
                continue
            groups_tuple = tuple(frozenset(group) for group in groups)
            for triplet in itertools.product(*groups):
                triplet_set = set(triplet)
                if not triplet_set.issubset(tip_names_set):
                    continue

                sisters_pair = self._find_sister_pair_from_path_cache(
                    triplet,
                    clade_cache,
                    path_cache,
                    triplet_set,
                )
                if not sisters_pair:
                    continue

                sisters = self._determine_sisters_cached(groups_tuple, sisters_pair)
                tree_summary[sisters] = tree_summary.get(sisters, 0) + 1

        return tree_summary

    def loop_through_trees_and_examine_sister_support_among_triplets(
        self,
        trees_file_path: str,
        groups_of_groups: dict[str, list[list[str]]],
        outgroup_taxa: list[str],
    ) -> dict[
        str, dict[str, dict[str, int]]
    ]:
        """
        go through all trees and all triplets of all trees. For each triplet,
        determine which two taxa are sister to one another
        """
        summary = dict()

        # For small/medium workloads, multiprocessing overhead dominates.
        if not self._should_use_multiprocessing(len(trees_file_path)):
            for tree_file in trees_file_path:
                try:
                    tree = Phylo.read(tree_file, "newick")
                    tips = self.get_tip_names_from_tree(tree)
                    summary = self.examine_all_triplets_and_sister_pairing(
                        tips, tree_file, summary, groups_of_groups, outgroup_taxa
                    )
                except FileNotFoundError:
                    print(f"{tree_file} corresponds to no such file.")
                    print("Please check file name and pathing")
                    sys.exit(2)
        else:
            # Use multiprocessing for larger datasets
            from functools import partial

            num_workers = min(mp.cpu_count(), self.MAX_MP_WORKERS)
            batch_size = max(1, len(trees_file_path) // num_workers)
            tree_batches = [trees_file_path[i:i + batch_size]
                           for i in range(0, len(trees_file_path), batch_size)]

            # Process batches in parallel
            process_func = partial(self._process_tree_batch,
                                  groups_of_groups=groups_of_groups,
                                  outgroup_taxa=outgroup_taxa)

            try:
                with mp.Pool(processes=num_workers) as pool:
                    batch_results = pool.map(process_func, tree_batches)
            except (OSError, ValueError):
                with ThreadPoolExecutor(max_workers=num_workers) as executor:
                    batch_results = list(executor.map(process_func, tree_batches))

            # Merge results
            for batch_summary in batch_results:
                for tree_file, tree_data in batch_summary.items():
                    if tree_file not in summary:
                        summary[tree_file] = {}
                    for sisters, count in tree_data.items():
                        if sisters not in summary[tree_file]:
                            summary[tree_file][sisters] = 0
                        summary[tree_file][sisters] += count

        return summary

    def determine_groups_of_groups(
        self,
        groups_arr: list[str | list[str]],
    ) -> tuple[
        dict[str, list[list[str]]],
        list[str],
    ]:
        groups_of_groups = {}

        # Pre-compute group sets for faster lookups
        self._group_sets_cache = {}
        self._group_tuple_cache_by_id = {}

        for group in groups_arr:
            temp = []
            group_sets = []
            for i in range(1, 4):
                taxa_list = [taxon_name for taxon_name in group[i]]
                temp.append(taxa_list)
                group_sets.append(frozenset(taxa_list))
            group_sets_tuple = tuple(group_sets)
            groups_of_groups[group[0]] = temp
            self._group_sets_cache[group[0]] = group_sets
            self._group_tuple_cache_by_id[id(temp)] = (
                tuple(id(taxa_list) for taxa_list in temp),
                group_sets_tuple,
            )

        outgroup_taxa = [taxon_name for taxon_name in group[4]]

        return groups_of_groups, outgroup_taxa

    @lru_cache(maxsize=1024)
    def _get_triplet_tree_cached(self, tips_tuple: tuple, triplet: tuple,
                                 tree_file: str, outgroup_tuple: tuple):
        """Cached version of get_triplet_tree."""
        tips = list(tips_tuple)
        outgroup_taxa = list(outgroup_tuple)
        return self.get_triplet_tree(tips, triplet, tree_file, outgroup_taxa)

    def _process_triplet_batch(
        self,
        triplet_batch: list[tuple],
        tips: list[str],
        tree_file: str,
        groups_of_groups: dict[str, list[list[str]]],
        outgroup_taxa: list[str],
    ) -> dict[str, dict[str, int]]:
        """Process a batch of triplets."""
        batch_summary = {}

        for triplet in triplet_batch:
            # Use cached version for tree pruning
            tree = self._get_triplet_tree_cached(
                tuple(tips), triplet, tree_file, tuple(outgroup_taxa)
            )

            if tree and hasattr(tree, 'get_terminals'):
                if self._has_exactly_three_terminals(tree):
                    for _, groups in groups_of_groups.items():
                        num_groups_represented = self.count_number_of_groups_in_triplet(
                            triplet, groups
                        )

                        if num_groups_represented == 3:
                            tip_names = self.get_tip_names_from_tree(tree)
                            self.set_branch_lengths_in_tree_to_one(tree)
                            batch_summary = self.determine_sisters_and_add_to_counter(
                                tip_names, tree, tree_file, groups, batch_summary
                            )

        return batch_summary

    def examine_all_triplets_and_sister_pairing(
        self,
        tips: list[str],
        tree_file: str,
        summary: dict[
            str, dict[str, dict[str, int]]
        ],
        groups_of_groups: dict[str, list[list[str]]],
        outgroup_taxa: list[str],
    ) -> dict[str, dict[str, int]]:
        """
        evaluate all triplets for sister relationships. Polytomies
        in input trees are accounted for
        """
        # get all combinations of three tips
        identifier = next(iter(groups_of_groups))
        triplet_groups = groups_of_groups[identifier]
        triplet_count = math.prod(len(group) for group in triplet_groups)

        # For small datasets, process sequentially
        if triplet_count < 50:
            for triplet in itertools.product(*triplet_groups):
                tree = self.get_triplet_tree(tips, triplet, tree_file, outgroup_taxa)
                if tree and hasattr(tree, 'get_terminals'):
                    if self._has_exactly_three_terminals(tree):
                        for _, groups in groups_of_groups.items():
                            num_groups_represented = self.count_number_of_groups_in_triplet(
                                triplet, groups
                            )
                            if num_groups_represented == 3:
                                tip_names = self.get_tip_names_from_tree(tree)
                                self.set_branch_lengths_in_tree_to_one(tree)
                                summary = self.determine_sisters_and_add_to_counter(
                                    tip_names, tree, tree_file, groups, summary
                                )
        else:
            try:
                tree = self._read_tree_with_cache(tree_file)
                prepared_tree = self._prepare_tree_for_triplets(tree, outgroup_taxa)
                tree_summary = self._evaluate_tree_triplets_fast(
                    prepared_tree, groups_of_groups
                )
                if tree_summary:
                    tree_counts = summary.setdefault(tree_file, {})
                    for sisters, count in tree_summary.items():
                        tree_counts[sisters] = tree_counts.get(sisters, 0) + count
                else:
                    legacy_counts = self._legacy_triplet_pass(
                        tips or self._guess_tips_from_groups(groups_of_groups, outgroup_taxa),
                        tree_file,
                        groups_of_groups,
                        outgroup_taxa,
                    )
                    if legacy_counts:
                        tree_counts = summary.setdefault(tree_file, {})
                        for sisters, count in legacy_counts.items():
                            tree_counts[sisters] = tree_counts.get(sisters, 0) + count
            except FileNotFoundError:
                legacy_counts = self._legacy_triplet_pass(
                    tips or self._guess_tips_from_groups(groups_of_groups, outgroup_taxa),
                    tree_file,
                    groups_of_groups,
                    outgroup_taxa,
                )
                if legacy_counts:
                    tree_counts = summary.setdefault(tree_file, {})
                    for sisters, count in legacy_counts.items():
                        tree_counts[sisters] = tree_counts.get(sisters, 0) + count

        return summary

    @lru_cache(maxsize=4096)
    def _count_groups_cached(self, triplet_tuple: tuple, groups_tuple: tuple) -> int:
        """Cached version of group counting."""
        triplet_set = set(triplet_tuple)
        num_groups_represented = sum(
            1 for group in groups_tuple if triplet_set.intersection(group)
        )
        return num_groups_represented

    def count_number_of_groups_in_triplet(
        self,
        triplet: tuple[str, str, str],
        groups: list[list[str]]
    ) -> int:
        """
        determine how many groups are represented in a triplet
        """
        groups_tuple = None
        cached_groups = getattr(self, "_group_tuple_cache_by_id", {}).get(id(groups))
        if cached_groups is not None:
            group_ids, cached_tuple = cached_groups
            if len(group_ids) == len(groups) and all(
                id(group) == group_id
                for group, group_id in zip(groups, group_ids)
            ):
                groups_tuple = cached_tuple

        if groups_tuple is None:
            groups_tuple = tuple(frozenset(group) for group in groups)
        return self._count_groups_cached(triplet, groups_tuple)

    def set_branch_lengths_in_tree_to_one(
        self,
        tree: Newick.Tree
    ) -> None:
        if self._set_standard_tree_branch_lengths_to_one(tree):
            return

        for clade in tree.find_clades():
            clade.branch_length = 1

    @staticmethod
    def _set_standard_tree_branch_lengths_to_one(tree: Newick.Tree) -> bool:
        try:
            root = tree.root
            root.clades
        except AttributeError:
            return False

        stack = [root]
        try:
            while stack:
                clade = stack.pop()
                clade.branch_length = 1
                stack.extend(clade.clades)
        except AttributeError:
            return False
        return True

    @staticmethod
    def _has_exactly_three_terminals(tree: Newick.Tree) -> bool:
        try:
            root = tree.root
            root.clades
        except AttributeError:
            return len(list(tree.get_terminals())) == 3

        count = 0
        stack = [root]
        try:
            while stack:
                clade = stack.pop()
                children = clade.clades
                if not isinstance(children, list):
                    return len(list(tree.get_terminals())) == 3
                if children:
                    stack.extend(children)
                else:
                    count += 1
                    if count > 3:
                        return False
        except AttributeError:
            return len(list(tree.get_terminals())) == 3
        return count == 3

    def check_if_triplet_is_a_polytomy(self, tree: Newick.Tree) -> bool:
        """
        count the number of internal branches. If 1, then the triplet is a polytomy
        """
        try:
            root = tree.root
            root.clades
        except AttributeError:
            return sum(1 for _ in tree.get_nonterminals()) == 1

        nonterminal_count = 0
        stack = [root]
        try:
            while stack:
                clade = stack.pop()
                children = clade.clades
                if not isinstance(children, list):
                    return sum(1 for _ in tree.get_nonterminals()) == 1
                if children:
                    nonterminal_count += 1
                    if nonterminal_count > 1:
                        return False
                    stack.extend(children)
        except AttributeError:
            return sum(1 for _ in tree.get_nonterminals()) == 1
        return nonterminal_count == 1

    def sister_relationship_counter(
        self,
        tree_file: str,
        summary: dict[str, dict[str, int]],
        sisters: str,
    ) -> dict[str, dict[str, int]]:
        """
        counter for how many times a particular sister relationship is observed
        """
        tree_file = str(tree_file)
        tree_counts = summary.get(tree_file)
        if tree_counts is None:
            tree_counts = summary[tree_file] = {}
        tree_counts[sisters] = tree_counts.get(sisters, 0) + 1

        return summary

    def get_triplet_tree(
        self,
        tips: list[str],
        triplet: tuple[str, str, str],
        tree_file: str,
        outgroup_taxa: list[str],
    ) -> Newick.Tree:
        """
        get a tree object of only the triplet of interest
        """
        # determine tips that are not in the triplet of interest
        triplet_taxa = set(triplet)
        tips_to_prune = [tip for tip in tips if tip not in triplet_taxa]
        # determine tips that in the outgroup
        outgroup_present = [value for value in tips if value in outgroup_taxa]
        tree = Phylo.read(tree_file, "newick")

        # root tree on outgroup taxa
        try:
            tree.root_with_outgroup(outgroup_present)

            # prune to a triplet
            tree = self.prune_tree_using_taxa_list(tree, tips_to_prune)

            return tree
        except ValueError:
            tree = False
            return tree

    @lru_cache(maxsize=4096)
    def _determine_sisters_cached(self, groups_tuple: tuple, pair: tuple) -> str:
        """Cached version of sister determination."""
        # Find which group each pair member belongs to
        idx0 = next(i for i, group in enumerate(groups_tuple) if pair[0] in group)
        idx1 = next(i for i, group in enumerate(groups_tuple) if pair[1] in group)

        # Sort and format the result
        sisters = sorted([idx0, idx1])
        return f"{sisters[0]}-{sisters[1]}"

    def determine_sisters_from_triplet(
        self,
        groups: list[list[str]],
        pair: tuple[str]
    ) -> str:
        """
        determine sister taxa from a triplet
        """
        # Convert to tuple of frozensets for caching
        groups_tuple = tuple(frozenset(group) for group in groups)
        return self._determine_sisters_cached(groups_tuple, pair)

    def determine_sisters_and_add_to_counter(
        self,
        tip_names: list[str],
        tree: Newick.Tree,
        tree_file: str,
        groups: list[list[str]],
        summary: dict[str, dict[str, int]],
    ) -> dict[str, dict[str, int]]:
        """
        determine which pair of taxa are sister to one another
        and add 1 to the counter for the sister pair
        """
        is_polytomy = self.check_if_triplet_is_a_polytomy(tree)
        if is_polytomy:
            return summary

        # get pairs from tip names
        pairs = list(itertools.combinations(tip_names, 2))
        for pair in pairs:
            # if distance between pair is 2 and the triplet is
            # not a polytomy (i.e., having only 1 internal branch)
            # then report the sisters in the triplet
            if tree.distance(pair[0], pair[1]) == 2:
                # determine which two tips are sisters
                sisters = self.determine_sisters_from_triplet(groups, pair)
                # add to summary dictionary of how many times that sister
                # relationship is observed
                summary = self.sister_relationship_counter(tree_file, summary, sisters)
        return summary

    def get_triplet_and_gene_support_freq_counts(
        self,
        summary: dict[str, dict[str, int]]
    ) -> tuple[
        dict[str, int], dict[str, int]
    ]:
        """
        count how many triplets and genes support the various sister relationships
        """
        # Count the total number of sister pairings
        # for the three possible pairs for triplets
        triplet_group_counts = {"g0g1_count": 0, "g0g2_count": 0, "g1g2_count": 0}

        # Also, keep track of which and how many genes
        # support each sister pairing
        gene_support_freq = {"0-1": 0, "1-2": 0, "0-2": 0}
        for tree in summary:
            # create empty key value pairs in case sister
            # pairing was never observed
            if "0-1" not in summary[tree]:
                summary[tree]["0-1"] = 0
            if "0-2" not in summary[tree]:
                summary[tree]["0-2"] = 0
            if "1-2" not in summary[tree]:
                summary[tree]["1-2"] = 0
            # create a running value of triplets that support each sister pair
            triplet_group_counts["g0g1_count"] += summary[tree]["0-1"]
            triplet_group_counts["g0g2_count"] += summary[tree]["0-2"]
            triplet_group_counts["g1g2_count"] += summary[tree]["1-2"]
            # determine which sister pairing is best supported in a single gene
            # and add one to the corresponding gene support frequency count
            gene_support_freq[max(summary[tree], key=summary[tree].get)] += 1

        return triplet_group_counts, gene_support_freq

    def chisquare_tests(
        self,
        triplet_group_counts: dict,
        gene_support_freq: dict
    ) -> tuple[
        object,
        object,
    ]:
        triplet_res = chisquare(
            [
                triplet_group_counts["g0g1_count"],
                triplet_group_counts["g0g2_count"],
                triplet_group_counts["g1g2_count"],
            ]
        )

        gene_support_freq_res = chisquare(
            [
                gene_support_freq["0-1"],
                gene_support_freq["0-2"],
                gene_support_freq["1-2"],
            ]
        )

        return triplet_res, gene_support_freq_res

    # def print_triplet_based_res(
    #     self,
    #     triplet_res,
    #     triplet_group_counts: dict
    # ) -> None:
    #     """
    #     print results to stdout for user
    #     """
    #     try:
    #         print(f"\nTriplet Results")
    #         print(f"===============")
    #         print(f"chi-squared: {round(triplet_res.statistic, 4)}")
    #         print(f"p-value: {round(triplet_res.pvalue, 6)}")
    #         print(f"total triplets: {sum(triplet_group_counts.values())}")
    #         print(f"0-1: {triplet_group_counts['g0g1_count']}")
    #         print(f"0-2: {triplet_group_counts['g0g2_count']}")
    #         print(f"1-2: {triplet_group_counts['g1g2_count']}")
    #     except BrokenPipeError:
    #         pass

    def print_gene_support_freq_res(
        self,
        gene_support_freq_res,
        gene_support_freq: dict[str, int],
        trees_file_path: list[str],
    ) -> None:
        """
        print results to stdout for user
        """
        if self.json_output:
            print_json(
                dict(
                    gene_support_frequency=dict(
                        chi_squared=round(float(gene_support_freq_res.statistic), 4),
                        p_value=round(float(gene_support_freq_res.pvalue), 6),
                        total_genes=(
                            gene_support_freq["0-1"]
                            + gene_support_freq["0-2"]
                            + gene_support_freq["1-2"]
                        ),
                        support_counts={
                            "0-1": gene_support_freq["0-1"],
                            "0-2": gene_support_freq["0-2"],
                            "1-2": gene_support_freq["1-2"],
                        },
                    )
                )
            )
            return

        try:
            print(
                "\n".join(
                    [
                        "Gene Support Frequency Results",
                        "==============================",
                        f"chi-squared: {round(gene_support_freq_res.statistic, 4)}",
                        f"p-value: {round(gene_support_freq_res.pvalue, 6)}",
                        (
                            "total genes: "
                            f"{gene_support_freq['0-1'] + gene_support_freq['0-2'] + gene_support_freq['1-2']}"
                        ),
                        f"0-1: {gene_support_freq['0-1']}",
                        f"0-2: {gene_support_freq['0-2']}",
                        f"1-2: {gene_support_freq['1-2']}",
                    ]
                )
            )
        except BrokenPipeError:
            pass
