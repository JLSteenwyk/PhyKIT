import sys
import itertools
import copy
from concurrent.futures import ThreadPoolExecutor
from scipy.stats import chisquare
from scipy.stats import _stats_py
from typing import Dict, List, Tuple, Union
import multiprocessing as mp
from functools import partial, lru_cache
import hashlib
import pickle
from unittest.mock import Mock

from Bio import Phylo
from Bio.Phylo import Newick
import numpy as np

from .base import Tree
from ...helpers.files import read_single_column_file_to_list


class PolytomyTest(Tree):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

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

    def process_args(self, args) -> Dict[str, str]:
        return dict(trees=args.trees, groups=args.groups)

    def _read_tree_with_cache(self, tree_path: str) -> Newick.Tree:
        if not hasattr(self, "_tree_cache"):
            self._tree_cache = {}
        if tree_path not in self._tree_cache:
            tree = Phylo.read(tree_path, self.tree_format)
            self._tree_cache[tree_path] = copy.deepcopy(tree)
        return copy.deepcopy(self._tree_cache[tree_path])

    def read_in_groups(
        self
    ) -> List[
        List[
            Union[str, List[str]]
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
        tree_files_batch: List[str],
        groups_of_groups: Dict[str, List[List[str]]],
        outgroup_taxa: List[str],
    ) -> Dict[str, Dict[str, Dict[str, int]]]:
        """Process a batch of trees in parallel."""
        batch_summary = {}
        if isinstance(self.examine_all_triplets_and_sister_pairing, Mock):
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

    def _prepare_tree_for_triplets(self, tree: Newick.Tree, outgroup_taxa: List[str]) -> Newick.Tree:
        prepared = copy.deepcopy(tree)
        if outgroup_taxa:
            try:
                prepared.root_with_outgroup(outgroup_taxa)
            except ValueError:
                pass
        return prepared

    @staticmethod
    def _build_clade_terminal_cache(tree: Newick.Tree) -> Dict[int, Tuple[str, ...]]:
        cache: Dict[int, Tuple[str, ...]] = {}
        for clade in tree.find_clades(order="postorder"):
            if clade.is_terminal():
                cache[id(clade)] = (clade.name,)
            else:
                names: List[str] = []
                for child in clade.clades:
                    names.extend(cache.get(id(child), ()))
                cache[id(clade)] = tuple(names)
        return cache

    def _find_sister_pair(
        self,
        tree: Newick.Tree,
        triplet: Tuple[str, str, str],
        clade_cache: Dict[int, Tuple[str, ...]],
    ) -> Union[Tuple[str, str], None]:
        triplet_set = set(triplet)
        try:
            lca = tree.common_ancestor(triplet)
        except ValueError:
            return None

        assignments: List[set] = []
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

    @staticmethod
    def _guess_tips_from_groups(
        groups_of_groups: Dict[str, List[List[str]]],
        outgroup_taxa: List[str],
    ) -> List[str]:
        tips = set(outgroup_taxa)
        for group_lists in groups_of_groups.values():
            for group in group_lists:
                tips.update(group)
        return list(tips)

    def _legacy_triplet_pass(
        self,
        tips: List[str],
        tree_file: str,
        groups_of_groups: Dict[str, List[List[str]]],
        outgroup_taxa: List[str],
    ) -> Dict[str, int]:
        identifier = list(groups_of_groups.keys())[0]
        triplet_tips = list(itertools.product(*groups_of_groups[identifier]))
        legacy_summary: Dict[str, int] = {}
        for triplet in triplet_tips:
            tree = self.get_triplet_tree(tips, triplet, tree_file, outgroup_taxa)
            if tree and hasattr(tree, "get_terminals"):
                terminal_count = len(list(tree.get_terminals()))
                if terminal_count == 3:
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
        groups_of_groups: Dict[str, List[List[str]]],
    ) -> Dict[str, int]:
        if not groups_of_groups:
            return {}

        tip_names_set = set(self.get_tip_names_from_tree(tree))
        clade_cache = self._build_clade_terminal_cache(tree)
        tree_summary: Dict[str, int] = {}

        for groups in groups_of_groups.values():
            if not groups:
                continue
            for triplet in itertools.product(*groups):
                if not set(triplet).issubset(tip_names_set):
                    continue

                sisters_pair = self._find_sister_pair(tree, triplet, clade_cache)
                if not sisters_pair:
                    continue

                sisters = self.determine_sisters_from_triplet(groups, sisters_pair)
                tree_summary[sisters] = tree_summary.get(sisters, 0) + 1

        return tree_summary

    def loop_through_trees_and_examine_sister_support_among_triplets(
        self,
        trees_file_path: str,
        groups_of_groups: Dict[str, List[List[str]]],
        outgroup_taxa: List[str],
    ) -> Dict[
        str, Dict[str, Dict[str, int]]
    ]:
        """
        go through all trees and all triplets of all trees. For each triplet,
        determine which two taxa are sister to one another
        """
        summary = dict()

        # For small datasets, process sequentially
        if len(trees_file_path) < 10:
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
            num_workers = min(mp.cpu_count(), 8)
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
        groups_arr: List[Union[str, List[str]]],
    ) -> Tuple[
        Dict[str, List[List[str]]],
        List[str],
    ]:
        groups_of_groups = {}

        # Pre-compute group sets for faster lookups
        self._group_sets_cache = {}

        for group in groups_arr:
            temp = []
            group_sets = []
            for i in range(1, 4):
                taxa_list = [taxon_name for taxon_name in group[i]]
                temp.append(taxa_list)
                group_sets.append(frozenset(taxa_list))
            groups_of_groups[group[0]] = temp
            self._group_sets_cache[group[0]] = group_sets

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
        triplet_batch: List[Tuple],
        tips: List[str],
        tree_file: str,
        groups_of_groups: Dict[str, List[List[str]]],
        outgroup_taxa: List[str],
    ) -> Dict[str, Dict[str, int]]:
        """Process a batch of triplets."""
        batch_summary = {}

        for triplet in triplet_batch:
            # Use cached version for tree pruning
            tree = self._get_triplet_tree_cached(
                tuple(tips), triplet, tree_file, tuple(outgroup_taxa)
            )

            if tree and hasattr(tree, 'get_terminals'):
                terminal_count = len(list(tree.get_terminals()))
                if terminal_count == 3:
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
        tips: List[str],
        tree_file: str,
        summary: Dict[
            str, Dict[str, Dict[str, int]]
        ],
        groups_of_groups: Dict[str, List[List[str]]],
        outgroup_taxa: List[str],
    ) -> Dict[str, Dict[str, int]]:
        """
        evaluate all triplets for sister relationships. Polytomies
        in input trees are accounted for
        """
        # get all combinations of three tips
        identifier = list(groups_of_groups.keys())[0]
        triplet_tips = list(itertools.product(*groups_of_groups[identifier]))

        # For small datasets, process sequentially
        if len(triplet_tips) < 50:
            for triplet in triplet_tips:
                tree = self.get_triplet_tree(tips, triplet, tree_file, outgroup_taxa)
                if tree and hasattr(tree, 'get_terminals'):
                    terminal_count = len(list(tree.get_terminals()))
                    if terminal_count == 3:
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
        triplet: Tuple[str, str, str],
        groups: List[List[str]]
    ) -> int:
        """
        determine how many groups are represented in a triplet
        """
        # Convert groups to tuple of frozensets for caching
        groups_tuple = tuple(frozenset(group) for group in groups)
        return self._count_groups_cached(triplet, groups_tuple)

    def set_branch_lengths_in_tree_to_one(
        self,
        tree: Newick.Tree
    ) -> None:
        # Single pass through all clades
        for clade in tree.find_clades():
            clade.branch_length = 1

    def check_if_triplet_is_a_polytomy(self, tree: Newick.Tree) -> bool:
        """
        count the number of internal branches. If 1, then the triplet is a polytomy
        """
        # Direct check without intermediate list creation
        nonterminal_count = sum(1 for _ in tree.get_nonterminals())
        return nonterminal_count == 1

    def sister_relationship_counter(
        self,
        tree_file: str,
        summary: Dict[str, Dict[str, int]],
        sisters: str,
    ) -> Dict[str, Dict[str, int]]:
        """
        counter for how many times a particular sister relationship is observed
        """
        # if tree is not in summary, create a key for it
        if tree_file not in summary.keys():
            summary[str(tree_file)] = {}
        # if the sister relationship is not in the tree file dict, create a key for it
        if sisters not in summary[str(tree_file)].keys():
            summary[str(tree_file)][sisters] = 1
        else:
            summary[str(tree_file)][sisters] += 1

        return summary

    def get_triplet_tree(
        self,
        tips: List[str],
        triplet: Tuple[str, str, str],
        tree_file: str,
        outgroup_taxa: List[str],
    ) -> Newick.Tree:
        """
        get a tree object of only the triplet of interest
        """
        # determine tips that are not in the triplet of interest
        tips_to_prune = list(set(tips) - set(list(triplet)))
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
        groups: List[List[str]],
        pair: Tuple[str]
    ) -> str:
        """
        determine sister taxa from a triplet
        """
        # Convert to tuple of frozensets for caching
        groups_tuple = tuple(frozenset(group) for group in groups)
        return self._determine_sisters_cached(groups_tuple, pair)

    def determine_sisters_and_add_to_counter(
        self,
        tip_names: List[str],
        tree: Newick.Tree,
        tree_file: str,
        groups: List[List[str]],
        summary: Dict[str, Dict[str, int]],
    ) -> Dict[str, Dict[str, int]]:
        """
        determine which pair of taxa are sister to one another
        and add 1 to the counter for the sister pair
        """
        # get pairs from tip names
        pairs = list(itertools.combinations(tip_names, 2))
        for pair in pairs:
            is_polytomy = self.check_if_triplet_is_a_polytomy(tree)
            # if distance between pair is 2 and the triplet is
            # not a polytomy (i.e., having only 1 internal branch)
            # then report the sisters in the triplet
            if tree.distance(pair[0], pair[1]) == 2 and not is_polytomy:
                # determine which two tips are sisters
                sisters = self.determine_sisters_from_triplet(groups, pair)
                # add to summary dictionary of how many times that sister
                # relationship is observed
                summary = self.sister_relationship_counter(tree_file, summary, sisters)
        return summary

    def get_triplet_and_gene_support_freq_counts(
        self,
        summary: Dict[str, Dict[str, int]]
    ) -> Tuple[
        Dict[str, int], Dict[str, int]
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
            if "0-1" not in summary[tree].keys():
                summary[tree]["0-1"] = 0
            if "0-2" not in summary[tree].keys():
                summary[tree]["0-2"] = 0
            if "1-2" not in summary[tree].keys():
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
    ) -> Tuple[
        _stats_py.Power_divergenceResult,
        _stats_py.Power_divergenceResult,
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
        gene_support_freq: Dict[str, int],
        trees_file_path: List[str],
    ) -> None:
        """
        print results to stdout for user
        """
        try:
            print(f"Gene Support Frequency Results")
            print(f"==============================")
            print(f"chi-squared: {round(gene_support_freq_res.statistic, 4)}")
            print(f"p-value: {round(gene_support_freq_res.pvalue, 6)}")
            print(
                f"total genes: {(gene_support_freq['0-1'] + gene_support_freq['0-2'] + gene_support_freq['1-2'])}"
            )
            print(f"0-1: {gene_support_freq['0-1']}")
            print(f"0-2: {gene_support_freq['0-2']}")
            print(f"1-2: {gene_support_freq['1-2']}")
        except BrokenPipeError:
            pass
