from __future__ import annotations

import math
import os
from io import StringIO
from pathlib import Path

from .base import Tree
from ...errors import PhykitUserError

_path_isabs = os.path.isabs


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


class _LazyNumpy:
    def __getattr__(self, name):
        import numpy as _np

        return getattr(_np, name)


class _LazyPhylo:
    def read(self, *args, **kwargs):
        from Bio import Phylo as _Phylo

        return _Phylo.read(*args, **kwargs)


np = _LazyNumpy()
Phylo = _LazyPhylo()


class _LazyPickle:
    def __getattr__(self, name):
        import pickle as _pickle

        return getattr(_pickle, name)

    def dumps(self, *args, **kwargs):
        import pickle as _pickle

        return _pickle.dumps(*args, **kwargs)

    def loads(self, *args, **kwargs):
        import pickle as _pickle

        return _pickle.loads(*args, **kwargs)


pickle = _LazyPickle()


def _mean_and_population_variance(values):
    n = len(values)
    mean_value = sum(values) / n
    variance = (
        sum((value - mean_value) * (value - mean_value) for value in values) / n
    )
    return float(mean_value), float(variance)


class ConcordanceAsr(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"])
        self.gene_trees_path = parsed["gene_trees_path"]
        self.trait_data_path = parsed["trait_data_path"]
        self.trait_column = parsed["trait_column"]
        self.method = parsed["method"]
        self.ci = parsed["ci"]
        self.plot_output = parsed["plot_output"]
        self.plot_uncertainty = parsed["plot_uncertainty"]
        self.missing_taxa = parsed["missing_taxa"]
        self.json_output = parsed["json_output"]
        self.plot_config = parsed["plot_config"]

    def process_args(self, args) -> dict:
        from ...helpers.plot_config import PlotConfig

        return dict(
            tree_file_path=args.tree,
            gene_trees_path=args.gene_trees,
            trait_data_path=args.trait_data,
            trait_column=getattr(args, "trait", None),
            method=getattr(args, "method", "weighted"),
            ci=getattr(args, "ci", False),
            plot_output=getattr(args, "plot", None),
            plot_uncertainty=getattr(args, "plot_uncertainty", None),
            missing_taxa=getattr(args, "missing_taxa", "shared"),
            json_output=getattr(args, "json", False),
            plot_config=PlotConfig.from_args(args),
        )

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
                append(children[-1])
                if len(children) == 2:
                    append(children[0])
                else:
                    for idx in range(len(children) - 2, -1, -1):
                        append(children[idx])

    def run(self) -> None:
        species_tree = self.read_tree_file()

        # Create a helper ASR instance for reusing methods
        from .ancestral_reconstruction import AncestralReconstruction

        self._asr = AncestralReconstruction.__new__(AncestralReconstruction)
        self._asr.ci = True  # always compute CIs internally

        self._asr.validate_tree(species_tree, min_tips=3, require_branch_lengths=True, context="concordance ASR")

        gene_trees = self._parse_gene_trees(self.gene_trees_path)

        if len(gene_trees) < 2:
            raise PhykitUserError(
                [
                    "At least 2 gene trees are required for concordance-aware ASR.",
                    f"Only {len(gene_trees)} gene tree(s) provided.",
                ],
                code=2,
            )

        # Validate gene tree branch lengths
        for i, gt in enumerate(gene_trees):
            for clade in gt.find_clades():
                if clade.branch_length is None and clade != gt.root:
                    raise PhykitUserError(
                        [
                            f"Gene tree {i + 1} has branches without lengths.",
                            "All gene trees must have branch lengths.",
                        ],
                        code=2,
                    )

        # Get taxa and handle mismatches
        species_tips = set(self.get_tip_names_from_tree(species_tree))
        gene_trees, all_taxa = self._normalize_taxa(
            species_tree, gene_trees, species_tips
        )

        # Parse trait data
        tree_tips = sorted(all_taxa)
        if self.trait_column is not None:
            trait_values = self._asr._parse_multi_trait_data(
                self.trait_data_path, tree_tips, self.trait_column
            )
        else:
            trait_values = self._asr._parse_single_trait_data(
                self.trait_data_path, tree_tips
            )

        # Prune species tree to shared taxa with trait data
        sp_tips = self.get_tip_names_from_tree(species_tree)
        tips_to_prune = self._tips_to_prune_for_ordered_mapping(
            sp_tips, trait_values
        )
        needs_species_copy = bool(tips_to_prune) or self.plot_config.ladderize
        species_copy = (
            self._fast_tree_copy(species_tree) if needs_species_copy else species_tree
        )
        if tips_to_prune:
            species_copy = self.prune_tree_using_taxa_list(
                species_copy, tips_to_prune
            )

        if self.plot_config.ladderize:
            species_copy.ladderize()

        if self.method == "distribution":
            result = self._run_distribution(
                species_copy, gene_trees, trait_values, all_taxa
            )
        else:
            result = self._run_weighted(
                species_copy, gene_trees, trait_values, all_taxa
            )

        if self.plot_output:
            self._plot_concordance_contmap(
                species_copy, result, self.plot_output
            )
            result["plot_output"] = self.plot_output

        if self.plot_uncertainty:
            self._plot_uncertainty(
                species_copy, result, self.plot_uncertainty
            )
            result["plot_uncertainty"] = self.plot_uncertainty

        if self.json_output:
            print_json(result)
        else:
            self._print_text_output(result)

    # ------------------------------------------------------------------
    # Gene tree parsing
    # ------------------------------------------------------------------

    def _parse_gene_trees(self, path: str) -> list:
        source = Path(path)
        trees = []
        trees_append = trees.append
        parent_str = str(source.parent)
        parent_prefix = "" if parent_str == "." else parent_str + os.sep
        try:
            with source.open() as handle:
                for line in handle:
                    stripped = line.strip()
                    if not stripped or stripped[0] == "#":
                        continue
                    if stripped[0] == "(":
                        trees_append(Phylo.read(StringIO(stripped), "newick"))
                    else:
                        tree_path = (
                            stripped
                            if _path_isabs(stripped)
                            else parent_prefix + stripped
                        )
                        trees_append(Phylo.read(tree_path, "newick"))
        except FileNotFoundError:
            raise PhykitUserError(
                [
                    f"{path} corresponds to no such file or directory.",
                    "Please check filename and pathing",
                ],
                code=2,
            )
        return trees

    # ------------------------------------------------------------------
    # Taxa normalization
    # ------------------------------------------------------------------

    def _normalize_taxa(self, species_tree, gene_trees, species_tips):
        gene_tip_sets = []
        for gt in gene_trees:
            gene_tip_sets.append(set(self.get_tip_names_from_tree(gt)))

        all_gene_taxa = set()
        for ts in gene_tip_sets:
            all_gene_taxa |= ts

        shared = species_tips & all_gene_taxa
        for ts in gene_tip_sets:
            shared &= ts

        if self.missing_taxa == "error":
            if shared != species_tips or any(ts != species_tips for ts in gene_tip_sets):
                raise PhykitUserError(
                    [
                        "Taxa mismatch between species tree and gene trees.",
                        "Use --missing-taxa shared to prune to shared taxa.",
                    ],
                    code=2,
                )

        if len(shared) < 3:
            raise PhykitUserError(
                [
                    f"Only {len(shared)} shared taxa across species and gene trees.",
                    "At least 3 shared taxa are required.",
                ],
                code=2,
            )

        # Prune species tree
        sp_tips_to_prune = [
            name for name in self.get_tip_names_from_tree(species_tree)
            if name not in shared
        ]
        if sp_tips_to_prune:
            species_tree = self.prune_tree_using_taxa_list(
                species_tree, sp_tips_to_prune
            )

        # Prune gene trees
        pruned_gene_trees = []
        for gt in gene_trees:
            gt_tips_to_prune = [
                name for name in self.get_tip_names_from_tree(gt)
                if name not in shared
            ]
            if gt_tips_to_prune:
                gt = self.prune_tree_using_taxa_list(gt, gt_tips_to_prune)
            pruned_gene_trees.append(gt)

        return pruned_gene_trees, shared

    # ------------------------------------------------------------------
    # gCF computation
    # ------------------------------------------------------------------

    @staticmethod
    def _canonical_split(taxa_side, all_taxa):
        taxa_side_len = len(taxa_side)
        all_taxa_len = len(all_taxa)
        if taxa_side_len * 2 < all_taxa_len:
            return frozenset(taxa_side)
        complement = all_taxa - taxa_side
        if taxa_side_len * 2 > all_taxa_len:
            return frozenset(complement)
        if not taxa_side:
            return frozenset(taxa_side)
        if min(taxa_side) <= min(complement):
            return frozenset(taxa_side)
        return frozenset(complement)

    @staticmethod
    def _collect_clade_tip_sets(tree):
        try:
            root = tree.root
            root.clades
        except AttributeError:
            root = None

        if root is not None:
            preorder = []
            stack = [root]
            pop = stack.pop
            append = stack.append
            append_preorder = preorder.append
            try:
                while stack:
                    clade = pop()
                    append_preorder(clade)
                    children = clade.clades
                    if children:
                        child_count = len(children)
                        if child_count == 2:
                            append(children[1])
                            append(children[0])
                        else:
                            for idx in range(child_count - 1, -1, -1):
                                append(children[idx])

                clade_tips = {}
                for clade in reversed(preorder):
                    children = clade.clades
                    if not children:
                        clade_tips[id(clade)] = frozenset((clade.name,))
                    elif len(children) == 2:
                        clade_tips[id(clade)] = (
                            clade_tips[id(children[0])]
                            | clade_tips[id(children[1])]
                        )
                    else:
                        tips = frozenset()
                        for child in children:
                            tips = tips | clade_tips.get(id(child), frozenset())
                        clade_tips[id(clade)] = tips
                return clade_tips
            except AttributeError:
                pass

        clade_tips = {}
        for clade in tree.find_clades(order="postorder"):
            if clade.is_terminal():
                clade_tips[id(clade)] = frozenset({clade.name})
            else:
                tips = frozenset()
                for child in clade.clades:
                    tips = tips | clade_tips.get(id(child), frozenset())
                clade_tips[id(clade)] = tips
        return clade_tips

    def _get_four_groups(
        self,
        tree,
        node,
        parent_map,
        all_taxa_fs,
        clade_tip_sets=None,
    ):
        """Identify the four subtree groups around an internal branch.

        For the branch connecting node to its parent:
          C1 = tips of node's first child
          C2 = tips of node's second child
          S  = tips of node's sibling under parent
          D  = remaining tips (everything above parent)

        Returns (C1, C2, S, D) as frozensets, or None if decomposition
        is not possible (e.g., node is root with <3 children, or leaf).
        """
        if node.is_terminal() or len(node.clades) < 2:
            return None

        if clade_tip_sets is None:
            clade_tip_sets = self._collect_clade_tip_sets(tree)

        empty = frozenset()
        children = node.clades
        C1 = clade_tip_sets.get(id(children[0]), empty)
        C2 = clade_tip_sets.get(id(children[1]), empty)
        # If node has >2 children (polytomy), merge extras into C2
        if len(children) > 2:
            C2 = C2.union(
                *(
                    clade_tip_sets.get(id(extra_child), empty)
                    for extra_child in children[2:]
                )
            )

        parent = parent_map.get(id(node))
        if parent is None:
            # node is root
            # For trifurcating root, can't define four groups for the root itself
            # (root has no branch above it).
            return None

        # Get siblings of node under parent
        siblings = [c for c in parent.clades if id(c) != id(node)]
        if not siblings:
            return None

        S = clade_tip_sets.get(id(siblings[0]), empty)
        # D = everything else (other siblings + above parent)
        D = all_taxa_fs - C1 - C2 - S

        return C1, C2, S, D

    def _compute_gcf_per_node(
        self, species_tree, gene_trees, all_taxa
    ) -> dict[int, tuple[float, float, float]]:
        """For each internal node, compute (gCF, gDF1, gDF2).

        Uses the four-group decomposition (C1, C2, S, D) around each
        internal branch to identify the concordant bipartition and the
        two NNI alternative bipartitions:
          Concordant: C1+C2 | S+D
          NNI alt 1 (swap C1<->S): S+C2 | C1+D
          NNI alt 2 (swap C2<->S): C1+S | C2+D
        """
        all_taxa_fs = frozenset(all_taxa)
        parent_map = self._asr._build_parent_map(species_tree)
        species_clade_tips = self._collect_clade_tip_sets(species_tree)

        # Extract bipartitions from all gene trees
        gene_tree_splits = []
        for gt in gene_trees:
            gt_clade_tips = self._collect_clade_tip_sets(gt)
            splits = set()
            try:
                root = gt.root
                root.clades
                clade_iter = self._iter_preorder(root)
                direct_clades = True
            except AttributeError:
                clade_iter = gt.get_nonterminals()
                direct_clades = False

            for clade in clade_iter:
                if direct_clades and not clade.clades:
                    continue
                tips = gt_clade_tips.get(id(clade), frozenset())
                if len(tips) <= 1 or tips == all_taxa_fs:
                    continue
                splits.add(self._canonical_split(tips, all_taxa_fs))
            gene_tree_splits.append(splits)

        result = {}
        try:
            root = species_tree.root
            root.clades
            species_clades = self._iter_preorder(root)
            direct_species_clades = True
        except AttributeError:
            species_clades = species_tree.find_clades(order="preorder")
            direct_species_clades = False

        for clade in species_clades:
            if direct_species_clades:
                if not clade.clades:
                    continue
            elif clade.is_terminal():
                continue

            groups = self._get_four_groups(
                species_tree,
                clade,
                parent_map,
                all_taxa_fs,
                species_clade_tips,
            )
            if groups is None:
                continue

            C1, C2, S, D = groups

            # The three quartet resolutions as bipartitions
            concordant_bp = self._canonical_split(C1 | C2, all_taxa_fs)
            nni_alt1_bp = self._canonical_split(S | C2, all_taxa_fs)
            nni_alt2_bp = self._canonical_split(C1 | S, all_taxa_fs)

            concordant, disc1, disc2 = self._count_gcf_topologies(
                gene_tree_splits,
                concordant_bp,
                nni_alt1_bp,
                nni_alt2_bp,
            )

            total = concordant + disc1 + disc2
            if total > 0:
                gcf = concordant / total
                gdf1 = disc1 / total
                gdf2 = disc2 / total
            else:
                gcf = 1.0
                gdf1 = 0.0
                gdf2 = 0.0

            result[id(clade)] = (gcf, gdf1, gdf2)

        return result

    @staticmethod
    def _count_gcf_topologies(
        gene_tree_splits,
        concordant_bp,
        nni_alt1_bp,
        nni_alt2_bp,
    ):
        concordant = 0
        disc1 = 0
        disc2 = 0
        for splits in gene_tree_splits:
            if concordant_bp in splits:
                concordant += 1
            if nni_alt1_bp in splits:
                disc1 += 1
            if nni_alt2_bp in splits:
                disc2 += 1
        return concordant, disc1, disc2

    # ------------------------------------------------------------------
    # NNI alternative tree building
    # ------------------------------------------------------------------

    @staticmethod
    def _fast_tree_copy(tree):
        return pickle.loads(pickle.dumps(tree, protocol=pickle.HIGHEST_PROTOCOL))

    def _build_nni_alternative_at_node(self, species_tree, target_node, parent_map):
        """Build two NNI alternative trees at a given internal branch.

        Returns list of (tree, expected_desc_frozenset) tuples, where
        expected_desc_frozenset is the descendant set of the target node
        position in the NNI alternative tree.

        For target_node with children C1, C2 and sibling S:
          swap_idx=0: swap C1 <-> S  =>  target gets {S, C2}
          swap_idx=1: swap C2 <-> S  =>  target gets {C1, S}
        """
        alternatives = []

        if id(target_node) not in parent_map:
            # target_node is root — no branch above it, skip
            return alternatives

        parent = parent_map[id(target_node)]
        siblings = [c for c in parent.clades if id(c) != id(target_node)]
        if not siblings or len(target_node.clades) < 2:
            return alternatives

        sibling = siblings[0]

        for swap_idx in range(min(2, len(target_node.clades))):
            tree_copy = self._fast_tree_copy(species_tree)
            copy_clades = list(tree_copy.find_clades(order="preorder"))
            orig_clades = list(species_tree.find_clades(order="preorder"))
            idx_map = {id(o): i for i, o in enumerate(orig_clades)}

            t_idx = idx_map[id(target_node)]
            p_idx = idx_map[id(parent)]
            s_idx = idx_map[id(sibling)]

            copy_target = copy_clades[t_idx]
            copy_parent = copy_clades[p_idx]
            copy_sibling = copy_clades[s_idx]

            # Swap target's child with sibling
            swapped_child = copy_target.clades[swap_idx]
            sib_pos = copy_parent.clades.index(copy_sibling)
            copy_parent.clades[sib_pos] = swapped_child
            copy_target.clades[swap_idx] = copy_sibling

            # Compute expected descendant set after the swap
            expected_desc = frozenset(
                t.name for t in copy_target.get_terminals()
            )
            alternatives.append((tree_copy, expected_desc))

        return alternatives

    # ------------------------------------------------------------------
    # ASR on a single tree
    # ------------------------------------------------------------------

    def _run_asr_on_tree(self, tree, trait_values):
        """Run _fast_anc on one tree, return results keyed by descendant frozensets."""
        tip_names = self.get_tip_names_from_tree(tree)
        tips_to_prune = self._tips_to_prune_for_ordered_mapping(
            tip_names, trait_values
        )
        tree_for_analysis = self._fast_tree_copy(tree) if tips_to_prune else tree
        if tips_to_prune:
            tree_for_analysis = self.prune_tree_using_taxa_list(
                tree_for_analysis, tips_to_prune
            )
            remaining_tips = self.get_tip_names_from_tree(tree_for_analysis)
        else:
            remaining_tips = tip_names

        if len(remaining_tips) < 3:
            return {}, {}, 0.0

        ordered_names = sorted(trait_values.keys() & set(remaining_tips))
        if len(ordered_names) < 3:
            return {}, {}, 0.0

        x = np.array([trait_values[name] for name in ordered_names])
        node_labels = self._asr._label_internal_nodes(tree_for_analysis)

        node_estimates, node_cis, sigma2, _ = self._asr._fast_anc(
            tree_for_analysis, x, ordered_names, node_labels
        )

        # Re-key by descendant frozensets
        estimates_by_desc = {}
        cis_by_desc = {}
        clade_tip_sets = self._collect_clade_tip_sets(tree_for_analysis)
        for clade in tree_for_analysis.find_clades(order="preorder"):
            if clade.is_terminal():
                continue
            if id(clade) not in node_labels:
                continue
            label = node_labels[id(clade)]
            if label in node_estimates:
                desc = clade_tip_sets.get(id(clade), frozenset())
                estimates_by_desc[desc] = node_estimates[label]
                if label in node_cis:
                    cis_by_desc[desc] = node_cis[label]

        return estimates_by_desc, cis_by_desc, sigma2

    # ------------------------------------------------------------------
    # Law of total variance
    # ------------------------------------------------------------------

    @staticmethod
    def _law_of_total_variance(weights, means, variances):
        """Compute law of total variance.

        Returns (total_variance, within_variance, between_variance).
        total = E[Var(X|T)] + Var(E[X|T])
        """
        total_weight = float(sum(weights))
        if total_weight == 0.0:
            return 0.0, 0.0, 0.0

        inv_total_weight = 1.0 / total_weight

        # Within-group variance: E[Var(X|T)]
        within_var = 0.0

        # Between-group variance: Var(E[X|T])
        weighted_mean = 0.0
        for weight, mean, variance in zip(weights, means, variances):
            normalized_weight = weight * inv_total_weight
            within_var += normalized_weight * variance
            weighted_mean += normalized_weight * mean

        between_var = 0.0
        for weight, mean in zip(weights, means):
            diff = mean - weighted_mean
            between_var += weight * inv_total_weight * diff * diff

        total_var = within_var + between_var
        return total_var, within_var, between_var

    # ------------------------------------------------------------------
    # Weighted method
    # ------------------------------------------------------------------

    def _run_weighted(self, species_tree, gene_trees, trait_values, all_taxa):
        """Concordance-weighted ASR on the species tree."""
        node_labels = self._asr._label_internal_nodes(species_tree)
        parent_map = self._asr._build_parent_map(species_tree)

        # Run ASR on species tree
        sp_estimates, sp_cis, sp_sigma2 = self._run_asr_on_tree(
            species_tree, trait_values
        )

        # Compute gCF
        gcf_per_node = self._compute_gcf_per_node(
            species_tree, gene_trees, all_taxa
        )
        species_clade_tips = self._collect_clade_tip_sets(species_tree)

        # Build NNI alternatives and run ASR on them
        # nni_estimates[node_id] = list of (estimates_dict, cis_dict, sigma2, expected_desc)
        nni_estimates = {}
        for clade in species_tree.find_clades(order="preorder"):
            if clade.is_terminal():
                continue
            if id(clade) not in gcf_per_node:
                continue

            alts = self._build_nni_alternative_at_node(
                species_tree, clade, parent_map
            )
            alt_results = []
            for alt_tree, expected_desc in alts:
                est, ci, sig = self._run_asr_on_tree(alt_tree, trait_values)
                alt_results.append((est, ci, sig, expected_desc))
            nni_estimates[id(clade)] = alt_results

        # Combine estimates
        ancestral_estimates = {}
        root_id = id(species_tree.root)

        for clade in species_tree.find_clades(order="preorder"):
            if clade.is_terminal():
                continue
            if id(clade) not in node_labels:
                continue

            label = node_labels[id(clade)]
            desc = species_clade_tips.get(id(clade), frozenset())
            descendants = sorted(desc)

            sp_est = sp_estimates.get(desc)
            if sp_est is None:
                continue

            gcf, gdf1, gdf2 = gcf_per_node.get(id(clade), (1.0, 0.0, 0.0))

            # Get sp CI variance
            sp_ci = sp_cis.get(desc)
            sp_var = 0.0
            if sp_ci:
                # CI = est +/- 1.96*se, so se = (ci_upper - ci_lower) / (2*1.96)
                sp_var = ((sp_ci[1] - sp_ci[0]) / (2 * 1.96)) ** 2

            # Collect estimates from NNI alternatives
            weights = [gcf]
            means = [sp_est]
            variances = [sp_var]

            alt_results = nni_estimates.get(id(clade), [])
            alt_weights = [gdf1, gdf2]
            for i, (alt_est, alt_ci, _, expected_desc) in enumerate(alt_results):
                if i >= 2:
                    break
                # Use the expected descendant set for the target node in the NNI tree
                alt_e = alt_est.get(expected_desc)
                if alt_e is not None:
                    weights.append(alt_weights[i])
                    means.append(alt_e)
                    alt_v = 0.0
                    alt_c = alt_ci.get(expected_desc)
                    if alt_c:
                        alt_v = ((alt_c[1] - alt_c[0]) / (2 * 1.96)) ** 2
                    variances.append(alt_v)

            # Normalize weights
            w_sum = sum(weights)
            if w_sum > 0:
                weights = [w / w_sum for w in weights]

            # Weighted estimate
            weighted_est = sum(w * m for w, m in zip(weights, means))

            # Law of total variance
            total_var, within_var, between_var = self._law_of_total_variance(
                weights, means, variances
            )

            entry = {
                "estimate": float(weighted_est),
                "descendants": descendants,
                "is_root": id(clade) == root_id,
                "gcf": float(gcf),
                "gdf1": float(gdf1),
                "gdf2": float(gdf2),
                "var_topology": float(between_var),
                "var_parameter": float(within_var),
                "source_estimates": [float(m) for m in means],
                "source_weights": [float(w) for w in weights],
            }

            if self.ci:
                se = math.sqrt(total_var) if total_var > 0 else 0.0
                entry["ci_lower"] = float(weighted_est - 1.96 * se)
                entry["ci_upper"] = float(weighted_est + 1.96 * se)

            ancestral_estimates[label] = entry

        return {
            "method": "weighted",
            "n_tips": len(trait_values),
            "n_gene_trees": len(gene_trees),
            "sigma2": float(sp_sigma2),
            "ancestral_estimates": ancestral_estimates,
        }

    # ------------------------------------------------------------------
    # Distribution method
    # ------------------------------------------------------------------

    def _run_distribution(self, species_tree, gene_trees, trait_values, all_taxa):
        """Reconstruct across gene tree distribution."""
        node_labels = self._asr._label_internal_nodes(species_tree)
        species_clade_tips = self._collect_clade_tip_sets(species_tree)

        # Compute gCF
        gcf_per_node = self._compute_gcf_per_node(
            species_tree, gene_trees, all_taxa
        )

        # Run ASR on each gene tree
        gene_tree_results = []
        for gt in gene_trees:
            est, ci, sig = self._run_asr_on_tree(gt, trait_values)
            gene_tree_results.append((est, ci, sig))

        # For each species tree node, collect estimates from gene trees
        ancestral_estimates = {}
        root_id = id(species_tree.root)

        for clade in species_tree.find_clades(order="preorder"):
            if clade.is_terminal():
                continue
            if id(clade) not in node_labels:
                continue

            label = node_labels[id(clade)]
            desc = species_clade_tips.get(id(clade), frozenset())
            descendants = sorted(desc)

            gcf, gdf1, gdf2 = gcf_per_node.get(id(clade), (1.0, 0.0, 0.0))

            # Collect estimates from gene trees that contain this bipartition
            estimates_list = []
            for est_dict, ci_dict, sig in gene_tree_results:
                if desc in est_dict:
                    estimates_list.append(est_dict[desc])

            if not estimates_list:
                continue

            mean_est, var_topology = _mean_and_population_variance(estimates_list)
            entry = {
                "estimate": mean_est,
                "descendants": descendants,
                "is_root": id(clade) == root_id,
                "gcf": float(gcf),
                "gdf1": float(gdf1),
                "gdf2": float(gdf2),
                "n_gene_trees_with_node": len(estimates_list),
                "var_topology": var_topology,
                "var_parameter": 0.0,
                "gene_tree_estimates": estimates_list,
            }

            if self.ci and len(estimates_list) >= 2:
                ci_lower = float(np.percentile(estimates_list, 2.5))
                ci_upper = float(np.percentile(estimates_list, 97.5))
                entry["ci_lower"] = ci_lower
                entry["ci_upper"] = ci_upper

            ancestral_estimates[label] = entry

        # Average sigma2 across gene trees
        sigma2_values = [sig for _, _, sig in gene_tree_results if sig > 0]
        avg_sigma2 = (
            float(sum(sigma2_values) / len(sigma2_values))
            if sigma2_values
            else 0.0
        )

        return {
            "method": "distribution",
            "n_tips": len(trait_values),
            "n_gene_trees": len(gene_trees),
            "sigma2": avg_sigma2,
            "ancestral_estimates": ancestral_estimates,
        }

    # ------------------------------------------------------------------
    # Output formatting
    # ------------------------------------------------------------------

    def _print_text_output(self, result) -> None:
        lines = [
            "Concordance-Aware Ancestral State Reconstruction",
            f"\nMethod: {result['method']}",
            f"Number of tips: {result['n_tips']}",
            f"Number of gene trees: {result['n_gene_trees']}",
            f"Sigma-squared (BM rate): {result['sigma2']:.6f}",
            "\nAncestral estimates:",
        ]
        estimates = result["ancestral_estimates"]

        has_ci = any("ci_lower" in e for e in estimates.values())

        if has_ci:
            lines.append(
                f"  {'Node':<12s}{'Desc':>6s}{'Estimate':>12s}"
                f"{'gCF':>8s}{'95% CI':>24s}"
                f"{'Var_topo':>12s}{'Var_param':>12s}"
            )
        else:
            lines.append(
                f"  {'Node':<12s}{'Desc':>6s}{'Estimate':>12s}"
                f"{'gCF':>8s}{'Var_topo':>12s}{'Var_param':>12s}"
            )

        for label, entry in estimates.items():
            n_desc = len(entry["descendants"])
            est = entry["estimate"]
            gcf = entry["gcf"]
            root_tag = " (root)" if entry.get("is_root", False) else ""
            var_t = entry.get("var_topology", 0.0)
            var_p = entry.get("var_parameter", 0.0)

            if has_ci and "ci_lower" in entry:
                ci_str = f"[{entry['ci_lower']:.4f}, {entry['ci_upper']:.4f}]"
                lines.append(
                    f"  {label + root_tag:<12s}{n_desc:>6d}{est:>12.4f}"
                    f"{gcf:>8.3f}{ci_str:>24s}"
                    f"{var_t:>12.6f}{var_p:>12.6f}"
                )
            elif has_ci:
                lines.append(
                    f"  {label + root_tag:<12s}{n_desc:>6d}{est:>12.4f}"
                    f"{gcf:>8.3f}{'':>24s}"
                    f"{var_t:>12.6f}{var_p:>12.6f}"
                )
            else:
                lines.append(
                    f"  {label + root_tag:<12s}{n_desc:>6d}{est:>12.4f}"
                    f"{gcf:>8.3f}{var_t:>12.6f}{var_p:>12.6f}"
                )
        print("\n".join(lines))

    # ------------------------------------------------------------------
    # Plotting
    # ------------------------------------------------------------------

    def _plot_concordance_contmap(self, tree, result, output_path) -> None:
        from ...helpers.plot_config import compute_node_positions
        from ...helpers.circular_layout import (
            compute_circular_coords,
            draw_circular_branches,
            draw_circular_tip_labels,
        )
        from ...helpers.color_annotations import (
            parse_color_file,
            resolve_mrca,
            draw_range_rect,
            draw_range_wedge,
            build_color_legend_handles,
            apply_label_colors,
        )

        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            from matplotlib.collections import LineCollection
            from matplotlib.colors import Normalize
        except ImportError:
            print(
                "matplotlib is required for concordance contMap plotting. "
                "Install matplotlib and retry."
            )
            raise SystemExit(2)

        estimates = result["ancestral_estimates"]
        preorder_clades = list(self._iter_preorder(tree.root))
        node_labels = {}
        counter = 1
        for clade in preorder_clades:
            if not clade.clades:
                continue
            if clade.name:
                node_labels[id(clade)] = clade.name
            else:
                node_labels[id(clade)] = f"N{counter}"
                counter += 1

        parent_map = {}
        for clade in preorder_clades:
            for child in clade.clades:
                parent_map[id(child)] = clade

        # Map estimates back to tree nodes
        label_to_est = {}
        label_to_gcf = {}
        for label, entry in estimates.items():
            label_to_est[label] = entry["estimate"]
            label_to_gcf[label] = entry.get("gcf", 1.0)

        # Build estimates dict keyed by id(clade)
        all_estimates = {}
        gcf_values = {}
        for clade in preorder_clades:
            if not clade.clades or id(clade) not in node_labels:
                continue
            label = node_labels[id(clade)]
            if label in label_to_est:
                all_estimates[id(clade)] = label_to_est[label]
                gcf_values[id(clade)] = label_to_gcf[label]

        tips = [clade for clade in preorder_clades if not clade.clades]
        root = tree.root
        node_x, node_y = compute_node_positions(
            tree,
            parent_map,
            cladogram=self.plot_config.cladogram,
            preorder_clades=preorder_clades,
        )

        config = self.plot_config
        config.resolve(n_rows=len(tips), n_cols=None)
        fig, ax = plt.subplots(figsize=(config.fig_width, config.fig_height))

        if self.plot_config.circular:
            # --- Circular mode ---
            coords = compute_circular_coords(
                tree,
                node_x,
                parent_map,
                preorder_clades=preorder_clades,
                terminal_clades=tips,
            )
            ax.set_aspect("equal")
            ax.axis("off")

            # Draw branches (gray)
            draw_circular_branches(ax, tree, coords, parent_map, color="gray", lw=2)

            # Draw gCF-sized dots at internal nodes
            gcf_x = []
            gcf_y = []
            gcf_sizes = []
            gcf_colors = []
            for clade in preorder_clades:
                if not clade.clades:
                    continue
                cid = id(clade)
                if cid in gcf_values and cid in coords:
                    gcf = gcf_values[cid]
                    size = 50 + 200 * gcf
                    color = plt.cm.RdYlGn(gcf)
                    gcf_x.append(coords[cid]["x"])
                    gcf_y.append(coords[cid]["y"])
                    gcf_sizes.append(size)
                    gcf_colors.append(color)
                    ax.annotate(
                        f"{gcf:.2f}",
                        (coords[cid]["x"], coords[cid]["y"]),
                        textcoords="offset points",
                        xytext=(5, 5),
                        fontsize=7,
                    )
            if gcf_x:
                ax.scatter(
                    gcf_x, gcf_y,
                    s=gcf_sizes, c=gcf_colors, zorder=5,
                    edgecolors="black", linewidths=0.5,
                )

            # Tip labels
            max_x = max(node_x.values()) if node_x else 1.0
            label_fontsize = config.ylabel_fontsize if config.ylabel_fontsize is not None else 9
            if label_fontsize > 0:
                draw_circular_tip_labels(ax, tree, coords, fontsize=label_fontsize, offset=max_x * 0.03)

            # Apply color annotations (range + label only; branches are trait-colored)
            if self.plot_config.color_file:
                color_data = parse_color_file(self.plot_config.color_file)
                for taxa_list, clr, lbl in color_data["ranges"]:
                    mrca = resolve_mrca(tree, taxa_list)
                    if mrca is not None:
                        draw_range_wedge(ax, tree, mrca, clr, coords)
                apply_label_colors(ax, color_data["labels"])
                color_legend = build_color_legend_handles(color_data)
                if color_legend:
                    ax.legend(handles=color_legend, loc="upper right", fontsize=8, frameon=True)

            # Colorbar for gCF
            sm = plt.cm.ScalarMappable(
                cmap=plt.cm.RdYlGn, norm=Normalize(vmin=0, vmax=1)
            )
            sm.set_array([])
            cbar = fig.colorbar(sm, ax=ax, pad=0.15)
            cbar.set_label("Gene Concordance Factor (gCF)")

            if config.show_title:
                ax.set_title(config.title or "Concordance-Aware ASR", fontsize=config.title_fontsize)
        else:
            # --- Rectangular mode ---
            # Draw branches
            horizontal_segments = []
            vertical_segments = []
            for clade in preorder_clades:
                if clade == root:
                    continue
                if id(clade) not in parent_map:
                    continue
                parent = parent_map[id(clade)]
                if id(parent) not in node_x or id(clade) not in node_x:
                    continue

                x0 = node_x[id(parent)]
                x1 = node_x[id(clade)]
                y0 = node_y.get(id(parent), 0)
                y1 = node_y.get(id(clade), 0)

                horizontal_segments.append(((x0, y1), (x1, y1)))
                vertical_segments.append(((x0, y0), (x0, y1)))

            if horizontal_segments:
                ax.add_collection(
                    LineCollection(
                        horizontal_segments,
                        colors="gray",
                        linewidths=2,
                        zorder=1,
                    ),
                    autolim=True,
                )
            if vertical_segments:
                ax.add_collection(
                    LineCollection(
                        vertical_segments,
                        colors="gray",
                        linewidths=2,
                        zorder=1,
                    ),
                    autolim=True,
                )
            if horizontal_segments or vertical_segments:
                ax.autoscale_view()

            # Draw gCF-sized dots at internal nodes
            gcf_x = []
            gcf_y = []
            gcf_sizes = []
            gcf_colors = []
            for clade in preorder_clades:
                if not clade.clades:
                    continue
                if id(clade) in gcf_values and id(clade) in node_x:
                    gcf = gcf_values[id(clade)]
                    size = 50 + 200 * gcf
                    color = plt.cm.RdYlGn(gcf)
                    gcf_x.append(node_x[id(clade)])
                    gcf_y.append(node_y.get(id(clade), 0))
                    gcf_sizes.append(size)
                    gcf_colors.append(color)
                    ax.annotate(
                        f"{gcf:.2f}",
                        (node_x[id(clade)], node_y.get(id(clade), 0)),
                        textcoords="offset points",
                        xytext=(5, 5),
                        fontsize=7,
                    )
            if gcf_x:
                ax.scatter(
                    gcf_x, gcf_y,
                    s=gcf_sizes, c=gcf_colors, zorder=5,
                    edgecolors="black", linewidths=0.5,
                )

            # Tip labels
            max_x = max(node_x.values()) if node_x else 0
            offset = max_x * 0.02
            label_fontsize = config.ylabel_fontsize if config.ylabel_fontsize is not None else 9
            if label_fontsize > 0:
                for tip in tips:
                    ax.text(
                        node_x[id(tip)] + offset, node_y[id(tip)],
                        tip.name, va="center", fontsize=label_fontsize,
                    )

            # Apply color annotations (range + label only; branches are trait-colored)
            if self.plot_config.color_file:
                color_data = parse_color_file(self.plot_config.color_file)
                for taxa_list, clr, lbl in color_data["ranges"]:
                    mrca = resolve_mrca(tree, taxa_list)
                    if mrca is not None:
                        draw_range_rect(ax, tree, mrca, clr, node_x, node_y)
                apply_label_colors(ax, color_data["labels"])
                color_legend = build_color_legend_handles(color_data)
                if color_legend:
                    ax.legend(handles=color_legend, loc="upper right", fontsize=8, frameon=True)

            # Colorbar for gCF
            sm = plt.cm.ScalarMappable(
                cmap=plt.cm.RdYlGn, norm=Normalize(vmin=0, vmax=1)
            )
            sm.set_array([])
            cbar = fig.colorbar(sm, ax=ax, pad=0.15)
            cbar.set_label("Gene Concordance Factor (gCF)")

            ax.set_xlabel("Branch length")
            ax.set_yticks([])
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.spines["left"].set_visible(False)
            if config.show_title:
                ax.set_title(config.title or "Concordance-Aware ASR", fontsize=config.title_fontsize)

        fig.tight_layout()
        fig.savefig(output_path, dpi=config.dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved concordance ASR plot: {output_path}")

    # ------------------------------------------------------------------
    # Uncertainty plot
    # ------------------------------------------------------------------

    def _collect_uncertainty_node_data(self, species_tree, result):
        """Collect uncertainty plot rows matched by descendant labels."""
        anc = result.get("ancestral_estimates", {})
        method = result.get("method")

        entries_by_desc = {}
        for entry in anc.values():
            desc = entry.get("descendants")
            if desc is not None:
                entries_by_desc.setdefault(tuple(desc), entry)

        node_data = []
        clade_tip_sets = self._collect_clade_tip_sets(species_tree)
        for clade in species_tree.find_clades(order="preorder"):
            if clade.is_terminal():
                continue
            desc = sorted(clade_tip_sets.get(id(clade), frozenset()))
            entry = entries_by_desc.get(tuple(desc))
            if entry is None:
                continue

            if method == "distribution":
                estimates = entry.get("gene_tree_estimates", [])
            else:
                estimates = entry.get("source_estimates", [])
            if len(estimates) < 2:
                continue

            if len(desc) <= 2:
                short = f"({', '.join(desc)})"
            else:
                short = f"({desc[0]}, ..., {desc[-1]})"
            node_data.append((
                short,
                estimates,
                entry.get("gcf", 1.0),
                entry.get("estimate", 0.0),
            ))

        return node_data

    def _plot_uncertainty(self, species_tree, result, output_path) -> None:
        """Violin + boxplot showing distribution of ancestral estimates.

        For the distribution method: per-gene-tree estimates.
        For the weighted method: concordant + discordant source estimates.
        """
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        except ImportError:
            print("matplotlib is required for plotting.")
            return

        anc = result.get("ancestral_estimates", {})
        if not anc:
            return

        config = self.plot_config

        # Collect data: (node_label, estimates_list, gcf, estimate)
        node_data = self._collect_uncertainty_node_data(species_tree, result)

        if not node_data:
            print("No nodes with enough estimates for uncertainty plot.")
            return

        # Sort by weighted/mean estimate
        node_data.sort(key=lambda x: x[3])

        labels = [d[0] for d in node_data]
        data = [d[1] for d in node_data]
        gcfs = [d[2] for d in node_data]

        n_nodes = len(node_data)
        fig_h = max(4, n_nodes * 0.5)
        fig_w = config.fig_width or 8
        if config.fig_height:
            fig_h = config.fig_height
        fig, ax = plt.subplots(figsize=(fig_w, fig_h))

        positions = list(range(n_nodes))

        # Color by gCF: high concordance = blue, low = red
        gcf_colors = []
        for g in gcfs:
            if g >= 0.7:
                gcf_colors.append("#2b8cbe")
            elif g >= 0.4:
                gcf_colors.append("#969696")
            else:
                gcf_colors.append("#d62728")

        # Violin plots (horizontal)
        vp = ax.violinplot(
            data, positions=positions, vert=False,
            showmeans=False, showmedians=False, showextrema=False,
        )
        for i, body in enumerate(vp["bodies"]):
            body.set_facecolor(gcf_colors[i])
            body.set_alpha(0.3)
            body.set_edgecolor(gcf_colors[i])

        # Boxplots
        bp = ax.boxplot(
            data, positions=positions, vert=False,
            widths=0.3, patch_artist=True, showfliers=True, zorder=3,
        )
        for i, box in enumerate(bp["boxes"]):
            box.set_facecolor(gcf_colors[i])
            box.set_alpha(0.6)
        for element in ["whiskers", "caps", "medians"]:
            for line in bp[element]:
                line.set_color("black")

        # Mean markers
        for i, d in enumerate(data):
            mean_val = sum(d) / len(d)
            ax.plot(
                mean_val, i, "D", color="white",
                markeredgecolor="black", markersize=5, zorder=4,
            )

        ax.set_yticks(positions)
        ax.set_yticklabels(labels, fontsize=7)
        xlabel = "Ancestral state estimate"
        ax.set_xlabel(xlabel)

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        # Legend for gCF colors
        from matplotlib.patches import Patch
        legend_handles = [
            Patch(facecolor="#2b8cbe", alpha=0.6, label="gCF >= 0.7"),
            Patch(facecolor="#969696", alpha=0.6, label="0.4 <= gCF < 0.7"),
            Patch(facecolor="#d62728", alpha=0.6, label="gCF < 0.4"),
        ]
        ax.legend(
            handles=legend_handles, loc="lower right",
            fontsize=7, title="Concordance", title_fontsize=8,
        )

        method_label = (
            "gene trees" if result["method"] == "distribution"
            else "concordance sources"
        )
        if config.show_title:
            ax.set_title(
                config.title
                or f"Ancestral state uncertainty across {method_label}",
                fontsize=config.title_fontsize,
            )

        fig.tight_layout()
        fig.savefig(output_path, dpi=config.dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved uncertainty plot: {output_path}")
