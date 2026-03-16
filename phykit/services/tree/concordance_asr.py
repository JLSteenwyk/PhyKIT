import copy
import math
import pickle
from io import StringIO
from pathlib import Path
from typing import Dict, Tuple

import numpy as np
from Bio import Phylo

from .base import Tree
from .ancestral_reconstruction import AncestralReconstruction
from ...helpers.json_output import print_json
from ...helpers.plot_config import PlotConfig, compute_node_x_cladogram
from ...errors import PhykitUserError


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
        self.missing_taxa = parsed["missing_taxa"]
        self.json_output = parsed["json_output"]
        self.plot_config = parsed["plot_config"]

    def process_args(self, args) -> Dict:
        return dict(
            tree_file_path=args.tree,
            gene_trees_path=args.gene_trees,
            trait_data_path=args.trait_data,
            trait_column=getattr(args, "trait", None),
            method=getattr(args, "method", "weighted"),
            ci=getattr(args, "ci", False),
            plot_output=getattr(args, "plot", None),
            missing_taxa=getattr(args, "missing_taxa", "shared"),
            json_output=getattr(args, "json", False),
            plot_config=PlotConfig.from_args(args),
        )

    def run(self) -> None:
        species_tree = self.read_tree_file()

        # Create a helper ASR instance for reusing methods
        self._asr = AncestralReconstruction.__new__(AncestralReconstruction)
        self._asr.ci = True  # always compute CIs internally

        self._asr._validate_tree(species_tree)

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
        species_tips = set(t.name for t in species_tree.get_terminals())
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
        species_copy = copy.deepcopy(species_tree)
        sp_tips = [t.name for t in species_copy.get_terminals()]
        tips_to_prune = [t for t in sp_tips if t not in trait_values]
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

        if self.json_output:
            print_json(result)
        else:
            self._print_text_output(result)

    # ------------------------------------------------------------------
    # Gene tree parsing
    # ------------------------------------------------------------------

    def _parse_gene_trees(self, path: str) -> list:
        try:
            lines = Path(path).read_text().splitlines()
        except FileNotFoundError:
            raise PhykitUserError(
                [
                    f"{path} corresponds to no such file or directory.",
                    "Please check filename and pathing",
                ],
                code=2,
            )

        cleaned = [l.strip() for l in lines if l.strip() and not l.strip().startswith("#")]
        trees = []
        for line in cleaned:
            if line.startswith("("):
                trees.append(Phylo.read(StringIO(line), "newick"))
            else:
                tree_path = Path(path).parent / line
                trees.append(Phylo.read(str(tree_path), "newick"))
        return trees

    # ------------------------------------------------------------------
    # Taxa normalization
    # ------------------------------------------------------------------

    def _normalize_taxa(self, species_tree, gene_trees, species_tips):
        gene_tip_sets = []
        for gt in gene_trees:
            gene_tip_sets.append(set(t.name for t in gt.get_terminals()))

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
            t.name for t in species_tree.get_terminals()
            if t.name not in shared
        ]
        if sp_tips_to_prune:
            species_tree = self.prune_tree_using_taxa_list(
                species_tree, sp_tips_to_prune
            )

        # Prune gene trees
        pruned_gene_trees = []
        for gt in gene_trees:
            gt_tips_to_prune = [
                t.name for t in gt.get_terminals()
                if t.name not in shared
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
        complement = all_taxa - taxa_side
        if len(taxa_side) < len(complement):
            return frozenset(taxa_side)
        elif len(taxa_side) > len(complement):
            return frozenset(complement)
        else:
            return min(frozenset(taxa_side), frozenset(complement),
                       key=lambda s: sorted(s))

    def _get_four_groups(self, tree, node, parent_map, all_taxa_fs):
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

        C1 = frozenset(t.name for t in node.clades[0].get_terminals())
        C2 = frozenset(t.name for t in node.clades[1].get_terminals())
        # If node has >2 children (polytomy), merge extras into C2
        for extra_child in node.clades[2:]:
            C2 = C2 | frozenset(t.name for t in extra_child.get_terminals())

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

        S = frozenset(t.name for t in siblings[0].get_terminals())
        # D = everything else (other siblings + above parent)
        D = all_taxa_fs - C1 - C2 - S

        return C1, C2, S, D

    def _compute_gcf_per_node(
        self, species_tree, gene_trees, all_taxa
    ) -> Dict[int, Tuple[float, float, float]]:
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

        # Extract bipartitions from all gene trees
        gene_tree_splits = []
        for gt in gene_trees:
            splits = set()
            for clade in gt.get_nonterminals():
                tips = frozenset(t.name for t in clade.get_terminals())
                if len(tips) <= 1 or tips == all_taxa_fs:
                    continue
                splits.add(self._canonical_split(tips, all_taxa_fs))
            gene_tree_splits.append(splits)

        result = {}
        for clade in species_tree.find_clades(order="preorder"):
            if clade.is_terminal():
                continue

            groups = self._get_four_groups(
                species_tree, clade, parent_map, all_taxa_fs
            )
            if groups is None:
                continue

            C1, C2, S, D = groups

            # The three quartet resolutions as bipartitions
            concordant_bp = self._canonical_split(C1 | C2, all_taxa_fs)
            nni_alt1_bp = self._canonical_split(S | C2, all_taxa_fs)
            nni_alt2_bp = self._canonical_split(C1 | S, all_taxa_fs)

            concordant = sum(
                1 for splits in gene_tree_splits if concordant_bp in splits
            )
            disc1 = sum(
                1 for splits in gene_tree_splits if nni_alt1_bp in splits
            )
            disc2 = sum(
                1 for splits in gene_tree_splits if nni_alt2_bp in splits
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
        tree_copy = copy.deepcopy(tree)

        # Prune to trait taxa
        tip_names = [t.name for t in tree_copy.get_terminals()]
        tips_to_prune = [t for t in tip_names if t not in trait_values]
        if tips_to_prune:
            tree_copy = self.prune_tree_using_taxa_list(tree_copy, tips_to_prune)

        remaining_tips = [t.name for t in tree_copy.get_terminals()]
        if len(remaining_tips) < 3:
            return {}, {}, 0.0

        ordered_names = sorted(trait_values.keys() & set(remaining_tips))
        if len(ordered_names) < 3:
            return {}, {}, 0.0

        x = np.array([trait_values[name] for name in ordered_names])
        node_labels = self._asr._label_internal_nodes(tree_copy)

        node_estimates, node_cis, sigma2, _ = self._asr._fast_anc(
            tree_copy, x, ordered_names, node_labels
        )

        # Re-key by descendant frozensets
        estimates_by_desc = {}
        cis_by_desc = {}
        for clade in tree_copy.find_clades(order="preorder"):
            if clade.is_terminal():
                continue
            if id(clade) not in node_labels:
                continue
            label = node_labels[id(clade)]
            if label in node_estimates:
                desc = frozenset(t.name for t in clade.get_terminals())
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
        weights = np.array(weights, dtype=float)
        means = np.array(means, dtype=float)
        variances = np.array(variances, dtype=float)

        if weights.sum() == 0:
            return 0.0, 0.0, 0.0

        weights = weights / weights.sum()

        # Within-group variance: E[Var(X|T)]
        within_var = float(np.sum(weights * variances))

        # Between-group variance: Var(E[X|T])
        weighted_mean = float(np.sum(weights * means))
        between_var = float(np.sum(weights * (means - weighted_mean) ** 2))

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
            desc = frozenset(t.name for t in clade.get_terminals())
            descendants = sorted(t.name for t in clade.get_terminals())

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
            desc = frozenset(t.name for t in clade.get_terminals())
            descendants = sorted(t.name for t in clade.get_terminals())

            gcf, gdf1, gdf2 = gcf_per_node.get(id(clade), (1.0, 0.0, 0.0))

            # Collect estimates from gene trees that contain this bipartition
            estimates_list = []
            for est_dict, ci_dict, sig in gene_tree_results:
                if desc in est_dict:
                    estimates_list.append(est_dict[desc])

            if not estimates_list:
                continue

            mean_est = float(np.mean(estimates_list))
            entry = {
                "estimate": mean_est,
                "descendants": descendants,
                "is_root": id(clade) == root_id,
                "gcf": float(gcf),
                "gdf1": float(gdf1),
                "gdf2": float(gdf2),
                "n_gene_trees_with_node": len(estimates_list),
                "var_topology": float(np.var(estimates_list)),
                "var_parameter": 0.0,
            }

            if self.ci and len(estimates_list) >= 2:
                ci_lower = float(np.percentile(estimates_list, 2.5))
                ci_upper = float(np.percentile(estimates_list, 97.5))
                entry["ci_lower"] = ci_lower
                entry["ci_upper"] = ci_upper

            ancestral_estimates[label] = entry

        # Average sigma2 across gene trees
        sigma2_values = [sig for _, _, sig in gene_tree_results if sig > 0]
        avg_sigma2 = float(np.mean(sigma2_values)) if sigma2_values else 0.0

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
        print("Concordance-Aware Ancestral State Reconstruction")
        print(f"\nMethod: {result['method']}")
        print(f"Number of tips: {result['n_tips']}")
        print(f"Number of gene trees: {result['n_gene_trees']}")
        print(f"Sigma-squared (BM rate): {result['sigma2']:.6f}")

        print("\nAncestral estimates:")
        estimates = result["ancestral_estimates"]

        has_ci = any("ci_lower" in e for e in estimates.values())

        if has_ci:
            print(
                f"  {'Node':<12s}{'Desc':>6s}{'Estimate':>12s}"
                f"{'gCF':>8s}{'95% CI':>24s}"
                f"{'Var_topo':>12s}{'Var_param':>12s}"
            )
        else:
            print(
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
                print(
                    f"  {label + root_tag:<12s}{n_desc:>6d}{est:>12.4f}"
                    f"{gcf:>8.3f}{ci_str:>24s}"
                    f"{var_t:>12.6f}{var_p:>12.6f}"
                )
            elif has_ci:
                print(
                    f"  {label + root_tag:<12s}{n_desc:>6d}{est:>12.4f}"
                    f"{gcf:>8.3f}{'':>24s}"
                    f"{var_t:>12.6f}{var_p:>12.6f}"
                )
            else:
                print(
                    f"  {label + root_tag:<12s}{n_desc:>6d}{est:>12.4f}"
                    f"{gcf:>8.3f}{var_t:>12.6f}{var_p:>12.6f}"
                )

    # ------------------------------------------------------------------
    # Plotting
    # ------------------------------------------------------------------

    def _plot_concordance_contmap(self, tree, result, output_path) -> None:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            from matplotlib.colors import Normalize
        except ImportError:
            print(
                "matplotlib is required for concordance contMap plotting. "
                "Install matplotlib and retry."
            )
            raise SystemExit(2)

        estimates = result["ancestral_estimates"]
        node_labels = self._asr._label_internal_nodes(tree)
        parent_map = self._asr._build_parent_map(tree)

        # Map estimates back to tree nodes
        label_to_est = {}
        label_to_gcf = {}
        for label, entry in estimates.items():
            label_to_est[label] = entry["estimate"]
            label_to_gcf[label] = entry.get("gcf", 1.0)

        # Build estimates dict keyed by id(clade)
        all_estimates = {}
        gcf_values = {}
        for clade in tree.find_clades(order="preorder"):
            if clade.is_terminal():
                # Use trait values if available in result
                pass
            else:
                if id(clade) in node_labels:
                    label = node_labels[id(clade)]
                    if label in label_to_est:
                        all_estimates[id(clade)] = label_to_est[label]
                        gcf_values[id(clade)] = label_to_gcf[label]

        tips = list(tree.get_terminals())
        node_x = {}
        node_y = {}

        for i, tip in enumerate(tips):
            node_y[id(tip)] = i

        root = tree.root
        if self.plot_config.cladogram:
            node_x = compute_node_x_cladogram(tree, parent_map)
        else:
            for clade in tree.find_clades(order="preorder"):
                if clade == root:
                    node_x[id(clade)] = 0.0
                else:
                    if id(clade) in parent_map:
                        parent = parent_map[id(clade)]
                        t = clade.branch_length if clade.branch_length else 0.0
                        node_x[id(clade)] = node_x.get(id(parent), 0.0) + t

        for clade in tree.find_clades(order="postorder"):
            if not clade.is_terminal() and id(clade) not in node_y:
                child_ys = [
                    node_y[id(c)] for c in clade.clades if id(c) in node_y
                ]
                if child_ys:
                    node_y[id(clade)] = np.mean(child_ys)
                else:
                    node_y[id(clade)] = 0.0

        config = self.plot_config
        config.resolve(n_rows=len(tips), n_cols=None)
        fig, ax = plt.subplots(figsize=(config.fig_width, config.fig_height))

        # Draw branches
        for clade in tree.find_clades(order="preorder"):
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

            ax.plot([x0, x1], [y1, y1], color="gray", lw=2)
            ax.plot([x0, x0], [y0, y1], color="gray", lw=2)

        # Draw gCF-sized dots at internal nodes
        for clade in tree.find_clades(order="preorder"):
            if clade.is_terminal():
                continue
            if id(clade) in gcf_values and id(clade) in node_x:
                gcf = gcf_values[id(clade)]
                size = 50 + 200 * gcf
                color = plt.cm.RdYlGn(gcf)
                ax.scatter(
                    node_x[id(clade)], node_y.get(id(clade), 0),
                    s=size, c=[color], zorder=5, edgecolors="black", linewidths=0.5
                )
                ax.annotate(
                    f"{gcf:.2f}",
                    (node_x[id(clade)], node_y.get(id(clade), 0)),
                    textcoords="offset points",
                    xytext=(5, 5),
                    fontsize=7,
                )

        # Tip labels
        max_x = max(node_x.values()) if node_x else 0
        offset = max_x * 0.02
        for tip in tips:
            ax.text(
                node_x[id(tip)] + offset, node_y[id(tip)],
                tip.name, va="center", fontsize=9,
            )

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
