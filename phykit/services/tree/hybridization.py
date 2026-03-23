from io import StringIO
from pathlib import Path
from typing import Dict, List

from Bio import Phylo
from scipy.stats import binomtest

from .base import Tree
from ...helpers.json_output import print_json
from ...helpers.plot_config import PlotConfig, compute_node_x_cladogram
from ...helpers.circular_layout import (
    compute_circular_coords,
    draw_circular_branches,
    draw_circular_tip_labels,
    draw_circular_colored_branch,
    draw_circular_colored_arc,
    circular_branch_points,
    radial_offset,
)
from ...helpers.color_annotations import (
    parse_color_file,
    resolve_mrca,
    draw_range_rect,
    draw_range_wedge,
    build_color_legend_handles,
)
from ...errors import PhykitUserError


class Hybridization(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"])
        self.gene_trees_path = parsed["gene_trees_path"]
        self.support_threshold = parsed["support_threshold"]
        self.alpha = parsed["alpha"]
        self.json_output = parsed["json_output"]
        self.plot_output = parsed["plot_output"]
        self.plot_config = parsed["plot_config"]

    def process_args(self, args) -> Dict:
        return dict(
            tree_file_path=args.tree,
            gene_trees_path=args.gene_trees,
            support_threshold=getattr(args, "support", None),
            alpha=getattr(args, "alpha", 0.05),
            json_output=getattr(args, "json", False),
            plot_output=getattr(args, "plot_output", None),
            plot_config=PlotConfig.from_args(args),
        )

    def run(self) -> None:
        species_tree = self.read_tree_file()
        gene_trees = self._parse_gene_trees(self.gene_trees_path)

        topology_counts, shared_taxa = self._count_topologies(
            species_tree, gene_trees, self.support_threshold
        )

        # Test each branch and collect results
        branch_results = []
        for branch_key in sorted(topology_counts.keys()):
            data = topology_counts[branch_key]
            total = data["n_concordant"] + data["n_alt1"] + data["n_alt2"]
            gcf = data["n_concordant"] / total if total > 0 else 1.0

            test_result = self._test_asymmetry(data["n_alt1"], data["n_alt2"])

            entry = dict(
                taxa=data["split"],
                n_concordant=data["n_concordant"],
                n_alt1=data["n_alt1"],
                n_alt2=data["n_alt2"],
                gcf=round(gcf, 4),
            )
            entry.update(test_result)
            branch_results.append(entry)

        # FDR correction across testable p-values
        testable_indices = []
        testable_pvals = []
        for i, entry in enumerate(branch_results):
            if entry["p_value"] is not None:
                testable_indices.append(i)
                testable_pvals.append(entry["p_value"])

        fdr_corrected = self._fdr(testable_pvals)
        for idx, fdr_p in zip(testable_indices, fdr_corrected):
            branch_results[idx]["fdr_p"] = fdr_p

        # Set fdr_p to None for untestable branches
        for entry in branch_results:
            if "fdr_p" not in entry:
                entry["fdr_p"] = None

        # Compute hybrid score and significance per branch
        alpha = self.alpha
        for entry in branch_results:
            if entry["fdr_p"] is not None and entry["fdr_p"] < alpha:
                entry["significant"] = True
                entry["hybrid_score"] = entry["asymmetry_ratio"]
            else:
                entry["significant"] = False
                entry["hybrid_score"] = 0.0

        n_reticulations = sum(1 for e in branch_results if e["significant"])

        # Summary
        summary = dict(
            n_branches=len(branch_results),
            n_gene_trees=len(gene_trees),
            support_threshold=self.support_threshold,
            alpha=alpha,
            n_reticulations=n_reticulations,
        )

        # Output
        if self.json_output:
            self._output_json(branch_results, summary)
        else:
            self._output_text(branch_results, summary)

        if self.plot_output:
            if self.plot_config.ladderize:
                species_tree.ladderize()
            self._plot(species_tree, branch_results, self.plot_output,
                       shared_taxa=shared_taxa)

    # ------------------------------------------------------------------
    # Output methods
    # ------------------------------------------------------------------

    def _output_text(self, branch_results, summary) -> None:
        try:
            print("Hybridization Analysis")
            print("======================")
            print(f"Species tree branches: {summary['n_branches']}")
            print(f"Gene trees: {summary['n_gene_trees']}")
            if summary["support_threshold"] is not None:
                print(f"Support threshold: {summary['support_threshold']}")
            print()
            print(f"Estimated reticulation events: {summary['n_reticulations']}")
            print()

            sig_branches = [e for e in branch_results if e["significant"]]
            nonsig_count = len(branch_results) - len(sig_branches)

            if sig_branches:
                print(f"Significant branches (FDR < {summary['alpha']}):")
                header = (
                    f"  {'Branch (taxa)':<35}"
                    f"{'gCF':>8}"
                    f"{'gDF1':>8}"
                    f"{'gDF2':>8}"
                    f"{'Ratio':>8}"
                    f"{'FDR_p':>10}"
                    f"{'Favored':>12}"
                )
                print(header)
                for entry in sig_branches:
                    branch_label = ",".join(entry["taxa"])
                    total = entry["n_concordant"] + entry["n_alt1"] + entry["n_alt2"]
                    gdf1 = entry["n_alt1"] / total if total > 0 else 0.0
                    gdf2 = entry["n_alt2"] / total if total > 0 else 0.0
                    ratio_str = f"{entry['asymmetry_ratio']:.2f}" if entry["asymmetry_ratio"] is not None else "NA"
                    fdr_str = f"{entry['fdr_p']:.4f}" if entry["fdr_p"] is not None else "NA"
                    favored = entry.get("favored_alt", "-") or "-"
                    print(
                        f"  {branch_label:<35}"
                        f"{entry['gcf']:>8.2f}"
                        f"{gdf1:>8.2f}"
                        f"{gdf2:>8.2f}"
                        f"{ratio_str:>8}"
                        f"{fdr_str:>10}"
                        f"{favored:>12}"
                    )
                print()

            print(f"Non-significant branches: {nonsig_count}")
        except BrokenPipeError:
            pass

    def _output_json(self, branch_results, summary) -> None:
        result = dict(
            n_branches=summary["n_branches"],
            n_gene_trees=summary["n_gene_trees"],
            support_threshold=summary["support_threshold"],
            alpha=summary["alpha"],
            n_reticulations=summary["n_reticulations"],
            branches=branch_results,
        )
        print_json(result)

    def _plot(self, species_tree, branch_results, output_path,
              shared_taxa=None) -> None:
        """Phylogram colored by hybridization score at each branch."""
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.colors import Normalize
        import numpy as np

        # Build lookup from taxa label -> branch result
        branch_lookup = {}
        for entry in branch_results:
            key = ",".join(entry["taxa"])
            branch_lookup[key] = entry

        parent_map = self._build_parent_map(species_tree)
        tips = list(species_tree.get_terminals())
        all_taxa_fs = (shared_taxa if shared_taxa is not None
                       else frozenset(t.name for t in tips))

        # Compute node positions
        node_x = {}
        node_y = {}

        for i, tip in enumerate(tips):
            node_y[id(tip)] = i

        root = species_tree.root
        if self.plot_config.cladogram:
            node_x = compute_node_x_cladogram(species_tree, parent_map)
        else:
            for clade in species_tree.find_clades(order="preorder"):
                if clade == root:
                    node_x[id(clade)] = 0.0
                else:
                    if id(clade) in parent_map:
                        parent = parent_map[id(clade)]
                        t = clade.branch_length if clade.branch_length else 0.0
                        node_x[id(clade)] = node_x.get(id(parent), 0.0) + t

        for clade in species_tree.find_clades(order="postorder"):
            if not clade.is_terminal() and id(clade) not in node_y:
                child_ys = [
                    node_y[id(c)] for c in clade.clades if id(c) in node_y
                ]
                if child_ys:
                    node_y[id(clade)] = np.mean(child_ys)
                else:
                    node_y[id(clade)] = 0.0

        # Map internal nodes to their branch result
        node_to_result = {}
        for clade in species_tree.find_clades(order="preorder"):
            if clade.is_terminal():
                continue
            node_tips = frozenset(t.name for t in clade.get_terminals()) & all_taxa_fs
            split_label = (
                sorted(node_tips)
                if len(node_tips) <= len(all_taxa_fs) - len(node_tips)
                else sorted(all_taxa_fs - node_tips)
            )
            key = ",".join(split_label)
            if key in branch_lookup:
                node_to_result[id(clade)] = branch_lookup[key]

        # Color setup: gray for 0, warm gradient for significant branches
        # Use OrRd colormap for significant branches (orange to red)
        cmap = plt.cm.OrRd
        norm = Normalize(vmin=0.5, vmax=1.0)

        config = self.plot_config
        config.resolve(n_rows=len(tips), n_cols=None)
        fig, ax = plt.subplots(figsize=(config.fig_width, config.fig_height))

        if self.plot_config.circular:
            self._plot_circular(ax, fig, species_tree, root, tips, parent_map,
                                node_x, node_y, node_to_result, cmap, norm, config)
        else:
            self._plot_rectangular(ax, fig, species_tree, root, tips, parent_map,
                                   node_x, node_y, node_to_result, cmap, norm, config)

        fig.tight_layout()
        fig.savefig(output_path, dpi=config.dpi, bbox_inches="tight")
        import matplotlib.pyplot as plt2
        plt.close(fig)

    def _plot_circular(self, ax, fig, species_tree, root, tips, parent_map,
                       node_x, node_y, node_to_result, cmap, norm, config):
        import matplotlib.pyplot as plt
        import numpy as np

        coords = compute_circular_coords(species_tree, node_x, parent_map)
        ax.set_aspect("equal")
        ax.axis("off")

        # Draw each branch individually with its hybridization color
        for clade in species_tree.find_clades(order="preorder"):
            if clade == root:
                continue
            cid = id(clade)
            if cid not in parent_map:
                continue
            parent = parent_map[cid]
            if id(parent) not in coords or cid not in coords:
                continue

            color = "gray"
            lw = 2
            if cid in node_to_result:
                entry = node_to_result[cid]
                if entry["hybrid_score"] > 0:
                    color = cmap(norm(entry["hybrid_score"]))
                    lw = 3

            draw_circular_colored_branch(ax, coords[id(parent)], coords[cid], color, lw=lw)

        # Draw arcs at internal nodes
        for clade in species_tree.find_clades(order="preorder"):
            if clade.is_terminal() or not clade.clades:
                continue
            cid = id(clade)
            pc = coords[cid]
            child_angles = [coords[id(ch)]["angle"] for ch in clade.clades]
            if len(child_angles) < 2:
                continue
            start_a = min(child_angles)
            end_a = max(child_angles)
            import math
            span = (end_a - start_a) % (2.0 * math.pi)
            if span > math.pi:
                start_a, end_a = end_a, start_a

            arc_color = "gray"
            if cid in node_to_result:
                entry = node_to_result[cid]
                if entry["hybrid_score"] > 0:
                    arc_color = cmap(norm(entry["hybrid_score"]))
            draw_circular_colored_arc(ax, 0.0, 0.0, pc["radius"], start_a, end_a, arc_color, lw=1.5)

        # Tip labels
        max_x = max(node_x.values()) if node_x else 1.0
        label_fontsize = config.ylabel_fontsize if config.ylabel_fontsize is not None else 9
        if label_fontsize > 0:
            draw_circular_tip_labels(ax, species_tree, coords, fontsize=label_fontsize, offset=max_x * 0.02)

        # Color annotations
        if self.plot_config.color_file:
            color_data = parse_color_file(self.plot_config.color_file)
            for taxa_list, clr, lbl in color_data["ranges"]:
                mrca = resolve_mrca(species_tree, taxa_list)
                if mrca is not None:
                    draw_range_wedge(ax, species_tree, mrca, clr, coords)
            for taxon, lbl_color in color_data["labels"].items():
                for text_obj in ax.texts:
                    if text_obj.get_text() == taxon:
                        text_obj.set_color(lbl_color)
                        break
            color_legend = build_color_legend_handles(color_data)
            if color_legend:
                ax.legend(handles=color_legend, loc="upper right", fontsize=8, frameon=True)

        # Mark significant branches with stars
        for clade in species_tree.find_clades(order="preorder"):
            if clade.is_terminal():
                continue
            cid = id(clade)
            if cid not in node_to_result:
                continue
            entry = node_to_result[cid]
            if cid not in coords:
                continue
            cx = coords[cid]["x"]
            cy = coords[cid]["y"]

            if entry["significant"]:
                ax.scatter(cx, cy, s=100, c="red", marker="*", zorder=5)

        # Colorbar
        legend_loc = config.legend_position or "upper right"
        if legend_loc != "none":
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            cbar = fig.colorbar(sm, ax=ax, pad=0.15)
            cbar_fontsize = config.axis_fontsize if config.axis_fontsize else 10
            cbar.set_label("Hybridization score", fontsize=cbar_fontsize)

        if config.show_title:
            ax.set_title(config.title or "Hybridization Analysis", fontsize=config.title_fontsize)

    def _plot_rectangular(self, ax, fig, species_tree, root, tips, parent_map,
                          node_x, node_y, node_to_result, cmap, norm, config):
        import matplotlib.pyplot as plt
        import numpy as np

        # Draw branches
        for clade in species_tree.find_clades(order="preorder"):
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

            color = "gray"
            lw = 2
            if id(clade) in node_to_result:
                entry = node_to_result[id(clade)]
                if entry["hybrid_score"] > 0:
                    color = cmap(norm(entry["hybrid_score"]))
                    lw = 3

            ax.plot([x0, x1], [y1, y1], color=color, lw=lw)
            ax.plot([x0, x0], [y0, y1], color="gray", lw=1.5)

        # Mark significant branches with stars
        for clade in species_tree.find_clades(order="preorder"):
            if clade.is_terminal():
                continue
            if id(clade) not in node_to_result:
                continue
            entry = node_to_result[id(clade)]
            x = node_x.get(id(clade), 0)
            y = node_y.get(id(clade), 0)

            if entry["significant"]:
                ax.scatter(x, y, s=100, c="red", marker="*", zorder=5)

        # Tip labels
        label_fontsize = config.ylabel_fontsize if config.ylabel_fontsize is not None else 9
        if label_fontsize > 0:
            max_x = max(node_x.values()) if node_x else 0
            offset = max_x * 0.02
            for tip in tips:
                ax.text(
                    node_x[id(tip)] + offset, node_y[id(tip)],
                    tip.name, va="center", fontsize=label_fontsize,
                )

        # Color annotations
        if self.plot_config.color_file:
            color_data = parse_color_file(self.plot_config.color_file)
            for taxa_list, clr, lbl in color_data["ranges"]:
                mrca = resolve_mrca(species_tree, taxa_list)
                if mrca is not None:
                    draw_range_rect(ax, species_tree, mrca, clr, node_x, node_y)
            for taxon, lbl_color in color_data["labels"].items():
                for text_obj in ax.texts:
                    if text_obj.get_text() == taxon:
                        text_obj.set_color(lbl_color)
                        break
            color_legend = build_color_legend_handles(color_data)
            if color_legend:
                ax.legend(handles=color_legend, loc="upper right", fontsize=8, frameon=True)

        # Colorbar
        legend_loc = config.legend_position or "upper right"
        if legend_loc != "none":
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            cbar = fig.colorbar(sm, ax=ax, pad=0.15)
            cbar_fontsize = config.axis_fontsize if config.axis_fontsize else 10
            cbar.set_label("Hybridization score", fontsize=cbar_fontsize)

        ax.set_xlabel("Branch length (subs/site)")
        ax.set_yticks([])
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)
        if config.show_title:
            ax.set_title(config.title or "Hybridization Analysis", fontsize=config.title_fontsize)
        if config.axis_fontsize:
            ax.xaxis.label.set_fontsize(config.axis_fontsize)

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
    # Bipartition extraction and topology counting
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

    @staticmethod
    def _build_parent_map(tree) -> Dict:
        parent_map = {}
        for clade in tree.find_clades(order="preorder"):
            for child in clade.clades:
                parent_map[id(child)] = clade
        return parent_map

    @staticmethod
    def _collect_taxa(trees) -> set:
        taxa = set()
        for tree in trees:
            stack = [tree.root]
            while stack:
                node = stack.pop()
                if not node.clades:
                    taxa.add(node.name)
                else:
                    stack.extend(node.clades)
        return taxa

    def _extract_splits(self, tree, all_taxa_fs, support_threshold=None) -> set:
        """Extract canonical bipartitions from a single tree.

        If support_threshold is provided, skip any internal node whose
        clade.confidence is not None and is below the threshold.
        """
        splits = set()
        tip_sets = {}
        stack = [(tree.root, False)]
        while stack:
            node, children_done = stack[-1]
            if not node.clades:
                stack.pop()
                name = node.name
                tip_sets[id(node)] = frozenset((name,)) if name in all_taxa_fs else frozenset()
            elif children_done:
                stack.pop()
                merged = frozenset().union(*(tip_sets[id(c)] for c in node.clades))
                tip_sets[id(node)] = merged
                if len(merged) > 1 and merged != all_taxa_fs:
                    # Check support threshold
                    if support_threshold is not None:
                        conf = getattr(node, 'confidence', None)
                        if conf is not None and conf < support_threshold:
                            # Skip this bipartition (collapse it)
                            continue
                    splits.add(self._canonical_split(merged, all_taxa_fs))
            else:
                stack[-1] = (node, True)
                for child in reversed(node.clades):
                    stack.append((child, False))
        return splits

    def _get_four_groups(self, tree, node, parent_map, all_taxa_fs):
        if node.is_terminal() or len(node.clades) < 2:
            return None

        C1 = frozenset(t.name for t in node.clades[0].get_terminals()) & all_taxa_fs
        C2 = frozenset(t.name for t in node.clades[1].get_terminals()) & all_taxa_fs
        for extra_child in node.clades[2:]:
            C2 = C2 | (frozenset(t.name for t in extra_child.get_terminals()) & all_taxa_fs)

        parent = parent_map.get(id(node))
        if parent is None:
            return None

        siblings = [c for c in parent.clades if id(c) != id(node)]
        if not siblings:
            return None

        chosen_sib = None
        S = frozenset()
        for sib in siblings:
            candidate = frozenset(t.name for t in sib.get_terminals()) & all_taxa_fs
            if candidate:
                S = candidate
                chosen_sib = sib
                break

        if not S:
            return None

        D = all_taxa_fs - C1 - C2 - S

        if not D:
            if chosen_sib.is_terminal() or len(chosen_sib.clades) < 2:
                return None
            S = frozenset(t.name for t in chosen_sib.clades[0].get_terminals()) & all_taxa_fs
            D = frozenset(t.name for t in chosen_sib.clades[1].get_terminals()) & all_taxa_fs
            for extra in chosen_sib.clades[2:]:
                D = D | (frozenset(t.name for t in extra.get_terminals()) & all_taxa_fs)

        if not C1 or not C2 or not S:
            return None

        return C1, C2, S, D

    def _count_topologies(self, species_tree, gene_trees,
                          support_threshold=None) -> Dict:
        """Count concordant and two NNI-alternative topologies for each
        internal branch of the species tree across gene trees.

        If support_threshold is provided, gene tree bipartitions with
        support below the threshold are collapsed before counting.
        """
        species_taxa = set(t.name for t in species_tree.get_terminals())
        gene_taxa = self._collect_taxa(gene_trees)
        all_taxa_fs = frozenset(sorted(species_taxa & gene_taxa))
        parent_map = self._build_parent_map(species_tree)

        gene_tree_splits = []
        for gt in gene_trees:
            splits = self._extract_splits(gt, all_taxa_fs, support_threshold)
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

            concordant_bp = self._canonical_split(C1 | C2, all_taxa_fs)
            nni_alt1_bp = self._canonical_split(S | C2, all_taxa_fs)
            nni_alt2_bp = self._canonical_split(C1 | S, all_taxa_fs)

            n_concordant = sum(1 for splits in gene_tree_splits if concordant_bp in splits)
            n_alt1 = sum(1 for splits in gene_tree_splits if nni_alt1_bp in splits)
            n_alt2 = sum(1 for splits in gene_tree_splits if nni_alt2_bp in splits)

            node_tips = frozenset(t.name for t in clade.get_terminals()) & all_taxa_fs
            split_label = (
                sorted(node_tips)
                if len(node_tips) <= len(all_taxa_fs) - len(node_tips)
                else sorted(all_taxa_fs - node_tips)
            )
            branch_key = ",".join(split_label)
            result[branch_key] = dict(
                split=split_label,
                n_concordant=n_concordant,
                n_alt1=n_alt1,
                n_alt2=n_alt2,
            )
        return result, all_taxa_fs

    # ------------------------------------------------------------------
    # Statistical testing
    # ------------------------------------------------------------------

    def _test_asymmetry(self, n_alt1: int, n_alt2: int) -> Dict:
        total = n_alt1 + n_alt2
        if total == 0:
            return dict(
                asymmetry_ratio=None,
                p_value=None,
                favored_alt=None,
            )

        asymmetry_ratio = max(n_alt1, n_alt2) / total
        result = binomtest(n_alt1, total, p=0.5, alternative='two-sided')
        p_value = result.pvalue

        if n_alt1 > n_alt2:
            favored_alt = "alt1"
        elif n_alt2 > n_alt1:
            favored_alt = "alt2"
        else:
            favored_alt = None

        return dict(
            asymmetry_ratio=asymmetry_ratio,
            p_value=p_value,
            favored_alt=favored_alt,
        )

    @staticmethod
    def _fdr(p_values: List[float]) -> List[float]:
        """Benjamini-Hochberg FDR correction."""
        n = len(p_values)
        if n == 0:
            return []
        indexed = sorted(enumerate(p_values), key=lambda x: x[1])
        corrected = [0.0] * n
        prev = 1.0
        for rank_minus_1 in range(n - 1, -1, -1):
            orig_idx, p = indexed[rank_minus_1]
            rank = rank_minus_1 + 1
            adjusted = min(p * n / rank, prev)
            adjusted = min(adjusted, 1.0)
            corrected[orig_idx] = adjusted
            prev = adjusted
        return corrected
