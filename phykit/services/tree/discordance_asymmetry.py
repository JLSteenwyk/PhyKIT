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


class DiscordanceAsymmetry(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"])
        self.gene_trees_path = parsed["gene_trees_path"]
        self.verbose = parsed["verbose"]
        self.annotate = parsed["annotate"]
        self.json_output = parsed["json_output"]
        self.plot_output = parsed["plot_output"]
        self.plot_config = parsed["plot_config"]

    def process_args(self, args) -> Dict:
        return dict(
            tree_file_path=args.tree,
            gene_trees_path=args.gene_trees,
            verbose=getattr(args, "verbose", False),
            annotate=getattr(args, "annotate", False),
            json_output=getattr(args, "json", False),
            plot_output=getattr(args, "plot_output", None),
            plot_config=PlotConfig.from_args(args),
        )

    def run(self) -> None:
        species_tree = self.read_tree_file()
        gene_trees = self._parse_gene_trees(self.gene_trees_path)

        topology_counts, shared_taxa = self._count_topologies(species_tree, gene_trees)

        # Test each branch and collect results
        branch_results = []
        for branch_key in sorted(topology_counts.keys()):
            data = topology_counts[branch_key]
            test_result = self._test_asymmetry(data["n_alt1"], data["n_alt2"])
            entry = dict(
                split=data["split"],
                n_concordant=data["n_concordant"],
                n_alt1=data["n_alt1"],
                n_alt2=data["n_alt2"],
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

        # Summary
        summary = dict(
            n_gene_trees=len(gene_trees),
            n_branches_tested=len(testable_indices),
            n_significant_fdr05=sum(
                1 for entry in branch_results
                if entry["fdr_p"] is not None and entry["fdr_p"] < 0.05
            ),
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
            header = (
                f"{'branch':<30}"
                f"{'n_conc':>8}"
                f"{'n_alt1':>8}"
                f"{'n_alt2':>8}"
                f"{'asym_ratio':>12}"
                f"{'binom_p':>12}"
                f"{'fdr_p':>12}"
                f"{'gene_flow':>12}"
            )
            print(header)
            print("-" * len(header))

            for entry in branch_results:
                branch_label = ",".join(entry["split"])
                asym = (
                    f"{entry['asymmetry_ratio']:.3f}"
                    if entry["asymmetry_ratio"] is not None
                    else "NA"
                )
                binom_p = (
                    f"{entry['p_value']:.4f}"
                    if entry["p_value"] is not None
                    else "NA"
                )
                fdr_p = (
                    f"{entry['fdr_p']:.4f}"
                    if entry["fdr_p"] is not None
                    else "NA"
                )
                gene_flow = "-"
                if (entry["fdr_p"] is not None and entry["fdr_p"] < 0.05
                        and entry["favored_alt"] is not None):
                    gene_flow = entry["favored_alt"]
                print(
                    f"{branch_label:<30}"
                    f"{entry['n_concordant']:>8}"
                    f"{entry['n_alt1']:>8}"
                    f"{entry['n_alt2']:>8}"
                    f"{asym:>12}"
                    f"{binom_p:>12}"
                    f"{fdr_p:>12}"
                    f"{gene_flow:>12}"
                )

            print("---")
            print(
                f"Summary: {summary['n_branches_tested']} branches tested, "
                f"{summary['n_significant_fdr05']} significant (FDR<0.05)"
            )

            if self.verbose:
                print()
                for entry in branch_results:
                    branch_label = ",".join(entry["split"])
                    total = entry["n_concordant"] + entry["n_alt1"] + entry["n_alt2"]
                    gcf = entry["n_concordant"] / total if total > 0 else 1.0
                    print(f"Branch: {branch_label}")
                    print(
                        f"  gCF={gcf:.3f}  gDF1={entry['n_alt1']}/{total}  "
                        f"gDF2={entry['n_alt2']}/{total}"
                    )
        except BrokenPipeError:
            pass

    def _output_json(self, branch_results, summary) -> None:
        result = dict(
            branches=branch_results,
            summary=summary,
        )
        print_json(result)

    def _plot(self, species_tree, branch_results, output_path,
              shared_taxa=None) -> None:
        """Phylogram colored by asymmetry ratio at each branch."""
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.colors import Normalize
        import numpy as np

        # Build lookup from split label -> branch result
        branch_lookup = {}
        for entry in branch_results:
            key = ",".join(entry["split"])
            branch_lookup[key] = entry

        parent_map = self._build_parent_map(species_tree)
        tips = list(species_tree.get_terminals())
        # Use shared taxa (intersection with gene trees) for split label
        # matching so labels are consistent with _count_topologies.
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

        # Color setup: diverging from blue (0.5 = symmetric) to red (1.0 = asymmetric)
        cmap = plt.cm.RdYlBu_r
        norm = Normalize(vmin=0.5, vmax=1.0)

        config = self.plot_config
        config.resolve(n_rows=len(tips), n_cols=None)
        fig, ax = plt.subplots(figsize=(config.fig_width, config.fig_height))

        if self.plot_config.circular:
            # --- Circular mode ---
            coords = compute_circular_coords(species_tree, node_x, parent_map)
            ax.set_aspect("equal")
            ax.axis("off")

            # Draw each branch individually with its asymmetry color
            for clade in species_tree.find_clades(order="preorder"):
                if clade == root:
                    continue
                cid = id(clade)
                if cid not in parent_map:
                    continue
                parent = parent_map[cid]
                if id(parent) not in coords or cid not in coords:
                    continue

                # Determine color
                color = "gray"
                lw = 2
                if cid in node_to_result:
                    entry = node_to_result[cid]
                    if entry["asymmetry_ratio"] is not None:
                        color = cmap(norm(entry["asymmetry_ratio"]))
                        lw = 3

                # Draw radial segment
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

                # Color arc by the node's own asymmetry if available
                arc_color = "gray"
                if cid in node_to_result:
                    entry = node_to_result[cid]
                    if entry["asymmetry_ratio"] is not None:
                        arc_color = cmap(norm(entry["asymmetry_ratio"]))
                draw_circular_colored_arc(ax, 0.0, 0.0, pc["radius"], start_a, end_a, arc_color, lw=1.5)

            # Tip labels
            max_x = max(node_x.values()) if node_x else 1.0
            label_fontsize = config.ylabel_fontsize if config.ylabel_fontsize is not None else 9
            if label_fontsize > 0:
                draw_circular_tip_labels(ax, species_tree, coords, fontsize=label_fontsize, offset=max_x * 0.02)

            # Apply color annotations (range + label only; branches are trait-colored)
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

            # Annotate internal nodes
            for clade in species_tree.find_clades(order="preorder"):
                if clade.is_terminal():
                    continue
                cid = id(clade)
                if cid not in node_to_result:
                    continue
                entry = node_to_result[cid]
                cx = coords[cid]["x"]
                cy = coords[cid]["y"]

                total = entry["n_concordant"] + entry["n_alt1"] + entry["n_alt2"]
                gcf = entry["n_concordant"] / total if total > 0 else 1.0

                if self.annotate:
                    ax.annotate(
                        f"{gcf:.2f}",
                        (cx, cy),
                        textcoords="offset points",
                        xytext=(5, 5),
                        fontsize=max(4, 7 - n_tips * 0.01),
                    )

                if (entry["fdr_p"] is not None and entry["fdr_p"] < 0.05
                        and entry["favored_alt"] is not None):
                    ax.scatter(cx, cy, s=100, c="red", marker="*", zorder=5)

            # Colorbar (hide if legend-position is none)
            legend_loc = config.legend_position or "upper right"
            if legend_loc != "none":
                sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
                sm.set_array([])
                cbar = fig.colorbar(sm, ax=ax, pad=0.15)
                cbar_fontsize = config.axis_fontsize if config.axis_fontsize else 10
                cbar.set_label("Asymmetry ratio", fontsize=cbar_fontsize)

            if config.show_title:
                ax.set_title(config.title or "Discordance Asymmetry", fontsize=config.title_fontsize)

            fig.tight_layout()
            fig.savefig(output_path, dpi=config.dpi, bbox_inches="tight")
            plt.close(fig)

        else:
            # --- Rectangular mode ---
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

                # Color the horizontal branch by asymmetry ratio if this is an internal node
                color = "gray"
                lw = 2
                if id(clade) in node_to_result:
                    entry = node_to_result[id(clade)]
                    if entry["asymmetry_ratio"] is not None:
                        color = cmap(norm(entry["asymmetry_ratio"]))
                        lw = 3

                ax.plot([x0, x1], [y1, y1], color=color, lw=lw)
                ax.plot([x0, x0], [y0, y1], color="gray", lw=1.5)

            # Annotate internal nodes
            for clade in species_tree.find_clades(order="preorder"):
                if clade.is_terminal():
                    continue
                if id(clade) not in node_to_result:
                    continue
                entry = node_to_result[id(clade)]
                x = node_x.get(id(clade), 0)
                y = node_y.get(id(clade), 0)

                # Show gCF value (only if --annotate)
                total = entry["n_concordant"] + entry["n_alt1"] + entry["n_alt2"]
                gcf = entry["n_concordant"] / total if total > 0 else 1.0

                if self.annotate:
                    ax.annotate(
                        f"{gcf:.2f}",
                        (x, y),
                        textcoords="offset points",
                        xytext=(5, 5),
                        fontsize=max(4, 7 - n_tips * 0.01),
                    )

                # Mark significant branches (FDR < 0.05)
                if (entry["fdr_p"] is not None and entry["fdr_p"] < 0.05
                        and entry["favored_alt"] is not None):
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

            # Apply color annotations (range + label only; branches are trait-colored)
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

            # Colorbar (hide if legend-position is none)
            legend_loc = config.legend_position or "upper right"
            if legend_loc != "none":
                sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
                sm.set_array([])
                cbar = fig.colorbar(sm, ax=ax, pad=0.15)
                cbar_fontsize = config.axis_fontsize if config.axis_fontsize else 10
                cbar.set_label("Asymmetry ratio", fontsize=cbar_fontsize)

            ax.set_xlabel("Branch length (subs/site)")
            ax.set_yticks([])
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.spines["left"].set_visible(False)
            if config.show_title:
                ax.set_title(config.title or "Discordance Asymmetry", fontsize=config.title_fontsize)
            if config.axis_fontsize:
                ax.xaxis.label.set_fontsize(config.axis_fontsize)
            fig.tight_layout()
            fig.savefig(output_path, dpi=config.dpi, bbox_inches="tight")
            plt.close(fig)

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
        """Normalize a bipartition to canonical form.

        Returns the smaller side as a frozenset; ties are broken
        lexicographically.
        """
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
        """Build a dict mapping child id -> parent clade."""
        parent_map = {}
        for clade in tree.find_clades(order="preorder"):
            for child in clade.clades:
                parent_map[id(child)] = clade
        return parent_map

    @staticmethod
    def _collect_taxa(trees) -> set:
        """Collect the union of all tip names across a list of trees.

        Uses a direct stack walk over .clades to avoid the overhead of
        Bio.Phylo's find_clades / match_attrs dispatch layer.
        """
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

    def _extract_splits(self, tree, all_taxa_fs) -> set:
        """Extract canonical bipartitions from a single tree.

        Performs a single postorder traversal, building tip sets
        bottom-up via .clades, and canonicalises each non-trivial split.
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
                    splits.add(self._canonical_split(merged, all_taxa_fs))
            else:
                stack[-1] = (node, True)
                for child in reversed(node.clades):
                    stack.append((child, False))
        return splits

    def _get_four_groups(self, tree, node, parent_map, all_taxa_fs):
        """Identify the four subtree groups around an internal branch.

        For the branch connecting *node* to its parent:
          C1 = tips of node's first child
          C2 = tips of node's second child (extra children merged for polytomies)
          S  = tips of node's sibling under parent
          D  = remaining tips (everything above parent)

        Returns (C1, C2, S, D) as frozensets, or None if decomposition
        is not possible (e.g., node is root, leaf, or has <2 children).
        """
        if node.is_terminal() or len(node.clades) < 2:
            return None

        C1 = frozenset(t.name for t in node.clades[0].get_terminals()) & all_taxa_fs
        C2 = frozenset(t.name for t in node.clades[1].get_terminals()) & all_taxa_fs
        # If node has >2 children (polytomy), merge extras into C2
        for extra_child in node.clades[2:]:
            C2 = C2 | (frozenset(t.name for t in extra_child.get_terminals()) & all_taxa_fs)

        parent = parent_map.get(id(node))
        if parent is None:
            # node is root — no branch above it
            return None

        # Get siblings of node under parent; pick the first sibling
        # whose shared taxa are non-empty (avoids skipping valid branches
        # when the first sibling's taxa are all absent from gene trees).
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

        # D = everything else (other siblings + above parent)
        D = all_taxa_fs - C1 - C2 - S

        # For the root branch (D empty), the sibling encompasses the
        # entire other side of the tree.  Decompose the sibling into its
        # children to obtain proper 4-subtree NNI groups.
        if not D:
            if chosen_sib.is_terminal() or len(chosen_sib.clades) < 2:
                return None
            S = frozenset(t.name for t in chosen_sib.clades[0].get_terminals()) & all_taxa_fs
            D = frozenset(t.name for t in chosen_sib.clades[1].get_terminals()) & all_taxa_fs
            for extra in chosen_sib.clades[2:]:
                D = D | (frozenset(t.name for t in extra.get_terminals()) & all_taxa_fs)

        # Skip if any group is empty (branch is degenerate after taxon filtering)
        if not C1 or not C2 or not S:
            return None

        return C1, C2, S, D

    def _count_topologies(self, species_tree, gene_trees) -> Dict:
        """Count concordant and two NNI-alternative topologies for each
        internal branch of the species tree across gene trees.

        Returns a dict keyed by branch label (comma-joined sorted taxa
        in the smaller partition side) with:
          split: list of sorted taxon names
          n_concordant: int
          n_alt1: int
          n_alt2: int
        """
        species_taxa = set(t.name for t in species_tree.get_terminals())
        gene_taxa = self._collect_taxa(gene_trees)
        all_taxa_fs = frozenset(sorted(species_taxa & gene_taxa))
        parent_map = self._build_parent_map(species_tree)

        # Extract bipartitions from all gene trees (topology only, no lengths).
        # Builds tip sets bottom-up in a single postorder pass per gene tree,
        # avoiding repeated get_terminals() calls (O(n) vs O(n²) per tree).
        gene_tree_splits = []
        for gt in gene_trees:
            splits = self._extract_splits(gt, all_taxa_fs)
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
        """Run a two-sided binomial test on the two NNI alternatives.

        Returns a dict with asymmetry_ratio, p_value, and favored_alt.
        """
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
