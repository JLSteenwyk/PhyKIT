import copy
import math
import sys
from typing import Dict, List, Tuple

import numpy as np

from .base import Tree
from ...helpers.json_output import print_json
from ...helpers.plot_config import PlotConfig, compute_node_x_cladogram
from ...errors import PhykitUserError


class ContMap(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"])
        self.trait_data_path = parsed["trait_data_path"]
        self.output_path = parsed["output_path"]
        self.json_output = parsed["json_output"]
        self.plot_config = parsed["plot_config"]

    def run(self) -> None:
        tree = self.read_tree_file()
        self._validate_tree(tree)

        tree_tips = self.get_tip_names_from_tree(tree)
        trait_values = self._parse_single_trait_data(
            self.trait_data_path, tree_tips
        )

        ordered_names = sorted(trait_values.keys())
        n = len(ordered_names)
        x = np.array([trait_values[name] for name in ordered_names])

        # Prune tree to shared taxa
        tree_copy = copy.deepcopy(tree)
        tip_names_in_tree = [t.name for t in tree_copy.get_terminals()]
        tips_to_prune = [t for t in tip_names_in_tree if t not in trait_values]
        if tips_to_prune:
            tree_copy = self.prune_tree_using_taxa_list(tree_copy, tips_to_prune)

        if self.plot_config.ladderize:
            tree_copy.ladderize()

        # Label internal nodes
        node_labels = self._label_internal_nodes(tree_copy)

        # Run fast ancestral reconstruction
        node_estimates, sigma2 = self._fast_anc(
            tree_copy, x, ordered_names, node_labels
        )

        # Always produce contMap plot
        self._plot_contmap(
            tree_copy, node_estimates, node_labels,
            trait_values, "trait", self.output_path
        )

        # Text output (always printed)
        print("Continuous Trait Map (contMap)")
        print(f"\nNumber of tips: {n}")
        print(f"Sigma-squared (BM rate): {sigma2:.4f}")
        print(f"Saved contMap plot: {self.output_path}")

        if self.json_output:
            result = {
                "n_tips": n,
                "sigma2": float(sigma2),
                "tip_values": {k: float(v) for k, v in sorted(trait_values.items())},
                "ancestral_estimates": {
                    label: float(node_estimates[label])
                    for label in sorted(node_estimates.keys())
                },
                "plot_output": self.output_path,
            }
            print_json(result)

    def process_args(self, args) -> Dict:
        return dict(
            tree_file_path=args.tree,
            trait_data_path=args.trait_data,
            output_path=args.output,
            json_output=getattr(args, "json", False),
            plot_config=PlotConfig.from_args(args),
        )

    def _validate_tree(self, tree) -> None:
        tips = list(tree.get_terminals())
        if len(tips) < 3:
            raise PhykitUserError(
                ["Tree must have at least 3 tips for contMap visualization."],
                code=2,
            )
        for clade in tree.find_clades():
            if clade.branch_length is None and clade != tree.root:
                raise PhykitUserError(
                    ["All branches in the tree must have lengths."],
                    code=2,
                )

    def _parse_single_trait_data(
        self, path: str, tree_tips: List[str]
    ) -> Dict[str, float]:
        try:
            with open(path) as f:
                lines = f.readlines()
        except FileNotFoundError:
            raise PhykitUserError(
                [
                    f"{path} corresponds to no such file or directory.",
                    "Please check filename and pathing",
                ],
                code=2,
            )

        traits = {}
        for line_num, line in enumerate(lines, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) != 2:
                raise PhykitUserError(
                    [
                        f"Line {line_num} in trait file has {len(parts)} columns; expected 2.",
                        "Each line should be: taxon_name<tab>trait_value",
                    ],
                    code=2,
                )
            taxon, value_str = parts
            try:
                traits[taxon] = float(value_str)
            except ValueError:
                raise PhykitUserError(
                    [
                        f"Non-numeric trait value '{value_str}' for taxon '{taxon}' on line {line_num}.",
                    ],
                    code=2,
                )

        tree_tip_set = set(tree_tips)
        trait_taxa_set = set(traits.keys())
        shared = tree_tip_set & trait_taxa_set

        tree_only = tree_tip_set - trait_taxa_set
        trait_only = trait_taxa_set - tree_tip_set

        if tree_only:
            print(
                f"Warning: {len(tree_only)} taxa in tree but not in trait file: "
                f"{', '.join(sorted(tree_only))}",
                file=sys.stderr,
            )
        if trait_only:
            print(
                f"Warning: {len(trait_only)} taxa in trait file but not in tree: "
                f"{', '.join(sorted(trait_only))}",
                file=sys.stderr,
            )

        if len(shared) < 3:
            raise PhykitUserError(
                [
                    f"Only {len(shared)} shared taxa between tree and trait file.",
                    "At least 3 shared taxa are required.",
                ],
                code=2,
            )

        return {taxon: traits[taxon] for taxon in shared}

    def _label_internal_nodes(self, tree) -> Dict:
        """Assign N1, N2, ... labels to unnamed internal nodes.

        Returns dict mapping id(clade) -> label for all internal nodes.
        """
        labels = {}
        counter = 1
        for clade in tree.find_clades(order="preorder"):
            if not clade.is_terminal():
                if clade.name:
                    labels[id(clade)] = clade.name
                else:
                    labels[id(clade)] = f"N{counter}"
                    counter += 1
        return labels

    def _build_parent_map(self, tree) -> Dict:
        parent_map = {}
        for clade in tree.find_clades(order="preorder"):
            for child in clade.clades:
                parent_map[id(child)] = clade
        return parent_map

    def _fast_anc(
        self, tree, x: np.ndarray, ordered_names: List[str],
        node_labels: Dict,
    ) -> Tuple[Dict, float]:
        """Two-pass ML ancestral reconstruction (Felsenstein's algorithm).

        Pass 1 (postorder): compute subtree estimates and variances.
        Pass 2 (preorder): incorporate information from the rest of the
        tree so that every internal node gets the full-tree ML estimate.

        Returns (node_estimates, sigma2).
        """
        name_to_idx = {name: i for i, name in enumerate(ordered_names)}
        parent_map = self._build_parent_map(tree)

        # --- Pass 1: postorder (subtree estimates) ---
        est_down = {}
        var_down = {}

        for clade in tree.find_clades(order="postorder"):
            if clade.is_terminal():
                name = clade.name
                if name in name_to_idx:
                    est_down[id(clade)] = x[name_to_idx[name]]
                    var_down[id(clade)] = 0.0
            else:
                children = clade.clades
                precisions = []
                weighted_vals = []
                for child in children:
                    child_id = id(child)
                    if child_id not in est_down:
                        continue
                    v_i = child.branch_length if child.branch_length else 0.0
                    child_var = var_down[child_id]
                    denom = v_i + child_var
                    if denom == 0:
                        denom = 1e-10
                    prec = 1.0 / denom
                    precisions.append(prec)
                    weighted_vals.append(prec * est_down[child_id])

                if precisions:
                    total_prec = sum(precisions)
                    est_down[id(clade)] = sum(weighted_vals) / total_prec
                    var_down[id(clade)] = 1.0 / total_prec

        # --- Pass 2: preorder (incorporate full-tree information) ---
        est_up = {}
        var_up = {}
        est_final = {}
        var_final = {}

        root = tree.root
        # Root has no information from above
        est_final[id(root)] = est_down[id(root)]
        var_final[id(root)] = var_down[id(root)]

        for clade in tree.find_clades(order="preorder"):
            if clade == root:
                continue
            if id(clade) not in parent_map:
                continue
            parent = parent_map[id(clade)]
            bl = clade.branch_length if clade.branch_length else 0.0

            # Gather precision from parent's "above" message
            prec_above = 0.0
            est_above_weighted = 0.0
            if id(parent) in est_up:
                v_up_parent = var_up[id(parent)]
                if v_up_parent > 0:
                    p = 1.0 / v_up_parent
                    prec_above += p
                    est_above_weighted += p * est_up[id(parent)]
                elif v_up_parent == 0 and id(parent) != id(root):
                    p = 1.0 / 1e-10
                    prec_above += p
                    est_above_weighted += p * est_up[id(parent)]

            # Gather precision from siblings
            for sibling in parent.clades:
                if id(sibling) == id(clade):
                    continue
                if id(sibling) not in est_down:
                    continue
                sib_bl = sibling.branch_length if sibling.branch_length else 0.0
                sib_var = var_down[id(sibling)] + sib_bl
                if sib_var == 0:
                    sib_var = 1e-10
                p = 1.0 / sib_var
                prec_above += p
                est_above_weighted += p * est_down[id(sibling)]

            # The "above" message for this node
            if prec_above > 0:
                est_up[id(clade)] = est_above_weighted / prec_above
                var_up[id(clade)] = bl + 1.0 / prec_above
            else:
                # No information from above (shouldn't happen in a valid tree)
                est_up[id(clade)] = est_down.get(id(clade), 0.0)
                var_up[id(clade)] = bl + 1e10

            # Combine down and up for final estimate
            if id(clade) in est_down:
                vd = var_down[id(clade)]
                vu = var_up[id(clade)]
                if vd == 0:
                    vd = 1e-10
                if vu == 0:
                    vu = 1e-10
                prec_d = 1.0 / vd
                prec_u = 1.0 / vu
                total = prec_d + prec_u
                est_final[id(clade)] = (
                    prec_d * est_down[id(clade)] + prec_u * est_up[id(clade)]
                ) / total
                var_final[id(clade)] = 1.0 / total

        # Compute sigma^2 from phylogenetic independent contrasts
        sigma2 = self._compute_sigma2_from_contrasts(
            tree, x, est_down, var_down, ordered_names
        )

        # Build output dicts keyed by node label
        node_estimates = {}
        for clade in tree.find_clades(order="preorder"):
            if not clade.is_terminal() and id(clade) in node_labels:
                label = node_labels[id(clade)]
                if id(clade) in est_final:
                    node_estimates[label] = est_final[id(clade)]

        return node_estimates, sigma2

    def _compute_sigma2_from_contrasts(
        self, tree, x: np.ndarray, node_estimates: Dict,
        node_variances: Dict, ordered_names: List[str],
    ) -> float:
        """Compute BM rate (sigma^2) from phylogenetic independent contrasts."""
        contrasts_sq = []
        for clade in tree.find_clades(order="postorder"):
            if clade.is_terminal():
                continue
            children = clade.clades
            valid_children = [
                c for c in children if id(c) in node_estimates
            ]
            if len(valid_children) < 2:
                continue

            # Process first two children
            c1, c2 = valid_children[0], valid_children[1]
            v1 = (c1.branch_length or 0.0) + node_variances[id(c1)]
            v2 = (c2.branch_length or 0.0) + node_variances[id(c2)]
            denom = v1 + v2
            if denom > 0:
                contrast = (
                    node_estimates[id(c1)] - node_estimates[id(c2)]
                ) / math.sqrt(denom)
                contrasts_sq.append(contrast ** 2)

            # Combined estimate and variance after first two children
            if denom > 0:
                combined_est = (
                    v2 * node_estimates[id(c1)] + v1 * node_estimates[id(c2)]
                ) / denom
                combined_var = (v1 * v2) / denom
            else:
                combined_est = node_estimates[id(c1)]
                combined_var = 0.0

            # Process additional children (for polytomies)
            for k in range(2, len(valid_children)):
                ck = valid_children[k]
                vk = (ck.branch_length or 0.0) + node_variances[id(ck)]
                vc = combined_var
                denom_k = vk + vc
                if denom_k > 0:
                    contrast_k = (
                        node_estimates[id(ck)] - combined_est
                    ) / math.sqrt(denom_k)
                    contrasts_sq.append(contrast_k ** 2)
                    combined_est = (
                        vc * node_estimates[id(ck)] + vk * combined_est
                    ) / denom_k
                    combined_var = (vk * vc) / denom_k

        if contrasts_sq:
            return sum(contrasts_sq) / len(contrasts_sq)
        return 0.0

    def _plot_contmap(
        self, tree, node_estimates, node_labels, trait_values,
        trait_name, output_path,
    ) -> None:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            from matplotlib.colors import Normalize
        except ImportError:
            print(
                "matplotlib is required for contMap plotting. "
                "Install matplotlib and retry."
            )
            raise SystemExit(2)

        parent_map = self._build_parent_map(tree)

        # Build estimates dict keyed by id(clade) for all nodes
        all_estimates = {}
        for clade in tree.find_clades(order="preorder"):
            if clade.is_terminal():
                if clade.name in trait_values:
                    all_estimates[id(clade)] = trait_values[clade.name]
            else:
                if id(clade) in node_labels:
                    label = node_labels[id(clade)]
                    if label in node_estimates:
                        all_estimates[id(clade)] = node_estimates[label]

        tips = list(tree.get_terminals())
        node_x = {}
        node_y = {}

        # Assign tip y-positions (evenly spaced)
        for i, tip in enumerate(tips):
            node_y[id(tip)] = i

        # Assign x-positions via preorder traversal
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
                        node_x[id(clade)] = node_x[id(parent)] + t

        # Assign internal y-positions (mean of children)
        for clade in tree.find_clades(order="postorder"):
            if not clade.is_terminal() and id(clade) not in node_y:
                child_ys = [
                    node_y[id(c)] for c in clade.clades if id(c) in node_y
                ]
                if child_ys:
                    node_y[id(clade)] = np.mean(child_ys)
                else:
                    node_y[id(clade)] = 0.0

        # Normalize trait values for colormap
        all_vals = list(all_estimates.values())
        if not all_vals:
            return
        vmin, vmax = min(all_vals), max(all_vals)
        if vmin == vmax:
            vmax = vmin + 1.0
        norm = Normalize(vmin=vmin, vmax=vmax)
        cmap = plt.get_cmap("coolwarm")

        config = self.plot_config
        config.resolve(n_rows=len(tips), n_cols=None)
        fig, ax = plt.subplots(figsize=(config.fig_width, config.fig_height))

        n_seg = 50

        for clade in tree.find_clades(order="preorder"):
            if clade == root:
                continue
            if id(clade) not in parent_map:
                continue

            parent = parent_map[id(clade)]
            if id(parent) not in node_x or id(clade) not in node_x:
                continue
            if id(parent) not in all_estimates or id(clade) not in all_estimates:
                continue

            x0 = node_x[id(parent)]
            x1 = node_x[id(clade)]
            y0 = node_y[id(parent)]
            y1 = node_y[id(clade)]
            val0 = all_estimates[id(parent)]
            val1 = all_estimates[id(clade)]

            # Horizontal segment (at child's y): colored gradient
            for s in range(n_seg):
                frac_s = s / n_seg
                frac_e = (s + 1) / n_seg
                xs = x0 + frac_s * (x1 - x0)
                xe = x0 + frac_e * (x1 - x0)
                v = val0 + ((frac_s + frac_e) / 2) * (val1 - val0)
                ax.plot(
                    [xs, xe], [y1, y1],
                    color=cmap(norm(v)), lw=3, solid_capstyle="butt",
                )

            # Vertical connector (at parent's x): parent's trait color
            ax.plot(
                [x0, x0], [y0, y1],
                color=cmap(norm(val0)), lw=3, solid_capstyle="butt",
            )

        # Tip labels
        max_x = max(node_x.values()) if node_x else 0
        offset = max_x * 0.02
        for tip in tips:
            ax.text(
                node_x[id(tip)] + offset, node_y[id(tip)],
                tip.name, va="center", fontsize=9,
            )

        # Colorbar
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax, pad=0.15)
        cbar.set_label(trait_name)

        ax.set_xlabel("Branch length")
        ax.set_yticks([])
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)
        if config.show_title:
            ax.set_title(config.title or "Continuous Trait Map (contMap)", fontsize=config.title_fontsize)
        fig.tight_layout()
        fig.savefig(output_path, dpi=config.dpi, bbox_inches="tight")
        plt.close(fig)
