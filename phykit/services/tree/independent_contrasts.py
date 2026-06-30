"""
Phylogenetically independent contrasts (PIC).

Computes Felsenstein's (1985) phylogenetically independent contrasts
for a continuous trait on a phylogeny. Each internal node yields one
contrast (standardized difference), producing n-1 contrasts for n tips.

Cross-validated against R's ape::pic().
"""
from __future__ import annotations

import math

from .base import Tree
from ...errors import PhykitUserError


def _json_dumps(*args, **kwargs):
    import json

    return json.dumps(*args, **kwargs)


class _LazyNumpy:
    def __getattr__(self, name):
        import numpy as _np

        return getattr(_np, name)


np = _LazyNumpy()


class IndependentContrasts(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"])
        self.trait_data_path = parsed["trait_data_path"]
        self.json_output = parsed["json_output"]

    def run(self) -> None:
        tree = self.read_tree_file_unmodified()
        copied_tree = False
        if self._needs_default_branch_lengths(tree):
            tree = self._fast_copy(tree)
            copied_tree = True
        self.validate_tree(
            tree,
            min_tips=3,
            assign_default_branch_length=1e-8,
            context="independent contrasts",
        )

        tree_tips = self.get_tip_names_from_tree(tree)
        tip_traits = self._parse_trait_data(self.trait_data_path, tree_tips)

        # Prune tree to shared taxa
        if len(tip_traits) != len(tree_tips):
            shared = set(tip_traits.keys())
            tips_to_prune = [t for t in tree_tips if t not in shared]
            if tips_to_prune:
                if not copied_tree:
                    tree = self._fast_copy(tree)
                    copied_tree = True
                tree = self.prune_tree_using_taxa_list(tree, tips_to_prune)

        # Resolve multifurcations (PIC requires fully dichotomous tree)
        has_polytomies = self._has_polytomies(tree)
        if has_polytomies is not False:
            if not copied_tree:
                tree = self._fast_copy(tree)
            self._resolve_polytomies(tree)

        if self.json_output:
            contrasts, node_labels = self._compute_pic(tree, tip_traits)
            self._print_json(contrasts, node_labels, tip_traits)
        elif self._should_use_text_summary_path(tree):
            contrasts, node_label_summaries = self._compute_pic_text_summaries(
                tree, tip_traits
            )
            self._print_text_summaries(contrasts, node_label_summaries)
        else:
            contrasts, node_labels = self._compute_pic(tree, tip_traits)
            self._print_text(contrasts, node_labels)

    def process_args(self, args) -> dict:
        return dict(
            tree_file_path=args.tree,
            trait_data_path=args.trait_data,
            json_output=getattr(args, "json", False),
        )

    @staticmethod
    def _needs_default_branch_lengths(tree) -> bool:
        try:
            root = tree.root
            root.clades
        except AttributeError:
            return True

        stack = [root]
        pop = stack.pop
        extend = stack.extend
        try:
            while stack:
                clade = pop()
                children = clade.clades
                if clade is not root and clade.branch_length is None:
                    return True
                if children:
                    extend(children)
        except AttributeError:
            return True

        return False

    @staticmethod
    def _has_polytomies(tree):
        try:
            root = tree.root
            root.clades
        except AttributeError:
            return None

        stack = [root]
        pop = stack.pop
        extend = stack.extend
        try:
            while stack:
                clade = pop()
                children = clade.clades
                if len(children) > 2:
                    return True
                if children:
                    extend(children)
        except AttributeError:
            return None

        return False

    def _parse_trait_data(
        self, path: str, tree_tips: list[str]
    ) -> dict[str, float]:
        """Parse two-column trait file (taxon<tab>value)."""
        try:
            with open(path) as f:
                traits = {}
                for line in f:
                    stripped = line.strip()
                    if not stripped or stripped.startswith("#"):
                        continue
                    parts = stripped.split("\t")
                    if len(parts) != 2:
                        continue
                    try:
                        traits[parts[0]] = float(parts[1])
                    except ValueError:
                        continue
        except FileNotFoundError:
            raise PhykitUserError(
                [f"{path} not found. Check filename and path."], code=2
            )

        if len(tree_tips) >= 3 and len(tree_tips) == len(traits):
            for tip, taxon in zip(tree_tips, traits):
                if tip != taxon:
                    break
            else:
                return traits

        tip_traits = {tip: traits[tip] for tip in tree_tips if tip in traits}
        if len(tip_traits) < 3:
            raise PhykitUserError(
                [
                    f"Only {len(tip_traits)} shared taxa between tree and trait file.",
                    "At least 3 shared taxa are required.",
                ],
                code=2,
            )

        return tip_traits

    def _resolve_polytomies(self, tree) -> None:
        """Resolve multifurcations by adding zero-length branches."""
        try:
            root = tree.root
            root.clades
        except AttributeError:
            newick_clade = None
            for clade in tree.find_clades(order="postorder"):
                while len(clade.clades) > 2:
                    if newick_clade is None:
                        from Bio.Phylo import Newick

                        newick_clade = Newick.Clade
                    child1 = clade.clades.pop()
                    child2 = clade.clades.pop()
                    new_internal = newick_clade(branch_length=0.0)
                    new_internal.clades = [child1, child2]
                    clade.clades.append(new_internal)
            return
        else:
            stack = [root]
            newick_clade = None
            try:
                pop = stack.pop
                extend = stack.extend
                while stack:
                    clade = pop()
                    children = clade.clades
                    if children:
                        extend(children)
                    while len(children) > 2:
                        # Take the last two children and merge them
                        if newick_clade is None:
                            from Bio.Phylo import Newick

                            newick_clade = Newick.Clade
                        child1 = children.pop()
                        child2 = children.pop()
                        new_internal = newick_clade(branch_length=0.0)
                        new_internal.clades = [child1, child2]
                        children.append(new_internal)
            except AttributeError:
                self._resolve_polytomies_fallback(tree)

    @staticmethod
    def _resolve_polytomies_fallback(tree) -> None:
        newick_clade = None
        for clade in tree.find_clades(order="postorder"):
            while len(clade.clades) > 2:
                if newick_clade is None:
                    from Bio.Phylo import Newick

                    newick_clade = Newick.Clade
                child1 = clade.clades.pop()
                child2 = clade.clades.pop()
                new_internal = newick_clade(branch_length=0.0)
                new_internal.clades = [child1, child2]
                clade.clades.append(new_internal)

    def _compute_pic(
        self, tree, tip_traits: dict[str, float]
    ) -> tuple[list[float], list[list[str]]]:
        """Compute phylogenetically independent contrasts.

        Felsenstein (1985) algorithm:
        - Postorder traversal through tree
        - At each internal node with children L, R:
          contrast = (x_L - x_R) / sqrt(v_L + v_R)
          x_node = (x_L/v_L + x_R/v_R) / (1/v_L + 1/v_R)
          v_node += v_L * v_R / (v_L + v_R)

        Returns (contrasts, node_tip_labels) where node_tip_labels[i]
        is the list of tip names descending from the node that produced
        contrast[i].
        """
        result = self._compute_pic_standard_tree(tree, tip_traits)
        if result is not None:
            return result

        # Store trait values and effective branch lengths per node
        node_val = {}
        node_bl = {}  # effective branch length (adjusted during traversal)
        node_labels_by_id = {}

        # Initialize tips
        for tip in tree.get_terminals():
            if tip.name in tip_traits:
                node_val[id(tip)] = tip_traits[tip.name]
                node_bl[id(tip)] = tip.branch_length if tip.branch_length else 1e-8
                node_labels_by_id[id(tip)] = [tip.name]

        contrasts = []
        node_labels = []
        merge_sorted_labels = self._merge_sorted_labels

        for clade in tree.find_clades(order="postorder"):
            if clade.is_terminal():
                continue
            if len(clade.clades) != 2:
                continue

            left, right = clade.clades
            lid, rid = id(left), id(right)

            if lid not in node_val or rid not in node_val:
                continue

            x_l = node_val[lid]
            x_r = node_val[rid]
            v_l = node_bl[lid]
            v_r = node_bl[rid]

            # Standardized contrast
            contrast = (x_l - x_r) / math.sqrt(v_l + v_r)
            contrasts.append(float(contrast))

            # Tip labels for this node
            tips_here = merge_sorted_labels(
                node_labels_by_id[lid],
                node_labels_by_id[rid],
            )
            node_labels.append(tips_here)

            # Weighted average trait value for this node
            node_val[id(clade)] = (x_l / v_l + x_r / v_r) / (1.0 / v_l + 1.0 / v_r)

            # Adjusted branch length
            parent_bl = clade.branch_length if clade.branch_length else 0.0
            node_bl[id(clade)] = parent_bl + (v_l * v_r) / (v_l + v_r)
            node_labels_by_id[id(clade)] = tips_here

        return contrasts, node_labels

    @staticmethod
    def _merge_sorted_labels(left: list[str], right: list[str]) -> list[str]:
        if not left:
            return list(right)
        if not right:
            return list(left)
        if left[-1] <= right[0]:
            return left + right
        if right[-1] <= left[0]:
            return right + left

        merged = []
        append = merged.append
        i = j = 0
        left_len = len(left)
        right_len = len(right)
        while i < left_len and j < right_len:
            if left[i] <= right[j]:
                append(left[i])
                i += 1
            else:
                append(right[j])
                j += 1
        if i < left_len:
            merged.extend(left[i:])
        if j < right_len:
            merged.extend(right[j:])
        return merged

    @staticmethod
    def _postorder_clades_fast(tree):
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

    def _compute_pic_standard_tree(
        self, tree, tip_traits: dict[str, float]
    ):
        clades = self._postorder_clades_fast(tree)
        if clades is None:
            return None

        node_val = {}
        node_bl = {}
        node_labels_by_id = {}
        contrasts = []
        node_labels = []
        merge_sorted_labels = self._merge_sorted_labels

        for clade in clades:
            clade_id = id(clade)
            children = clade.clades
            if not children:
                if clade.name in tip_traits:
                    node_val[clade_id] = tip_traits[clade.name]
                    node_bl[clade_id] = (
                        clade.branch_length if clade.branch_length else 1e-8
                    )
                    node_labels_by_id[clade_id] = [clade.name]
                continue

            if len(children) != 2:
                continue

            left, right = children
            lid, rid = id(left), id(right)

            if lid not in node_val or rid not in node_val:
                continue

            x_l = node_val[lid]
            x_r = node_val[rid]
            v_l = node_bl[lid]
            v_r = node_bl[rid]

            contrast = (x_l - x_r) / math.sqrt(v_l + v_r)
            contrasts.append(float(contrast))

            tips_here = merge_sorted_labels(
                node_labels_by_id[lid],
                node_labels_by_id[rid],
            )
            node_labels.append(tips_here)

            node_val[clade_id] = (x_l / v_l + x_r / v_r) / (
                1.0 / v_l + 1.0 / v_r
            )

            parent_bl = clade.branch_length if clade.branch_length else 0.0
            node_bl[clade_id] = parent_bl + (v_l * v_r) / (v_l + v_r)
            node_labels_by_id[clade_id] = tips_here

        return contrasts, node_labels

    def _compute_pic_text_summaries(
        self, tree, tip_traits: dict[str, float]
    ) -> tuple[list[float], list[tuple[list[str], int]]]:
        result = self._compute_pic_text_summaries_standard_tree(tree, tip_traits)
        if result is not None:
            return result

        contrasts, node_labels = self._compute_pic(tree, tip_traits)
        return contrasts, [(tips[:3], len(tips)) for tips in node_labels]

    def _compute_pic_text_summaries_standard_tree(
        self, tree, tip_traits: dict[str, float]
    ):
        clades = self._postorder_clades_fast(tree)
        if clades is None:
            return None

        node_val = {}
        node_bl = {}
        node_label_summaries_by_id = {}
        contrasts = []
        node_label_summaries = []
        merge_limited_sorted_labels = self._merge_limited_sorted_labels

        for clade in clades:
            clade_id = id(clade)
            children = clade.clades
            if not children:
                if clade.name in tip_traits:
                    node_val[clade_id] = tip_traits[clade.name]
                    node_bl[clade_id] = (
                        clade.branch_length if clade.branch_length else 1e-8
                    )
                    node_label_summaries_by_id[clade_id] = ([clade.name], 1)
                continue

            if len(children) != 2:
                continue

            left, right = children
            lid, rid = id(left), id(right)

            if lid not in node_val or rid not in node_val:
                continue

            x_l = node_val[lid]
            x_r = node_val[rid]
            v_l = node_bl[lid]
            v_r = node_bl[rid]

            contrast = (x_l - x_r) / math.sqrt(v_l + v_r)
            contrasts.append(float(contrast))

            left_labels, left_count = node_label_summaries_by_id[lid]
            right_labels, right_count = node_label_summaries_by_id[rid]
            tips_here = merge_limited_sorted_labels(
                left_labels,
                left_count,
                right_labels,
                right_count,
                3,
            )
            node_label_summaries.append(tips_here)

            node_val[clade_id] = (x_l / v_l + x_r / v_r) / (
                1.0 / v_l + 1.0 / v_r
            )

            parent_bl = clade.branch_length if clade.branch_length else 0.0
            node_bl[clade_id] = parent_bl + (v_l * v_r) / (v_l + v_r)
            node_label_summaries_by_id[clade_id] = tips_here

        return contrasts, node_label_summaries

    @staticmethod
    def _should_use_text_summary_path(tree) -> bool:
        try:
            root = tree.root
            root.clades
        except AttributeError:
            return False

        max_nodes_to_probe = 4096
        stack = [(root, 0)]
        nodes_seen = 0
        try:
            while stack and nodes_seen < max_nodes_to_probe:
                clade, depth = stack.pop()
                nodes_seen += 1
                if depth > 64:
                    return True
                children = clade.clades
                child_depth = depth + 1
                for child in children:
                    stack.append((child, child_depth))
        except AttributeError:
            return False

        return False

    @staticmethod
    def _merge_limited_sorted_labels(
        left: list[str],
        left_count: int,
        right: list[str],
        right_count: int,
        limit: int,
    ) -> tuple[list[str], int]:
        total_count = left_count + right_count
        if limit <= 0:
            return [], total_count
        if not left:
            return list(right[:limit]), total_count
        if not right:
            return list(left[:limit]), total_count
        if left[-1] <= right[0]:
            if len(left) >= limit:
                return list(left[:limit]), total_count
            return left + right[: limit - len(left)], total_count
        if right[-1] <= left[0]:
            if len(right) >= limit:
                return list(right[:limit]), total_count
            return right + left[: limit - len(right)], total_count

        merged = []
        append = merged.append
        i = j = 0
        left_len = len(left)
        right_len = len(right)
        while len(merged) < limit and i < left_len and j < right_len:
            if left[i] <= right[j]:
                append(left[i])
                i += 1
            else:
                append(right[j])
                j += 1
        while len(merged) < limit and i < left_len:
            append(left[i])
            i += 1
        while len(merged) < limit and j < right_len:
            append(right[j])
            j += 1
        return merged, total_count

    def _print_text(self, contrasts, node_labels):
        lines = [
            f"Number of contrasts: {len(contrasts)}",
            "",
            f"{'Node':<6}{'Contrast':>12}  Tips",
            "-" * 60,
        ]
        for i, (c, tips) in enumerate(zip(contrasts, node_labels), 1):
            tips_str = ", ".join(tips[:3])
            if len(tips) > 3:
                tips_str += f", ... ({len(tips)} total)"
            lines.append(f"{i:<6}{c:>12.6f}  {tips_str}")
        lines.extend(
            [
                "",
                f"Mean absolute contrast: {np.mean(np.abs(contrasts)):.6f}",
                f"Variance of contrasts:  {np.var(contrasts, ddof=1):.6f}",
            ]
        )
        print("\n".join(lines))

    def _print_text_summaries(self, contrasts, node_label_summaries):
        lines = [
            f"Number of contrasts: {len(contrasts)}",
            "",
            f"{'Node':<6}{'Contrast':>12}  Tips",
            "-" * 60,
        ]
        for i, (c, (tip_preview, tip_count)) in enumerate(
            zip(contrasts, node_label_summaries), 1
        ):
            tips_str = ", ".join(tip_preview)
            if tip_count > 3:
                tips_str += f", ... ({tip_count} total)"
            lines.append(f"{i:<6}{c:>12.6f}  {tips_str}")
        lines.extend(
            [
                "",
                f"Mean absolute contrast: {np.mean(np.abs(contrasts)):.6f}",
                f"Variance of contrasts:  {np.var(contrasts, ddof=1):.6f}",
            ]
        )
        print("\n".join(lines))

    def _print_json(self, contrasts, node_labels, tip_traits):
        contrast_values = np.asarray(contrasts, dtype=float)
        nodes = [
            {
                "node": i + 1,
                "contrast": round(c, 6),
                "tips": tips,
            }
            for i, (c, tips) in enumerate(zip(contrasts, node_labels))
        ]
        payload = {
            "n_taxa": len(tip_traits),
            "n_contrasts": len(contrasts),
            "contrasts": nodes,
            "mean_absolute_contrast": round(float(np.mean(np.abs(contrast_values))), 6),
            "variance_of_contrasts": round(float(np.var(contrast_values, ddof=1)), 6),
        }
        try:
            print(_json_dumps(payload, sort_keys=True))
        except BrokenPipeError:
            pass
