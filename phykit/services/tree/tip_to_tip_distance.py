from __future__ import annotations

import sys
import itertools
import math


from .base import Tree


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)

class _LazyNumpy:
    _module = None

    def __getattr__(self, name):
        module = self._module
        if module is None:
            import numpy as _np

            module = _np
            self._module = module
        value = getattr(module, name)
        setattr(self, name, value)
        return value


class _LazyTreeMixin:
    def find_any(self, *args, **kwargs):
        from Bio.Phylo.BaseTree import TreeMixin as _TreeMixin

        find_any = _TreeMixin.find_any
        self.find_any = find_any
        return find_any(*args, **kwargs)

    def distance(self, *args, **kwargs):
        from Bio.Phylo.BaseTree import TreeMixin as _TreeMixin

        distance = _TreeMixin.distance
        self.distance = distance
        return distance(*args, **kwargs)


np = _LazyNumpy()
TreeMixin = _LazyTreeMixin()


class TipToTipDistance(Tree):
    _MATRIX_FAST_FILL_MIN_ROWS = 200_000
    _JSON_NO_COMBO_MIN_TIPS = 512
    _TEXT_NO_COMBO_MIN_TIPS = 1024

    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(
            tree_file_path=parsed["tree_file_path"],
            tip_1=parsed["tip_1"],
            tip_2=parsed["tip_2"],
        )
        self.json_output = parsed["json_output"]
        self.all_pairs = parsed["all_pairs"]
        self.plot = parsed["plot"]
        self.plot_output = parsed["plot_output"]
        self.plot_config = parsed["plot_config"]

    def run(self):
        tree_zero = self.read_tree_file_unmodified()

        if self.plot and not self.all_pairs:
            print("--plot requires --all-pairs for tip_to_tip_distance.")
            sys.exit(2)

        if self.all_pairs:
            if not self.json_output and not self.plot:
                output = self._format_all_pairwise_distances_fast(tree_zero)
                if output is not None:
                    if output:
                        print(output)
                    return

            rows = self.calculate_all_pairwise_distances(tree_zero)
            if self.plot:
                self._plot_tip_distance_heatmap(rows)

            if self.json_output:
                payload = dict(all_pairs=True, rows=rows, pairs=rows)
                if self.plot:
                    payload["plot_output"] = self.plot_output
                print_json(payload)
                return

            output = "\n".join(
                f"{row['taxon_a']}\t{row['taxon_b']}\t{row['tip_to_tip_distance']}"
                for row in rows
            )
            if output:
                print(output)
            if self.plot:
                print(f"Saved tip-to-tip heatmap: {self.plot_output}")
            return

        if not self.tip_1 or not self.tip_2:
            print("tip_1 and tip_2 are required unless --all-pairs is provided.\nExiting...")
            sys.exit(2)

        distance = round(
            self.calculate_tip_to_tip_distance(tree_zero, self.tip_1, self.tip_2),
            4,
        )
        if self.json_output:
            print_json(
                dict(
                    taxon_a=self.tip_1,
                    taxon_b=self.tip_2,
                    tip_to_tip_distance=distance,
                )
            )
            return
        print(distance)

    def check_leaves(
        self,
        tree_zero: Newick.Tree,
        tip_1: str,
        tip_2: str,
    ) -> None:
        leaf1 = TreeMixin.find_any(tree_zero, tip_1)
        if not bool(leaf1):
            print(tip_1, "not on tree\nExiting...")
            sys.exit(2)
        leaf2 = TreeMixin.find_any(tree_zero, tip_2)
        if not bool(leaf2):
            print(tip_2, "not on tree\nExiting...")
            sys.exit(2)

    def calculate_tip_to_tip_distance(
        self,
        tree_zero: Newick.Tree,
        tip_1: str,
        tip_2: str,
    ) -> float:
        distance = self._calculate_tip_to_tip_distance_fast(tree_zero, tip_1, tip_2)
        if distance is not None:
            return distance

        self.check_leaves(tree_zero, tip_1, tip_2)
        return TreeMixin.distance(tree_zero, tip_1, tip_2)

    @staticmethod
    def _calculate_tip_to_tip_distance_fast(
        tree_zero: Newick.Tree,
        tip_1: str,
        tip_2: str,
    ):
        try:
            root = tree_zero.root
            if tip_1 == tip_2:
                return TipToTipDistance._calculate_same_tip_distance_fast(
                    root,
                    tip_1,
                )

            parent_map = {}
            depth_by_node = {root: 0.0}
            found = {}
            targets = {tip_1, tip_2}
            stack = [root]
            pop = stack.pop
            append = stack.append

            while stack and len(found) < len(targets):
                clade = pop()
                children = clade.clades
                if children:
                    parent_depth = depth_by_node[clade]
                    child_count = len(children)
                    if child_count == 2:
                        child = children[1]
                        parent_map[child] = clade
                        depth_by_node[child] = parent_depth + (
                            child.branch_length or 0.0
                        )
                        append(child)
                        child = children[0]
                        parent_map[child] = clade
                        depth_by_node[child] = parent_depth + (
                            child.branch_length or 0.0
                        )
                        append(child)
                    else:
                        for child in reversed(children):
                            parent_map[child] = clade
                            depth_by_node[child] = parent_depth + (
                                child.branch_length or 0.0
                            )
                            append(child)
                elif clade.name in targets:
                    found.setdefault(clade.name, clade)

            leaf1 = found.get(tip_1)
            if leaf1 is None:
                print(tip_1, "not on tree\nExiting...")
                sys.exit(2)

            leaf2 = found.get(tip_2)
            if leaf2 is None:
                print(tip_2, "not on tree\nExiting...")
                sys.exit(2)

            if leaf1 is leaf2:
                return 0.0

            ancestors = set()
            current = leaf1
            while True:
                ancestors.add(current)
                if current is root:
                    break
                current = parent_map[current]

            current = leaf2
            while current not in ancestors:
                current = parent_map[current]

            lca = current
            return depth_by_node[leaf1] + depth_by_node[leaf2] - 2.0 * depth_by_node[lca]
        except (AttributeError, KeyError, TypeError):
            return None

    @staticmethod
    def _calculate_same_tip_distance_fast(root, tip_name: str):
        stack = [root]
        pop = stack.pop
        append = stack.append
        try:
            while stack:
                clade = pop()
                children = clade.clades
                if children:
                    child_count = len(children)
                    if child_count == 2:
                        append(children[1])
                        append(children[0])
                    else:
                        for index in range(child_count - 1, -1, -1):
                            append(children[index])
                elif clade.name == tip_name:
                    return 0.0
        except AttributeError:
            return None

        print(tip_name, "not on tree\nExiting...")
        sys.exit(2)

    def process_args(self, args) -> dict[str, str]:
        from ...helpers.plot_config import PlotConfig

        return dict(
            tree_file_path=args.tree_zero,
            tip_1=getattr(args, "tip_1", None),
            tip_2=getattr(args, "tip_2", None),
            json_output=getattr(args, "json", False),
            all_pairs=getattr(args, "all_pairs", False),
            plot=getattr(args, "plot", False),
            plot_output=getattr(args, "plot_output", "tip_to_tip_distance_heatmap.png"),
            plot_config=PlotConfig.from_args(args),
        )

    def calculate_all_pairwise_distances(self, tree_zero: Newick.Tree) -> list[dict[str, object]]:
        tips = self.calculate_terminal_names_fast(tree_zero)
        if tips is None:
            tips = [tip.name for tip in tree_zero.get_terminals()]

        if len(tips) >= self._JSON_NO_COMBO_MIN_TIPS:
            fast_result = self.calculate_pairwise_tip_distances_fast(
                tree_zero,
                tips,
                include_combos=False,
            )
            if fast_result is not None:
                _, distances = fast_result
                return self._format_pairwise_rows_without_combos(tips, distances)

        fast_result = self.calculate_pairwise_tip_distances_fast(tree_zero, tips)
        if fast_result is not None:
            combos, distances = fast_result
            return [
                {
                    "taxon_a": taxon_a,
                    "taxon_b": taxon_b,
                    "tip_to_tip_distance": round(float(distance), 4),
                }
                for (taxon_a, taxon_b), distance in zip(combos, distances)
            ]

        rows = []
        append = rows.append
        distance = TreeMixin.distance
        round_ = round
        float_ = float
        for taxon_a, taxon_b in itertools.combinations(tips, 2):
            append(
                {
                    "taxon_a": taxon_a,
                    "taxon_b": taxon_b,
                    "tip_to_tip_distance": round_(float_(distance(tree_zero, taxon_a, taxon_b)), 4),
                }
            )
        return rows

    def _format_pairwise_rows_without_combos(
        self,
        tips: list[str],
        distances: list[float],
    ) -> list[dict[str, object]]:
        if len(tips) < self._TEXT_NO_COMBO_MIN_TIPS:
            return [
                {
                    "taxon_a": taxon_a,
                    "taxon_b": taxon_b,
                    "tip_to_tip_distance": round(float(distance), 4),
                }
                for (taxon_a, taxon_b), distance in zip(
                    itertools.combinations(tips, 2),
                    distances,
                )
            ]

        rows = []
        append = rows.append
        round_ = round
        float_ = float
        idx = 0
        tip_count = len(tips)
        for i in range(tip_count - 1):
            taxon_a = tips[i]
            for j in range(i + 1, tip_count):
                append(
                    {
                        "taxon_a": taxon_a,
                        "taxon_b": tips[j],
                        "tip_to_tip_distance": round_(float_(distances[idx]), 4),
                    }
                )
                idx += 1
        return rows

    def _format_all_pairwise_distances_fast(self, tree_zero: Newick.Tree) -> str | None:
        tips = self.calculate_terminal_names_fast(tree_zero)
        if tips is None:
            return None

        if len(tips) >= self._TEXT_NO_COMBO_MIN_TIPS:
            fast_result = self.calculate_pairwise_tip_distances_fast(
                tree_zero,
                tips,
                include_combos=False,
            )
            if fast_result is None:
                return None

            _, distances = fast_result
            return self._format_pairwise_distances_without_combos(tips, distances)

        fast_result = self.calculate_pairwise_tip_distances_fast(tree_zero, tips)
        if fast_result is None:
            return None

        combos, distances = fast_result
        return "\n".join(
            f"{taxon_a}\t{taxon_b}\t{round(float(distance), 4)}"
            for (taxon_a, taxon_b), distance in zip(combos, distances)
        )

    @staticmethod
    def _format_pairwise_distances_without_combos(
        tips: list[str],
        distances: list[float],
    ) -> str:
        lines = []
        append = lines.append
        round_ = round
        float_ = float
        idx = 0
        tip_count = len(tips)
        for i in range(tip_count - 1):
            taxon_a = tips[i]
            for j in range(i + 1, tip_count):
                append(
                    f"{taxon_a}\t{tips[j]}\t{round_(float_(distances[idx]), 4)}"
                )
                idx += 1
        return "\n".join(lines)

    def _build_distance_matrix(
        self,
        rows: list[dict[str, object]],
    ) -> tuple[list[str], np.ndarray]:
        if (
            len(rows) >= self._MATRIX_FAST_FILL_MIN_ROWS
        ):
            taxa = self._sorted_upper_triangle_taxa_from_rows(rows)
            if taxa is not None and self._rows_are_sorted_upper_triangle(taxa, rows):
                n_taxa = len(taxa)
                matrix = np.zeros((n_taxa, n_taxa), dtype=np.float64)
                distances = np.fromiter(
                    (float(row["tip_to_tip_distance"]) for row in rows),
                    dtype=np.float64,
                    count=len(rows),
                )
                upper_i, upper_j = np.triu_indices(n_taxa, 1)
                matrix[upper_i, upper_j] = distances
                matrix[upper_j, upper_i] = distances
                return taxa, matrix

        taxa = sorted(
            {str(row["taxon_a"]) for row in rows}
            | {str(row["taxon_b"]) for row in rows}
        )
        n_taxa = len(taxa)

        taxon_to_index = {taxon: idx for idx, taxon in enumerate(taxa)}
        matrix = np.zeros((n_taxa, n_taxa), dtype=np.float64)

        for row in rows:
            i = taxon_to_index[str(row["taxon_a"])]
            j = taxon_to_index[str(row["taxon_b"])]
            distance = float(row["tip_to_tip_distance"])
            matrix[i, j] = distance
            matrix[j, i] = distance
        return taxa, matrix

    @staticmethod
    def _sorted_upper_triangle_taxa_from_rows(
        rows: list[dict[str, object]],
    ) -> list[str] | None:
        row_count = len(rows)
        if row_count == 0:
            return []

        discriminant = 1 + 8 * row_count
        root = math.isqrt(discriminant)
        if root * root != discriminant:
            return None

        n_taxa = (1 + root) // 2
        if n_taxa * (n_taxa - 1) // 2 != row_count:
            return None

        try:
            first_row = rows[0]
            taxa = [str(first_row["taxon_a"])]
            taxa.extend(str(rows[idx]["taxon_b"]) for idx in range(n_taxa - 1))
        except (IndexError, KeyError, TypeError):
            return None
        return taxa

    @staticmethod
    def _rows_are_sorted_upper_triangle(
        taxa: list[str],
        rows: list[dict[str, object]],
    ) -> bool:
        n_taxa = len(taxa)
        if len(rows) != n_taxa * (n_taxa - 1) // 2:
            return False

        cursor = 0
        for idx, taxon_a in enumerate(taxa[:-1]):
            for taxon_b in taxa[idx + 1:]:
                row = rows[cursor]
                row_a = row["taxon_a"]
                row_b = row["taxon_b"]
                if (
                    (row_a != taxon_a and str(row_a) != taxon_a)
                    or (row_b != taxon_b and str(row_b) != taxon_b)
                ):
                    return False
                cursor += 1
        return True

    def _plot_tip_distance_heatmap(self, rows: list[dict[str, object]]) -> None:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        except ImportError:
            print("matplotlib is required for --plot in tip_to_tip_distance. Install matplotlib and retry.")
            raise SystemExit(2)

        if not rows:
            return

        from scipy.cluster.hierarchy import linkage, leaves_list
        from scipy.spatial.distance import squareform

        taxa, matrix = self._build_distance_matrix(rows)
        n_taxa = len(taxa)
        if n_taxa >= 3:
            condensed = squareform(matrix, checks=False)
            order = leaves_list(linkage(condensed, method="average"))
        else:
            order = np.arange(n_taxa)

        ordered_matrix = matrix[np.ix_(order, order)]
        ordered_taxa = [taxa[idx] for idx in order]

        config = self.plot_config
        config.resolve(n_rows=n_taxa, n_cols=n_taxa)

        fig_w = config.fig_width or max(6, min(20, n_taxa * 0.35))
        fig_h = config.fig_height or fig_w
        fig, ax = plt.subplots(figsize=(fig_w, fig_h))
        image = ax.imshow(ordered_matrix, cmap="viridis", interpolation="nearest")

        if config.ylabel_fontsize and config.ylabel_fontsize > 0:
            ax.set_xticks(np.arange(n_taxa))
            ax.set_yticks(np.arange(n_taxa))
            ax.set_xticklabels(ordered_taxa, rotation=90, fontsize=config.xlabel_fontsize or config.ylabel_fontsize)
            ax.set_yticklabels(ordered_taxa, fontsize=config.ylabel_fontsize)
        else:
            ax.set_xticks([])
            ax.set_yticks([])

        ax.set_xlabel("Taxa (clustered)")
        ax.set_ylabel("Taxa (clustered)")
        colorbar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
        colorbar.set_label("Distance")

        if config.show_title:
            ax.set_title(config.title or "Tip-to-Tip Distance Heatmap", fontsize=config.title_fontsize)
        if config.axis_fontsize:
            ax.xaxis.label.set_fontsize(config.axis_fontsize)
            ax.yaxis.label.set_fontsize(config.axis_fontsize)

        fig.savefig(self.plot_output, dpi=config.dpi, bbox_inches="tight")
        plt.close(fig)
