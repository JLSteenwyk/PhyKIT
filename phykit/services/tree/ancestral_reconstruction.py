from __future__ import annotations

import math
import sys

from .base import Tree
from ...helpers.tree_paths import build_root_path_map
from ...errors import PhykitUserError


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


def _value_range(values):
    iterator = iter(values)
    try:
        vmin = vmax = next(iterator)
    except StopIteration:
        return None

    for value in iterator:
        if value < vmin:
            vmin = value
        elif value > vmax:
            vmax = value
    return vmin, vmax


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


np = _LazyNumpy()
_CHO_FACTOR = None
_CHO_SOLVE = None
_DISCRETE_ASR_PLOT_CACHE = {}
_DISCRETE_ASR_PLOT_CACHE_MAX = 4
_DISCRETE_ASR_PLOT_CACHE_MAX_NODES = 5000


class _LazyPickle:
    _module = None

    def _load(self):
        module = self._module
        if module is None:
            import pickle as module

            self._module = module
        return module

    def __getattr__(self, name):
        value = getattr(self._load(), name)
        setattr(self, name, value)
        return value

    def dumps(self, *args, **kwargs):
        module = self._load()
        dumps = module.dumps
        self.dumps = dumps
        if "loads" not in self.__dict__:
            self.loads = module.loads

        return dumps(*args, **kwargs)

    def loads(self, *args, **kwargs):
        module = self._load()
        loads = module.loads
        self.loads = loads
        if "dumps" not in self.__dict__:
            self.dumps = module.dumps

        return loads(*args, **kwargs)


pickle = _LazyPickle()


class _LazyHeapq:
    def merge(self, *args, **kwargs):
        import heapq as _heapq

        return _heapq.merge(*args, **kwargs)


heapq = _LazyHeapq()


def cho_factor(*args, **kwargs):
    global _CHO_FACTOR

    if _CHO_FACTOR is None:
        from scipy.linalg import cho_factor as _cho_factor

        _CHO_FACTOR = _cho_factor

    return _CHO_FACTOR(*args, **kwargs)


def cho_solve(*args, **kwargs):
    global _CHO_SOLVE

    if _CHO_SOLVE is None:
        from scipy.linalg import cho_solve as _cho_solve

        _CHO_SOLVE = _cho_solve

    return _CHO_SOLVE(*args, **kwargs)


def build_vcv_matrix(*args, **kwargs):
    from .vcv_utils import build_vcv_matrix as _build_vcv_matrix

    return _build_vcv_matrix(*args, **kwargs)


def compute_node_positions(*args, **kwargs):
    from ...helpers.plot_config import compute_node_positions as _compute_node_positions

    return _compute_node_positions(*args, **kwargs)


def build_q_matrix(*args, **kwargs):
    from ...helpers.discrete_models import build_q_matrix as _build_q_matrix

    return _build_q_matrix(*args, **kwargs)


def matrix_exp(*args, **kwargs):
    from ...helpers.discrete_models import matrix_exp as _matrix_exp

    return _matrix_exp(*args, **kwargs)


def felsenstein_pruning(*args, **kwargs):
    from ...helpers.discrete_models import (
        felsenstein_pruning as _felsenstein_pruning,
    )

    return _felsenstein_pruning(*args, **kwargs)


def fit_q_matrix(*args, **kwargs):
    from ...helpers.discrete_models import fit_q_matrix as _fit_q_matrix

    return _fit_q_matrix(*args, **kwargs)


def parse_discrete_traits(*args, **kwargs):
    from ...helpers.discrete_models import parse_discrete_traits as _parse_discrete_traits

    return _parse_discrete_traits(*args, **kwargs)


def compute_circular_coords(*args, **kwargs):
    from ...helpers.circular_layout import compute_circular_coords as _compute_circular_coords

    return _compute_circular_coords(*args, **kwargs)


def draw_circular_branches(*args, **kwargs):
    from ...helpers.circular_layout import draw_circular_branches as _draw_circular_branches

    return _draw_circular_branches(*args, **kwargs)


def draw_circular_tip_labels(*args, **kwargs):
    from ...helpers.circular_layout import draw_circular_tip_labels as _draw_circular_tip_labels

    return _draw_circular_tip_labels(*args, **kwargs)


def draw_circular_gradient_branches(*args, **kwargs):
    from ...helpers.circular_layout import (
        draw_circular_gradient_branches as _draw_circular_gradient_branches,
    )

    return _draw_circular_gradient_branches(*args, **kwargs)


def draw_circular_colored_arcs(*args, **kwargs):
    from ...helpers.circular_layout import (
        draw_circular_colored_arcs as _draw_circular_colored_arcs,
    )

    return _draw_circular_colored_arcs(*args, **kwargs)


def draw_circular_scalar_arcs(*args, **kwargs):
    from ...helpers.circular_layout import (
        draw_circular_scalar_arcs as _draw_circular_scalar_arcs,
    )

    return _draw_circular_scalar_arcs(*args, **kwargs)


def parse_color_file(*args, **kwargs):
    from ...helpers.color_annotations import parse_color_file as _parse_color_file

    return _parse_color_file(*args, **kwargs)


def resolve_mrca(*args, **kwargs):
    from ...helpers.color_annotations import resolve_mrca as _resolve_mrca

    return _resolve_mrca(*args, **kwargs)


def draw_range_rect(*args, **kwargs):
    from ...helpers.color_annotations import draw_range_rect as _draw_range_rect

    return _draw_range_rect(*args, **kwargs)


def draw_range_wedge(*args, **kwargs):
    from ...helpers.color_annotations import draw_range_wedge as _draw_range_wedge

    return _draw_range_wedge(*args, **kwargs)


def get_clade_branch_ids(*args, **kwargs):
    from ...helpers.color_annotations import get_clade_branch_ids as _get_clade_branch_ids

    return _get_clade_branch_ids(*args, **kwargs)


def build_color_legend_handles(*args, **kwargs):
    from ...helpers.color_annotations import (
        build_color_legend_handles as _build_color_legend_handles,
    )

    return _build_color_legend_handles(*args, **kwargs)


def apply_label_colors(*args, **kwargs):
    from ...helpers.color_annotations import apply_label_colors as _apply_label_colors

    return _apply_label_colors(*args, **kwargs)


class AncestralReconstruction(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"])
        self.trait_data_path = parsed["trait_data_path"]
        self.trait_column = parsed["trait_column"]
        self.method = parsed["method"]
        self.ci = parsed["ci"]
        self.plot_output = parsed["plot_output"]
        self.json_output = parsed["json_output"]
        self.trait_type = parsed["trait_type"]
        self.model = parsed["model"]
        self.plot_ci = parsed["plot_ci"]
        self.ci_size = parsed["ci_size"]
        self.plot_config = parsed["plot_config"]

    def run(self) -> None:
        tree = self.read_tree_file_unmodified()
        self.validate_tree(tree, min_tips=3, require_branch_lengths=True, context="ancestral reconstruction")
        tree_tips = self.get_tip_names_from_tree(tree)
        if self.trait_type == "discrete":
            self._run_discrete(tree, tree_tips)
        else:
            self._run_continuous(tree, tree_tips)

    def _run_continuous(self, tree, tree_tips) -> None:
        if self.trait_column is not None:
            trait_values = self._parse_multi_trait_data(
                self.trait_data_path, tree_tips, self.trait_column
            )
            trait_name = self.trait_column
        else:
            trait_values = self._parse_single_trait_data(
                self.trait_data_path, tree_tips
            )
            trait_name = "trait"

        ordered_names = sorted(trait_values.keys())
        n = len(ordered_names)
        x = np.array([trait_values[name] for name in ordered_names])

        tips_to_prune = self._tips_to_prune_for_ordered_mapping(
            tree_tips, trait_values
        )
        needs_working_copy = bool(tips_to_prune) or self.plot_config.ladderize
        tree_for_analysis = self._fast_copy(tree) if needs_working_copy else tree
        if tips_to_prune:
            tree_for_analysis = self.prune_tree_using_taxa_list(
                tree_for_analysis, tips_to_prune
            )

        if self.plot_config.ladderize:
            tree_for_analysis.ladderize()

        # Label internal nodes
        node_labels = self._label_internal_nodes(tree_for_analysis)

        if self.method == "ml":
            (
                node_estimates, node_cis, sigma2, log_likelihood
            ) = self._anc_ml(tree_for_analysis, x, ordered_names, node_labels)
        else:
            (
                node_estimates, node_cis, sigma2, log_likelihood
            ) = self._fast_anc(tree_for_analysis, x, ordered_names, node_labels)

        # Build result
        result = self._format_result(
            method=self.method,
            trait_name=trait_name,
            n_tips=n,
            sigma2=sigma2,
            log_likelihood=log_likelihood,
            node_estimates=node_estimates,
            node_cis=node_cis,
            node_labels=node_labels,
            tree=tree_for_analysis,
            trait_values=trait_values,
        )

        if self.plot_output:
            ci_data = node_cis if (self.plot_ci and self.ci and node_cis) else None
            self._plot_contmap(
                tree_for_analysis, node_estimates, node_labels,
                trait_values, trait_name, self.plot_output,
                node_cis=ci_data,
            )
            result["plot_output"] = self.plot_output

        if self.json_output:
            print_json(result)
        else:
            self._print_text_output(
                method=self.method,
                trait_name=trait_name,
                n_tips=n,
                sigma2=sigma2,
                log_likelihood=log_likelihood,
                node_estimates=node_estimates,
                node_cis=node_cis,
                node_labels=node_labels,
                tree=tree_for_analysis,
            )

    def process_args(self, args) -> Dict:
        from ...helpers.plot_config import PlotConfig

        return dict(
            tree_file_path=args.tree,
            trait_data_path=args.trait_data,
            trait_column=getattr(args, "trait", None),
            method=getattr(args, "method", "fast"),
            ci=getattr(args, "ci", False),
            plot_output=getattr(args, "plot", None),
            json_output=getattr(args, "json", False),
            trait_type=getattr(args, "type", "continuous"),
            model=getattr(args, "model", "ER"),
            plot_ci=getattr(args, "plot_ci", False),
            ci_size=getattr(args, "ci_size", 1.0),
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

    def _parse_single_trait_data(
        self, path: str, tree_tips: List[str]
    ) -> Dict[str, float]:
        try:
            traits = {}
            with open(path) as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    if not line or line[0] == "#":
                        continue
                    taxon, sep, value_str = line.partition("\t")
                    if not sep or "\t" in value_str:
                        column_count = line.count("\t") + 1
                        raise PhykitUserError(
                            [
                                f"Line {line_num} in trait file has {column_count} columns; expected 2.",
                                "Each line should be: taxon_name<tab>trait_value",
                            ],
                            code=2,
                        )
                    try:
                        traits[taxon] = float(value_str)
                    except ValueError:
                        raise PhykitUserError(
                            [
                                f"Non-numeric trait value '{value_str}' for taxon '{taxon}' on line {line_num}.",
                            ],
                            code=2,
                        )
        except FileNotFoundError:
            raise PhykitUserError(
                [
                    f"{path} corresponds to no such file or directory.",
                    "Please check filename and pathing",
                ],
                code=2,
            )

        if (
            len(tree_tips) >= 3
            and len(tree_tips) == len(traits)
            and next(iter(traits)) == tree_tips[0]
            and next(reversed(traits)) == tree_tips[-1]
            and list(traits) == tree_tips
        ):
            return traits

        tree_tip_set = set(tree_tips)
        if (
            len(tree_tip_set) >= 3
            and len(tree_tip_set) == len(traits)
            and tree_tip_set == traits.keys()
        ):
            return traits

        trait_taxa_set = set(traits)
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

    def _parse_multi_trait_data(
        self, path: str, tree_tips: List[str], trait_column: str
    ) -> Dict[str, float]:
        try:
            with open(path) as f:
                header_parts = None
                for line in f:
                    stripped = line.strip()
                    if not stripped or stripped.startswith("#"):
                        continue
                    header_parts = stripped.split("\t")
                    break

                if header_parts is None:
                    raise PhykitUserError(
                        [
                            "Multi-trait file must have a header row and at least one data row.",
                        ],
                        code=2,
                    )

                if len(header_parts) < 2:
                    raise PhykitUserError(
                        [
                            "Header must have at least 2 columns (taxon + at least 1 trait).",
                        ],
                        code=2,
                    )
                trait_names = header_parts[1:]

                if trait_column not in trait_names:
                    raise PhykitUserError(
                        [
                            f"Column '{trait_column}' not found in trait file.",
                            f"Available columns: {', '.join(trait_names)}",
                        ],
                        code=2,
                    )

                col_idx = trait_names.index(trait_column)

                traits = {}
                saw_data = False
                data_line_idx = 2
                for line in f:
                    stripped = line.strip()
                    if not stripped or stripped.startswith("#"):
                        continue
                    saw_data = True
                    parts = stripped.split("\t")
                    if len(parts) != len(header_parts):
                        raise PhykitUserError(
                            [
                                f"Line {data_line_idx} has {len(parts)} columns; expected {len(header_parts)}.",
                            ],
                            code=2,
                        )
                    taxon = parts[0]
                    val_str = parts[1 + col_idx]
                    try:
                        traits[taxon] = float(val_str)
                    except ValueError:
                        raise PhykitUserError(
                            [
                                f"Non-numeric trait value '{val_str}' for taxon '{taxon}' "
                                f"(trait '{trait_column}') on line {data_line_idx}.",
                            ],
                            code=2,
                        )
                    data_line_idx += 1

                if not saw_data:
                    raise PhykitUserError(
                        [
                            "Multi-trait file must have a header row and at least one data row.",
                        ],
                        code=2,
                    )
        except FileNotFoundError:
            raise PhykitUserError(
                [
                    f"{path} corresponds to no such file or directory.",
                    "Please check filename and pathing",
                ],
                code=2,
            )

        if (
            len(tree_tips) >= 3
            and len(tree_tips) == len(traits)
            and next(iter(traits)) == tree_tips[0]
            and next(reversed(traits)) == tree_tips[-1]
            and list(traits) == tree_tips
        ):
            return traits

        tree_tip_set = set(tree_tips)
        if (
            len(tree_tip_set) >= 3
            and len(tree_tip_set) == len(traits)
            and tree_tip_set == traits.keys()
        ):
            return traits

        trait_taxa_set = set(traits)
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
        try:
            root = tree.root
            root.clades
        except AttributeError:
            root = None

        if root is not None:
            stack = [root]
            try:
                while stack:
                    clade = stack.pop()
                    children = clade.clades
                    if children:
                        if clade.name:
                            labels[id(clade)] = clade.name
                        else:
                            labels[id(clade)] = f"N{counter}"
                            counter += 1
                        stack.extend(reversed(children))
                return labels
            except AttributeError:
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

    def _get_descendant_tips(self, tree, node) -> List[str]:
        tips = Tree.calculate_terminal_names_fast(node)
        if tips is not None:
            return sorted(tips)
        return sorted(tip.name for tip in node.get_terminals())

    def _collect_descendant_tip_counts(self, tree):
        try:
            preorder = list(self._iter_preorder(tree.root))
        except AttributeError:
            return None

        counts = {}
        try:
            for clade in reversed(preorder):
                children = clade.clades
                if children:
                    counts[id(clade)] = sum(counts[id(child)] for child in children)
                else:
                    counts[id(clade)] = 1
        except (AttributeError, KeyError):
            return None
        return counts

    def _collect_descendant_tip_names(self, tree):
        try:
            preorder = list(self._iter_preorder(tree.root))
        except AttributeError:
            return None

        names = {}
        try:
            for clade in reversed(preorder):
                children = clade.clades
                if children:
                    if len(children) == 2:
                        left = names[id(children[0])]
                        right = names[id(children[1])]
                        if left and right and left[-1] <= right[0]:
                            names[id(clade)] = left + right
                        elif left and right and right[-1] <= left[0]:
                            names[id(clade)] = right + left
                        else:
                            names[id(clade)] = tuple(heapq.merge(left, right))
                    else:
                        child_names = [names[id(child)] for child in children]
                        ordered = True
                        previous = child_names[0]
                        for current in child_names[1:]:
                            if previous and current and previous[-1] > current[0]:
                                ordered = False
                                break
                            previous = current
                        if ordered:
                            merged_names = []
                            for current in child_names:
                                merged_names.extend(current)
                            names[id(clade)] = tuple(merged_names)
                        else:
                            names[id(clade)] = tuple(heapq.merge(*child_names))
                else:
                    names[id(clade)] = (clade.name,)
        except (AttributeError, KeyError, TypeError):
            return None
        return names

    def _build_parent_map(self, tree) -> Dict:
        parent_map = {}
        try:
            root = tree.root
            root.clades
        except AttributeError:
            root = None

        if root is not None:
            stack = [root]
            pop = stack.pop
            extend = stack.extend
            try:
                while stack:
                    clade = pop()
                    children = clade.clades
                    for child in children:
                        parent_map[id(child)] = clade
                    if children:
                        extend(children)
                return parent_map
            except AttributeError:
                parent_map = {}

        for clade in tree.find_clades(order="preorder"):
            for child in clade.clades:
                parent_map[id(child)] = clade
        return parent_map

    def _fast_anc(
        self, tree, x: np.ndarray, ordered_names: List[str],
        node_labels: Dict,
    ) -> Tuple[Dict, Dict, float, float]:
        """Two-pass ML ancestral reconstruction (Felsenstein's algorithm).

        Pass 1 (postorder): compute subtree estimates and variances.
        Pass 2 (preorder): incorporate information from the rest of the
        tree so that every internal node gets the full-tree ML estimate.
        """
        name_to_idx = {name: i for i, name in enumerate(ordered_names)}
        try:
            root = tree.root
            root.clades
            preorder_clades = list(self._iter_preorder(root))
            parent_map = {}
            for clade in preorder_clades:
                for child in clade.clades:
                    parent_map[id(child)] = clade
        except AttributeError:
            root = tree.root
            preorder_clades = None
            parent_map = self._build_parent_map(tree)

        # --- Pass 1: postorder (subtree estimates) ---
        est_down = {}
        var_down = {}

        postorder_clades = (
            reversed(preorder_clades)
            if preorder_clades is not None
            else tree.find_clades(order="postorder")
        )
        for clade in postorder_clades:
            children = clade.clades
            clade_id = id(clade)
            if not children:
                idx = name_to_idx.get(clade.name)
                if idx is not None:
                    est_down[clade_id] = x[idx]
                    var_down[clade_id] = 0.0
            else:
                total_prec = 0.0
                weighted_sum = 0.0
                for child in children:
                    child_id = id(child)
                    child_estimate = est_down.get(child_id)
                    if child_estimate is None:
                        continue
                    v_i = child.branch_length if child.branch_length else 0.0
                    child_var = var_down[child_id]
                    denom = v_i + child_var
                    if denom == 0:
                        denom = 1e-10
                    prec = 1.0 / denom
                    total_prec += prec
                    weighted_sum += prec * child_estimate

                if total_prec:
                    est_down[clade_id] = weighted_sum / total_prec
                    var_down[clade_id] = 1.0 / total_prec

        # --- Pass 2: preorder (incorporate full-tree information) ---
        # For each node, compute the "message from above": an estimate
        # and variance representing all information NOT in the subtree
        # below that node.
        est_up = {}
        var_up = {}
        est_final = {}
        var_final = {}

        # Root has no information from above
        root_id = id(root)
        est_final[root_id] = est_down[root_id]
        var_final[root_id] = var_down[root_id]

        preorder_iter = (
            preorder_clades
            if preorder_clades is not None
            else tree.find_clades(order="preorder")
        )
        for clade in preorder_iter:
            clade_id = id(clade)
            if clade_id == root_id:
                continue
            parent = parent_map.get(clade_id)
            if parent is None:
                continue
            parent_id = id(parent)
            bl = clade.branch_length if clade.branch_length else 0.0

            # Gather precision from parent's "above" message
            prec_above = 0.0
            est_above_weighted = 0.0
            parent_estimate = est_up.get(parent_id)
            if parent_estimate is not None:
                v_up_parent = var_up[parent_id]
                if v_up_parent > 0:
                    p = 1.0 / v_up_parent
                    prec_above += p
                    est_above_weighted += p * parent_estimate
                elif v_up_parent == 0 and parent_id != root_id:
                    p = 1.0 / 1e-10
                    prec_above += p
                    est_above_weighted += p * parent_estimate

            # Gather precision from siblings
            for sibling in parent.clades:
                sibling_id = id(sibling)
                if sibling_id == clade_id:
                    continue
                sibling_estimate = est_down.get(sibling_id)
                if sibling_estimate is None:
                    continue
                sib_bl = sibling.branch_length if sibling.branch_length else 0.0
                sib_var = var_down[sibling_id] + sib_bl
                if sib_var == 0:
                    sib_var = 1e-10
                p = 1.0 / sib_var
                prec_above += p
                est_above_weighted += p * sibling_estimate

            # The "above" message for this node
            if prec_above > 0:
                up_estimate = est_above_weighted / prec_above
                est_up[clade_id] = up_estimate
                up_variance = bl + 1.0 / prec_above
                var_up[clade_id] = up_variance
            else:
                # No information from above (shouldn't happen in a valid tree)
                up_estimate = est_down.get(clade_id, 0.0)
                est_up[clade_id] = up_estimate
                up_variance = bl + 1e10
                var_up[clade_id] = up_variance

            # Combine down and up for final estimate
            down_estimate = est_down.get(clade_id)
            if down_estimate is not None:
                vd = var_down[clade_id]
                vu = up_variance
                if vd == 0:
                    vd = 1e-10
                if vu == 0:
                    vu = 1e-10
                prec_d = 1.0 / vd
                prec_u = 1.0 / vu
                total = prec_d + prec_u
                est_final[clade_id] = (
                    prec_d * down_estimate + prec_u * up_estimate
                ) / total
                var_final[clade_id] = 1.0 / total

        # Compute sigma^2 from phylogenetic independent contrasts
        sigma2 = self._compute_sigma2_from_contrasts(
            tree,
            x,
            est_down,
            var_down,
            ordered_names,
            postorder_clades=(
                reversed(preorder_clades)
                if preorder_clades is not None
                else None
            ),
        )

        # Compute log-likelihood
        log_likelihood = self._compute_log_likelihood(
            tree, x, ordered_names, sigma2
        )

        # Build output dicts keyed by node label
        node_estimates = {}
        node_cis = {}
        preorder_iter = (
            preorder_clades
            if preorder_clades is not None
            else tree.find_clades(order="preorder")
        )
        for clade in preorder_iter:
            if clade.clades:
                clade_id = id(clade)
                label = node_labels.get(clade_id)
                if label is None:
                    continue
                est = est_final.get(clade_id)
                if est is None:
                    continue
                node_estimates[label] = est
                if self.ci:
                    var = var_final.get(clade_id)
                    if var is not None:
                        se = math.sqrt(sigma2 * var) if sigma2 * var > 0 else 0.0
                        node_cis[label] = (est - 1.96 * se, est + 1.96 * se)

        return node_estimates, node_cis, sigma2, log_likelihood

    def _compute_sigma2_from_contrasts(
        self, tree, x: np.ndarray, node_estimates: Dict,
        node_variances: Dict, ordered_names: List[str], postorder_clades=None,
    ) -> float:
        """Compute BM rate (sigma^2) from phylogenetic independent contrasts.

        For nodes with >2 children (polytomies), iteratively compute
        pairwise contrasts and update the combined estimate/variance,
        producing one contrast per additional child (n_children - 1
        contrasts per node). This ensures the total number of contrasts
        equals n_tips - 1 for a fully resolved tree.
        """
        contrasts_sum = 0.0
        contrast_count = 0
        if postorder_clades is None:
            postorder_clades = tree.find_clades(order="postorder")
        for clade in postorder_clades:
            if not clade.clades:
                continue
            children = clade.clades
            first_child = None
            combined_est = 0.0
            combined_var = 0.0
            valid_count = 0

            for child in children:
                child_id = id(child)
                child_estimate = node_estimates.get(child_id)
                if child_estimate is None:
                    continue

                if valid_count == 0:
                    first_child = child
                    valid_count = 1
                    continue

                if valid_count == 1:
                    first_id = id(first_child)
                    v1 = (
                        first_child.branch_length or 0.0
                    ) + node_variances[first_id]
                    v2 = (child.branch_length or 0.0) + node_variances[child_id]
                    denom = v1 + v2
                    if denom > 0:
                        diff = node_estimates[first_id] - child_estimate
                        contrasts_sum += (diff * diff) / denom
                        contrast_count += 1
                        combined_est = (
                            v2 * node_estimates[first_id] + v1 * child_estimate
                        ) / denom
                        combined_var = (v1 * v2) / denom
                    else:
                        combined_est = node_estimates[first_id]
                        combined_var = 0.0
                    valid_count = 2
                    continue

                vk = (child.branch_length or 0.0) + node_variances[child_id]
                vc = combined_var
                denom_k = vk + vc
                if denom_k > 0:
                    diff = child_estimate - combined_est
                    contrasts_sum += (diff * diff) / denom_k
                    contrast_count += 1
                    combined_est = (
                        vc * child_estimate + vk * combined_est
                    ) / denom_k
                    combined_var = (vk * vc) / denom_k

        if contrast_count:
            return contrasts_sum / contrast_count
        return 0.0

    def _compute_log_likelihood(
        self, tree, x: np.ndarray, ordered_names: List[str], sigma2: float,
    ) -> float:
        """Compute log-likelihood under BM model."""
        n = len(ordered_names)
        if sigma2 <= 0 or n == 0:
            return float("-inf")

        vcv = self._build_vcv_matrix(tree, ordered_names)
        C = sigma2 * vcv

        try:
            factor = cho_factor(C, lower=True, check_finite=False)
            rhs = np.empty((n, 2), dtype=np.result_type(C, x))
            rhs[:, 0] = 1.0
            rhs[:, 1] = x
            solved = cho_solve(factor, rhs, check_finite=False)
            C_inv_ones = solved[:, 0]
            C_inv_x = solved[:, 1]
            logdet = 2.0 * float(np.log(np.diagonal(factor[0])).sum())
        except (np.linalg.LinAlgError, FloatingPointError, ValueError):
            try:
                sign, logdet = np.linalg.slogdet(C)
                if sign <= 0:
                    return float("-inf")
                C_inv = np.linalg.inv(C)
            except np.linalg.LinAlgError:
                return float("-inf")

            ones = np.ones(n)
            C_inv_ones = C_inv @ ones
            C_inv_x = C_inv @ x

        ones = np.ones(n)
        a_hat = float(ones @ C_inv_x) / float(ones @ C_inv_ones)
        residuals = x - a_hat
        C_inv_residuals = C_inv_x - a_hat * C_inv_ones

        ll = -0.5 * (
            n * np.log(2 * np.pi)
            + logdet
            + float(residuals @ C_inv_residuals)
        )
        return float(ll)

    def _build_vcv_matrix(
        self, tree, ordered_names: List[str]
    ) -> np.ndarray:
        return build_vcv_matrix(tree, ordered_names)

    def _anc_ml(
        self, tree, x: np.ndarray, ordered_names: List[str],
        node_labels: Dict,
    ) -> Tuple[Dict, Dict, float, float]:
        """VCV-based ML ancestral reconstruction with exact conditional CIs."""
        n = len(ordered_names)
        vcv = self._build_vcv_matrix(tree, ordered_names)

        try:
            C_inv = np.linalg.inv(vcv)
        except np.linalg.LinAlgError:
            raise PhykitUserError(
                [
                    "Singular VCV matrix: cannot invert.",
                    "Check that the tree has valid branch lengths.",
                ],
                code=2,
            )

        ones = np.ones(n)
        C_inv_x = C_inv @ x
        C_inv_ones = C_inv @ ones
        root_precision = float(ones @ C_inv_ones)
        # ML root estimate
        a_hat = float(ones @ C_inv_x) / root_precision
        # ML rate (sigma^2)
        residuals = x - a_hat
        C_inv_residuals = C_inv_x - C_inv_ones * a_hat
        residual_ss = float(residuals @ C_inv_residuals)
        sigma2_ml = residual_ss / n
        # REML rate for CIs (matches R's fastAnc)
        sigma2_reml = residual_ss / (n - 1)

        # Log-likelihood (uses ML sigma^2)
        sign, logdet = np.linalg.slogdet(vcv)
        ll = -0.5 * (
            n * np.log(2 * np.pi) + n * np.log(sigma2_ml) + logdet + n
        )

        # Variance of the root estimate
        root_var = 1.0 / root_precision

        # Compute cross-covariance and conditional estimates for each internal node
        node_estimates = {}
        node_cis = {}
        direct_cross_covariance = None
        try:
            depths = tree.depths()
            root = tree.root
            root_depth = depths[root]
            direct_cross_covariance = self._prepare_ml_cross_covariances_direct(
                tree,
                ordered_names,
                node_labels,
                depths,
                root_depth,
            )
        except (AttributeError, KeyError, TypeError):
            depths = None
            root = tree.root
            root_depth = None

        if direct_cross_covariance is not None:
            (
                estimate_clades,
                estimate_labels,
                d_values,
                cross_covariance_matrix,
            ) = direct_cross_covariance
        else:
            cross_covariances = []
            d_values = []
            estimate_clades = []
            estimate_labels = []
            try:
                if depths is None:
                    raise AttributeError
                terminal_by_name = {
                    terminal.name: terminal for terminal in tree.get_terminals()
                }
                root_paths = build_root_path_map(tree, root)
                if not root_paths:
                    raise KeyError
                tip_paths = {
                    name: root_paths[terminal_by_name[name]]
                    for name in ordered_names
                }
                fast_cross_covariance = True
            except (AttributeError, KeyError, TypeError):
                terminals = tree.get_terminals()
                fast_cross_covariance = False

            for clade in tree.find_clades(order="preorder"):
                if clade.is_terminal():
                    continue
                if id(clade) not in node_labels:
                    continue

                label = node_labels[id(clade)]
                estimate_clades.append(clade)
                estimate_labels.append(label)
                if fast_cross_covariance:
                    d_i = depths[clade] - root_depth
                else:
                    d_i = tree.distance(tree.root, clade)
                d_values.append(d_i)

                # Cross-covariance: c_i[j] = distance from root to MRCA(node_i, tip_j)
                c_i = np.zeros(n)
                if fast_cross_covariance:
                    clade_ancestors = set(root_paths[clade])
                    for j, tip_name in enumerate(ordered_names):
                        mrca = root
                        for ancestor in reversed(tip_paths[tip_name]):
                            if ancestor in clade_ancestors:
                                mrca = ancestor
                                break
                        c_i[j] = depths[mrca] - root_depth
                else:
                    for j, tip_name in enumerate(ordered_names):
                        tip_clade = None
                        for terminal in terminals:
                            if terminal.name == tip_name:
                                tip_clade = terminal
                                break
                        mrca = tree.common_ancestor(clade, tip_clade)
                        c_i[j] = tree.distance(tree.root, mrca)
                cross_covariances.append(c_i)

            if not cross_covariances:
                return node_estimates, node_cis, sigma2_reml, float(ll)
            cross_covariance_matrix = np.vstack(cross_covariances)
            d_values = np.asarray(d_values, dtype=float)

        if cross_covariance_matrix.size == 0:
            return node_estimates, node_cis, sigma2_reml, float(ll)

        cond_means = a_hat + cross_covariance_matrix @ C_inv_residuals
        if self.ci:
            # Batch the inverse products for all internal nodes to avoid one
            # dense vector-matrix multiply per node.
            cross_covariance_inv = cross_covariance_matrix @ C_inv
            covariance_terms = np.einsum(
                "ij,ij->i", cross_covariance_inv, cross_covariance_matrix
            )
        else:
            covariance_terms = None

        for idx, label in enumerate(estimate_labels):
            cond_mean = float(cond_means[idx])
            node_estimates[label] = cond_mean

            if self.ci:
                clade = estimate_clades[idx]
                if clade == tree.root:
                    # Root variance from GLS estimation
                    cond_var = sigma2_reml * root_var
                else:
                    # Conditional variance using REML sigma^2
                    cond_var = sigma2_reml * (
                        d_values[idx] - float(covariance_terms[idx])
                    )
                if cond_var < 0:
                    cond_var = 0.0
                se = math.sqrt(cond_var)
                node_cis[label] = (cond_mean - 1.96 * se, cond_mean + 1.96 * se)

        return node_estimates, node_cis, sigma2_reml, float(ll)

    def _prepare_ml_cross_covariances_direct(
        self,
        tree,
        ordered_names: List[str],
        node_labels: Dict,
        depths: Dict,
        root_depth: float,
    ):
        try:
            root = tree.root
            preorder = list(self._iter_preorder(root))
        except AttributeError:
            return None

        n = len(ordered_names)
        name_to_index = {name: idx for idx, name in enumerate(ordered_names)}
        if len(name_to_index) != n:
            return None

        descendant_indices = {}
        seen_tips = set()
        try:
            for clade in reversed(preorder):
                children = clade.clades
                if children:
                    indices = []
                    for child in children:
                        indices.extend(descendant_indices[id(child)])
                else:
                    idx = name_to_index.get(clade.name)
                    if idx is None:
                        indices = []
                    elif idx in seen_tips:
                        return None
                    else:
                        seen_tips.add(idx)
                        indices = [idx]
                descendant_indices[id(clade)] = np.asarray(indices, dtype=np.intp)
        except (AttributeError, KeyError):
            return None

        if len(seen_tips) != n:
            return None

        cov_by_id = {id(root): np.zeros(n)}
        estimate_clades = []
        estimate_labels = []
        d_values = []
        cross_covariances = []

        try:
            for clade in preorder:
                children = clade.clades
                if not children:
                    continue

                clade_id = id(clade)
                c_i = cov_by_id[clade_id]
                if clade_id in node_labels:
                    estimate_clades.append(clade)
                    estimate_labels.append(node_labels[clade_id])
                    d_values.append(depths[clade] - root_depth)
                    cross_covariances.append(c_i)

                for child in children:
                    child_vec = c_i.copy()
                    child_indices = descendant_indices[id(child)]
                    if child_indices.size:
                        child_vec[child_indices] = depths[child] - root_depth
                    cov_by_id[id(child)] = child_vec
        except (AttributeError, KeyError):
            return None

        if not cross_covariances:
            return estimate_clades, estimate_labels, np.array([]), np.empty((0, n))

        return (
            estimate_clades,
            estimate_labels,
            np.asarray(d_values, dtype=float),
            np.vstack(cross_covariances),
        )

    def _format_result(
        self, *, method, trait_name, n_tips, sigma2, log_likelihood,
        node_estimates, node_cis, node_labels, tree, trait_values,
    ) -> Dict:
        ancestral_estimates = {}
        root_id = id(tree.root)
        descendant_names = self._collect_descendant_tip_names(tree)

        clade_iter = (
            self._iter_preorder(tree.root)
            if descendant_names is not None
            else tree.find_clades(order="preorder")
        )
        for clade in clade_iter:
            if clade.is_terminal():
                continue
            if id(clade) not in node_labels:
                continue
            label = node_labels[id(clade)]
            if label not in node_estimates:
                continue

            entry = {
                "estimate": float(node_estimates[label]),
                "descendants": (
                    list(descendant_names[id(clade)])
                    if descendant_names is not None
                    else self._get_descendant_tips(tree, clade)
                ),
                "is_root": id(clade) == root_id,
            }
            if label in node_cis:
                entry["ci_lower"] = float(node_cis[label][0])
                entry["ci_upper"] = float(node_cis[label][1])
            ancestral_estimates[label] = entry

        result = {
            "method": method,
            "trait": trait_name,
            "n_tips": n_tips,
            "log_likelihood": float(log_likelihood),
            "sigma2": float(sigma2),
            "ancestral_estimates": ancestral_estimates,
            "tip_values": {k: float(v) for k, v in sorted(trait_values.items())},
        }
        return result

    def _print_text_output(
        self, *, method, trait_name, n_tips, sigma2, log_likelihood,
        node_estimates, node_cis, node_labels, tree,
    ) -> None:
        method_desc = (
            "fast (Felsenstein's contrasts)" if method == "fast"
            else "ml (maximum likelihood, VCV-based)"
        )
        lines = [
            "Ancestral State Reconstruction",
            f"\nMethod: {method_desc}",
            f"Trait: {trait_name}",
            f"Number of tips: {n_tips}",
            f"\nLog-likelihood: {log_likelihood:.4f}",
            f"Sigma-squared (BM rate): {sigma2:.6f}",
            "\nAncestral estimates:",
        ]
        append = lines.append

        root_id = id(tree.root)
        descendant_counts = self._collect_descendant_tip_counts(tree)

        if node_cis:
            append(
                f"  {'Node':<12s}{'Descendants':>12s}{'Estimate':>12s}{'95% CI':>24s}"
            )
            clade_iter = (
                self._iter_preorder(tree.root)
                if descendant_counts is not None
                else tree.find_clades(order="preorder")
            )
            for clade in clade_iter:
                if clade.is_terminal():
                    continue
                if id(clade) not in node_labels:
                    continue
                label = node_labels[id(clade)]
                if label not in node_estimates:
                    continue
                n_desc = (
                    descendant_counts[id(clade)]
                    if descendant_counts is not None
                    else len(self._get_descendant_tips(tree, clade))
                )
                est = node_estimates[label]
                root_tag = " (root)" if id(clade) == root_id else ""
                if label in node_cis:
                    ci_lo, ci_hi = node_cis[label]
                    ci_str = f"[{ci_lo:.4f}, {ci_hi:.4f}]"
                else:
                    ci_str = ""
                append(
                    f"  {label + root_tag:<12s}{n_desc:>12d}{est:>12.4f}{ci_str:>24s}"
                )
        else:
            append(
                f"  {'Node':<12s}{'Descendants':>12s}{'Estimate':>12s}"
            )
            clade_iter = (
                self._iter_preorder(tree.root)
                if descendant_counts is not None
                else tree.find_clades(order="preorder")
            )
            for clade in clade_iter:
                if clade.is_terminal():
                    continue
                if id(clade) not in node_labels:
                    continue
                label = node_labels[id(clade)]
                if label not in node_estimates:
                    continue
                n_desc = (
                    descendant_counts[id(clade)]
                    if descendant_counts is not None
                    else len(self._get_descendant_tips(tree, clade))
                )
                est = node_estimates[label]
                root_tag = " (root)" if id(clade) == root_id else ""
                append(
                    f"  {label + root_tag:<12s}{n_desc:>12d}{est:>12.4f}"
                )
        print("\n".join(lines))

    def _plot_contmap(
        self, tree, node_estimates, node_labels, trait_values,
        trait_name, output_path, node_cis=None,
    ) -> None:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            from matplotlib.collections import LineCollection
            from matplotlib.colors import Normalize
        except ImportError:
            print(
                "matplotlib is required for contMap plotting. "
                "Install matplotlib and retry."
            )
            raise SystemExit(2)

        preorder_clades = list(self._iter_preorder(tree.root))
        parent_map = {}
        for clade in preorder_clades:
            for child in clade.clades:
                parent_map[id(child)] = clade

        # Build estimates dict keyed by id(clade) for all nodes
        all_estimates = {}
        node_labels_map = {}  # id(clade) → label string (for CI lookup)
        for clade in preorder_clades:
            if not clade.clades:
                if clade.name in trait_values:
                    all_estimates[id(clade)] = trait_values[clade.name]
            else:
                if id(clade) in node_labels:
                    label = node_labels[id(clade)]
                    node_labels_map[id(clade)] = label
                    if label in node_estimates:
                        all_estimates[id(clade)] = node_estimates[label]

        tips = [clade for clade in preorder_clades if not clade.clades]
        root = tree.root
        node_x, node_y = compute_node_positions(
            tree,
            parent_map,
            cladogram=self.plot_config.cladogram,
            preorder_clades=preorder_clades,
        )

        # Normalize trait values for colormap
        value_range = _value_range(all_estimates.values())
        if value_range is None:
            return
        vmin, vmax = value_range
        if vmin == vmax:
            vmax = vmin + 1.0
        norm = Normalize(vmin=vmin, vmax=vmax)
        cmap = plt.get_cmap("coolwarm")

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

            gradient_branches = []
            for clade in preorder_clades:
                if clade == root:
                    continue
                if id(clade) not in parent_map:
                    continue
                parent = parent_map[id(clade)]
                pid = id(parent)
                cid = id(clade)
                if pid not in coords or cid not in coords:
                    continue
                if pid not in all_estimates or cid not in all_estimates:
                    continue

                parent_val = all_estimates[pid]
                child_val = all_estimates[cid]

                gradient_branches.append(
                    (coords[pid], coords[cid], parent_val, child_val)
                )
            draw_circular_gradient_branches(
                ax, gradient_branches, cmap, vmin, vmax, lw=3,
            )

            # Arcs at internal nodes colored by trait value
            scalar_arcs = []
            for clade in preorder_clades:
                if not clade.clades:
                    continue
                cid = id(clade)
                if cid not in coords or cid not in all_estimates:
                    continue
                child_angles = [coords[id(ch)]["angle"] for ch in clade.clades if id(ch) in coords]
                if len(child_angles) < 2:
                    continue
                node_val = all_estimates[cid]
                min_a = min(child_angles)
                max_a = max(child_angles)
                scalar_arcs.append(
                    (0, 0, coords[cid]["radius"], min_a, max_a, node_val)
                )
            draw_circular_scalar_arcs(ax, scalar_arcs, cmap, norm, lw=3)

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

            # Colorbar
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            cbar = fig.colorbar(sm, ax=ax, pad=0.15)
            cbar.set_label(trait_name)

            if config.show_title:
                ax.set_title(
                    config.title or "Continuous Trait Map (contMap)",
                    fontsize=config.title_fontsize,
                )
        else:
            # --- Rectangular mode ---
            n_seg = 50
            segment_start_fracs = [s / n_seg for s in range(n_seg)]
            segment_end_fracs = [(s + 1) / n_seg for s in range(n_seg)]
            horizontal_segments = []
            horizontal_values = []
            vertical_segments = []
            vertical_values = []

            for clade in preorder_clades:
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
                dx = x1 - x0
                dv = val1 - val0
                for frac_s, frac_e in zip(segment_start_fracs, segment_end_fracs):
                    xs = x0 + frac_s * dx
                    xe = x0 + frac_e * dx
                    v = val0 + ((frac_s + frac_e) / 2) * dv
                    horizontal_segments.append(((xs, y1), (xe, y1)))
                    horizontal_values.append(v)

                # Vertical connector (at parent's x): parent's trait color
                vertical_segments.append(((x0, y0), (x0, y1)))
                vertical_values.append(val0)

            if horizontal_segments:
                horizontal_collection = LineCollection(
                    horizontal_segments,
                    cmap=cmap,
                    norm=norm,
                    linewidths=3,
                    capstyle="butt",
                )
                horizontal_collection.set_array(
                    np.asarray(horizontal_values, dtype=float)
                )
                ax.add_collection(horizontal_collection, autolim=True)
            if vertical_segments:
                vertical_collection = LineCollection(
                    vertical_segments,
                    cmap=cmap,
                    norm=norm,
                    linewidths=3,
                    capstyle="butt",
                )
                vertical_collection.set_array(
                    np.asarray(vertical_values, dtype=float)
                )
                ax.add_collection(vertical_collection, autolim=True)
            ax.autoscale_view()

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
                ax.set_title(
                    config.title or "Continuous Trait Map (contMap)",
                    fontsize=config.title_fontsize,
                )

        # Draw CI bars at internal nodes if requested
        if node_cis:
            ci_scale = self.ci_size
            # Build label→id mapping for looking up CI values
            label_to_id = {
                label: cid for cid, label in node_labels_map.items()
            }

            if config.circular:
                # Circular: use data coordinates from coords dict
                ci_segments = []
                ci_x = []
                ci_y = []
                for label, (ci_lo, ci_hi) in node_cis.items():
                    cid = label_to_id.get(label)
                    if cid is None or cid not in coords:
                        continue
                    cx = coords[cid]["x"]
                    cy = coords[cid]["y"]
                    angle = coords[cid]["angle"]
                    # CI bar perpendicular to the radius (tangential)
                    bar_len = (ci_hi - ci_lo) * ci_scale * 0.3
                    dx = -np.sin(angle) * bar_len / 2
                    dy = np.cos(angle) * bar_len / 2
                    ci_segments.append(((cx - dx, cy - dy), (cx + dx, cy + dy)))
                    ci_x.append(cx)
                    ci_y.append(cy)
                if ci_segments:
                    ax.add_collection(
                        LineCollection(
                            ci_segments,
                            colors="black",
                            linewidths=1.5 * ci_scale,
                            zorder=8,
                        ),
                        autolim=True,
                    )
                    ax.autoscale_view()
                if ci_x:
                    ax.scatter(ci_x, ci_y, s=15 * ci_scale, c="black", zorder=9)
            else:
                # Rectangular: vertical bars at node positions
                max_x_val = max(node_x.values()) if node_x else 1.0
                ci_vertical_segments = []
                ci_cap_segments = []
                ci_x = []
                ci_y = []
                for label, (ci_lo, ci_hi) in node_cis.items():
                    cid = label_to_id.get(label)
                    if cid is None or cid not in node_x or cid not in node_y:
                        continue
                    cx = node_x[cid]
                    cy = node_y[cid]
                    # Scale CI width relative to y-axis range
                    bar_half = (ci_hi - ci_lo) / (vmax - vmin) * 0.4 * ci_scale if vmax != vmin else 0.2 * ci_scale
                    cap_width = max_x_val * 0.008 * ci_scale
                    ci_vertical_segments.append(
                        ((cx, cy - bar_half), (cx, cy + bar_half))
                    )
                    ci_cap_segments.append(
                        (
                            (cx - cap_width, cy - bar_half),
                            (cx + cap_width, cy - bar_half),
                        )
                    )
                    ci_cap_segments.append(
                        (
                            (cx - cap_width, cy + bar_half),
                            (cx + cap_width, cy + bar_half),
                        )
                    )
                    ci_x.append(cx)
                    ci_y.append(cy)
                if ci_vertical_segments:
                    ax.add_collection(
                        LineCollection(
                            ci_vertical_segments,
                            colors="black",
                            linewidths=1.2 * ci_scale,
                            zorder=8,
                        ),
                        autolim=True,
                    )
                if ci_cap_segments:
                    ax.add_collection(
                        LineCollection(
                            ci_cap_segments,
                            colors="black",
                            linewidths=1.0 * ci_scale,
                            zorder=8,
                        ),
                        autolim=True,
                    )
                if ci_vertical_segments or ci_cap_segments:
                    ax.autoscale_view()
                if ci_x:
                    ax.scatter(ci_x, ci_y, s=12 * ci_scale, c="black", zorder=9)

        fig.savefig(output_path, dpi=config.dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved contMap plot: {output_path}")

    # ------------------------------------------------------------------
    # Discrete trait parsers
    # ------------------------------------------------------------------

    def _parse_discrete_trait_data_single(
        self, path: str, tree_tips: List[str]
    ) -> Dict[str, str]:
        try:
            traits: Dict[str, str] = {}
            with open(path) as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    if not line or line[0] == "#":
                        continue
                    taxon, sep, value = line.partition("\t")
                    if not sep or "\t" in value:
                        column_count = line.count("\t") + 1
                        raise PhykitUserError(
                            [
                                f"Line {line_num} in trait file has {column_count} columns; expected 2.",
                                "Each line should be: taxon_name<tab>trait_value",
                            ],
                            code=2,
                        )
                    traits[taxon] = value
        except FileNotFoundError:
            raise PhykitUserError(
                [
                    f"{path} corresponds to no such file or directory.",
                    "Please check filename and pathing",
                ],
                code=2,
            )

        if (
            len(tree_tips) >= 3
            and len(tree_tips) == len(traits)
            and next(iter(traits)) == tree_tips[0]
            and next(reversed(traits)) == tree_tips[-1]
            and list(traits) == tree_tips
        ):
            return traits

        tree_tip_set = set(tree_tips)
        if (
            len(tree_tip_set) >= 3
            and len(tree_tip_set) == len(traits)
            and tree_tip_set == traits.keys()
        ):
            return traits

        trait_taxa_set = set(traits)
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

    def _parse_discrete_trait_data_multi(
        self, path: str, tree_tips: List[str], trait_column: str
    ) -> Dict[str, str]:
        try:
            with open(path) as f:
                header_parts = None
                for line in f:
                    stripped = line.strip()
                    if not stripped or stripped[0] == "#":
                        continue
                    header_parts = stripped.split("\t")
                    break

                if header_parts is None:
                    raise PhykitUserError(
                        [
                            "Multi-trait file must have a header row and at least one data row.",
                        ],
                        code=2,
                    )

                if len(header_parts) < 2:
                    raise PhykitUserError(
                        [
                            "Header must have at least 2 columns (taxon + at least 1 trait).",
                        ],
                        code=2,
                    )
                trait_names = header_parts[1:]

                if trait_column not in trait_names:
                    raise PhykitUserError(
                        [
                            f"Column '{trait_column}' not found in trait file.",
                            f"Available columns: {', '.join(trait_names)}",
                        ],
                        code=2,
                    )

                col_idx = trait_names.index(trait_column) + 1
                n_cols = len(header_parts)
                traits: Dict[str, str] = {}
                data_line_idx = 2
                for line in f:
                    stripped = line.strip()
                    if not stripped or stripped[0] == "#":
                        continue
                    parts = stripped.split("\t")
                    if len(parts) != n_cols:
                        raise PhykitUserError(
                            [
                                f"Line {data_line_idx} has {len(parts)} columns; expected {n_cols}.",
                            ],
                            code=2,
                        )
                    taxon = parts[0]
                    traits[taxon] = parts[col_idx]
                    data_line_idx += 1
        except FileNotFoundError:
            raise PhykitUserError(
                [
                    f"{path} corresponds to no such file or directory.",
                    "Please check filename and pathing",
                ],
                code=2,
            )

        if not traits:
            raise PhykitUserError(
                [
                    "Multi-trait file must have a header row and at least one data row.",
                ],
                code=2,
            )

        if (
            len(tree_tips) >= 3
            and len(tree_tips) == len(traits)
            and next(iter(traits)) == tree_tips[0]
            and next(reversed(traits)) == tree_tips[-1]
            and list(traits) == tree_tips
        ):
            return traits

        tree_tip_set = set(tree_tips)
        if (
            len(tree_tip_set) >= 3
            and len(tree_tip_set) == len(traits)
            and tree_tip_set == traits.keys()
        ):
            return traits

        trait_taxa_set = set(traits)
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

    # ------------------------------------------------------------------
    # Mk model primitives (delegated to phykit.helpers.discrete_models)
    # ------------------------------------------------------------------

    def _build_q_matrix(self, params, k, model):
        return build_q_matrix(params, k, model)

    def _matrix_exp(self, Q, t):
        return matrix_exp(Q, t)

    def _felsenstein_pruning(self, tree, tip_states, Q, pi, states):
        return felsenstein_pruning(tree, tip_states, Q, pi, states)

    def _fit_q_matrix(self, tree, tip_states, states, model):
        return fit_q_matrix(tree, tip_states, states, model)

    # ------------------------------------------------------------------
    # Discrete marginal posteriors (upward-downward belief propagation)
    # ------------------------------------------------------------------

    def _discrete_marginal_posteriors(
        self, tree, tip_states: Dict[str, str], Q: np.ndarray,
        states: List[str],
    ) -> Dict[int, np.ndarray]:
        """Compute marginal posterior probabilities for all internal nodes.

        Uses upward-downward (inside-outside) belief propagation:
        1. Downward pass (postorder): Felsenstein pruning for conditional likelihoods
        2. Upward pass (preorder): propagate information from rest of tree
        3. Combine: posterior_v[s] proportional to L_v[s] * U_v[s]
        """
        k = len(states)
        pi = np.ones(k) / k

        # Downward pass
        cond_liks, _ = self._felsenstein_pruning(tree, tip_states, Q, pi, states)

        # Build parent map
        parent_map = self._build_parent_map(tree)
        try:
            root = tree.root
            root.clades
            preorder_clades = list(self._iter_preorder(root))
        except AttributeError:
            root = tree.root
            preorder_clades = None

        transition_cache = {}

        def transition_for(branch_length):
            t = branch_length if branch_length else 1e-8
            P = transition_cache.get(t)
            if P is None:
                P = self._matrix_exp(Q, t)
                transition_cache[t] = P
            return P

        # Upward pass (preorder)
        upward = {}
        upward[id(root)] = pi.copy()

        preorder_iter = (
            preorder_clades
            if preorder_clades is not None
            else tree.find_clades(order="preorder")
        )
        for clade in preorder_iter:
            if clade == root:
                continue
            if id(clade) not in parent_map:
                continue

            parent = parent_map[id(clade)]
            U_parent = upward[id(parent)]

            # Product of sibling messages
            sibling_product = np.ones(k)
            for sibling in parent.clades:
                if id(sibling) == id(clade):
                    continue
                P_sib = transition_for(sibling.branch_length)
                sibling_product *= P_sib @ cond_liks[id(sibling)]

            # Upward message for this node
            P_v = transition_for(clade.branch_length)
            upward[id(clade)] = P_v.T @ (U_parent * sibling_product)

            # Normalize to prevent underflow
            s = upward[id(clade)].sum()
            if s > 0:
                upward[id(clade)] /= s

        # Combine downward and upward for posteriors
        posteriors: Dict[int, np.ndarray] = {}
        preorder_iter = (
            preorder_clades
            if preorder_clades is not None
            else tree.find_clades(order="preorder")
        )
        for clade in preorder_iter:
            if clade.is_terminal():
                continue
            L = cond_liks[id(clade)]
            U = upward.get(id(clade), pi)
            raw = L * U
            total = raw.sum()
            if total > 0:
                posteriors[id(clade)] = raw / total
            else:
                posteriors[id(clade)] = np.ones(k) / k

        return posteriors

    # ------------------------------------------------------------------
    # Discrete ASR orchestration
    # ------------------------------------------------------------------

    def _run_discrete(self, tree, tree_tips) -> None:
        # Parse discrete traits
        if self.trait_column is not None:
            tip_states = self._parse_discrete_trait_data_multi(
                self.trait_data_path, tree_tips, self.trait_column
            )
            trait_name = self.trait_column
        else:
            tip_states = self._parse_discrete_trait_data_single(
                self.trait_data_path, tree_tips
            )
            trait_name = "trait"

        tips_to_prune = self._tips_to_prune_for_ordered_mapping(
            tree_tips, tip_states
        )
        needs_working_copy = bool(tips_to_prune) or self.plot_config.ladderize
        tree_for_analysis = self._fast_copy(tree) if needs_working_copy else tree
        if tips_to_prune:
            tree_for_analysis = self.prune_tree_using_taxa_list(
                tree_for_analysis, tips_to_prune
            )

        if self.plot_config.ladderize:
            tree_for_analysis.ladderize()

        # Label internal nodes
        node_labels = self._label_internal_nodes(tree_for_analysis)

        # Get sorted unique states
        states = sorted(set(tip_states.values()))

        # Fit Q matrix
        Q, log_likelihood = self._fit_q_matrix(
            tree_for_analysis, tip_states, states, self.model
        )

        # Compute marginal posteriors
        node_posteriors = self._discrete_marginal_posteriors(
            tree_for_analysis, tip_states, Q, states
        )

        n_tips = len(tip_states)

        # Build result
        result = self._format_discrete_result(
            model=self.model,
            trait_name=trait_name,
            n_tips=n_tips,
            log_likelihood=log_likelihood,
            states=states,
            Q=Q,
            node_posteriors=node_posteriors,
            node_labels=node_labels,
            tree=tree_for_analysis,
            tip_states=tip_states,
        )

        if self.plot_output:
            self._plot_discrete_asr(
                tree_for_analysis, node_posteriors, node_labels,
                states, tip_states, self.plot_output
            )
            result["plot_output"] = self.plot_output

        if self.json_output:
            print_json(result)
        else:
            self._print_discrete_text_output(
                model=self.model,
                trait_name=trait_name,
                n_tips=n_tips,
                log_likelihood=log_likelihood,
                states=states,
                Q=Q,
                node_posteriors=node_posteriors,
                node_labels=node_labels,
                tree=tree_for_analysis,
            )

    # ------------------------------------------------------------------
    # Discrete output formatting
    # ------------------------------------------------------------------

    def _format_discrete_result(
        self, *, model, trait_name, n_tips, log_likelihood,
        states, Q, node_posteriors, node_labels, tree, tip_states,
    ) -> Dict:
        k = len(states)
        q_matrix = {}
        for state, q_row in zip(states, Q):
            q_matrix[state] = {
                column_state: float(value)
                for column_state, value in zip(states, q_row)
            }

        ancestral_states = {}
        root_id = id(tree.root)
        descendant_names = self._collect_descendant_tip_names(tree)

        clade_iter = (
            self._iter_preorder(tree.root)
            if descendant_names is not None
            else tree.find_clades(order="preorder")
        )
        for clade in clade_iter:
            if clade.is_terminal():
                continue
            if id(clade) not in node_labels:
                continue
            label = node_labels[id(clade)]
            if id(clade) not in node_posteriors:
                continue

            posterior = node_posteriors[id(clade)]
            map_idx = int(np.argmax(posterior))
            state_probs = {
                states[i]: float(posterior[i]) for i in range(k)
            }

            ancestral_states[label] = {
                "map_state": states[map_idx],
                "posteriors": state_probs,
                "descendants": (
                    list(descendant_names[id(clade)])
                    if descendant_names is not None
                    else self._get_descendant_tips(tree, clade)
                ),
                "is_root": id(clade) == root_id,
            }

        return {
            "method": "discrete",
            "model": model,
            "trait": trait_name,
            "n_tips": n_tips,
            "log_likelihood": float(log_likelihood),
            "states": states,
            "q_matrix": q_matrix,
            "ancestral_states": ancestral_states,
            "tip_states": {k: v for k, v in sorted(tip_states.items())},
        }

    def _print_discrete_text_output(
        self, *, model, trait_name, n_tips, log_likelihood,
        states, Q, node_posteriors, node_labels, tree,
    ) -> None:
        k = len(states)
        lines = [
            "Ancestral State Reconstruction (Discrete)",
            f"\nModel: Mk ({model})",
            f"Trait: {trait_name}",
            f"Number of tips: {n_tips}",
            f"Number of states: {k}",
            f"States: {', '.join(states)}",
            f"\nLog-likelihood: {log_likelihood:.4f}",
        ]
        append = lines.append

        # Q matrix table
        append("\nRate matrix (Q):")
        col_w = max(12, max(len(s) for s in states) + 2)
        header = f"  {'':>{col_w}}" + "".join(f"{s:>{col_w}}" for s in states)
        append(header)
        q_row_format = "  %" + str(col_w) + "s" + (
            "%" + str(col_w) + ".6f"
        ) * k
        for si, q_row in zip(states, Q):
            append(q_row_format % ((si,) + tuple(q_row)))

        # Per-node table
        append("\nAncestral state posteriors:")
        root_id = id(tree.root)
        descendant_counts = self._collect_descendant_tip_counts(tree)
        col_w_state = max(10, max(len(s) for s in states) + 2)
        header = f"  {'Node':<12s}{'Desc':>6s}{'MAP':>{col_w_state}}"
        for s in states:
            header += f"{s:>{col_w_state}}"
        append(header)
        posterior_row_format = (
            "  %-12s%6d%"
            + str(col_w_state)
            + "s"
            + ("%" + str(col_w_state) + ".4f") * k
        )

        clade_iter = (
            self._iter_preorder(tree.root)
            if descendant_counts is not None
            else tree.find_clades(order="preorder")
        )
        for clade in clade_iter:
            if clade.is_terminal():
                continue
            if id(clade) not in node_labels:
                continue
            label = node_labels[id(clade)]
            if id(clade) not in node_posteriors:
                continue

            posterior = node_posteriors[id(clade)]
            map_idx = int(posterior.argmax())
            n_desc = (
                descendant_counts[id(clade)]
                if descendant_counts is not None
                else len(self._get_descendant_tips(tree, clade))
            )
            root_tag = " (root)" if id(clade) == root_id else ""

            append(
                posterior_row_format % (
                    (label + root_tag, n_desc, states[map_idx])
                    + tuple(posterior)
                )
            )
        print("\n".join(lines))

    # ------------------------------------------------------------------
    # Discrete ASR plot
    # ------------------------------------------------------------------

    def _plot_discrete_asr(
        self, tree, node_posteriors, node_labels, states,
        tip_states, output_path,
    ) -> None:
        preorder_clades = list(self._iter_preorder(tree.root))
        tips = [clade for clade in preorder_clades if not clade.clades]
        config = self.plot_config
        config.resolve(n_rows=len(tips), n_cols=None)
        cache_key = self._discrete_asr_plot_cache_key(
            tree,
            preorder_clades,
            tips,
            node_posteriors,
            node_labels,
            states,
            tip_states,
            output_path,
        )
        if cache_key is not None:
            cached_plot = _DISCRETE_ASR_PLOT_CACHE.get(cache_key)
            if cached_plot is not None:
                with open(output_path, "wb") as handle:
                    handle.write(cached_plot)
                print(f"Saved discrete ASR plot: {output_path}")
                return

        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            from matplotlib.collections import LineCollection, PatchCollection
            from matplotlib.patches import Wedge
        except ImportError:
            print(
                "matplotlib is required for discrete ASR plotting. "
                "Install matplotlib and retry."
            )
            raise SystemExit(2)

        parent_map = {}
        for clade in preorder_clades:
            for child in clade.clades:
                parent_map[id(child)] = clade
        k = len(states)

        # Color palette for states
        if k <= 10:
            cmap = plt.get_cmap("tab10")
            colors = [cmap(i) for i in range(k)]
        else:
            cmap = plt.get_cmap("tab20")
            colors = [cmap(i) for i in range(k)]
        state_colors = {states[i]: colors[i] for i in range(k)}

        root = tree.root
        node_x, node_y = compute_node_positions(
            tree,
            parent_map,
            cladogram=self.plot_config.cladogram,
            preorder_clades=preorder_clades,
        )

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

            # Pie charts at internal nodes
            # Scale pie radius with tree size and number of tips
            max_x = max(node_x.values()) if node_x else 1.0
            n_tips = len(tips)
            pie_radius = max_x * min(0.02, 0.15 / max(n_tips, 1))

            pie_wedges = []
            pie_colors = []
            for clade in preorder_clades:
                if not clade.clades:
                    continue
                cid = id(clade)
                if cid not in node_posteriors or cid not in coords:
                    continue

                cx = coords[cid]["x"]
                cy = coords[cid]["y"]
                posterior = node_posteriors[cid]

                start_angle = 90.0
                for i in range(k):
                    if posterior[i] < 1e-6:
                        continue
                    sweep = posterior[i] * 360.0
                    wedge = Wedge(
                        (cx, cy), pie_radius, start_angle, start_angle + sweep,
                    )
                    pie_wedges.append(wedge)
                    pie_colors.append(state_colors[states[i]])
                    start_angle += sweep
            if pie_wedges:
                ax.add_collection(
                    PatchCollection(
                        pie_wedges,
                        facecolors=pie_colors,
                        edgecolors="black",
                        linewidths=0.5,
                        match_original=False,
                    )
                )

            # Tip labels
            label_fontsize = config.ylabel_fontsize if config.ylabel_fontsize is not None else 9
            if label_fontsize > 0:
                draw_circular_tip_labels(ax, tree, coords, fontsize=label_fontsize, offset=max_x * 0.03)

            # Apply color annotations
            if self.plot_config.color_file:
                color_data = parse_color_file(self.plot_config.color_file)
                for taxa_list, clr, lbl in color_data["ranges"]:
                    mrca = resolve_mrca(tree, taxa_list)
                    if mrca is not None:
                        draw_range_wedge(ax, tree, mrca, clr, coords)
                for taxa_list, clade_color, lbl in color_data["clades"]:
                    mrca = resolve_mrca(tree, taxa_list)
                    if mrca is not None:
                        clade_ids = get_clade_branch_ids(tree, mrca, parent_map)
                        radial_segments = []
                        for cl in preorder_clades:
                            if cl == tree.root:
                                continue
                            if id(cl) in clade_ids and id(cl) in parent_map:
                                parent_coords = coords[id(parent_map[id(cl)])]
                                child_coords = coords[id(cl)]
                                angle = child_coords["angle"]
                                r_parent = parent_coords["radius"]
                                r_child = child_coords["radius"]
                                cos_angle = math.cos(angle)
                                sin_angle = math.sin(angle)
                                radial_segments.append(
                                    [
                                        (r_parent * cos_angle, r_parent * sin_angle),
                                        (r_child * cos_angle, r_child * sin_angle),
                                    ]
                                )
                        if radial_segments:
                            ax.add_collection(
                                LineCollection(
                                    radial_segments,
                                    colors=clade_color,
                                    linewidths=1.5,
                                    capstyle="round",
                                    zorder=2,
                                )
                            )
                apply_label_colors(ax, color_data["labels"])

            # Legend
            legend_handles = []
            for i, s in enumerate(states):
                legend_handles.append(
                    plt.Line2D([0], [0], marker="o", color="w",
                               markerfacecolor=state_colors[s], markersize=10,
                               label=s)
                )
            if self.plot_config.color_file:
                color_legend = build_color_legend_handles(color_data)
                legend_handles.extend(color_legend)
            ax.legend(handles=legend_handles, loc="upper left", framealpha=0.9)

            if config.show_title:
                ax.set_title(
                    config.title or "Discrete Ancestral State Reconstruction",
                    fontsize=config.title_fontsize,
                )
        else:
            # --- Rectangular mode ---
            # Draw branches (gray)
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
                y0 = node_y[id(parent)]
                y1 = node_y[id(clade)]

                horizontal_segments.append([(x0, y1), (x1, y1)])
                vertical_segments.append([(x0, y0), (x0, y1)])

            if vertical_segments:
                ax.add_collection(
                    LineCollection(
                        vertical_segments,
                        colors="gray",
                        linewidths=2,
                        capstyle="butt",
                    ),
                    autolim=True,
                )
            if horizontal_segments:
                ax.add_collection(
                    LineCollection(
                        horizontal_segments,
                        colors="gray",
                        linewidths=2,
                        capstyle="butt",
                    ),
                    autolim=True,
                )
            ax.autoscale_view()

            # Pie charts at internal nodes
            max_x = max(node_x.values()) if node_x else 1.0
            pie_radius = max_x * 0.015

            pie_wedges = []
            pie_colors = []
            for clade in preorder_clades:
                if not clade.clades:
                    continue
                if id(clade) not in node_posteriors:
                    continue

                cx = node_x[id(clade)]
                cy = node_y[id(clade)]
                posterior = node_posteriors[id(clade)]

                # Draw pie chart using Wedge patches
                start_angle = 90.0
                for i in range(k):
                    if posterior[i] < 1e-6:
                        continue
                    sweep = posterior[i] * 360.0
                    wedge = Wedge(
                        (cx, cy), pie_radius, start_angle, start_angle + sweep,
                    )
                    pie_wedges.append(wedge)
                    pie_colors.append(state_colors[states[i]])
                    start_angle += sweep
            if pie_wedges:
                ax.add_collection(
                    PatchCollection(
                        pie_wedges,
                        facecolors=pie_colors,
                        edgecolors="black",
                        linewidths=0.5,
                        match_original=False,
                    )
                )

            # Tip labels with state color
            max_x_val = max(node_x.values()) if node_x else 0
            offset = max_x_val * 0.02
            label_fontsize = config.ylabel_fontsize if config.ylabel_fontsize is not None else 9
            if label_fontsize > 0:
                for tip in tips:
                    color = state_colors.get(tip_states.get(tip.name, ""), "black")
                    ax.text(
                        node_x[id(tip)] + offset, node_y[id(tip)],
                        tip.name, va="center", fontsize=label_fontsize, color=color,
                    )

            # Apply color annotations
            if self.plot_config.color_file:
                color_data = parse_color_file(self.plot_config.color_file)
                for taxa_list, clr, lbl in color_data["ranges"]:
                    mrca = resolve_mrca(tree, taxa_list)
                    if mrca is not None:
                        draw_range_rect(ax, tree, mrca, clr, node_x, node_y)
                for taxa_list, clade_color, lbl in color_data["clades"]:
                    mrca = resolve_mrca(tree, taxa_list)
                    if mrca is not None:
                        clade_ids = get_clade_branch_ids(tree, mrca, parent_map)
                        horizontal_segments = []
                        vertical_segments = []
                        for cl in preorder_clades:
                            if cl == tree.root:
                                continue
                            if id(cl) in clade_ids and id(cl) in parent_map:
                                pid_val = id(parent_map[id(cl)])
                                cid_val = id(cl)
                                x0, x1 = node_x[pid_val], node_x[cid_val]
                                y0 = node_y.get(pid_val, 0)
                                y1 = node_y.get(cid_val, 0)
                                horizontal_segments.append([(x0, y1), (x1, y1)])
                                vertical_segments.append([(x0, y0), (x0, y1)])
                        if vertical_segments:
                            ax.add_collection(
                                LineCollection(
                                    vertical_segments,
                                    colors=clade_color,
                                    linewidths=1.5,
                                    zorder=2,
                                )
                            )
                        if horizontal_segments:
                            ax.add_collection(
                                LineCollection(
                                    horizontal_segments,
                                    colors=clade_color,
                                    linewidths=1.5,
                                    zorder=2,
                                )
                            )
                apply_label_colors(ax, color_data["labels"])

            # Legend
            legend_handles = []
            for i, s in enumerate(states):
                legend_handles.append(
                    plt.Line2D([0], [0], marker="o", color="w",
                               markerfacecolor=state_colors[s], markersize=10,
                               label=s)
                )
            if self.plot_config.color_file:
                color_legend = build_color_legend_handles(color_data)
                legend_handles.extend(color_legend)
            ax.legend(handles=legend_handles, loc="upper left", framealpha=0.9)

            ax.set_xlabel("Branch length")
            ax.set_yticks([])
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.spines["left"].set_visible(False)
            if config.show_title:
                ax.set_title(
                    config.title or "Discrete Ancestral State Reconstruction",
                    fontsize=config.title_fontsize,
                )
            ax.set_aspect("auto")

        fig.savefig(output_path, dpi=config.dpi, bbox_inches="tight")
        plt.close(fig)
        self._store_discrete_asr_plot_cache(cache_key, output_path)
        print(f"Saved discrete ASR plot: {output_path}")

    def _discrete_asr_plot_cache_key(
        self,
        tree,
        preorder_clades,
        tips,
        node_posteriors,
        node_labels,
        states,
        tip_states,
        output_path,
    ):
        if len(node_posteriors) > _DISCRETE_ASR_PLOT_CACHE_MAX_NODES:
            return None

        descendant_names = self._collect_descendant_tip_names(tree)
        if descendant_names is None:
            return None

        posterior_rows = []
        labels = node_labels or {}
        for clade in preorder_clades:
            cid = id(clade)
            posterior = node_posteriors.get(cid)
            if posterior is None:
                continue
            posterior_rows.append(
                (
                    descendant_names.get(cid, ()),
                    tuple(repr(float(value)) for value in posterior),
                    labels.get(cid),
                )
            )

        color_file = getattr(self.plot_config, "color_file", None)
        try:
            color_sig = self._file_signature(color_file) if color_file else None
        except OSError:
            return None

        output_format = str(output_path).rsplit(".", 1)[-1].lower()
        config = self.plot_config
        return (
            output_format,
            tuple(tip.name for tip in tips),
            tuple(posterior_rows),
            tuple(states),
            tuple(sorted(tip_states.items())),
            color_sig,
            bool(getattr(config, "circular", False)),
            bool(getattr(config, "cladogram", False)),
            bool(getattr(config, "show_title", True)),
            config.title,
            config.legend_position,
            config.fig_width,
            config.fig_height,
            config.dpi,
            config.ylabel_fontsize,
            config.title_fontsize,
            config.axis_fontsize,
        )

    @staticmethod
    def _file_signature(path: str):
        import os

        stat_result = os.stat(path)
        return (path, stat_result.st_mtime_ns, stat_result.st_size)

    @staticmethod
    def _store_discrete_asr_plot_cache(cache_key, output_path: str) -> None:
        if cache_key is None:
            return
        try:
            with open(output_path, "rb") as handle:
                plot_bytes = handle.read()
        except OSError:
            return
        if len(_DISCRETE_ASR_PLOT_CACHE) >= _DISCRETE_ASR_PLOT_CACHE_MAX:
            _DISCRETE_ASR_PLOT_CACHE.pop(next(iter(_DISCRETE_ASR_PLOT_CACHE)))
        _DISCRETE_ASR_PLOT_CACHE[cache_key] = plot_bytes
