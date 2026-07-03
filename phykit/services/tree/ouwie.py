from __future__ import annotations

import math
import sys

from .base import Tree
from ...errors import PhykitUserError

ALL_MODELS = ["BM1", "BMS", "OU1", "OUM", "OUMV", "OUMA", "OUMVA"]


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


def _same_ordered_keys(left, right) -> bool:
    if len(left) != len(right):
        return False
    for left_key, right_key in zip(left, right):
        if left_key != right_key:
            return False
    return True


class _LazyPickle:
    def dumps(self, *args, **kwargs):
        import pickle as _pickle

        return _pickle.dumps(*args, **kwargs)

    def loads(self, *args, **kwargs):
        import pickle as _pickle

        return _pickle.loads(*args, **kwargs)

    def __getattr__(self, name):
        import pickle as _pickle

        return getattr(_pickle, name)


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


pickle = _LazyPickle()
np = _LazyNumpy()
_CHO_FACTOR = None
_CHO_SOLVE = None
_MINIMIZE = None
_MINIMIZE_SCALAR = None


def cho_factor(*args, **kwargs):
    global _CHO_FACTOR
    if _CHO_FACTOR is None:
        from scipy.linalg import cho_factor as _CHO_FACTOR

    return _CHO_FACTOR(*args, **kwargs)


def cho_solve(*args, **kwargs):
    global _CHO_SOLVE
    if _CHO_SOLVE is None:
        from scipy.linalg import cho_solve as _CHO_SOLVE

    return _CHO_SOLVE(*args, **kwargs)


def minimize(*args, **kwargs):
    global _MINIMIZE
    if _MINIMIZE is None:
        from scipy.optimize import minimize as _MINIMIZE

    return _MINIMIZE(*args, **kwargs)


def minimize_scalar(*args, **kwargs):
    global _MINIMIZE_SCALAR
    if _MINIMIZE_SCALAR is None:
        from scipy.optimize import minimize_scalar as _MINIMIZE_SCALAR

    return _MINIMIZE_SCALAR(*args, **kwargs)


def _merge_nonempty_child_state_sets(state_sets, children):
    intersection = None
    union = None
    for child in children:
        child_set = state_sets.get(id(child))
        if not child_set:
            continue
        if intersection is None:
            intersection = child_set
            union = child_set
        else:
            intersection = intersection & child_set
            union = union | child_set

    if intersection is None:
        return set()
    return intersection if intersection else union


class OUwie(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"])
        self.trait_data_path = parsed["trait_data_path"]
        self.regime_data_path = parsed["regime_data_path"]
        self.selected_models = parsed["models"]
        self.json_output = parsed["json_output"]

    def run(self) -> None:
        tree = self.read_tree_file_unmodified()
        self.validate_tree(tree, min_tips=3, require_branch_lengths=True, context="OUwie model fitting")

        tree_tips = self.get_tip_names_from_tree(tree)
        traits = self._parse_trait_file(self.trait_data_path, tree_tips)
        regime_assignments = self._parse_regime_file(
            self.regime_data_path, tree_tips
        )

        (
            traits,
            regime_assignments,
            tips_to_prune,
            ordered_names,
            regimes,
        ) = self._prepare_shared_trait_regime_data(
            tree_tips,
            traits,
            regime_assignments,
        )
        tree_for_analysis = self._fast_copy(tree) if tips_to_prune else tree
        if tips_to_prune:
            tree_for_analysis = self.prune_tree_using_taxa_list(
                tree_for_analysis, tips_to_prune
            )

        n = len(ordered_names)
        x = np.array([traits[name] for name in ordered_names])

        n_regimes = len(regimes)

        if n_regimes < 2:
            raise PhykitUserError(
                ["At least 2 distinct regimes are required for OUwie models."],
                code=2,
            )

        parent_map = self._build_parent_map(tree_for_analysis)
        branch_regimes, root_regime = self._assign_branch_regimes_and_root(
            tree_for_analysis, regime_assignments, parent_map
        )

        lineage_info = self._build_lineage_info(
            tree_for_analysis, ordered_names, parent_map, branch_regimes
        )

        per_regime_vcv = self._build_per_regime_vcv(
            tree_for_analysis,
            ordered_names,
            regimes,
            branch_regimes,
            parent_map,
            lineage_info=lineage_info,
        )

        vcv_total = sum(per_regime_vcv.values())
        tree_height = float(np.diagonal(vcv_total).max())

        results = []
        for model_name in self.selected_models:
            res = self._fit_model(
                model_name, x, vcv_total, per_regime_vcv,
                ordered_names, regimes, lineage_info,
                root_regime, tree_height, n_regimes,
            )
            results.append(res)

        results = self._compute_model_comparison(
            results, n, x, vcv_total, regime_assignments, ordered_names
        )

        if self.json_output:
            self._print_json_output(results, n, regimes)
        else:
            self._print_text_output(results, n, regimes)

    @staticmethod
    def _shared_trait_regime_taxa(
        traits: dict[str, float],
        regime_assignments: dict[str, str],
    ) -> set[str]:
        if len(regime_assignments) < len(traits):
            return {name for name in regime_assignments if name in traits}
        return traits.keys() & regime_assignments.keys()

    @staticmethod
    def _prepare_shared_trait_regime_data(
        tree_tips,
        traits,
        regime_assignments,
    ):
        if _same_ordered_keys(traits, regime_assignments):
            if len(traits) < 3:
                raise PhykitUserError(
                    [
                        "Only "
                        f"{len(traits)} shared taxa among tree, trait, and regime files.",
                        "At least 3 shared taxa are required.",
                    ],
                    code=2,
                )
            tips_to_prune = Tree._tips_to_prune_for_ordered_mapping(
                tree_tips, traits
            )
            ordered_names = sorted(traits)
            regimes = sorted(set(regime_assignments.values()))
            return (
                traits,
                regime_assignments,
                tips_to_prune,
                ordered_names,
                regimes,
            )

        shared = OUwie._shared_trait_regime_taxa(traits, regime_assignments)
        if len(shared) < 3:
            raise PhykitUserError(
                [
                    f"Only {len(shared)} shared taxa among tree, trait, and regime files.",
                    "At least 3 shared taxa are required.",
                ],
                code=2,
            )

        if len(shared) == len(traits) == len(regime_assignments):
            shared_traits = traits
            shared_regime_assignments = regime_assignments
        else:
            shared_traits = {name: traits[name] for name in shared}
            shared_regime_assignments = {
                name: regime_assignments[name]
                for name in shared
            }

        tips_to_prune = [name for name in tree_tips if name not in shared]
        ordered_names = sorted(shared)
        regimes = sorted(set(shared_regime_assignments.values()))
        return (
            shared_traits,
            shared_regime_assignments,
            tips_to_prune,
            ordered_names,
            regimes,
        )

    def process_args(self, args) -> dict:
        models = ALL_MODELS[:]
        if hasattr(args, "models") and args.models:
            requested = [m.strip() for m in args.models.split(",")]
            valid = {m.lower(): m for m in ALL_MODELS}
            models = []
            for m in requested:
                canonical = valid.get(m.lower())
                if canonical is None:
                    raise PhykitUserError(
                        [
                            f"Unknown model '{m}'.",
                            f"Available models: {', '.join(ALL_MODELS)}",
                        ],
                        code=2,
                    )
                models.append(canonical)
        return dict(
            tree_file_path=args.tree,
            trait_data_path=args.trait_data,
            regime_data_path=args.regime_data,
            models=models,
            json_output=getattr(args, "json", False),
        )

    # ── Tree & data parsing ──────────────────────────────────────────

    def _parse_trait_file(
        self, path: str, tree_tips: list[str]
    ) -> dict[str, float]:
        try:
            traits = {}
            with open(path) as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    if not line or line[0] == "#":
                        continue
                    parts = line.split("\t", 2)
                    if len(parts) != 2:
                        column_count = line.count("\t") + 1
                        raise PhykitUserError(
                            [
                                f"Line {line_num} in trait file has {column_count} columns; expected 2.",
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
        except FileNotFoundError:
            raise PhykitUserError(
                [
                    f"{path} corresponds to no such file or directory.",
                    "Please check filename and pathing",
                ],
                code=2,
            )

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

    def _parse_regime_file(
        self, path: str, tree_tips: list[str]
    ) -> dict[str, str]:
        try:
            regimes = {}
            with open(path) as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    if not line or line[0] == "#":
                        continue
                    parts = line.split("\t", 2)
                    if len(parts) != 2:
                        column_count = line.count("\t") + 1
                        raise PhykitUserError(
                            [
                                f"Line {line_num} in regime file has {column_count} columns; expected 2.",
                                "Each line should be: taxon_name<tab>regime_label",
                            ],
                            code=2,
                        )
                    taxon, regime = parts
                    regimes[taxon] = regime
        except FileNotFoundError:
            raise PhykitUserError(
                [
                    f"{path} corresponds to no such file or directory.",
                    "Please check filename and pathing",
                ],
                code=2,
            )

        tree_tip_set = set(tree_tips)
        if (
            len(tree_tip_set) >= 3
            and len(tree_tip_set) == len(regimes)
            and tree_tip_set == regimes.keys()
        ):
            return regimes

        regime_taxa_set = set(regimes)
        shared = tree_tip_set & regime_taxa_set

        tree_only = tree_tip_set - regime_taxa_set
        regime_only = regime_taxa_set - tree_tip_set

        if tree_only:
            print(
                f"Warning: {len(tree_only)} taxa in tree but not in regime file: "
                f"{', '.join(sorted(tree_only))}",
                file=sys.stderr,
            )
        if regime_only:
            print(
                f"Warning: {len(regime_only)} taxa in regime file but not in tree: "
                f"{', '.join(sorted(regime_only))}",
                file=sys.stderr,
            )

        if len(shared) < 3:
            raise PhykitUserError(
                [
                    f"Only {len(shared)} shared taxa between tree and regime file.",
                    "At least 3 shared taxa are required.",
                ],
                code=2,
            )

        return {taxon: regimes[taxon] for taxon in shared}

    # ── Tree structure helpers ────────────────────────────────────────

    def _build_parent_map(self, tree) -> dict:
        try:
            stack = [tree.root]
            parent_map = {}
            pop = stack.pop
            extend = stack.extend
            while stack:
                clade = pop()
                children = clade.clades
                for child in children:
                    parent_map[id(child)] = clade
                if children:
                    extend(children)
        except AttributeError:
            parent_map = {}
            for clade in tree.find_clades(order="preorder"):
                for child in clade.clades:
                    parent_map[id(child)] = clade
        return parent_map

    def _assign_branch_regimes(
        self, tree, tip_regimes: dict[str, str], parent_map: dict
    ) -> dict:
        branch_regimes, _ = self._assign_branch_regimes_and_root(
            tree, tip_regimes, parent_map
        )
        return branch_regimes

    def _assign_branch_regimes_and_root(
        self, tree, tip_regimes: dict[str, str], parent_map: dict
    ) -> tuple[dict, str]:
        state_sets = {}
        preorder = []
        postorder = []
        try:
            root = tree.root
            stack = [(root, False)]
            while stack:
                clade, visited = stack.pop()
                if visited:
                    postorder.append(clade)
                    continue
                preorder.append(clade)
                stack.append((clade, True))
                for child in reversed(clade.clades):
                    stack.append((child, False))
        except AttributeError:
            root = tree.root
            postorder = list(tree.find_clades(order="postorder"))
            preorder = list(tree.find_clades(order="preorder"))

        for clade in postorder:
            if clade.is_terminal():
                if clade.name in tip_regimes:
                    state_sets[id(clade)] = {tip_regimes[clade.name]}
                else:
                    state_sets[id(clade)] = set()
            else:
                children = clade.clades
                if len(children) == 2:
                    left = state_sets.get(id(children[0]))
                    right = state_sets.get(id(children[1]))
                    if left and right:
                        intersection = left & right
                        state_sets[id(clade)] = (
                            intersection if intersection else left | right
                        )
                    elif left:
                        state_sets[id(clade)] = left
                    elif right:
                        state_sets[id(clade)] = right
                    else:
                        state_sets[id(clade)] = set()
                else:
                    state_sets[id(clade)] = _merge_nonempty_child_state_sets(
                        state_sets,
                        children,
                    )

        root_set = state_sets.get(id(root), set())
        all_regimes = sorted(set(tip_regimes.values()))
        if root_set:
            root_regime = min(root_set)
        else:
            root_regime = all_regimes[0] if all_regimes else "unknown"

        node_regimes = {id(root): root_regime}
        for clade in preorder:
            if clade == root:
                continue
            clade_set = state_sets.get(id(clade), set())
            if not clade_set:
                parent = parent_map.get(id(clade))
                if parent and id(parent) in node_regimes:
                    node_regimes[id(clade)] = node_regimes[id(parent)]
                continue

            parent = parent_map.get(id(clade))
            if parent and id(parent) in node_regimes:
                parent_regime = node_regimes[id(parent)]
                if parent_regime in clade_set:
                    node_regimes[id(clade)] = parent_regime
                else:
                    node_regimes[id(clade)] = min(clade_set)
            else:
                node_regimes[id(clade)] = min(clade_set)

        branch_regimes = {}
        for clade in preorder:
            if clade == root:
                continue
            branch_regimes[id(clade)] = node_regimes.get(id(clade), "unknown")

        return branch_regimes, root_regime

    def _get_root_regime(self, tree, tip_regimes, parent_map):
        _, root_regime = self._assign_branch_regimes_and_root(
            tree, tip_regimes, parent_map
        )
        return root_regime

    def _build_root_to_tip_paths(
        self, tree, ordered_names: list[str], parent_map: dict
    ) -> dict[str, list]:
        tip_map = self._terminal_map_for_ordered_names(tree, ordered_names)

        paths = {}
        for name in ordered_names:
            tip = tip_map[name]
            path = []
            current = tip
            while id(current) in parent_map:
                bl = current.branch_length if current.branch_length else 0.0
                path.append((id(current), bl))
                current = parent_map[id(current)]
            path.reverse()
            paths[name] = path

        return paths

    def _build_per_regime_vcv(
        self, tree, ordered_names: list[str], regimes: list[str],
        branch_regimes: dict, parent_map: dict,
        lineage_info: dict[str, list[tuple]] | None = None,
    ) -> dict[str, np.ndarray]:
        if lineage_info is not None:
            return self._build_per_regime_vcv_from_lineage_info(
                ordered_names, regimes, lineage_info
            )

        n = len(ordered_names)
        paths = self._build_root_to_tip_paths(tree, ordered_names, parent_map)

        per_regime_vcv = {r: np.zeros((n, n)) for r in regimes}
        branch_to_indices = {}
        branch_lengths = {}
        for idx, name in enumerate(ordered_names):
            for clade_id, bl in paths[name]:
                branch_to_indices.setdefault(clade_id, []).append(idx)
                branch_lengths[clade_id] = bl

        for clade_id, indices in branch_to_indices.items():
            regime = branch_regimes.get(clade_id, regimes[0])
            if regime not in per_regime_vcv:
                continue
            branch_length = branch_lengths[clade_id]
            if len(indices) == 1:
                per_regime_vcv[regime][indices[0], indices[0]] += branch_length
                continue
            idx = np.asarray(indices, dtype=np.intp)
            per_regime_vcv[regime][np.ix_(idx, idx)] += branch_length

        return per_regime_vcv

    def _build_per_regime_vcv_from_lineage_info(
        self,
        ordered_names: list[str],
        regimes: list[str],
        lineage_info: dict[str, list[tuple]],
    ) -> dict[str, np.ndarray]:
        n = len(ordered_names)
        per_regime_vcv = {r: np.zeros((n, n)) for r in regimes}
        default_regime = regimes[0] if regimes else None
        branch_to_indices = {}
        branch_info = {}

        for idx, name in enumerate(ordered_names):
            for clade_id, bl, regime, _d_start, _d_end in lineage_info[name]:
                branch_to_indices.setdefault(clade_id, []).append(idx)
                branch_info[clade_id] = (bl, regime)

        for clade_id, indices in branch_to_indices.items():
            branch_length, regime = branch_info[clade_id]
            matrix = per_regime_vcv.get(regime)
            if matrix is None:
                if regime == "unknown" and default_regime is not None:
                    matrix = per_regime_vcv[default_regime]
                else:
                    continue
            if len(indices) == 1:
                matrix[indices[0], indices[0]] += branch_length
                continue
            idx = np.asarray(indices, dtype=np.intp)
            matrix[np.ix_(idx, idx)] += branch_length

        return per_regime_vcv

    # ── Lineage info for OU weight/VCV matrices ──────────────────────

    def _build_lineage_info(
        self, tree, ordered_names: list[str], parent_map: dict,
        branch_regimes: dict,
    ) -> dict[str, list[tuple]]:
        """Build root-to-tip lineage info for each tip.

        Returns dict mapping tip_name -> list of tuples:
            (clade_id, branch_length, regime, dist_from_root_start, dist_from_root_end)
        """
        tip_map = self._terminal_map_for_ordered_names(tree, ordered_names)

        lineage_info = {}
        for name in ordered_names:
            tip = tip_map[name]
            path = []
            current = tip
            while id(current) in parent_map:
                bl = current.branch_length if current.branch_length else 0.0
                regime = branch_regimes.get(id(current), "unknown")
                path.append((id(current), bl, regime))
                current = parent_map[id(current)]
            path.reverse()

            # Add cumulative distances
            enriched = []
            cum = 0.0
            for clade_id, bl, regime in path:
                d_start = cum
                d_end = cum + bl
                enriched.append((clade_id, bl, regime, d_start, d_end))
                cum += bl

            lineage_info[name] = enriched

        return lineage_info

    @staticmethod
    def _terminal_map_for_ordered_names(tree, ordered_names):
        ordered_name_set = set(ordered_names)
        tip_map = {}

        try:
            root = tree.root
            root.clades
        except AttributeError:
            for tip in tree.get_terminals():
                if tip.name in ordered_name_set:
                    tip_map[tip.name] = tip
            return tip_map

        stack = [root]
        try:
            pop = stack.pop
            append = stack.append
            contains = ordered_name_set.__contains__
            while stack:
                clade = pop()
                children = clade.clades
                if children:
                    child_count = len(children)
                    if child_count == 2:
                        append(children[1])
                        append(children[0])
                    else:
                        for idx in range(child_count - 1, -1, -1):
                            append(children[idx])
                elif contains(clade.name):
                    tip_map[clade.name] = clade
        except AttributeError:
            tip_map = {}
            for tip in tree.get_terminals():
                if tip.name in ordered_name_set:
                    tip_map[tip.name] = tip

        return tip_map

    # ── Weight matrix construction ───────────────────────────────────

    def _build_weight_matrix_single_alpha(
        self, ordered_names: list[str], lineage_info: dict,
        regimes: list[str], alpha: float, root_regime: str,
    ) -> np.ndarray:
        """Build W matrix for OUM/OUMV: single alpha, multiple optima.

        W[i,r] = contribution of regime r's optimum to tip i's expected value.
        With root.station=TRUE, rows sum to 1.
        """
        n = len(ordered_names)
        R = len(regimes)
        regime_idx = {r: j for j, r in enumerate(regimes)}
        W = np.zeros((n, R))

        for i, name in enumerate(ordered_names):
            path = lineage_info[name]
            T_i = sum(bl for _, bl, _, _, _ in path)

            for _, bl, regime, d_start, d_end in path:
                if regime not in regime_idx:
                    continue
                j = regime_idx[regime]
                # Weight = (1 - exp(-alpha*bl)) * exp(-alpha*(T_i - d_end))
                if alpha < 1e-10:
                    w = bl / T_i if T_i > 0 else 0.0
                else:
                    w = (1.0 - np.exp(-alpha * bl)) * np.exp(-alpha * (T_i - d_end))
                W[i, j] += w

            # root.station=TRUE: add exp(-alpha*T_i) to root regime column
            if root_regime in regime_idx:
                j_root = regime_idx[root_regime]
                if alpha < 1e-10:
                    W[i, j_root] += 0.0  # BM limit: root contributes nothing extra
                else:
                    W[i, j_root] += np.exp(-alpha * T_i)

        return W

    def _prepare_single_alpha_weight_cache(
        self, ordered_names: list[str], lineage_info: dict,
        regimes: list[str], root_regime: str,
    ) -> dict[str, np.ndarray]:
        """Cache alpha-invariant arrays for single-alpha OUM weight matrices."""
        regime_idx = {r: j for j, r in enumerate(regimes)}
        n = len(ordered_names)
        rows = []
        cols = []
        branch_lengths = []
        tail_lengths = []
        bm_weights = []
        tip_heights = np.empty(n, dtype=float)

        for i, name in enumerate(ordered_names):
            path = lineage_info[name]
            tip_height = sum(bl for _, bl, _, _, _ in path)
            tip_heights[i] = tip_height
            inv_tip_height = 1.0 / tip_height if tip_height > 0 else 0.0

            for _, bl, regime, _d_start, d_end in path:
                col = regime_idx.get(regime)
                if col is None:
                    continue
                rows.append(i)
                cols.append(col)
                branch_lengths.append(bl)
                tail_lengths.append(tip_height - d_end)
                bm_weights.append(bl * inv_tip_height)

        return {
            "rows": np.asarray(rows, dtype=np.intp),
            "cols": np.asarray(cols, dtype=np.intp),
            "branch_lengths": np.asarray(branch_lengths, dtype=float),
            "tail_lengths": np.asarray(tail_lengths, dtype=float),
            "bm_weights": np.asarray(bm_weights, dtype=float),
            "tip_heights": tip_heights,
            "root_col": regime_idx.get(root_regime),
            "shape": (n, len(regimes)),
        }

    def _build_weight_matrix_single_alpha_from_cache(
        self, cache: dict[str, np.ndarray], alpha: float,
    ) -> np.ndarray:
        W = np.zeros(cache["shape"])
        if alpha < 1e-10:
            weights = cache["bm_weights"]
        else:
            weights = (
                1.0 - np.exp(-alpha * cache["branch_lengths"])
            ) * np.exp(-alpha * cache["tail_lengths"])

        np.add.at(W, (cache["rows"], cache["cols"]), weights)

        root_col = cache["root_col"]
        if root_col is not None and alpha >= 1e-10:
            W[:, root_col] += np.exp(-alpha * cache["tip_heights"])

        return W

    def _build_weight_matrix_multi_alpha(
        self, ordered_names: list[str], lineage_info: dict,
        regimes: list[str], alphas_dict: dict[str, float],
        root_regime: str,
    ) -> np.ndarray:
        """Build W matrix for OUMA/OUMVA: per-regime alpha values."""
        n = len(ordered_names)
        R = len(regimes)
        regime_idx = {r: j for j, r in enumerate(regimes)}
        root_col = regime_idx.get(root_regime)
        W = np.zeros((n, R))
        exp = math.exp

        for i, name in enumerate(ordered_names):
            path = lineage_info[name]

            # Compute cumulative decay from each branch end to tip
            # For multi-alpha, decay from branch b to tip uses alphas
            # of all downstream branches
            n_branches = len(path)
            # Precompute cumulative decay factor from branch end to tip
            # decay_to_tip[b] = product of exp(-alpha_r * bl_r) for all branches after b
            decay_to_tip = [1.0] * n_branches
            for b in range(n_branches - 2, -1, -1):
                _, bl_next, regime_next, _, _ = path[b + 1]
                alpha_next = alphas_dict.get(regime_next, 0.01)
                decay_to_tip[b] = decay_to_tip[b + 1] * exp(
                    -alpha_next * bl_next
                )

            for b_idx, (_, bl, regime, _d_start, _d_end) in enumerate(path):
                if regime not in regime_idx:
                    continue
                j = regime_idx[regime]
                alpha_r = alphas_dict.get(regime, 0.01)

                # Weight from this branch segment
                if alpha_r < 1e-10:
                    contribution = bl
                else:
                    contribution = 1.0 - exp(-alpha_r * bl)

                # Decay to tip from end of this branch
                if b_idx < n_branches - 1:
                    contribution *= decay_to_tip[b_idx]

                W[i, j] += contribution

            # root.station: add decay from root to tip for root regime
            if root_col is not None:
                total_decay = 1.0
                if n_branches:
                    _, bl_first, regime_first, _, _ = path[0]
                    alpha_first = alphas_dict.get(regime_first, 0.01)
                    total_decay = decay_to_tip[0] * exp(-alpha_first * bl_first)
                W[i, root_col] += total_decay

        return W

    # ── OU VCV construction ──────────────────────────────────────────

    def _build_ou_vcv_single_alpha(
        self, ordered_names: list[str], lineage_info: dict,
        alpha: float, sigma2: float,
    ) -> np.ndarray:
        """Build OU VCV matrix (conditional on root value).

        V_ij = sigma2 / (2*alpha) * exp(-alpha*(d_i + d_j)) * (1 - exp(-2*alpha*s_ij))

        where d_i = T_i - s_ij (unique path for tip i), s_ij = shared path.
        For a single global alpha and sigma2. Uses Martins & Hansen (1997).
        """
        cache = self._prepare_single_alpha_ou_cache(ordered_names, lineage_info)
        return self._build_ou_vcv_single_alpha_from_cache(cache, alpha, sigma2)

    def _prepare_single_alpha_ou_cache(
        self, ordered_names: list[str], lineage_info: dict,
    ) -> dict[str, np.ndarray]:
        """Cache shared path lengths used by single-alpha OU VCV builders."""
        n = len(ordered_names)
        shared_paths = np.zeros((n, n))
        paths = [lineage_info[name] for name in ordered_names]
        tip_heights = np.array(
            [sum(bl for _, bl, _, _, _ in path) for path in paths],
            dtype=float,
        )

        branch_to_indices = {}
        branch_lengths = {}
        for idx, path in enumerate(paths):
            for clade_id, bl, _regime, _d_start, _d_end in path:
                branch_to_indices.setdefault(clade_id, []).append(idx)
                branch_lengths[clade_id] = bl

        for clade_id, indices in branch_to_indices.items():
            branch_length = branch_lengths[clade_id]
            if len(indices) == 1:
                shared_paths[indices[0], indices[0]] += branch_length
                continue
            idx = np.asarray(indices, dtype=np.intp)
            shared_paths[np.ix_(idx, idx)] += branch_length

        unique_paths = (
            tip_heights[:, None] + tip_heights[None, :] - 2.0 * shared_paths
        )
        return {
            "tip_heights": tip_heights,
            "shared_paths": shared_paths,
            "unique_paths": unique_paths,
        }

    def _build_ou_vcv_single_alpha_from_cache(
        self, cache: dict[str, np.ndarray], alpha: float, sigma2: float,
    ) -> np.ndarray:
        shared_paths = cache["shared_paths"]
        if alpha < 1e-10:
            return sigma2 * shared_paths

        return (
            sigma2 / (2.0 * alpha)
        ) * np.exp(
            -alpha * cache["unique_paths"]
        ) * (
            1.0 - np.exp(-2.0 * alpha * shared_paths)
        )

    def _build_ou_H_matrices(
        self, ordered_names: list[str], lineage_info: dict,
        regimes: list[str], alpha: float,
    ) -> dict[str, np.ndarray]:
        """Build per-regime H_r(alpha) matrices for OUMV.

        V = sum_r sigma^2_r * H_r(alpha)

        H_r[i,j] captures the OU-transformed contribution of regime r's
        portion of the shared path to the covariance between tips i and j.
        """
        n = len(ordered_names)
        H = {r: np.zeros((n, n)) for r in regimes}
        tip_heights = np.array(
            [
                sum(bl for _, bl, _, _, _ in lineage_info[name])
                for name in ordered_names
            ],
            dtype=float,
        )
        branch_to_indices = {}
        branch_info = {}
        for idx, name in enumerate(ordered_names):
            for clade_id, bl, regime, _d_start, d_end in lineage_info[name]:
                branch_to_indices.setdefault(clade_id, []).append(idx)
                branch_info[clade_id] = (bl, regime, d_end)

        for clade_id, indices in branch_to_indices.items():
            bl, regime, d_end = branch_info[clade_id]
            if regime not in H:
                continue
            if len(indices) == 1:
                tip_idx = indices[0]
                if alpha < 1e-10:
                    H[regime][tip_idx, tip_idx] += bl
                else:
                    decay = np.exp(-alpha * (tip_heights[tip_idx] - d_end))
                    coeff = (1.0 - np.exp(-2.0 * alpha * bl)) / (2.0 * alpha)
                    H[regime][tip_idx, tip_idx] += coeff * decay * decay
                continue
            idx = np.asarray(indices, dtype=np.intp)
            if alpha < 1e-10:
                H[regime][idx[:, None], idx] += bl
            else:
                decay = np.exp(-alpha * (tip_heights[idx] - d_end))
                coeff = (1.0 - np.exp(-2.0 * alpha * bl)) / (2.0 * alpha)
                H[regime][idx[:, None], idx] += coeff * np.outer(decay, decay)

        return H

    def _build_ou_vcv_multi_alpha(
        self, ordered_names: list[str], lineage_info: dict,
        regimes: list[str], alphas_dict: dict[str, float],
        sigmas_dict: dict[str, float],
    ) -> np.ndarray:
        """Build OU VCV for OUMA/OUMVA with per-regime alpha and sigma^2.

        Uses per-branch computation with regime-specific decay factors.
        """
        n = len(ordered_names)
        V = np.zeros((n, n))
        branch_to_indices = {}
        branch_info = {}
        branch_decays = {}
        exp = math.exp

        for idx, name in enumerate(ordered_names):
            path = lineage_info[name]
            decay_to_tip = [1.0] * len(path)
            for b in range(len(path) - 2, -1, -1):
                _, bl_next, r_next, _, _ = path[b + 1]
                alpha_next = alphas_dict.get(r_next, 0.01)
                decay_to_tip[b] = decay_to_tip[b + 1] * exp(-alpha_next * bl_next)

            for branch_idx, (
                clade_id, bl, regime, _d_start, _d_end,
            ) in enumerate(path):
                branch_to_indices.setdefault(clade_id, []).append(idx)
                branch_info[clade_id] = (bl, regime)
                branch_decays.setdefault(clade_id, []).append(
                    decay_to_tip[branch_idx]
                )

        for clade_id, indices in branch_to_indices.items():
            bl, regime = branch_info[clade_id]
            alpha_r = alphas_dict.get(regime, 0.01)
            sigma2_r = sigmas_dict.get(regime, 0.01)
            if alpha_r < 1e-10:
                coeff = sigma2_r * bl
            else:
                coeff = sigma2_r / (2.0 * alpha_r) * (
                    1.0 - exp(-2.0 * alpha_r * bl)
                )
            if len(indices) == 1:
                decay = branch_decays[clade_id][0]
                V[indices[0], indices[0]] += coeff * decay * decay
                continue
            idx = np.asarray(indices, dtype=np.intp)
            decay = np.asarray(branch_decays[clade_id], dtype=float)
            V[np.ix_(idx, idx)] += coeff * np.outer(decay, decay)

        return V

    # ── GLS estimation ───────────────────────────────────────────────

    def _gls_theta_hat(
        self, x: np.ndarray, W: np.ndarray, V_inv: np.ndarray
    ) -> np.ndarray:
        """GLS estimate of theta: (W'V^{-1}W)^{-1} W'V^{-1}x"""
        WtVi = W.T @ V_inv
        try:
            theta = np.linalg.solve(WtVi @ W, WtVi @ x)
        except np.linalg.LinAlgError:
            theta = np.linalg.lstsq(WtVi @ W, WtVi @ x, rcond=None)[0]
        return theta

    def _ou_log_likelihood(
        self, x: np.ndarray, W: np.ndarray, V: np.ndarray,
        theta: np.ndarray,
    ) -> float:
        """Multivariate normal log-likelihood for OU model."""
        n = len(x)
        e = x - W @ theta

        try:
            sign, logdet = np.linalg.slogdet(V)
            if sign <= 0:
                return float("-inf")
            V_inv = np.linalg.inv(V)
        except np.linalg.LinAlgError:
            return float("-inf")

        quad = float(e @ V_inv @ e)
        ll = -0.5 * (n * np.log(2 * np.pi) + logdet + quad)
        return float(ll)

    def _concentrated_ll_bm(
        self, x: np.ndarray, C: np.ndarray
    ) -> tuple[float, float, float]:
        """Concentrated log-likelihood for BM (sigma^2 profiled out)."""
        try:
            return self._concentrated_ll_bm_cholesky(x, C)
        except (np.linalg.LinAlgError, FloatingPointError, ValueError):
            return self._concentrated_ll_bm_inverse(x, C)

    def _concentrated_ll_bm_cholesky(
        self, x: np.ndarray, C: np.ndarray
    ) -> tuple[float, float, float]:
        """Cholesky-backed concentrated BM likelihood."""
        n = len(x)
        ones = np.ones(n)
        factor = cho_factor(C, lower=True, check_finite=False)
        solve_rhs = np.empty((n, 2), dtype=np.result_type(x, C))
        solve_rhs[:, 0] = 1.0
        solve_rhs[:, 1] = x
        solved = cho_solve(factor, solve_rhs, check_finite=False)
        C_inv_ones = solved[:, 0]
        C_inv_x = solved[:, 1]

        denom = float(ones @ C_inv_ones)
        if abs(denom) < 1e-300:
            return float("-inf"), 0.0, 0.0
        z0 = float(ones @ C_inv_x) / denom

        e = x - z0
        C_inv_e = C_inv_x - C_inv_ones * z0
        sig2 = float(e @ C_inv_e) / n
        if sig2 <= 0:
            return float("-inf"), 0.0, z0

        logdet_C = 2.0 * float(np.log(np.diagonal(factor[0])).sum())
        ll = -0.5 * (n * np.log(2 * np.pi) + n * np.log(sig2) + logdet_C + n)
        return float(ll), float(sig2), float(z0)

    def _concentrated_ll_bm_inverse(
        self, x: np.ndarray, C: np.ndarray
    ) -> tuple[float, float, float]:
        """Inverse-based BM likelihood retained as a fallback."""
        n = len(x)
        ones = np.ones(n)

        try:
            C_inv = np.linalg.inv(C)
        except np.linalg.LinAlgError:
            return float("-inf"), 0.0, 0.0

        denom = float(ones @ C_inv @ ones)
        if abs(denom) < 1e-300:
            return float("-inf"), 0.0, 0.0
        z0 = float(ones @ C_inv @ x) / denom

        e = x - z0
        sig2 = float(e @ C_inv @ e) / n

        if sig2 <= 0:
            return float("-inf"), 0.0, z0

        sign, logdet_C = np.linalg.slogdet(C)
        if sign <= 0:
            return float("-inf"), sig2, z0

        ll = -0.5 * (n * np.log(2 * np.pi) + n * np.log(sig2) + logdet_C + n)
        return float(ll), float(sig2), float(z0)

    def _bm_log_likelihood_from_vcv(
        self, x: np.ndarray, V: np.ndarray, ones: np.ndarray
    ) -> tuple[float, float]:
        try:
            return self._bm_log_likelihood_from_vcv_cholesky(x, V, ones)
        except (np.linalg.LinAlgError, FloatingPointError, ValueError):
            return self._bm_log_likelihood_from_vcv_inverse(x, V, ones)

    def _bm_log_likelihood_from_vcv_cholesky(
        self, x: np.ndarray, V: np.ndarray, ones: np.ndarray
    ) -> tuple[float, float]:
        n = len(x)
        factor = cho_factor(V, lower=True, check_finite=False)
        solve_rhs = np.empty((n, 2), dtype=np.result_type(x, V))
        solve_rhs[:, 0] = ones
        solve_rhs[:, 1] = x
        solved = cho_solve(factor, solve_rhs, check_finite=False)
        V_inv_ones = solved[:, 0]
        V_inv_x = solved[:, 1]
        denom = float(ones @ V_inv_ones)
        if denom == 0:
            return 0.0, float("-inf")

        z0 = float(ones @ V_inv_x) / denom
        e = x - z0
        V_inv_e = V_inv_x - V_inv_ones * z0
        quad = float(e @ V_inv_e)
        logdet = 2.0 * float(np.log(np.diagonal(factor[0])).sum())
        ll = -0.5 * (n * np.log(2 * np.pi) + logdet + quad)
        return z0, float(ll)

    def _bm_log_likelihood_from_vcv_inverse(
        self, x: np.ndarray, V: np.ndarray, ones: np.ndarray
    ) -> tuple[float, float]:
        n = len(x)
        try:
            sign, logdet = np.linalg.slogdet(V)
            if sign <= 0:
                return 0.0, float("-inf")
            V_inv = np.linalg.inv(V)
        except np.linalg.LinAlgError:
            return 0.0, float("-inf")

        denom = float(ones @ V_inv @ ones)
        if denom == 0:
            return 0.0, float("-inf")
        z0 = float(ones @ V_inv @ x) / denom
        e = x - z0
        quad = float(e @ V_inv @ e)
        ll = -0.5 * (n * np.log(2 * np.pi) + logdet + quad)
        return z0, float(ll)

    def _ou_profile_likelihood(
        self, x: np.ndarray, W: np.ndarray, V: np.ndarray
    ) -> tuple[np.ndarray, float, float]:
        try:
            return self._ou_profile_likelihood_cholesky(x, W, V)
        except (np.linalg.LinAlgError, FloatingPointError, ValueError):
            return self._ou_profile_likelihood_inverse(x, W, V)

    def _ou_profile_likelihood_cholesky(
        self, x: np.ndarray, W: np.ndarray, V: np.ndarray
    ) -> tuple[np.ndarray, float, float]:
        n = len(x)
        factor = cho_factor(V, lower=True, check_finite=False)
        n_predictors = W.shape[1]
        solve_rhs = np.empty(
            (n, n_predictors + 1),
            dtype=np.result_type(W, x, V),
        )
        solve_rhs[:, :n_predictors] = W
        solve_rhs[:, n_predictors] = x
        solved = cho_solve(factor, solve_rhs, check_finite=False)
        V_inv_W = solved[:, :n_predictors]
        V_inv_x = solved[:, n_predictors]

        info = W.T @ V_inv_W
        rhs = W.T @ V_inv_x
        try:
            theta = np.linalg.solve(info, rhs)
        except np.linalg.LinAlgError:
            theta = np.linalg.lstsq(info, rhs, rcond=None)[0]

        e = x - W @ theta
        V_inv_e = V_inv_x - V_inv_W @ theta
        sig2 = float(e @ V_inv_e) / n
        if sig2 <= 0:
            return theta, sig2, float("-inf")

        logdet = 2.0 * float(np.log(np.diagonal(factor[0])).sum())
        ll = -0.5 * (n * np.log(2 * np.pi) + n * np.log(sig2) + logdet + n)
        return theta, float(sig2), float(ll)

    def _ou_profile_likelihood_inverse(
        self, x: np.ndarray, W: np.ndarray, V: np.ndarray
    ) -> tuple[np.ndarray, float, float]:
        n = len(x)
        try:
            sign, logdet = np.linalg.slogdet(V)
            if sign <= 0:
                return np.zeros(W.shape[1]), 0.0, float("-inf")
            V_inv = np.linalg.inv(V)
        except np.linalg.LinAlgError:
            return np.zeros(W.shape[1]), 0.0, float("-inf")

        theta = self._gls_theta_hat(x, W, V_inv)
        e = x - W @ theta
        sig2 = float(e @ V_inv @ e) / n
        if sig2 <= 0:
            return theta, sig2, float("-inf")

        ll = -0.5 * (n * np.log(2 * np.pi) + n * np.log(sig2) + logdet + n)
        return theta, float(sig2), float(ll)

    def _ou_gls_log_likelihood(
        self, x: np.ndarray, W: np.ndarray, V: np.ndarray
    ) -> tuple[np.ndarray, float]:
        try:
            return self._ou_gls_log_likelihood_cholesky(x, W, V)
        except (np.linalg.LinAlgError, FloatingPointError, ValueError):
            return self._ou_gls_log_likelihood_inverse(x, W, V)

    def _ou_gls_log_likelihood_cholesky(
        self, x: np.ndarray, W: np.ndarray, V: np.ndarray
    ) -> tuple[np.ndarray, float]:
        n = len(x)
        factor = cho_factor(V, lower=True, check_finite=False)
        n_predictors = W.shape[1]
        solve_rhs = np.empty(
            (n, n_predictors + 1),
            dtype=np.result_type(W, x, V),
        )
        solve_rhs[:, :n_predictors] = W
        solve_rhs[:, n_predictors] = x
        solved = cho_solve(factor, solve_rhs, check_finite=False)
        V_inv_W = solved[:, :n_predictors]
        V_inv_x = solved[:, n_predictors]

        info = W.T @ V_inv_W
        rhs = W.T @ V_inv_x
        try:
            theta = np.linalg.solve(info, rhs)
        except np.linalg.LinAlgError:
            theta = np.linalg.lstsq(info, rhs, rcond=None)[0]

        e = x - W @ theta
        V_inv_e = V_inv_x - V_inv_W @ theta
        quad = float(e @ V_inv_e)
        logdet = 2.0 * float(np.log(np.diagonal(factor[0])).sum())
        ll = -0.5 * (n * np.log(2 * np.pi) + logdet + quad)
        return theta, float(ll)

    def _ou_gls_log_likelihood_inverse(
        self, x: np.ndarray, W: np.ndarray, V: np.ndarray
    ) -> tuple[np.ndarray, float]:
        n = len(x)
        try:
            sign, logdet = np.linalg.slogdet(V)
            if sign <= 0:
                return np.zeros(W.shape[1]), float("-inf")
            V_inv = np.linalg.inv(V)
        except np.linalg.LinAlgError:
            return np.zeros(W.shape[1]), float("-inf")

        theta = self._gls_theta_hat(x, W, V_inv)
        e = x - W @ theta
        quad = float(e @ V_inv @ e)
        ll = -0.5 * (n * np.log(2 * np.pi) + logdet + quad)
        return theta, float(ll)

    # ── Model fitting ────────────────────────────────────────────────

    def _fit_model(
        self, name: str, x: np.ndarray, vcv_total: np.ndarray,
        per_regime_vcv: dict, ordered_names: list[str],
        regimes: list[str], lineage_info: dict,
        root_regime: str, tree_height: float, n_regimes: int,
    ) -> dict:
        if name == "BM1":
            return self._fit_bm1(x, vcv_total)
        elif name == "BMS":
            return self._fit_bms(x, per_regime_vcv, regimes)
        elif name == "OU1":
            return self._fit_ou1(x, ordered_names, lineage_info, tree_height)
        elif name == "OUM":
            return self._fit_oum(
                x, ordered_names, lineage_info, regimes,
                root_regime, tree_height,
            )
        elif name == "OUMV":
            return self._fit_oumv(
                x, ordered_names, lineage_info, regimes,
                root_regime, tree_height,
            )
        elif name == "OUMA":
            return self._fit_ouma(
                x, ordered_names, lineage_info, regimes,
                root_regime, tree_height,
            )
        elif name == "OUMVA":
            return self._fit_oumva(
                x, ordered_names, lineage_info, regimes,
                root_regime, tree_height,
            )
        else:
            raise PhykitUserError([f"Unknown model: {name}"], code=2)

    def _fit_bm1(self, x: np.ndarray, C: np.ndarray) -> dict:
        """BM1: single-rate Brownian motion. k=2 (z0, sigma2)."""
        ll, sig2, z0 = self._concentrated_ll_bm(x, C)
        return dict(
            model="BM1", log_likelihood=ll, k_params=2,
            params=dict(z0=z0, sigma2=sig2),
        )

    def _fit_bms(
        self, x: np.ndarray, per_regime_vcv: dict[str, np.ndarray],
        regimes: list[str],
    ) -> dict:
        """BMS: multi-rate BM. k = R+1 (z0, sigma2_1..R)."""
        n = len(x)
        k = len(regimes)
        regime_matrices = [per_regime_vcv[r] for r in regimes]
        ones = np.ones(n)

        def neg_log_likelihood(log_sigma2s):
            sigma2s = np.exp(log_sigma2s)
            V = np.zeros((n, n))
            for i_r in range(k):
                V += sigma2s[i_r] * regime_matrices[i_r]

            _, ll = self._bm_log_likelihood_from_vcv(x, V, ones)
            if not np.isfinite(ll):
                return 1e20
            return -ll

        starting_scales = [0.001, 0.01, 0.1, 1.0, 10.0]
        best_negll = np.inf
        best_params = np.zeros(k)

        for sv in starting_scales:
            x0 = np.log(np.ones(k) * sv)
            try:
                result = minimize(
                    neg_log_likelihood, x0, method="L-BFGS-B",
                    options={"maxiter": 5000},
                )
                if result.fun < best_negll:
                    best_negll = result.fun
                    best_params = result.x
            except (ValueError, np.linalg.LinAlgError):
                continue

        try:
            result = minimize(
                neg_log_likelihood, best_params, method="Nelder-Mead",
                options={"maxiter": 10000, "xatol": 1e-10, "fatol": 1e-10},
            )
            if result.fun < best_negll:
                best_negll = result.fun
                best_params = result.x
        except (ValueError, np.linalg.LinAlgError):
            pass

        sigma2s = np.exp(best_params)

        V = np.zeros((n, n))
        for i_r in range(k):
            V += sigma2s[i_r] * regime_matrices[i_r]

        z0, ll = self._bm_log_likelihood_from_vcv(x, V, ones)
        if not np.isfinite(ll):
            return dict(
                model="BMS", log_likelihood=float("-inf"),
                k_params=k + 1,
                params=dict(z0=0.0, sigma2={r: float(sigma2s[i]) for i, r in enumerate(regimes)}),
            )

        sigma2_dict = {r: float(sigma2s[i]) for i, r in enumerate(regimes)}
        return dict(
            model="BMS", log_likelihood=float(ll), k_params=k + 1,
            params=dict(z0=float(z0), sigma2=sigma2_dict),
        )

    def _fit_ou1(
        self, x: np.ndarray, ordered_names: list[str],
        lineage_info: dict, tree_height: float,
    ) -> dict:
        """OU1: single-regime OU. k=3 (theta, alpha, sigma2)."""
        n = len(x)
        upper = 100.0 / tree_height if tree_height > 0 else 100.0
        ou_cache = self._prepare_single_alpha_ou_cache(
            ordered_names, lineage_info
        )
        ones = np.ones(n)[:, None]

        def neg_ll(alpha):
            V = self._build_ou_vcv_single_alpha_from_cache(
                ou_cache, alpha, 1.0
            )
            _, _, ll = self._ou_profile_likelihood(x, ones, V)
            if not np.isfinite(ll):
                return 1e20
            return -ll

        alpha_hat, _ = self._optimize_parameter(neg_ll, (1e-8, upper))

        V = self._build_ou_vcv_single_alpha_from_cache(
            ou_cache, alpha_hat, 1.0
        )
        theta, sig2, ll = self._ou_profile_likelihood(x, ones, V)
        if not np.isfinite(ll):
            return dict(
                model="OU1", log_likelihood=float("-inf"), k_params=3,
                params=dict(theta=0.0, alpha=float(alpha_hat), sigma2=0.0),
            )

        return dict(
            model="OU1", log_likelihood=float(ll), k_params=3,
            params=dict(theta=float(theta[0]), alpha=float(alpha_hat), sigma2=float(sig2)),
        )

    def _fit_oum(
        self, x: np.ndarray, ordered_names: list[str],
        lineage_info: dict, regimes: list[str],
        root_regime: str, tree_height: float,
    ) -> dict:
        """OUM: multi-optimum OU. k = R+2 (theta_1..R, alpha, sigma2)."""
        R = len(regimes)
        upper = 100.0 / tree_height if tree_height > 0 else 100.0
        ou_cache = self._prepare_single_alpha_ou_cache(
            ordered_names, lineage_info
        )
        weight_cache = self._prepare_single_alpha_weight_cache(
            ordered_names, lineage_info, regimes, root_regime
        )

        def neg_ll(alpha):
            V = self._build_ou_vcv_single_alpha_from_cache(
                ou_cache, alpha, 1.0
            )
            W = self._build_weight_matrix_single_alpha_from_cache(
                weight_cache, alpha
            )

            _, _, ll = self._ou_profile_likelihood(x, W, V)
            if not np.isfinite(ll):
                return 1e20
            return -ll

        alpha_hat, _ = self._optimize_parameter(neg_ll, (1e-8, upper))

        V = self._build_ou_vcv_single_alpha_from_cache(
            ou_cache, alpha_hat, 1.0
        )
        W = self._build_weight_matrix_single_alpha_from_cache(
            weight_cache, alpha_hat
        )

        theta, sig2, ll = self._ou_profile_likelihood(x, W, V)
        if not np.isfinite(ll):
            return dict(
                model="OUM", log_likelihood=float("-inf"), k_params=R + 2,
                params=dict(
                    theta={r: 0.0 for r in regimes},
                    alpha=float(alpha_hat), sigma2=0.0,
                ),
            )

        theta_dict = {regimes[j]: float(theta[j]) for j in range(R)}
        return dict(
            model="OUM", log_likelihood=float(ll), k_params=R + 2,
            params=dict(
                theta=theta_dict,
                alpha=float(alpha_hat), sigma2=float(sig2),
            ),
        )

    def _fit_oumv(
        self, x: np.ndarray, ordered_names: list[str],
        lineage_info: dict, regimes: list[str],
        root_regime: str, tree_height: float,
    ) -> dict:
        """OUMV: multi-optimum + per-regime sigma^2. k = 2R+1."""
        n = len(x)
        R = len(regimes)
        weight_cache = self._prepare_single_alpha_weight_cache(
            ordered_names, lineage_info, regimes, root_regime
        )

        def neg_ll(params):
            log_alpha = params[0]
            log_sigma2s = params[1:]
            alpha = np.exp(log_alpha)
            sigma2s = np.exp(log_sigma2s)

            H = self._build_ou_H_matrices(
                ordered_names, lineage_info, regimes, alpha,
            )
            V = np.zeros((n, n))
            for r_idx, r in enumerate(regimes):
                V += sigma2s[r_idx] * H[r]

            W = self._build_weight_matrix_single_alpha_from_cache(
                weight_cache, alpha
            )

            _, ll = self._ou_gls_log_likelihood(x, W, V)
            if not np.isfinite(ll):
                return 1e20
            return -ll

        starting_scales = [0.001, 0.01, 0.1, 1.0, 10.0]
        best_negll = np.inf
        best_params = np.zeros(1 + R)

        for sv in starting_scales:
            x0 = np.array([np.log(0.1 / tree_height if tree_height > 0 else 0.1)]
                          + [np.log(sv)] * R)
            try:
                result = minimize(
                    neg_ll, x0, method="L-BFGS-B",
                    options={"maxiter": 5000},
                )
                if result.fun < best_negll:
                    best_negll = result.fun
                    best_params = result.x
            except (ValueError, np.linalg.LinAlgError):
                continue

        try:
            result = minimize(
                neg_ll, best_params, method="Nelder-Mead",
                options={"maxiter": 10000, "xatol": 1e-10, "fatol": 1e-10},
            )
            if result.fun < best_negll:
                best_negll = result.fun
                best_params = result.x
        except (ValueError, np.linalg.LinAlgError):
            pass

        alpha = np.exp(best_params[0])
        sigma2s = np.exp(best_params[1:])

        H = self._build_ou_H_matrices(
            ordered_names, lineage_info, regimes, alpha,
        )
        V = np.zeros((n, n))
        for r_idx, r in enumerate(regimes):
            V += sigma2s[r_idx] * H[r]

        W = self._build_weight_matrix_single_alpha_from_cache(
            weight_cache, alpha
        )

        theta, ll = self._ou_gls_log_likelihood(x, W, V)
        if not np.isfinite(ll):
            sigma2_dict = {r: float(sigma2s[i]) for i, r in enumerate(regimes)}
            return dict(
                model="OUMV", log_likelihood=float("-inf"), k_params=2 * R + 1,
                params=dict(
                    theta={r: 0.0 for r in regimes},
                    alpha=float(alpha), sigma2=sigma2_dict,
                ),
            )

        theta_dict = {regimes[j]: float(theta[j]) for j in range(R)}
        sigma2_dict = {r: float(sigma2s[i]) for i, r in enumerate(regimes)}
        return dict(
            model="OUMV", log_likelihood=float(ll), k_params=2 * R + 1,
            params=dict(
                theta=theta_dict,
                alpha=float(alpha), sigma2=sigma2_dict,
            ),
        )

    def _fit_ouma(
        self, x: np.ndarray, ordered_names: list[str],
        lineage_info: dict, regimes: list[str],
        root_regime: str, tree_height: float,
    ) -> dict:
        """OUMA: multi-optimum + per-regime alpha. k = 2R+1."""
        n = len(x)
        R = len(regimes)

        def neg_ll(log_alphas):
            alphas = np.exp(log_alphas)
            alphas_dict = {regimes[i]: float(alphas[i]) for i in range(R)}

            sigmas_dict = {r: 1.0 for r in regimes}
            V = self._build_ou_vcv_multi_alpha(
                ordered_names, lineage_info, regimes, alphas_dict, sigmas_dict,
            )

            W = self._build_weight_matrix_multi_alpha(
                ordered_names, lineage_info, regimes, alphas_dict, root_regime,
            )

            _, _, ll = self._ou_profile_likelihood(x, W, V)
            if not np.isfinite(ll):
                return 1e20
            return -ll

        starting_scales = [0.001, 0.01, 0.1, 1.0, 10.0]
        best_negll = np.inf
        best_params = np.zeros(R)

        for sv in starting_scales:
            x0 = np.log(np.ones(R) * sv / (tree_height if tree_height > 0 else 1.0))
            try:
                result = minimize(
                    neg_ll, x0, method="L-BFGS-B",
                    options={"maxiter": 5000},
                )
                if result.fun < best_negll:
                    best_negll = result.fun
                    best_params = result.x
            except (ValueError, np.linalg.LinAlgError):
                continue

        try:
            result = minimize(
                neg_ll, best_params, method="Nelder-Mead",
                options={"maxiter": 10000, "xatol": 1e-10, "fatol": 1e-10},
            )
            if result.fun < best_negll:
                best_negll = result.fun
                best_params = result.x
        except (ValueError, np.linalg.LinAlgError):
            pass

        alphas = np.exp(best_params)
        alphas_dict = {regimes[i]: float(alphas[i]) for i in range(R)}

        sigmas_dict = {r: 1.0 for r in regimes}
        V = self._build_ou_vcv_multi_alpha(
            ordered_names, lineage_info, regimes, alphas_dict, sigmas_dict,
        )
        W = self._build_weight_matrix_multi_alpha(
            ordered_names, lineage_info, regimes, alphas_dict, root_regime,
        )

        theta, sig2, ll = self._ou_profile_likelihood(x, W, V)
        if not np.isfinite(ll):
            alpha_dict = {r: float(alphas[i]) for i, r in enumerate(regimes)}
            return dict(
                model="OUMA", log_likelihood=float("-inf"), k_params=2 * R + 1,
                params=dict(
                    theta={r: 0.0 for r in regimes},
                    alpha=alpha_dict, sigma2=0.0,
                ),
            )

        theta_dict = {regimes[j]: float(theta[j]) for j in range(R)}
        alpha_dict = {r: float(alphas[i]) for i, r in enumerate(regimes)}
        return dict(
            model="OUMA", log_likelihood=float(ll), k_params=2 * R + 1,
            params=dict(
                theta=theta_dict,
                alpha=alpha_dict, sigma2=float(sig2),
            ),
        )

    def _fit_oumva(
        self, x: np.ndarray, ordered_names: list[str],
        lineage_info: dict, regimes: list[str],
        root_regime: str, tree_height: float,
    ) -> dict:
        """OUMVA: all regime-specific. k = 3R."""
        n = len(x)
        R = len(regimes)

        def neg_ll(params):
            log_alphas = params[:R]
            log_sigma2s = params[R:]
            alphas = np.exp(log_alphas)
            sigma2s = np.exp(log_sigma2s)

            alphas_dict = {regimes[i]: float(alphas[i]) for i in range(R)}
            sigmas_dict = {regimes[i]: float(sigma2s[i]) for i in range(R)}

            V = self._build_ou_vcv_multi_alpha(
                ordered_names, lineage_info, regimes, alphas_dict, sigmas_dict,
            )
            W = self._build_weight_matrix_multi_alpha(
                ordered_names, lineage_info, regimes, alphas_dict, root_regime,
            )

            _, ll = self._ou_gls_log_likelihood(x, W, V)
            if not np.isfinite(ll):
                return 1e20
            return -ll

        starting_scales = [0.001, 0.01, 0.1, 1.0, 10.0]
        best_negll = np.inf
        best_params = np.zeros(2 * R)

        for sv in starting_scales:
            x0 = np.array(
                [np.log(sv / (tree_height if tree_height > 0 else 1.0))] * R
                + [np.log(sv)] * R
            )
            try:
                result = minimize(
                    neg_ll, x0, method="L-BFGS-B",
                    options={"maxiter": 5000},
                )
                if result.fun < best_negll:
                    best_negll = result.fun
                    best_params = result.x
            except (ValueError, np.linalg.LinAlgError):
                continue

        try:
            result = minimize(
                neg_ll, best_params, method="Nelder-Mead",
                options={"maxiter": 10000, "xatol": 1e-10, "fatol": 1e-10},
            )
            if result.fun < best_negll:
                best_negll = result.fun
                best_params = result.x
        except (ValueError, np.linalg.LinAlgError):
            pass

        alphas = np.exp(best_params[:R])
        sigma2s = np.exp(best_params[R:])

        alphas_dict = {regimes[i]: float(alphas[i]) for i in range(R)}
        sigmas_dict = {regimes[i]: float(sigma2s[i]) for i in range(R)}

        V = self._build_ou_vcv_multi_alpha(
            ordered_names, lineage_info, regimes, alphas_dict, sigmas_dict,
        )
        W = self._build_weight_matrix_multi_alpha(
            ordered_names, lineage_info, regimes, alphas_dict, root_regime,
        )

        theta, ll = self._ou_gls_log_likelihood(x, W, V)
        if not np.isfinite(ll):
            alpha_dict = {r: float(alphas[i]) for i, r in enumerate(regimes)}
            sigma2_dict = {r: float(sigma2s[i]) for i, r in enumerate(regimes)}
            return dict(
                model="OUMVA", log_likelihood=float("-inf"), k_params=3 * R,
                params=dict(
                    theta={r: 0.0 for r in regimes},
                    alpha=alpha_dict, sigma2=sigma2_dict,
                ),
            )

        theta_dict = {regimes[j]: float(theta[j]) for j in range(R)}
        alpha_dict = {r: float(alphas[i]) for i, r in enumerate(regimes)}
        sigma2_dict = {r: float(sigma2s[i]) for i, r in enumerate(regimes)}
        return dict(
            model="OUMVA", log_likelihood=float(ll), k_params=3 * R,
            params=dict(
                theta=theta_dict,
                alpha=alpha_dict, sigma2=sigma2_dict,
            ),
        )

    # ── Optimization helper ──────────────────────────────────────────

    def _optimize_parameter(self, neg_ll_func, bounds, niter=10):
        lo_bound, hi_bound = bounds
        bounds_lo = np.linspace(
            lo_bound, hi_bound - (hi_bound - lo_bound) / niter, niter
        )
        bounds_hi = np.linspace(
            lo_bound + (hi_bound - lo_bound) / niter, hi_bound, niter
        )

        best_ll = -np.inf
        best_param = (lo_bound + hi_bound) / 2.0

        for lo, hi in zip(bounds_lo, bounds_hi):
            if lo >= hi:
                continue
            try:
                res = minimize_scalar(
                    neg_ll_func, bounds=(lo, hi), method="bounded"
                )
                ll_val = -res.fun
                if ll_val > best_ll:
                    best_ll = ll_val
                    best_param = res.x
            except (ValueError, RuntimeError):
                continue

        return best_param, best_ll

    # ── Model comparison ─────────────────────────────────────────────

    def _compute_model_comparison(
        self, results: list[dict], n: int,
        x: np.ndarray = None, vcv_total: np.ndarray = None,
        regime_assignments: dict = None, ordered_names: list[str] = None,
    ) -> list[dict]:
        for r in results:
            k = r["k_params"]
            ll = r["log_likelihood"]
            r["aic"] = -2.0 * ll + 2.0 * k
            if n - k - 1 > 0:
                r["aicc"] = r["aic"] + 2.0 * k * (k + 1) / (n - k - 1)
            else:
                r["aicc"] = float("inf")
            r["bic"] = -2.0 * ll + k * np.log(n)

        # Sort by AICc
        results.sort(key=lambda r: r["aicc"])

        min_aic = min(r["aic"] for r in results)
        min_aicc = results[0]["aicc"]
        min_bic = min(r["bic"] for r in results)

        for r in results:
            r["delta_aic"] = r["aic"] - min_aic
            r["delta_aicc"] = r["aicc"] - min_aicc
            r["delta_bic"] = r["bic"] - min_bic

        # AICc weights
        raw_weights = [math.exp(-0.5 * r["delta_aicc"]) for r in results]
        total = sum(raw_weights)
        if total > 0.0:
            for r, w in zip(results, raw_weights):
                r["aicc_weight"] = w / total
        else:
            for r in results:
                r["aicc_weight"] = 0.0

        # Compute R² = 1 - (σ²_model / σ²_BM1)
        bm1_sig2 = None
        for r in results:
            if r["model"] == "BM1":
                bm1_sig2 = r["params"]["sigma2"]
                break

        # If BM1 not in selected models, fit silently as baseline
        if bm1_sig2 is None and x is not None and vcv_total is not None:
            bm1_result = self._fit_bm1(x, vcv_total)
            bm1_sig2 = bm1_result["params"]["sigma2"]

        # Build tip-count-per-regime map for weighted averaging
        regime_tip_counts = {}
        if regime_assignments and ordered_names:
            for name in ordered_names:
                r_name = regime_assignments.get(name)
                if r_name is not None:
                    regime_tip_counts[r_name] = regime_tip_counts.get(r_name, 0) + 1

        for r in results:
            params = r["params"]
            sig2 = params.get("sigma2")
            if isinstance(sig2, dict):
                # Multi-regime: weighted average by tips per regime
                if regime_tip_counts:
                    sig2_model = sum(
                        (regime_tip_counts.get(rname, 0) / n) * sv
                        for rname, sv in sig2.items()
                    )
                else:
                    sig2_model = sum(sig2.values()) / len(sig2) if sig2 else 0
            elif sig2 is not None:
                sig2_model = sig2
            else:
                sig2_model = 0

            if bm1_sig2 is not None and bm1_sig2 > 0:
                r["r_squared"] = 1.0 - sig2_model / bm1_sig2
            else:
                r["r_squared"] = float("nan")

        return results

    # ── Output ───────────────────────────────────────────────────────

    def _print_text_output(
        self, results: list[dict], n: int, regimes: list[str]
    ) -> None:
        lines = [
            "OUwie Model Comparison (Multi-Regime OU)",
            f"\nNumber of tips: {n}",
            f"Regimes: {len(regimes)} ({', '.join(regimes)})\n",
        ]

        header = (
            f"{'Model':<8}{'k':<5}{'LL':<12}"
            f"{'AIC':<10}{'AICc':<10}{'dAICc':<9}{'AICcW':<9}"
            f"{'BIC':<10}{'dBIC':<9}{'R2':<7}"
        )
        lines.append(header)

        for r in results:
            lines.append(
                f"{r['model']:<8}{r['k_params']:<5}{r['log_likelihood']:<12.3f}"
                f"{r['aic']:<10.2f}{r['aicc']:<10.2f}"
                f"{r['delta_aicc']:<9.2f}{r['aicc_weight']:<9.3f}"
                f"{r['bic']:<10.2f}{r['delta_bic']:<9.2f}"
                f"{r['r_squared']:<7.3f}"
            )

        lines.append(f"\nBest model (AICc): {results[0]['model']}")
        best_bic = min(results, key=lambda r: r["bic"])["model"]
        lines.append(f"Best model (BIC):  {best_bic}")

        # Print parameter details for best model
        best = results[0]
        lines.append(f"\nParameter estimates ({best['model']}):")
        params = best["params"]

        if isinstance(params.get("theta"), dict):
            lines.append("  Theta (trait optima):")
            lines.extend(f"    {r_name}: {val:.4f}" for r_name, val in params["theta"].items())
        elif "theta" in params:
            lines.append(f"  Theta: {params['theta']:.4f}")
        if "z0" in params:
            lines.append(f"  z0: {params['z0']:.4f}")

        if isinstance(params.get("alpha"), dict):
            lines.append("  Alpha (selection strength):")
            lines.extend(f"    {r_name}: {val:.4f}" for r_name, val in params["alpha"].items())
        elif "alpha" in params:
            lines.append(f"  Alpha: {params['alpha']:.4f}")

        if isinstance(params.get("sigma2"), dict):
            lines.append("  Sigma-squared (diffusion rate):")
            lines.extend(f"    {r_name}: {val:.4f}" for r_name, val in params["sigma2"].items())
        elif "sigma2" in params:
            lines.append(f"  Sigma-squared: {params['sigma2']:.4f}")
        print("\n".join(lines))

    def _print_json_output(
        self, results: list[dict], n: int, regimes: list[str]
    ) -> None:
        models = {}
        for r in results:
            models[r["model"]] = dict(
                log_likelihood=r["log_likelihood"],
                k_params=r["k_params"],
                aic=r["aic"],
                aicc=r["aicc"],
                delta_aic=r["delta_aic"],
                delta_aicc=r["delta_aicc"],
                aicc_weight=r["aicc_weight"],
                bic=r["bic"],
                delta_bic=r["delta_bic"],
                r_squared=r["r_squared"],
                params=r["params"],
            )

        best_aicc = results[0]["model"]
        best_bic = min(results, key=lambda r: r["bic"])["model"]

        payload = dict(
            n_tips=n,
            regimes=regimes,
            models=models,
            best_model_aicc=best_aicc,
            best_model_bic=best_bic,
        )
        print_json(payload, sort_keys=False)
