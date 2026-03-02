import copy
import sys
from typing import Dict, List, Tuple

import numpy as np
from scipy.optimize import minimize, minimize_scalar

from .base import Tree
from ...helpers.json_output import print_json
from ...errors import PhykitUserError

ALL_MODELS = ["BM1", "BMS", "OU1", "OUM", "OUMV", "OUMA", "OUMVA"]


class OUwie(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"])
        self.trait_data_path = parsed["trait_data_path"]
        self.regime_data_path = parsed["regime_data_path"]
        self.selected_models = parsed["models"]
        self.json_output = parsed["json_output"]

    def run(self) -> None:
        tree = self.read_tree_file()
        self._validate_tree(tree)

        tree_tips = self.get_tip_names_from_tree(tree)
        traits = self._parse_trait_file(self.trait_data_path, tree_tips)
        regime_assignments = self._parse_regime_file(
            self.regime_data_path, tree_tips
        )

        shared = set(traits.keys()) & set(regime_assignments.keys())
        if len(shared) < 3:
            raise PhykitUserError(
                [
                    f"Only {len(shared)} shared taxa among tree, trait, and regime files.",
                    "At least 3 shared taxa are required.",
                ],
                code=2,
            )

        traits = {k: traits[k] for k in shared}
        regime_assignments = {k: regime_assignments[k] for k in shared}

        tree_copy = copy.deepcopy(tree)
        tip_names_in_tree = [t.name for t in tree_copy.get_terminals()]
        tips_to_prune = [t for t in tip_names_in_tree if t not in shared]
        if tips_to_prune:
            tree_copy = self.prune_tree_using_taxa_list(tree_copy, tips_to_prune)

        ordered_names = sorted(traits.keys())
        n = len(ordered_names)
        x = np.array([traits[name] for name in ordered_names])

        regimes = sorted(set(regime_assignments.values()))
        n_regimes = len(regimes)

        if n_regimes < 2:
            raise PhykitUserError(
                ["At least 2 distinct regimes are required for OUwie models."],
                code=2,
            )

        parent_map = self._build_parent_map(tree_copy)
        branch_regimes = self._assign_branch_regimes(
            tree_copy, regime_assignments, parent_map
        )

        # Determine root regime
        root_regime = branch_regimes.get(
            id(tree_copy.root), regimes[0]
        )
        # Get root regime from Fitch parsimony node regimes
        root_regime = self._get_root_regime(tree_copy, regime_assignments, parent_map)

        lineage_info = self._build_lineage_info(
            tree_copy, ordered_names, parent_map, branch_regimes
        )

        per_regime_vcv = self._build_per_regime_vcv(
            tree_copy, ordered_names, regimes, branch_regimes, parent_map
        )

        vcv_total = sum(per_regime_vcv.values())
        tree_height = float(np.max(np.diag(vcv_total)))

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

    def process_args(self, args) -> Dict:
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

    def _validate_tree(self, tree) -> None:
        tips = list(tree.get_terminals())
        if len(tips) < 3:
            raise PhykitUserError(
                ["Tree must have at least 3 tips for OUwie model fitting."],
                code=2,
            )
        for clade in tree.find_clades():
            if clade.branch_length is None and clade != tree.root:
                raise PhykitUserError(
                    ["All branches in the tree must have lengths."],
                    code=2,
                )

    def _parse_trait_file(
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

    def _parse_regime_file(
        self, path: str, tree_tips: List[str]
    ) -> Dict[str, str]:
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

        regimes = {}
        for line_num, line in enumerate(lines, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) != 2:
                raise PhykitUserError(
                    [
                        f"Line {line_num} in regime file has {len(parts)} columns; expected 2.",
                        "Each line should be: taxon_name<tab>regime_label",
                    ],
                    code=2,
                )
            taxon, regime = parts
            regimes[taxon] = regime

        tree_tip_set = set(tree_tips)
        regime_taxa_set = set(regimes.keys())
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

    def _build_parent_map(self, tree) -> Dict:
        parent_map = {}
        for clade in tree.find_clades(order="preorder"):
            for child in clade.clades:
                parent_map[id(child)] = clade
        return parent_map

    def _assign_branch_regimes(
        self, tree, tip_regimes: Dict[str, str], parent_map: Dict
    ) -> Dict:
        state_sets = {}
        for clade in tree.find_clades(order="postorder"):
            if clade.is_terminal():
                if clade.name in tip_regimes:
                    state_sets[id(clade)] = {tip_regimes[clade.name]}
                else:
                    state_sets[id(clade)] = set()
            else:
                child_sets = [
                    state_sets[id(c)] for c in clade.clades
                    if id(c) in state_sets and state_sets[id(c)]
                ]
                if not child_sets:
                    state_sets[id(clade)] = set()
                else:
                    intersection = child_sets[0]
                    for cs in child_sets[1:]:
                        intersection = intersection & cs
                    if intersection:
                        state_sets[id(clade)] = intersection
                    else:
                        union = set()
                        for cs in child_sets:
                            union = union | cs
                        state_sets[id(clade)] = union

        node_regimes = {}
        root = tree.root
        root_set = state_sets.get(id(root), set())
        if root_set:
            node_regimes[id(root)] = sorted(root_set)[0]
        else:
            node_regimes[id(root)] = sorted(set(tip_regimes.values()))[0]

        for clade in tree.find_clades(order="preorder"):
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
                    node_regimes[id(clade)] = sorted(clade_set)[0]
            else:
                node_regimes[id(clade)] = sorted(clade_set)[0]

        branch_regimes = {}
        for clade in tree.find_clades(order="preorder"):
            if clade == root:
                continue
            branch_regimes[id(clade)] = node_regimes.get(id(clade), "unknown")

        return branch_regimes

    def _get_root_regime(self, tree, tip_regimes, parent_map):
        state_sets = {}
        for clade in tree.find_clades(order="postorder"):
            if clade.is_terminal():
                if clade.name in tip_regimes:
                    state_sets[id(clade)] = {tip_regimes[clade.name]}
                else:
                    state_sets[id(clade)] = set()
            else:
                child_sets = [
                    state_sets[id(c)] for c in clade.clades
                    if id(c) in state_sets and state_sets[id(c)]
                ]
                if not child_sets:
                    state_sets[id(clade)] = set()
                else:
                    intersection = child_sets[0]
                    for cs in child_sets[1:]:
                        intersection = intersection & cs
                    if intersection:
                        state_sets[id(clade)] = intersection
                    else:
                        union = set()
                        for cs in child_sets:
                            union = union | cs
                        state_sets[id(clade)] = union

        root = tree.root
        root_set = state_sets.get(id(root), set())
        if root_set:
            return sorted(root_set)[0]
        return sorted(set(tip_regimes.values()))[0]

    def _build_root_to_tip_paths(
        self, tree, ordered_names: List[str], parent_map: Dict
    ) -> Dict[str, List]:
        tip_map = {}
        for tip in tree.get_terminals():
            if tip.name in ordered_names:
                tip_map[tip.name] = tip

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
        self, tree, ordered_names: List[str], regimes: List[str],
        branch_regimes: Dict, parent_map: Dict,
    ) -> Dict[str, np.ndarray]:
        n = len(ordered_names)
        paths = self._build_root_to_tip_paths(tree, ordered_names, parent_map)

        per_regime_vcv = {r: np.zeros((n, n)) for r in regimes}

        for i in range(n):
            path_i = paths[ordered_names[i]]
            for clade_id, bl in path_i:
                r = branch_regimes.get(clade_id, regimes[0])
                per_regime_vcv[r][i, i] += bl

            for j in range(i + 1, n):
                path_j = paths[ordered_names[j]]
                min_len = min(len(path_i), len(path_j))
                for s in range(min_len):
                    if path_i[s][0] == path_j[s][0]:
                        clade_id = path_i[s][0]
                        bl = path_i[s][1]
                        r = branch_regimes.get(clade_id, regimes[0])
                        per_regime_vcv[r][i, j] += bl
                        per_regime_vcv[r][j, i] += bl
                    else:
                        break

        return per_regime_vcv

    # ── Lineage info for OU weight/VCV matrices ──────────────────────

    def _build_lineage_info(
        self, tree, ordered_names: List[str], parent_map: Dict,
        branch_regimes: Dict,
    ) -> Dict[str, List[Tuple]]:
        """Build root-to-tip lineage info for each tip.

        Returns dict mapping tip_name -> list of tuples:
            (clade_id, branch_length, regime, dist_from_root_start, dist_from_root_end)
        """
        tip_map = {}
        for tip in tree.get_terminals():
            if tip.name in ordered_names:
                tip_map[tip.name] = tip

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

    # ── Weight matrix construction ───────────────────────────────────

    def _build_weight_matrix_single_alpha(
        self, ordered_names: List[str], lineage_info: Dict,
        regimes: List[str], alpha: float, root_regime: str,
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

    def _build_weight_matrix_multi_alpha(
        self, ordered_names: List[str], lineage_info: Dict,
        regimes: List[str], alphas_dict: Dict[str, float],
        root_regime: str,
    ) -> np.ndarray:
        """Build W matrix for OUMA/OUMVA: per-regime alpha values."""
        n = len(ordered_names)
        R = len(regimes)
        regime_idx = {r: j for j, r in enumerate(regimes)}
        W = np.zeros((n, R))

        for i, name in enumerate(ordered_names):
            path = lineage_info[name]

            # Compute cumulative decay from each branch end to tip
            # For multi-alpha, decay from branch b to tip uses alphas
            # of all downstream branches
            n_branches = len(path)
            # Precompute cumulative decay factor from branch end to tip
            # decay_to_tip[b] = product of exp(-alpha_r * bl_r) for all branches after b
            decay_to_tip = np.ones(n_branches)
            for b in range(n_branches - 2, -1, -1):
                _, bl_next, regime_next, _, _ = path[b + 1]
                alpha_next = alphas_dict.get(regime_next, 0.01)
                decay_to_tip[b] = decay_to_tip[b + 1] * np.exp(-alpha_next * bl_next)

            for b_idx, (_, bl, regime, d_start, d_end) in enumerate(path):
                if regime not in regime_idx:
                    continue
                j = regime_idx[regime]
                alpha_r = alphas_dict.get(regime, 0.01)

                # Weight from this branch segment
                if alpha_r < 1e-10:
                    contribution = bl
                else:
                    contribution = (1.0 - np.exp(-alpha_r * bl))

                # Decay to tip from end of this branch
                if b_idx < n_branches - 1:
                    contribution *= decay_to_tip[b_idx]

                W[i, j] += contribution

            # root.station: add decay from root to tip for root regime
            if root_regime in regime_idx:
                j_root = regime_idx[root_regime]
                # Total decay from root to tip
                total_decay = 1.0
                for _, bl, regime, _, _ in path:
                    alpha_r = alphas_dict.get(regime, 0.01)
                    total_decay *= np.exp(-alpha_r * bl)
                W[i, j_root] += total_decay

        return W

    # ── OU VCV construction ──────────────────────────────────────────

    def _build_ou_vcv_single_alpha(
        self, ordered_names: List[str], lineage_info: Dict,
        alpha: float, sigma2: float,
    ) -> np.ndarray:
        """Build OU VCV matrix (conditional on root value).

        V_ij = sigma2 / (2*alpha) * exp(-alpha*(d_i + d_j)) * (1 - exp(-2*alpha*s_ij))

        where d_i = T_i - s_ij (unique path for tip i), s_ij = shared path.
        For a single global alpha and sigma2. Uses Martins & Hansen (1997).
        """
        n = len(ordered_names)
        V = np.zeros((n, n))

        tip_heights = {}
        for name in ordered_names:
            tip_heights[name] = sum(bl for _, bl, _, _, _ in lineage_info[name])

        for i in range(n):
            for j in range(i, n):
                T_i = tip_heights[ordered_names[i]]
                T_j = tip_heights[ordered_names[j]]

                if i == j:
                    s_ij = T_i
                else:
                    path_i = lineage_info[ordered_names[i]]
                    path_j = lineage_info[ordered_names[j]]
                    s_ij = 0.0
                    min_len = min(len(path_i), len(path_j))
                    for s in range(min_len):
                        if path_i[s][0] == path_j[s][0]:
                            s_ij += path_i[s][1]
                        else:
                            break

                if alpha < 1e-10:
                    val = sigma2 * s_ij
                else:
                    d_i = T_i - s_ij
                    d_j = T_j - s_ij
                    val = (sigma2 / (2.0 * alpha)) * np.exp(
                        -alpha * (d_i + d_j)
                    ) * (1.0 - np.exp(-2.0 * alpha * s_ij))

                V[i, j] = val
                V[j, i] = val

        return V

    def _build_ou_H_matrices(
        self, ordered_names: List[str], lineage_info: Dict,
        regimes: List[str], alpha: float,
    ) -> Dict[str, np.ndarray]:
        """Build per-regime H_r(alpha) matrices for OUMV.

        V = sum_r sigma^2_r * H_r(alpha)

        H_r[i,j] captures the OU-transformed contribution of regime r's
        portion of the shared path to the covariance between tips i and j.
        """
        n = len(ordered_names)
        H = {r: np.zeros((n, n)) for r in regimes}

        tip_heights = {}
        for name in ordered_names:
            tip_heights[name] = sum(bl for _, bl, _, _, _ in lineage_info[name])

        for i in range(n):
            path_i = lineage_info[ordered_names[i]]
            T_i = tip_heights[ordered_names[i]]

            # Diagonal: contributions from each branch along this lineage
            for _, bl, regime, d_start, d_end in path_i:
                if regime not in H:
                    continue
                if alpha < 1e-10:
                    H[regime][i, i] += bl
                else:
                    # OU contribution of this branch segment to variance
                    # = exp(-2*alpha*(T_i - d_end)) * (1 - exp(-2*alpha*bl)) / (2*alpha)
                    contrib = np.exp(-2.0 * alpha * (T_i - d_end)) * (
                        1.0 - np.exp(-2.0 * alpha * bl)
                    ) / (2.0 * alpha)
                    H[regime][i, i] += contrib

            for j in range(i + 1, n):
                path_j = lineage_info[ordered_names[j]]
                T_j = tip_heights[ordered_names[j]]

                # Shared branches
                min_len = min(len(path_i), len(path_j))
                for s in range(min_len):
                    if path_i[s][0] == path_j[s][0]:
                        _, bl, regime, d_start, d_end = path_i[s]
                        if regime not in H:
                            continue
                        if alpha < 1e-10:
                            H[regime][i, j] += bl
                            H[regime][j, i] += bl
                        else:
                            # OU covariance contribution from shared branch in regime r
                            contrib = np.exp(-alpha * (T_i + T_j - 2.0 * d_end)) * (
                                1.0 - np.exp(-2.0 * alpha * bl)
                            ) / (2.0 * alpha)
                            H[regime][i, j] += contrib
                            H[regime][j, i] += contrib
                    else:
                        break

        return H

    def _build_ou_vcv_multi_alpha(
        self, ordered_names: List[str], lineage_info: Dict,
        regimes: List[str], alphas_dict: Dict[str, float],
        sigmas_dict: Dict[str, float],
    ) -> np.ndarray:
        """Build OU VCV for OUMA/OUMVA with per-regime alpha and sigma^2.

        Uses per-branch computation with regime-specific decay factors.
        """
        n = len(ordered_names)
        V = np.zeros((n, n))

        for i in range(n):
            path_i = lineage_info[ordered_names[i]]

            # Precompute cumulative decay factors from each branch end to tip i
            n_i = len(path_i)
            decay_i = np.ones(n_i)  # decay from end of branch b to tip
            for b in range(n_i - 2, -1, -1):
                _, bl_next, r_next, _, _ = path_i[b + 1]
                alpha_next = alphas_dict.get(r_next, 0.01)
                decay_i[b] = decay_i[b + 1] * np.exp(-alpha_next * bl_next)

            # Diagonal: variance for tip i
            for b_idx, (_, bl, regime, d_start, d_end) in enumerate(path_i):
                alpha_r = alphas_dict.get(regime, 0.01)
                sigma2_r = sigmas_dict.get(regime, 0.01)
                d_b = decay_i[b_idx]  # decay from end of branch b to tip

                if alpha_r < 1e-10:
                    contrib = sigma2_r * bl * d_b * d_b
                else:
                    contrib = sigma2_r / (2.0 * alpha_r) * (
                        1.0 - np.exp(-2.0 * alpha_r * bl)
                    ) * d_b * d_b
                V[i, i] += contrib

            for j in range(i + 1, n):
                path_j = lineage_info[ordered_names[j]]

                # Precompute decay factors for tip j
                n_j = len(path_j)
                decay_j = np.ones(n_j)
                for b in range(n_j - 2, -1, -1):
                    _, bl_next, r_next, _, _ = path_j[b + 1]
                    alpha_next = alphas_dict.get(r_next, 0.01)
                    decay_j[b] = decay_j[b + 1] * np.exp(-alpha_next * bl_next)

                # Covariance: only shared branches contribute
                min_len = min(n_i, n_j)
                cov = 0.0
                for s in range(min_len):
                    if path_i[s][0] == path_j[s][0]:
                        _, bl, regime, d_start, d_end = path_i[s]
                        alpha_r = alphas_dict.get(regime, 0.01)
                        sigma2_r = sigmas_dict.get(regime, 0.01)
                        d_i_b = decay_i[s]
                        d_j_b = decay_j[s]

                        if alpha_r < 1e-10:
                            contrib = sigma2_r * bl * d_i_b * d_j_b
                        else:
                            contrib = sigma2_r / (2.0 * alpha_r) * (
                                1.0 - np.exp(-2.0 * alpha_r * bl)
                            ) * d_i_b * d_j_b
                        cov += contrib
                    else:
                        break

                V[i, j] = cov
                V[j, i] = cov

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
    ) -> Tuple[float, float, float]:
        """Concentrated log-likelihood for BM (sigma^2 profiled out)."""
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

    # ── Model fitting ────────────────────────────────────────────────

    def _fit_model(
        self, name: str, x: np.ndarray, vcv_total: np.ndarray,
        per_regime_vcv: Dict, ordered_names: List[str],
        regimes: List[str], lineage_info: Dict,
        root_regime: str, tree_height: float, n_regimes: int,
    ) -> Dict:
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

    def _fit_bm1(self, x: np.ndarray, C: np.ndarray) -> Dict:
        """BM1: single-rate Brownian motion. k=2 (z0, sigma2)."""
        ll, sig2, z0 = self._concentrated_ll_bm(x, C)
        return dict(
            model="BM1", log_likelihood=ll, k_params=2,
            params=dict(z0=z0, sigma2=sig2),
        )

    def _fit_bms(
        self, x: np.ndarray, per_regime_vcv: Dict[str, np.ndarray],
        regimes: List[str],
    ) -> Dict:
        """BMS: multi-rate BM. k = R+1 (z0, sigma2_1..R)."""
        n = len(x)
        k = len(regimes)
        regime_matrices = [per_regime_vcv[r] for r in regimes]

        def neg_log_likelihood(log_sigma2s):
            sigma2s = np.exp(log_sigma2s)
            V = np.zeros((n, n))
            for i_r in range(k):
                V += sigma2s[i_r] * regime_matrices[i_r]

            try:
                sign, logdet = np.linalg.slogdet(V)
                if sign <= 0:
                    return 1e20
                V_inv = np.linalg.inv(V)
            except np.linalg.LinAlgError:
                return 1e20

            ones = np.ones(n)
            denom = float(ones @ V_inv @ ones)
            if denom == 0:
                return 1e20
            z0 = float(ones @ V_inv @ x) / denom
            e = x - z0
            quad = float(e @ V_inv @ e)

            ll = -0.5 * (n * np.log(2 * np.pi) + logdet + quad)
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

        try:
            V_inv = np.linalg.inv(V)
        except np.linalg.LinAlgError:
            return dict(
                model="BMS", log_likelihood=float("-inf"),
                k_params=k + 1,
                params=dict(z0=0.0, sigma2={r: float(sigma2s[i]) for i, r in enumerate(regimes)}),
            )

        ones = np.ones(n)
        z0 = float(ones @ V_inv @ x) / float(ones @ V_inv @ ones)
        e = x - z0
        sign, logdet = np.linalg.slogdet(V)
        quad = float(e @ V_inv @ e)
        ll = -0.5 * (n * np.log(2 * np.pi) + logdet + quad)

        sigma2_dict = {r: float(sigma2s[i]) for i, r in enumerate(regimes)}
        return dict(
            model="BMS", log_likelihood=float(ll), k_params=k + 1,
            params=dict(z0=float(z0), sigma2=sigma2_dict),
        )

    def _fit_ou1(
        self, x: np.ndarray, ordered_names: List[str],
        lineage_info: Dict, tree_height: float,
    ) -> Dict:
        """OU1: single-regime OU. k=3 (theta, alpha, sigma2)."""
        n = len(x)
        upper = 100.0 / tree_height if tree_height > 0 else 100.0

        def neg_ll(alpha):
            V = self._build_ou_vcv_single_alpha(
                ordered_names, lineage_info, alpha, 1.0
            )
            # Profile sigma2 and theta (=z0 for single regime)
            try:
                V_inv = np.linalg.inv(V)
            except np.linalg.LinAlgError:
                return 1e20

            ones = np.ones(n)
            denom = float(ones @ V_inv @ ones)
            if abs(denom) < 1e-300:
                return 1e20
            theta = float(ones @ V_inv @ x) / denom
            e = x - theta
            sig2 = float(e @ V_inv @ e) / n
            if sig2 <= 0:
                return 1e20

            sign, logdet = np.linalg.slogdet(V)
            if sign <= 0:
                return 1e20

            ll = -0.5 * (n * np.log(2 * np.pi) + n * np.log(sig2) + logdet + n)
            return -ll

        alpha_hat, _ = self._optimize_parameter(neg_ll, (1e-8, upper))

        V = self._build_ou_vcv_single_alpha(
            ordered_names, lineage_info, alpha_hat, 1.0
        )
        try:
            V_inv = np.linalg.inv(V)
        except np.linalg.LinAlgError:
            return dict(
                model="OU1", log_likelihood=float("-inf"), k_params=3,
                params=dict(theta=0.0, alpha=float(alpha_hat), sigma2=0.0),
            )

        ones = np.ones(n)
        theta = float(ones @ V_inv @ x) / float(ones @ V_inv @ ones)
        e = x - theta
        sig2 = float(e @ V_inv @ e) / n

        sign, logdet = np.linalg.slogdet(V)
        ll = -0.5 * (n * np.log(2 * np.pi) + n * np.log(sig2) + logdet + n)

        return dict(
            model="OU1", log_likelihood=float(ll), k_params=3,
            params=dict(theta=float(theta), alpha=float(alpha_hat), sigma2=float(sig2)),
        )

    def _fit_oum(
        self, x: np.ndarray, ordered_names: List[str],
        lineage_info: Dict, regimes: List[str],
        root_regime: str, tree_height: float,
    ) -> Dict:
        """OUM: multi-optimum OU. k = R+2 (theta_1..R, alpha, sigma2)."""
        n = len(x)
        R = len(regimes)
        upper = 100.0 / tree_height if tree_height > 0 else 100.0

        def neg_ll(alpha):
            V = self._build_ou_vcv_single_alpha(
                ordered_names, lineage_info, alpha, 1.0
            )
            W = self._build_weight_matrix_single_alpha(
                ordered_names, lineage_info, regimes, alpha, root_regime,
            )

            try:
                sign, logdet = np.linalg.slogdet(V)
                if sign <= 0:
                    return 1e20
                V_inv = np.linalg.inv(V)
            except np.linalg.LinAlgError:
                return 1e20

            theta = self._gls_theta_hat(x, W, V_inv)
            e = x - W @ theta
            sig2 = float(e @ V_inv @ e) / n
            if sig2 <= 0:
                return 1e20

            ll = -0.5 * (n * np.log(2 * np.pi) + n * np.log(sig2) + logdet + n)
            return -ll

        alpha_hat, _ = self._optimize_parameter(neg_ll, (1e-8, upper))

        V = self._build_ou_vcv_single_alpha(
            ordered_names, lineage_info, alpha_hat, 1.0
        )
        W = self._build_weight_matrix_single_alpha(
            ordered_names, lineage_info, regimes, alpha_hat, root_regime,
        )

        try:
            V_inv = np.linalg.inv(V)
        except np.linalg.LinAlgError:
            return dict(
                model="OUM", log_likelihood=float("-inf"), k_params=R + 2,
                params=dict(
                    theta={r: 0.0 for r in regimes},
                    alpha=float(alpha_hat), sigma2=0.0,
                ),
            )

        theta = self._gls_theta_hat(x, W, V_inv)
        e = x - W @ theta
        sig2 = float(e @ V_inv @ e) / n

        sign, logdet = np.linalg.slogdet(V)
        ll = -0.5 * (n * np.log(2 * np.pi) + n * np.log(sig2) + logdet + n)

        theta_dict = {regimes[j]: float(theta[j]) for j in range(R)}
        return dict(
            model="OUM", log_likelihood=float(ll), k_params=R + 2,
            params=dict(
                theta=theta_dict,
                alpha=float(alpha_hat), sigma2=float(sig2),
            ),
        )

    def _fit_oumv(
        self, x: np.ndarray, ordered_names: List[str],
        lineage_info: Dict, regimes: List[str],
        root_regime: str, tree_height: float,
    ) -> Dict:
        """OUMV: multi-optimum + per-regime sigma^2. k = 2R+1."""
        n = len(x)
        R = len(regimes)

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

            W = self._build_weight_matrix_single_alpha(
                ordered_names, lineage_info, regimes, alpha, root_regime,
            )

            try:
                sign, logdet = np.linalg.slogdet(V)
                if sign <= 0:
                    return 1e20
                V_inv = np.linalg.inv(V)
            except np.linalg.LinAlgError:
                return 1e20

            theta = self._gls_theta_hat(x, W, V_inv)
            e = x - W @ theta
            quad = float(e @ V_inv @ e)

            ll = -0.5 * (n * np.log(2 * np.pi) + logdet + quad)
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

        W = self._build_weight_matrix_single_alpha(
            ordered_names, lineage_info, regimes, alpha, root_regime,
        )

        try:
            V_inv = np.linalg.inv(V)
        except np.linalg.LinAlgError:
            sigma2_dict = {r: float(sigma2s[i]) for i, r in enumerate(regimes)}
            return dict(
                model="OUMV", log_likelihood=float("-inf"), k_params=2 * R + 1,
                params=dict(
                    theta={r: 0.0 for r in regimes},
                    alpha=float(alpha), sigma2=sigma2_dict,
                ),
            )

        theta = self._gls_theta_hat(x, W, V_inv)
        e = x - W @ theta
        sign, logdet = np.linalg.slogdet(V)
        quad = float(e @ V_inv @ e)
        ll = -0.5 * (n * np.log(2 * np.pi) + logdet + quad)

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
        self, x: np.ndarray, ordered_names: List[str],
        lineage_info: Dict, regimes: List[str],
        root_regime: str, tree_height: float,
    ) -> Dict:
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

            try:
                sign, logdet = np.linalg.slogdet(V)
                if sign <= 0:
                    return 1e20
                V_inv = np.linalg.inv(V)
            except np.linalg.LinAlgError:
                return 1e20

            theta = self._gls_theta_hat(x, W, V_inv)
            e = x - W @ theta
            sig2 = float(e @ V_inv @ e) / n
            if sig2 <= 0:
                return 1e20

            ll = -0.5 * (n * np.log(2 * np.pi) + n * np.log(sig2) + logdet + n)
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

        try:
            V_inv = np.linalg.inv(V)
        except np.linalg.LinAlgError:
            alpha_dict = {r: float(alphas[i]) for i, r in enumerate(regimes)}
            return dict(
                model="OUMA", log_likelihood=float("-inf"), k_params=2 * R + 1,
                params=dict(
                    theta={r: 0.0 for r in regimes},
                    alpha=alpha_dict, sigma2=0.0,
                ),
            )

        theta = self._gls_theta_hat(x, W, V_inv)
        e = x - W @ theta
        sig2 = float(e @ V_inv @ e) / n

        sign, logdet = np.linalg.slogdet(V)
        ll = -0.5 * (n * np.log(2 * np.pi) + n * np.log(sig2) + logdet + n)

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
        self, x: np.ndarray, ordered_names: List[str],
        lineage_info: Dict, regimes: List[str],
        root_regime: str, tree_height: float,
    ) -> Dict:
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

            try:
                sign, logdet = np.linalg.slogdet(V)
                if sign <= 0:
                    return 1e20
                V_inv = np.linalg.inv(V)
            except np.linalg.LinAlgError:
                return 1e20

            theta = self._gls_theta_hat(x, W, V_inv)
            e = x - W @ theta
            quad = float(e @ V_inv @ e)

            ll = -0.5 * (n * np.log(2 * np.pi) + logdet + quad)
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

        try:
            V_inv = np.linalg.inv(V)
        except np.linalg.LinAlgError:
            alpha_dict = {r: float(alphas[i]) for i, r in enumerate(regimes)}
            sigma2_dict = {r: float(sigma2s[i]) for i, r in enumerate(regimes)}
            return dict(
                model="OUMVA", log_likelihood=float("-inf"), k_params=3 * R,
                params=dict(
                    theta={r: 0.0 for r in regimes},
                    alpha=alpha_dict, sigma2=sigma2_dict,
                ),
            )

        theta = self._gls_theta_hat(x, W, V_inv)
        e = x - W @ theta
        sign, logdet = np.linalg.slogdet(V)
        quad = float(e @ V_inv @ e)
        ll = -0.5 * (n * np.log(2 * np.pi) + logdet + quad)

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
        self, results: List[Dict], n: int,
        x: np.ndarray = None, vcv_total: np.ndarray = None,
        regime_assignments: Dict = None, ordered_names: List[str] = None,
    ) -> List[Dict]:
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
        delta_aiccs = np.array([r["delta_aicc"] for r in results])
        raw_weights = np.exp(-0.5 * delta_aiccs)
        total = raw_weights.sum()
        aicc_weights = raw_weights / total if total > 0 else raw_weights

        for r, w in zip(results, aicc_weights):
            r["aicc_weight"] = float(w)

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
                    sig2_vals = list(sig2.values())
                    sig2_model = sum(sig2_vals) / len(sig2_vals) if sig2_vals else 0
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
        self, results: List[Dict], n: int, regimes: List[str]
    ) -> None:
        print("OUwie Model Comparison (Multi-Regime OU)")
        print(f"\nNumber of tips: {n}")
        print(f"Regimes: {len(regimes)} ({', '.join(regimes)})\n")

        header = (
            f"{'Model':<8}{'k':<5}{'LL':<12}"
            f"{'AIC':<10}{'AICc':<10}{'dAICc':<9}{'AICcW':<9}"
            f"{'BIC':<10}{'dBIC':<9}{'R2':<7}"
        )
        print(header)

        for r in results:
            print(
                f"{r['model']:<8}{r['k_params']:<5}{r['log_likelihood']:<12.3f}"
                f"{r['aic']:<10.2f}{r['aicc']:<10.2f}"
                f"{r['delta_aicc']:<9.2f}{r['aicc_weight']:<9.3f}"
                f"{r['bic']:<10.2f}{r['delta_bic']:<9.2f}"
                f"{r['r_squared']:<7.3f}"
            )

        print(f"\nBest model (AICc): {results[0]['model']}")
        best_bic = min(results, key=lambda r: r["bic"])["model"]
        print(f"Best model (BIC):  {best_bic}")

        # Print parameter details for best model
        best = results[0]
        print(f"\nParameter estimates ({best['model']}):")
        params = best["params"]

        if isinstance(params.get("theta"), dict):
            print("  Theta (trait optima):")
            for r_name, val in params["theta"].items():
                print(f"    {r_name}: {val:.4f}")
        elif "theta" in params:
            print(f"  Theta: {params['theta']:.4f}")
        if "z0" in params:
            print(f"  z0: {params['z0']:.4f}")

        if isinstance(params.get("alpha"), dict):
            print("  Alpha (selection strength):")
            for r_name, val in params["alpha"].items():
                print(f"    {r_name}: {val:.4f}")
        elif "alpha" in params:
            print(f"  Alpha: {params['alpha']:.4f}")

        if isinstance(params.get("sigma2"), dict):
            print("  Sigma-squared (diffusion rate):")
            for r_name, val in params["sigma2"].items():
                print(f"    {r_name}: {val:.4f}")
        elif "sigma2" in params:
            print(f"  Sigma-squared: {params['sigma2']:.4f}")

    def _print_json_output(
        self, results: List[Dict], n: int, regimes: List[str]
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
