from __future__ import annotations

import math
import sys
from collections import deque

from .base import Tree
from ...errors import PhykitUserError


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


def _json_load(*args, **kwargs):
    import json

    return json.load(*args, **kwargs)


class _LazyNumpy:
    def __init__(self):
        self._module = None

    def __getattr__(self, name):
        module = self._module
        if module is None:
            import numpy as _np

            module = _np
            self._module = module

        attr = getattr(module, name)
        setattr(self, name, attr)
        return attr


np = _LazyNumpy()
_CHO_FACTOR = None
_CHO_SOLVE = None
_MINIMIZE_SCALAR = None
_COLUMN_SUM_METHOD_MIN_SIZE = 180


def _permutation_p_value_ge(permutations: np.ndarray, observed: float) -> float:
    if permutations.size == 0:
        return float("nan")
    return float(np.count_nonzero(permutations >= observed) / permutations.size)


def _matrix_column_sums(matrix: np.ndarray) -> np.ndarray:
    if matrix.shape[0] >= _COLUMN_SUM_METHOD_MIN_SIZE:
        return matrix.sum(axis=0)
    return np.sum(matrix, axis=0)


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


def minimize_scalar(*args, **kwargs):
    global _MINIMIZE_SCALAR

    if _MINIMIZE_SCALAR is None:
        from scipy.optimize import minimize_scalar as _minimize_scalar

        _MINIMIZE_SCALAR = _minimize_scalar

    return _MINIMIZE_SCALAR(*args, **kwargs)


def _chi2_sf_df1(chi2_stat: float) -> float:
    return math.erfc(math.sqrt(chi2_stat / 2.0))


def _inverse_from_spd_matrix(matrix):
    try:
        factor = cho_factor(matrix, lower=True, check_finite=False)
        return cho_solve(
            factor,
            np.eye(matrix.shape[0]),
            check_finite=False,
        )
    except (np.linalg.LinAlgError, FloatingPointError, ValueError):
        return np.linalg.inv(matrix)


class NetworkSignal(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"])
        self.trait_data_path = parsed["trait_data_path"]
        self.method = parsed["method"]
        self.permutations = parsed["permutations"]
        self.hybrid_edges = parsed["hybrid_edges"]
        self.quartet_json = parsed["quartet_json"]
        self.json_output = parsed["json_output"]

    def process_args(self, args) -> dict[str, str]:
        return dict(
            tree_file_path=args.tree,
            trait_data_path=args.trait_data,
            method=getattr(args, "method", "both"),
            permutations=getattr(args, "permutations", 1000),
            hybrid_edges=getattr(args, "hybrid", None),
            quartet_json=getattr(args, "quartet_json", None),
            json_output=getattr(args, "json", False),
        )

    # ------------------------------------------------------------------
    # Network DAG construction
    # ------------------------------------------------------------------

    @staticmethod
    def _tree_to_dag(tree):
        """Convert a Biopython tree to an internal DAG representation.

        Polytomous nodes (>2 children) are faithfully represented as star
        topologies in the DAG where all children share the same parent node.
        The VCV matrix will show equal covariance among all children of a
        polytomy, which is the correct representation of an unresolved
        relationship.

        Returns:
            nodes: list of integer node IDs in topological order (root first,
                   from level-order traversal)
            parents: dict mapping node_id -> list of (parent_id, edge_length, gamma_weight)
            tip_map: dict mapping node_id -> taxon_name (for terminal nodes)
        """
        nodes = []
        parents = {}
        tip_map = {}

        # Assign integer IDs via level-order traversal (BFS)
        clade_to_id = {}
        node_id = 0
        queue = deque()
        root_clade = tree.root
        clade_to_id[id(root_clade)] = node_id
        nodes.append(node_id)
        parents[node_id] = []
        node_id += 1
        queue.append(root_clade)

        while queue:
            clade = queue.popleft()
            parent_id = clade_to_id[id(clade)]

            for child in clade.clades:
                child_id = node_id
                clade_to_id[id(child)] = child_id
                nodes.append(child_id)
                edge_len = child.branch_length if child.branch_length is not None else 0.0
                parents[child_id] = [(parent_id, edge_len, 1.0)]
                node_id += 1

                if child.is_terminal():
                    tip_map[child_id] = child.name
                else:
                    queue.append(child)

        return nodes, parents, tip_map

    @staticmethod
    def _add_hybrid_edge(nodes, parents, tip_map, donor, recipient, gamma):
        """Make the recipient a hybrid node by adding a second parent edge
        from the donor's tree parent.

        The recipient gets two parents:
        - Original tree parent with weight (1 - gamma)
        - Donor's tree parent with weight gamma and donor's edge length

        Validates: donor/recipient exist, 0 < gamma < 0.5
        """
        if gamma <= 0 or gamma >= 0.5:
            raise PhykitUserError(
                [
                    f"gamma must be strictly between 0 and 0.5, got {gamma}.",
                    "gamma represents the minor inheritance proportion.",
                ],
                code=2,
            )

        # Build name -> node_id map
        name_to_id = {v: k for k, v in tip_map.items()}

        if donor not in name_to_id:
            raise PhykitUserError(
                [
                    f"Donor taxon '{donor}' not found in the tree.",
                    f"Available taxa: {', '.join(sorted(name_to_id.keys()))}",
                ],
                code=2,
            )
        if recipient not in name_to_id:
            raise PhykitUserError(
                [
                    f"Recipient taxon '{recipient}' not found in the tree.",
                    f"Available taxa: {', '.join(sorted(name_to_id.keys()))}",
                ],
                code=2,
            )

        donor_id = name_to_id[donor]
        recipient_id = name_to_id[recipient]

        # The donor's tree parent
        if not parents.get(donor_id):
            raise PhykitUserError(
                [f"Donor '{donor}' has no parent (is root?)."],
                code=2,
            )
        donor_parent_id = parents[donor_id][0][0]
        donor_edge_len = parents[donor_id][0][1]

        # Update recipient's original parent to weight (1 - gamma)
        orig_entries = parents[recipient_id]
        if len(orig_entries) != 1:
            raise PhykitUserError(
                [
                    f"Recipient '{recipient}' already has {len(orig_entries)} parents.",
                    "Only single-parent (tree) nodes can be converted to hybrids.",
                ],
                code=2,
            )
        orig_parent_id, orig_edge_len, _ = orig_entries[0]
        parents[recipient_id] = [
            (orig_parent_id, orig_edge_len, 1.0 - gamma),
            (donor_parent_id, donor_edge_len, gamma),
        ]

    # ------------------------------------------------------------------
    # Network VCV (Bastide et al. 2018)
    # ------------------------------------------------------------------

    @staticmethod
    def _compute_network_vcv(nodes, parents, tip_map, ordered_tips):
        """Compute the VCV matrix on a phylogenetic network using the
        Bastide et al. 2018 recursive algorithm.

        Process nodes in topological order (root first).  Maintain a full
        node x node V matrix.

        Root:      V[root, root] = 0
        Tree node c (parent p, length l):
                   V[c, c] = V[p, p] + l
                   V[c, j] = V[p, j]  for all j processed before c
        Hybrid node h (parents p1 with gamma, p2 with 1-gamma, lengths l1, l2):
                   V[h, h] = gamma^2*(V[p1,p1]+l1) + (1-gamma)^2*(V[p2,p2]+l2)
                              + 2*gamma*(1-gamma)*V[p1,p2]
                   V[h, j] = gamma*V[p1,j] + (1-gamma)*V[p2,j]

        Finally extract the tip x tip submatrix.
        """
        n_nodes = len(nodes)
        node_to_idx = {nid: i for i, nid in enumerate(nodes)}

        V = np.zeros((n_nodes, n_nodes))

        for nid in nodes:
            idx = node_to_idx[nid]
            parent_list = parents.get(nid, [])

            if len(parent_list) == 0:
                # Root node
                V[idx, idx] = 0.0
            elif len(parent_list) == 1:
                # Tree node
                p_id, edge_len, _ = parent_list[0]
                p_idx = node_to_idx[p_id]
                V[idx, idx] = V[p_idx, p_idx] + edge_len
                # Copy covariances from parent (only preceding nodes)
                row = V[p_idx, :idx]
                V[idx, :idx] = row
                V[:idx, idx] = row
            elif len(parent_list) == 2:
                # Hybrid node
                p1_id, l1, g1 = parent_list[0]
                p2_id, l2, g2 = parent_list[1]
                p1_idx = node_to_idx[p1_id]
                p2_idx = node_to_idx[p2_id]

                # Variance
                V[idx, idx] = (
                    g1 * g1 * (V[p1_idx, p1_idx] + l1)
                    + g2 * g2 * (V[p2_idx, p2_idx] + l2)
                    + 2 * g1 * g2 * V[p1_idx, p2_idx]
                )
                # Covariances (only preceding nodes)
                cov = g1 * V[p1_idx, :idx] + g2 * V[p2_idx, :idx]
                V[idx, :idx] = cov
                V[:idx, idx] = cov
            else:
                raise PhykitUserError(
                    [f"Node {nid} has {len(parent_list)} parents; expected 0, 1, or 2."],
                    code=2,
                )

        # Extract tip x tip submatrix
        name_to_id = {v: k for k, v in tip_map.items()}
        tip_indices = np.fromiter(
            (node_to_idx[name_to_id[tip_name]] for tip_name in ordered_tips),
            dtype=np.intp,
            count=len(ordered_tips),
        )
        return V[np.ix_(tip_indices, tip_indices)]

    # ------------------------------------------------------------------
    # Quartet JSON inference
    # ------------------------------------------------------------------

    @staticmethod
    def _parse_topology_string(topo_str):
        """Parse '{A, B} | {C, D}' to (frozenset({'A','B'}), frozenset({'C','D'}))."""
        parts = topo_str.split("|")
        if len(parts) != 2:
            raise PhykitUserError(
                [f"Cannot parse topology string: '{topo_str}'"],
                code=2,
            )

        def parse_set(s):
            s = s.strip().strip("{}")
            return frozenset(t.strip() for t in s.split(","))

        return (parse_set(parts[0]), parse_set(parts[1]))

    @staticmethod
    def _identify_swap_pair(dominant, minor):
        """Identify the 2 taxa that switch sides between dominant and minor topologies.

        Given dominant (left_dom, right_dom) and minor (left_minor, right_minor),
        the taxa that 'moved' are: (left_dom - left_minor) union (left_minor - left_dom).
        """
        left_dom, right_dom = dominant
        left_minor, right_minor = minor

        from_left = left_dom - left_minor
        to_left = left_minor - left_dom
        swap = from_left | to_left
        return swap

    @staticmethod
    def _top_two_topology_indices(counts):
        """Return the first two topology indices by descending count."""
        c0, c1, c2 = counts
        if c0 >= c1:
            if c1 >= c2:
                return 0, 1
            if c0 >= c2:
                return 0, 2
            return 2, 0
        if c0 >= c2:
            return 1, 0
        if c1 >= c2:
            return 1, 2
        return 2, 1

    @staticmethod
    def _infer_hybrid_edges(quartet_data):
        """Infer hybrid edges from quartet_network JSON data.

        For each 'hybrid' classified quartet, identify the swap pair between
        the dominant and the next-most-elevated topology to determine donor
        and recipient.

        Aggregates evidence across multiple quartets: counts how many quartets
        support each swap pair, averages gamma estimates, and returns edges
        ordered by descending frequency.

        Returns a list of dicts with keys: donor, recipient, gamma, n_quartets.
        """
        pair_stats = {}
        parsed_topology_cache = {}
        parse_topology = NetworkSignal._parse_topology_string
        top_two_indices = NetworkSignal._top_two_topology_indices
        pair_stats_get = pair_stats.get
        topology_cache_get = parsed_topology_cache.get

        for q in quartet_data.get("quartets", []):
            if q["classification"] != "hybrid":
                continue

            taxa = q["taxa"]
            counts = q["counts"]
            dominant_topo_str = q["dominant_topology"]

            # Preserve topology-string validation without reparsing repeated JSON rows.
            if topology_cache_get(dominant_topo_str) is None:
                parsed_topology_cache[dominant_topo_str] = parse_topology(
                    dominant_topo_str
                )

            # The three topologies for quartet (a,b,c,d):
            # topo 0: {a,b}|{c,d}, topo 1: {a,c}|{b,d}, topo 2: {a,d}|{b,c}
            a, b, c, d = taxa

            # Find the dominant and second-most elevated topologies.
            dominant_idx, minor_idx = top_two_indices(counts)
            if dominant_idx == 0:
                taxon_1, taxon_2 = (b, c) if minor_idx == 1 else (b, d)
            elif dominant_idx == 1:
                taxon_1, taxon_2 = (b, c) if minor_idx == 0 else (c, d)
            else:
                taxon_1, taxon_2 = (b, d) if minor_idx == 0 else (c, d)

            # Preserve the sorted pair orientation returned by sorted(swap).
            pair_key = (
                (taxon_1, taxon_2)
                if taxon_1 <= taxon_2
                else (taxon_2, taxon_1)
            )
            total = counts[0] + counts[1] + counts[2]
            gamma_est = counts[minor_idx] / total if total > 0 else 0.1
            if gamma_est > 0.499:
                gamma_est = 0.499
            elif gamma_est < 0.001:
                gamma_est = 0.001

            stats = pair_stats_get(pair_key)
            if stats is None:
                pair_stats[pair_key] = [1, gamma_est]
            else:
                stats[0] += 1
                stats[1] += gamma_est

        return [
            {
                "donor": donor,
                "recipient": recipient,
                "gamma": round(gamma_sum / n_quartets, 4),
                "n_quartets": n_quartets,
            }
            for (donor, recipient), (n_quartets, gamma_sum) in sorted(
                pair_stats.items(),
                key=lambda item: item[1][0],
                reverse=True,
            )
        ]

    # ------------------------------------------------------------------
    # K and Lambda (same formulas as phylogenetic_signal.py)
    # ------------------------------------------------------------------

    @staticmethod
    def _log_likelihood(x, C):
        """Concentrated log-likelihood with sigma^2 profiled out.

        Matches phytools::phylosig(method='lambda') internals.
        """
        try:
            return NetworkSignal._log_likelihood_cholesky(x, C)
        except (np.linalg.LinAlgError, FloatingPointError, ValueError):
            return NetworkSignal._log_likelihood_inverse(x, C)

    @staticmethod
    def _log_likelihood_cholesky(x, C):
        """Cholesky-backed likelihood for positive-definite VCV matrices."""
        n = len(x)
        ones = np.ones(n)
        factor = cho_factor(C, lower=True, check_finite=False)

        solve_rhs = np.empty((n, 2), dtype=np.result_type(x, C))
        solve_rhs[:, 0] = 1.0
        solve_rhs[:, 1] = x
        solved = cho_solve(factor, solve_rhs, check_finite=False)
        C_inv_ones = solved[:, 0]
        C_inv_x = solved[:, 1]

        # GLS estimate of phylogenetic mean
        a_hat = float((ones @ C_inv_x) / (ones @ C_inv_ones))

        # Residuals
        e = x - a_hat
        C_inv_e = C_inv_x - a_hat * C_inv_ones

        # MLE of sigma^2 (profiled out)
        sig2 = float(e @ C_inv_e / n)
        if sig2 <= 0:
            return NetworkSignal._log_likelihood_inverse(x, C)

        # Concentrated log-likelihood
        logdet_C = 2.0 * float(np.log(np.diagonal(factor[0])).sum())
        logdet_sig2C = n * np.log(sig2) + logdet_C
        ll = -0.5 * (n * np.log(2 * np.pi) + logdet_sig2C + n)

        return ll, a_hat

    @staticmethod
    def _log_likelihood_inverse(x, C):
        """Inverse-based likelihood implementation retained as a fallback."""
        n = len(x)
        ones = np.ones(n)

        C_inv = np.linalg.inv(C)

        # GLS estimate of phylogenetic mean
        a_hat = float((ones @ C_inv @ x) / (ones @ C_inv @ ones))

        # Residuals
        e = x - a_hat

        # MLE of sigma^2 (profiled out)
        sig2 = float(e @ C_inv @ e / n)

        # Concentrated log-likelihood
        sign, logdet_C = np.linalg.slogdet(C)
        logdet_sig2C = n * np.log(sig2) + logdet_C
        ll = -0.5 * (n * np.log(2 * np.pi) + logdet_sig2C + n)

        return ll, a_hat

    @staticmethod
    def _blombergs_k(x, vcv, n_perm):
        """Blomberg's K following phytools::phylosig(method='K').

        K = (e'e / e'C_inv'e) / ((trace(C) - n/sum(C_inv)) / (n-1))
        """
        n = len(x)
        ones = np.ones(n)
        C = vcv.copy()

        C_inv = _inverse_from_spd_matrix(C)
        sum_C_inv = float(ones @ C_inv @ ones)

        # GLS estimate of phylogenetic mean
        a_hat = float((ones @ C_inv @ x) / sum_C_inv)

        # Residuals
        e = x - a_hat

        # Observed ratio
        observed_ratio = float(e @ e) / float(e @ C_inv @ e)

        # Expected ratio under BM
        expected_ratio = (np.trace(C) - n / sum_C_inv) / (n - 1)

        # K statistic
        K = observed_ratio / expected_ratio

        # Permutation test
        k_perm = NetworkSignal._blombergs_k_permutations(
            x,
            C_inv,
            sum_C_inv,
            expected_ratio,
            n_perm,
        )

        p_value = _permutation_p_value_ge(k_perm, K)

        return dict(
            K=float(K),
            p_value=p_value,
            permutations=n_perm,
        )

    @staticmethod
    def _blombergs_k_permutations(
        x,
        C_inv,
        sum_C_inv,
        expected_ratio,
        n_perm,
        batch_size=128,
    ):
        rng = np.random.default_rng(seed=42)
        n = len(x)
        weights = _matrix_column_sums(C_inv)
        x_sum = float(x.sum())
        x_sumsq = float(x @ x)
        k_perm = np.empty(n_perm, dtype=np.float64)

        for start in range(0, n_perm, batch_size):
            stop = min(start + batch_size, n_perm)
            batch_len = stop - start
            perms = np.empty((batch_len, n), dtype=np.float64)
            for row in range(batch_len):
                perms[row] = rng.permutation(x)

            weighted_sums = perms @ weights
            a_perm = weighted_sums / sum_C_inv
            numerator = x_sumsq - (2.0 * a_perm * x_sum) + (n * a_perm * a_perm)
            C_inv_perm = perms @ C_inv
            perm_quadratics = np.einsum("ij,ij->i", perms, C_inv_perm)
            denominator = perm_quadratics - (a_perm * a_perm * sum_C_inv)
            k_perm[start:stop] = (numerator / denominator) / expected_ratio

        return k_perm

    @staticmethod
    def _pagels_lambda(x, vcv, max_lambda=1.0):
        """Pagel's lambda with multi-interval optimization."""
        n = len(x)
        diag_vals = vcv.diagonal().copy()
        diag_step = vcv.shape[0] + 1
        niter = 10

        def neg_ll(lam):
            C_lam = vcv * lam
            C_lam.ravel()[::diag_step] = diag_vals
            try:
                ll, _ = NetworkSignal._log_likelihood(x, C_lam)
                return -ll
            except (np.linalg.LinAlgError, FloatingPointError, ValueError):
                return 1e10

        bounds_lo = np.linspace(0, max_lambda - max_lambda / niter, niter)
        bounds_hi = np.linspace(max_lambda / niter, max_lambda, niter)

        best_ll = -np.inf
        lambda_hat = 0.0
        for lo, hi in zip(bounds_lo, bounds_hi):
            res = minimize_scalar(neg_ll, bounds=(lo, hi), method="bounded")
            ll_val = -res.fun
            if ll_val > best_ll:
                best_ll = ll_val
                lambda_hat = res.x

        # Log-likelihood at fitted lambda
        C_fitted = vcv * lambda_hat
        C_fitted.ravel()[::diag_step] = diag_vals
        ll_fitted, _ = NetworkSignal._log_likelihood(x, C_fitted)

        # Log-likelihood at lambda = 0
        C_zero = vcv * 0.0
        C_zero.ravel()[::diag_step] = diag_vals
        ll_zero, _ = NetworkSignal._log_likelihood(x, C_zero)

        # Likelihood ratio test (1 df)
        lr_stat = 2.0 * (ll_fitted - ll_zero)
        lr_stat = max(lr_stat, 0.0)
        p_value = _chi2_sf_df1(lr_stat)

        return dict(
            **{"lambda": float(lambda_hat)},
            log_likelihood=float(ll_fitted),
            p_value=p_value,
        )

    # ------------------------------------------------------------------
    # Trait parsing (same as phylogenetic_signal.py)
    # ------------------------------------------------------------------

    def _parse_trait_file(self, path, tree_tips):
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

    # ------------------------------------------------------------------
    # Main entry points
    # ------------------------------------------------------------------

    def _validate_tree(self, tree) -> None:
        root = tree.root
        tip_count = 0
        polytomy_count = 0
        missing_branch_lengths = False
        stack = [root]
        pop = stack.pop
        extend = stack.extend

        while stack:
            clade = pop()
            children = clade.clades

            if clade is not root and clade.branch_length is None:
                missing_branch_lengths = True

            if children:
                child_count = len(children)
                if child_count > 2 and not (clade is root and child_count == 3):
                    polytomy_count += 1
                extend(children)
            else:
                tip_count += 1

        if tip_count < 3:
            raise PhykitUserError(
                ["Tree must have at least 3 tips for phylogenetic signal analysis."],
                code=2,
            )
        if missing_branch_lengths:
            raise PhykitUserError(
                ["All branches in the tree must have lengths."],
                code=2,
            )
        if polytomy_count > 0:
            print(
                f"Warning: tree contains {polytomy_count} polytomous node(s). "
                "These are treated as star topologies in the network VCV.",
                file=sys.stderr,
            )

    def run(self):
        tree = self.read_tree_file_unmodified()
        self._validate_tree(tree)

        tree_tips = self.get_tip_names_from_tree(tree)
        traits = self._parse_trait_file(self.trait_data_path, tree_tips)

        ordered_names = sorted(traits.keys())
        x = np.array([traits[name] for name in ordered_names])

        # Build network DAG
        nodes, parents, tip_map = self._tree_to_dag(tree)

        # Add hybrid edges
        if self.quartet_json:
            try:
                with open(self.quartet_json) as f:
                    quartet_data = _json_load(f)
            except FileNotFoundError:
                raise PhykitUserError(
                    [
                        f"{self.quartet_json} corresponds to no such file or directory.",
                        "Please check filename and pathing",
                    ],
                    code=2,
                )
            edges = self._infer_hybrid_edges(quartet_data)
            for edge in edges:
                try:
                    self._add_hybrid_edge(
                        nodes, parents, tip_map,
                        edge["donor"], edge["recipient"], edge["gamma"],
                    )
                except PhykitUserError as e:
                    print(f"Warning: skipping hybrid edge {edge['donor']} -> {edge['recipient']}: {e.messages[0]}", file=sys.stderr)

        if self.hybrid_edges:
            for edge_spec in self.hybrid_edges:
                parts = edge_spec.split(":")
                if len(parts) != 3:
                    raise PhykitUserError(
                        [
                            f"Invalid hybrid edge spec: '{edge_spec}'.",
                            "Expected format: donor:recipient:gamma",
                        ],
                        code=2,
                    )
                donor, recipient, gamma_str = parts
                try:
                    gamma = float(gamma_str)
                except ValueError:
                    raise PhykitUserError(
                        [f"Invalid gamma value: '{gamma_str}'."],
                        code=2,
                    )
                self._add_hybrid_edge(nodes, parents, tip_map, donor, recipient, gamma)

        # Compute VCV
        vcv = self._compute_network_vcv(nodes, parents, tip_map, ordered_names)

        # Compute signal
        results = {}
        if self.method == "blombergs_k":
            results = self._blombergs_k(x, vcv, self.permutations)
        elif self.method == "lambda":
            results = self._pagels_lambda(x, vcv, max_lambda=1.0)
        elif self.method == "both":
            k_results = self._blombergs_k(x, vcv, self.permutations)
            l_results = self._pagels_lambda(x, vcv, max_lambda=1.0)
            results = {"blombergs_k": k_results, "pagels_lambda": l_results}

        self._output(results, ordered_names, nodes, parents, tip_map)

    def _output(self, results, ordered_names=None, nodes=None, parents=None, tip_map=None):
        if self.json_output:
            print_json(results)
            return

        try:
            lines = []
            # Print hybrid edge info
            if tip_map and parents:
                children_by_parent = {}
                for child_id, child_name in tip_map.items():
                    child_parents = parents.get(child_id)
                    if child_parents:
                        parent_id = child_parents[0][0]
                        children_by_parent.setdefault(parent_id, []).append(
                            (child_id, child_name)
                        )

                for nid, name in tip_map.items():
                    if len(parents.get(nid, [])) == 2:
                        gamma = parents[nid][1][2]
                        # Find the donor name from the second parent
                        donor_parent_id = parents[nid][1][0]
                        donor_name = None
                        for child_id, child_name in children_by_parent.get(
                            donor_parent_id, []
                        ):
                            if child_id != nid:
                                donor_name = child_name
                                break
                        if donor_name:
                            lines.append(
                                f"Hybrid edge: {donor_name} -> {name} (gamma={gamma:.4f})"
                            )

            if ordered_names:
                lines.append(f"Network taxa: {len(ordered_names)}")
                lines.append("---")

            if self.method == "blombergs_k":
                lines.append(
                    f"Blomberg's K: {results['K']:.4f}    "
                    f"p-value: {results['p_value']:.4f}"
                )
            elif self.method == "lambda":
                lines.append(
                    f"Pagel's lambda: {results['lambda']:.4f}    "
                    f"log-likelihood: {results['log_likelihood']:.4f}    "
                    f"p-value: {results['p_value']:.4f}"
                )
            elif self.method == "both":
                k = results["blombergs_k"]
                l = results["pagels_lambda"]
                lines.append(
                    f"Blomberg's K: {k['K']:.4f}    "
                    f"p-value: {k['p_value']:.4f}"
                )
                lines.append(
                    f"Pagel's lambda: {l['lambda']:.4f}    "
                    f"log-likelihood: {l['log_likelihood']:.4f}    "
                    f"p-value: {l['p_value']:.4f}"
                )
            if lines:
                print("\n".join(lines))
        except BrokenPipeError:
            pass
