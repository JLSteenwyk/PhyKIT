import json
import sys
from collections import Counter, deque
from typing import Dict, List, Tuple

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.stats import chi2

from .base import Tree
from ...helpers.json_output import print_json
from ...errors import PhykitUserError


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

    def process_args(self, args) -> Dict[str, str]:
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
                for j in range(idx):
                    V[idx, j] = V[p_idx, j]
                    V[j, idx] = V[p_idx, j]
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
                for j in range(idx):
                    cov = g1 * V[p1_idx, j] + g2 * V[p2_idx, j]
                    V[idx, j] = cov
                    V[j, idx] = cov
            else:
                raise PhykitUserError(
                    [f"Node {nid} has {len(parent_list)} parents; expected 0, 1, or 2."],
                    code=2,
                )

        # Extract tip x tip submatrix
        name_to_id = {v: k for k, v in tip_map.items()}
        tip_ids = [name_to_id[tip_name] for tip_name in ordered_tips]

        n_tips = len(ordered_tips)
        vcv = np.zeros((n_tips, n_tips))
        for i, tid_i in enumerate(tip_ids):
            for j, tid_j in enumerate(tip_ids):
                vcv[i, j] = V[node_to_idx[tid_i], node_to_idx[tid_j]]

        return vcv

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
    def _infer_hybrid_edges(quartet_data):
        """Infer hybrid edges from quartet_network JSON data.

        For each 'hybrid' classified quartet, identify the swap pair between
        the dominant and the next-most-elevated topology to determine donor
        and recipient.

        Aggregates evidence across multiple quartets: counts how many quartets
        support each swap pair, averages gamma estimates, and returns edges
        ordered by frequency via Counter.most_common().

        Returns a list of dicts with keys: donor, recipient, gamma, n_quartets.
        """
        pair_counts = Counter()
        pair_gammas = {}  # (donor, recipient) -> list of gamma values

        for q in quartet_data.get("quartets", []):
            if q["classification"] != "hybrid":
                continue

            taxa = q["taxa"]
            counts = q["counts"]
            dominant_topo_str = q["dominant_topology"]

            # Parse dominant topology
            dominant = NetworkSignal._parse_topology_string(dominant_topo_str)

            # The three topologies for quartet (a,b,c,d):
            # topo 0: {a,b}|{c,d}, topo 1: {a,c}|{b,d}, topo 2: {a,d}|{b,c}
            a, b, c, d = taxa
            topos = [
                (frozenset({a, b}), frozenset({c, d})),
                (frozenset({a, c}), frozenset({b, d})),
                (frozenset({a, d}), frozenset({b, c})),
            ]

            # Find the dominant index
            dominant_idx = counts.index(max(counts))

            # Find the second-most elevated (minor) topology
            sorted_with_idx = sorted(enumerate(counts), key=lambda p: -p[1])
            minor_idx = sorted_with_idx[1][0]

            minor_topo = topos[minor_idx]
            dominant_topo = topos[dominant_idx]

            swap = NetworkSignal._identify_swap_pair(dominant_topo, minor_topo)
            swap_list = sorted(swap)

            if len(swap_list) == 2:
                # Use the taxon from the minor side as donor and
                # the other as recipient; gamma estimated from count ratio
                total = sum(counts)
                minor_count = counts[minor_idx]
                gamma_est = minor_count / total if total > 0 else 0.1
                gamma_est = min(gamma_est, 0.499)
                gamma_est = max(gamma_est, 0.001)

                pair_key = (swap_list[0], swap_list[1])
                pair_counts[pair_key] += 1
                pair_gammas.setdefault(pair_key, []).append(gamma_est)

        edges = []
        for (donor, recipient), n_quartets in pair_counts.most_common():
            avg_gamma = sum(pair_gammas[(donor, recipient)]) / n_quartets
            avg_gamma = round(avg_gamma, 4)
            edges.append({
                "donor": donor,
                "recipient": recipient,
                "gamma": avg_gamma,
                "n_quartets": n_quartets,
            })

        return edges

    # ------------------------------------------------------------------
    # K and Lambda (same formulas as phylogenetic_signal.py)
    # ------------------------------------------------------------------

    @staticmethod
    def _log_likelihood(x, C):
        """Concentrated log-likelihood with sigma^2 profiled out.

        Matches phytools::phylosig(method='lambda') internals.
        """
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

        C_inv = np.linalg.inv(C)
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
        rng = np.random.default_rng(seed=42)
        k_perm = np.empty(n_perm)
        for perm_i in range(n_perm):
            x_perm = rng.permutation(x)
            a_p = float((ones @ C_inv @ x_perm) / sum_C_inv)
            e_p = x_perm - a_p
            obs_p = float(e_p @ e_p) / float(e_p @ C_inv @ e_p)
            k_perm[perm_i] = obs_p / expected_ratio

        p_value = float(np.mean(k_perm >= K))

        return dict(
            K=float(K),
            p_value=p_value,
            permutations=n_perm,
        )

    @staticmethod
    def _pagels_lambda(x, vcv, max_lambda=1.0):
        """Pagel's lambda with multi-interval optimization."""
        n = len(x)
        diag_vals = np.diag(vcv).copy()
        niter = 10

        def neg_ll(lam):
            C_lam = vcv * lam
            np.fill_diagonal(C_lam, diag_vals)
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
        np.fill_diagonal(C_fitted, diag_vals)
        ll_fitted, _ = NetworkSignal._log_likelihood(x, C_fitted)

        # Log-likelihood at lambda = 0
        C_zero = vcv * 0.0
        np.fill_diagonal(C_zero, diag_vals)
        ll_zero, _ = NetworkSignal._log_likelihood(x, C_zero)

        # Likelihood ratio test (1 df)
        lr_stat = 2.0 * (ll_fitted - ll_zero)
        lr_stat = max(lr_stat, 0.0)
        p_value = float(chi2.sf(lr_stat, df=1))

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

    # ------------------------------------------------------------------
    # Main entry points
    # ------------------------------------------------------------------

    def _validate_tree(self, tree) -> None:
        tips = list(tree.get_terminals())
        if len(tips) < 3:
            raise PhykitUserError(
                ["Tree must have at least 3 tips for phylogenetic signal analysis."],
                code=2,
            )
        for clade in tree.find_clades():
            if clade.branch_length is None and clade != tree.root:
                raise PhykitUserError(
                    ["All branches in the tree must have lengths."],
                    code=2,
                )

    def run(self):
        tree = self.read_tree_file()
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
                    quartet_data = json.load(f)
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
            # Print hybrid edge info
            if tip_map and parents:
                for nid, name in tip_map.items():
                    if len(parents.get(nid, [])) == 2:
                        gamma = parents[nid][1][2]
                        # Find the donor name from the second parent
                        donor_parent_id = parents[nid][1][0]
                        donor_name = None
                        for child_id, child_name in tip_map.items():
                            if child_id != nid and parents.get(child_id):
                                if parents[child_id][0][0] == donor_parent_id:
                                    donor_name = child_name
                                    break
                        if donor_name:
                            print(f"Hybrid edge: {donor_name} -> {name} (gamma={gamma:.4f})")

            if ordered_names:
                print(f"Network taxa: {len(ordered_names)}")
                print("---")

            if self.method == "blombergs_k":
                print(
                    f"Blomberg's K: {results['K']:.4f}    "
                    f"p-value: {results['p_value']:.4f}"
                )
            elif self.method == "lambda":
                print(
                    f"Pagel's lambda: {results['lambda']:.4f}    "
                    f"log-likelihood: {results['log_likelihood']:.4f}    "
                    f"p-value: {results['p_value']:.4f}"
                )
            elif self.method == "both":
                k = results["blombergs_k"]
                l = results["pagels_lambda"]
                print(
                    f"Blomberg's K: {k['K']:.4f}    "
                    f"p-value: {k['p_value']:.4f}"
                )
                print(
                    f"Pagel's lambda: {l['lambda']:.4f}    "
                    f"log-likelihood: {l['log_likelihood']:.4f}    "
                    f"p-value: {l['p_value']:.4f}"
                )
        except BrokenPipeError:
            pass
