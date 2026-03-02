import sys
from typing import Dict, List, Tuple

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.stats import chi2

from .base import Tree
from ...helpers.json_output import print_json
from ...errors import PhykitUserError


class PhylogeneticSignal(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"])
        self.trait_data_path = parsed["trait_data_path"]
        self.method = parsed["method"]
        self.permutations = parsed["permutations"]
        self.json_output = parsed["json_output"]
        self.gene_trees_path = parsed["gene_trees_path"]

    def run(self) -> None:
        from .vcv_utils import build_vcv_matrix, build_discordance_vcv, parse_gene_trees

        tree = self.read_tree_file()
        self._validate_tree(tree)

        tree_tips = self.get_tip_names_from_tree(tree)
        traits = self._parse_trait_file(self.trait_data_path, tree_tips)

        ordered_names = sorted(traits.keys())

        if self.gene_trees_path:
            gene_trees = parse_gene_trees(self.gene_trees_path)
            vcv, vcv_meta = build_discordance_vcv(tree, gene_trees, ordered_names)
            # Subset traits to shared taxa if needed
            shared = vcv_meta["shared_taxa"]
            if set(shared) != set(ordered_names):
                traits = {k: traits[k] for k in shared}
                ordered_names = shared
        else:
            vcv = build_vcv_matrix(tree, ordered_names)
            vcv_meta = None

        x = np.array([traits[name] for name in ordered_names])

        r2_phylo = self._compute_r2_phylo(x, vcv)

        if self.method == "blombergs_k":
            result = self._blombergs_k(x, vcv, self.permutations)
            result["r_squared_phylo"] = r2_phylo
            if vcv_meta is not None:
                result["vcv_metadata"] = vcv_meta
            if self.json_output:
                print_json(result)
                return
            print(
                f"{round(result['K'], 4)}\t"
                f"{round(result['p_value'], 4)}\t"
                f"{round(result['r_squared_phylo'], 4)}"
            )
        elif self.method == "lambda":
            max_lam = self._max_lambda(tree)
            result = self._pagels_lambda(x, vcv, max_lam)
            result["r_squared_phylo"] = r2_phylo
            if vcv_meta is not None:
                result["vcv_metadata"] = vcv_meta
            if self.json_output:
                print_json(result)
                return
            print(
                f"{round(result['lambda'], 4)}\t"
                f"{round(result['log_likelihood'], 4)}\t"
                f"{round(result['p_value'], 4)}\t"
                f"{round(result['r_squared_phylo'], 4)}"
            )

    def process_args(self, args) -> Dict[str, str]:
        return dict(
            tree_file_path=args.tree,
            trait_data_path=args.trait_data,
            method=getattr(args, "method", "blombergs_k"),
            permutations=getattr(args, "permutations", 1000),
            json_output=getattr(args, "json", False),
            gene_trees_path=getattr(args, "gene_trees", None),
        )

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

    def _build_vcv_matrix(
        self, tree, ordered_names: List[str]
    ) -> np.ndarray:
        from .vcv_utils import build_vcv_matrix
        return build_vcv_matrix(tree, ordered_names)

    def _log_likelihood(
        self, x: np.ndarray, C: np.ndarray
    ) -> Tuple[float, float]:
        """Concentrated log-likelihood with sigma^2 profiled out.

        Matches phytools::phylosig(method="lambda") internals.
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
        # logdet(sig2 * C) = n * log(sig2) + logdet(C)
        logdet_sig2C = n * np.log(sig2) + logdet_C
        ll = -0.5 * (n * np.log(2 * np.pi) + logdet_sig2C + n)

        return ll, a_hat

    def _compute_r2_phylo(self, x: np.ndarray, vcv: np.ndarray) -> float:
        """Compute R²_phylo = 1 - (σ²_BM / σ²_WN).

        σ²_BM is the MLE of trait variance under Brownian Motion.
        σ²_WN is the MLE of trait variance under White Noise (= sample variance).
        """
        n = len(x)
        C_inv = np.linalg.inv(vcv)
        ones = np.ones(n)
        a_hat = float((ones @ C_inv @ x) / (ones @ C_inv @ ones))
        e = x - a_hat
        sig2_bm = float(e @ C_inv @ e) / n
        sig2_wn = float(np.var(x))  # numpy uses /n by default
        if sig2_wn == 0:
            return float("nan")
        return 1.0 - sig2_bm / sig2_wn

    def _blombergs_k(
        self, x: np.ndarray, vcv: np.ndarray, n_perm: int
    ) -> Dict:
        """Blomberg's K following phytools::phylosig(method="K").

        K = (e'e / e'C_inv'e) / ((trace(C) - n/sum(C_inv)) / (n-1))
        Uses raw sums (not MSE/MSE0 with sample-size normalization).
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

        # Observed ratio (raw sums, matching phytools)
        observed_ratio = float(e @ e) / float(e @ C_inv @ e)

        # Expected ratio under BM (Blomberg et al. 2003)
        expected_ratio = (np.trace(C) - n / sum_C_inv) / (n - 1)

        # K statistic
        K = observed_ratio / expected_ratio

        # Permutation test (matching phytools: shuffle x, recompute K)
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

    def _max_lambda(self, tree) -> float:
        """Compute the upper bound for lambda, matching phytools::maxLambda.

        For ultrametric trees: max_tip_height / max_parent_node_height.
        For non-ultrametric trees: 1.0.
        """
        tips = tree.get_terminals()
        root = tree.root
        tip_heights = [tree.distance(root, tip) for tip in tips]
        max_tip_height = max(tip_heights)
        min_tip_height = min(tip_heights)

        # Check ultrametric: all tips equidistant from root (within tolerance)
        is_ultrametric = (max_tip_height - min_tip_height) / max_tip_height < 1e-6

        if not is_ultrametric:
            return 1.0

        # For ultrametric trees: max_tip_height / max_parent_node_height
        # Parent node height = distance from root to the start of each edge
        max_parent_height = 0.0
        for clade in tree.find_clades(order="level"):
            if clade == root:
                continue
            # Parent height = node height - branch length
            node_height = tree.distance(root, clade)
            parent_height = node_height - (clade.branch_length or 0.0)
            if parent_height > max_parent_height:
                max_parent_height = parent_height

        if max_parent_height == 0.0:
            return 1.0

        return max_tip_height / max_parent_height

    def _pagels_lambda(
        self, x: np.ndarray, vcv: np.ndarray, max_lambda: float = 1.0
    ) -> Dict:
        n = len(x)
        diag_vals = np.diag(vcv).copy()
        niter = 10

        def neg_ll(lam):
            C_lam = vcv * lam
            np.fill_diagonal(C_lam, diag_vals)
            try:
                ll, _ = self._log_likelihood(x, C_lam)
                return -ll
            except (np.linalg.LinAlgError, FloatingPointError, ValueError):
                return 1e10

        # Multi-interval optimization matching phytools (niter=10 sub-intervals)
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
        ll_fitted, _ = self._log_likelihood(x, C_fitted)

        # Log-likelihood at lambda = 0 (no phylogenetic signal)
        C_zero = vcv * 0.0
        np.fill_diagonal(C_zero, diag_vals)
        ll_zero, _ = self._log_likelihood(x, C_zero)

        # Likelihood ratio test (1 df)
        lr_stat = 2.0 * (ll_fitted - ll_zero)
        lr_stat = max(lr_stat, 0.0)
        p_value = float(chi2.sf(lr_stat, df=1))

        return dict(
            **{"lambda": float(lambda_hat)},
            log_likelihood=float(ll_fitted),
            p_value=p_value,
        )
