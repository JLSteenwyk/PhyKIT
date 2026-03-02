import sys
from typing import Dict, List, Tuple

import numpy as np
from scipy.optimize import minimize
from scipy.stats import norm as norm_dist

from .base import Tree
from ...helpers.json_output import print_json
from ...errors import PhykitUserError


class PhylogeneticGLM(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"])
        self.trait_data_path = parsed["trait_data_path"]
        self.response = parsed["response"]
        self.predictors = parsed["predictors"]
        self.family = parsed["family"]
        self.method = parsed["method"]
        self.btol = parsed["btol"]
        self.log_alpha_bound = parsed["log_alpha_bound"]
        self.json_output = parsed["json_output"]
        self.gene_trees_path = parsed["gene_trees_path"]

    def run(self) -> None:
        tree = self.read_tree_file()
        self._validate_tree(tree)

        tree_tips = self.get_tip_names_from_tree(tree)
        trait_names, traits = self._parse_multi_trait_file(
            self.trait_data_path, tree_tips
        )

        # Validate column names
        all_columns = [self.response] + list(self.predictors)
        for col in all_columns:
            if col not in trait_names:
                raise PhykitUserError(
                    [
                        f"Column '{col}' not found in trait file.",
                        f"Available columns: {', '.join(trait_names)}",
                    ],
                    code=2,
                )

        if self.response in self.predictors:
            raise PhykitUserError(
                [
                    f"Response variable '{self.response}' must not also be a predictor.",
                ],
                code=2,
            )

        ordered_names = sorted(traits.keys())

        # Pre-compute discordance VCV if gene trees provided
        vcv_meta = None
        if self.gene_trees_path:
            from .vcv_utils import build_discordance_vcv, parse_gene_trees
            gene_trees = parse_gene_trees(self.gene_trees_path)
            _, vcv_meta = build_discordance_vcv(tree, gene_trees, ordered_names)
            shared = vcv_meta["shared_taxa"]
            if set(shared) != set(ordered_names):
                traits = {k: traits[k] for k in shared}
                ordered_names = shared

        n = len(ordered_names)
        k = len(self.predictors)

        if n <= k + 1:
            raise PhykitUserError(
                [
                    f"Insufficient observations: n={n} but need at least {k + 2} "
                    f"(more than k+1={k + 1} parameters).",
                ],
                code=2,
            )

        # Build response vector y and design matrix X (with intercept)
        resp_idx = trait_names.index(self.response)
        pred_indices = [trait_names.index(p) for p in self.predictors]

        y = np.array([traits[name][resp_idx] for name in ordered_names])
        X_pred = np.array(
            [[traits[name][j] for j in pred_indices] for name in ordered_names]
        )
        X = np.column_stack([np.ones(n), X_pred])

        # Validate response variable
        if self.family == "binomial":
            unique_vals = set(y)
            if not unique_vals.issubset({0.0, 1.0}):
                raise PhykitUserError(
                    [
                        f"For binomial family, response variable '{self.response}' "
                        f"must contain only 0 and 1 values.",
                        f"Found values: {sorted(unique_vals)}",
                    ],
                    code=2,
                )
        elif self.family == "poisson":
            if np.any(y < 0):
                raise PhykitUserError(
                    [
                        f"For poisson family, response variable '{self.response}' "
                        f"must contain non-negative values.",
                    ],
                    code=2,
                )
            if not np.allclose(y, np.round(y)):
                raise PhykitUserError(
                    [
                        f"For poisson family, response variable '{self.response}' "
                        f"must contain integer values.",
                    ],
                    code=2,
                )

        # Pre-compute VCV for Poisson GEE if gene trees provided
        self._precomputed_vcv = None
        if self.gene_trees_path:
            from .vcv_utils import build_discordance_vcv, parse_gene_trees
            gene_trees = parse_gene_trees(self.gene_trees_path)
            self._precomputed_vcv, _ = build_discordance_vcv(
                tree, gene_trees, ordered_names
            )

        if self.family == "binomial":
            result = self._fit_logistic_mple(tree, y, X, ordered_names)
        else:
            result = self._fit_poisson_gee(tree, y, X, ordered_names)

        if vcv_meta is not None:
            result["vcv_metadata"] = vcv_meta

        if self.json_output:
            print_json(result)
        else:
            self._print_text_output(result)

    def process_args(self, args) -> Dict[str, str]:
        family = getattr(args, "family", "binomial")
        method = getattr(args, "method", None)
        if method is None:
            method = "logistic_MPLE" if family == "binomial" else "poisson_GEE"
        return dict(
            tree_file_path=args.tree,
            trait_data_path=args.trait_data,
            response=args.response,
            predictors=args.predictors,
            family=family,
            method=method,
            btol=getattr(args, "btol", 10),
            log_alpha_bound=getattr(args, "log_alpha_bound", 4),
            json_output=getattr(args, "json", False),
            gene_trees_path=getattr(args, "gene_trees", None),
        )

    def _validate_tree(self, tree) -> None:
        tips = list(tree.get_terminals())
        if len(tips) < 3:
            raise PhykitUserError(
                ["Tree must have at least 3 tips for phylogenetic GLM."],
                code=2,
            )
        for clade in tree.find_clades():
            if clade.branch_length is None and clade != tree.root:
                raise PhykitUserError(
                    ["All branches in the tree must have lengths."],
                    code=2,
                )

    def _parse_multi_trait_file(
        self, path: str, tree_tips: List[str]
    ) -> Tuple[List[str], Dict[str, List[float]]]:
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

        data_lines = []
        for line in lines:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            data_lines.append(stripped)

        if len(data_lines) < 2:
            raise PhykitUserError(
                [
                    "Multi-trait file must have a header row and at least one data row.",
                ],
                code=2,
            )

        header_parts = data_lines[0].split("\t")
        n_cols = len(header_parts)
        if n_cols < 2:
            raise PhykitUserError(
                [
                    "Header must have at least 2 columns (taxon + at least 1 trait).",
                ],
                code=2,
            )
        trait_names = header_parts[1:]

        traits = {}
        for line_idx, line in enumerate(data_lines[1:], 2):
            parts = line.split("\t")
            if len(parts) != n_cols:
                raise PhykitUserError(
                    [
                        f"Line {line_idx} has {len(parts)} columns; expected {n_cols}.",
                        f"Each line should have: taxon_name<tab>"
                        f"{'<tab>'.join(['trait'] * len(trait_names))}",
                    ],
                    code=2,
                )
            taxon = parts[0]
            values = []
            for i, val_str in enumerate(parts[1:]):
                try:
                    values.append(float(val_str))
                except ValueError:
                    raise PhykitUserError(
                        [
                            f"Non-numeric trait value '{val_str}' for taxon '{taxon}' "
                            f"(trait '{trait_names[i]}') on line {line_idx}.",
                        ],
                        code=2,
                    )
            traits[taxon] = values

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

        filtered = {taxon: traits[taxon] for taxon in shared}
        return trait_names, filtered

    def _build_vcv_matrix(
        self, tree, ordered_names: List[str]
    ) -> np.ndarray:
        from .vcv_utils import build_vcv_matrix
        return build_vcv_matrix(tree, ordered_names)

    def _signif_code(self, p: float) -> str:
        if p < 0.001:
            return "***"
        elif p < 0.01:
            return "**"
        elif p < 0.05:
            return "*"
        elif p < 0.1:
            return "."
        else:
            return " "

    # ----------------------------------------------------------------
    # Logistic MPLE (Ives & Garland 2010)
    # ----------------------------------------------------------------

    def _make_ultrametric(
        self, tree, ordered_names: List[str]
    ) -> Tuple[np.ndarray, float, float]:
        """Compute ultrametric correction D, Tmax, and mean tip height.

        D_i = max(root-to-tip) - root_to_tip[i]
        Tmax = max root-to-tip distance
        mean_height = mean root-to-tip distance (used for initial alpha)
        """
        root_to_tip = {}
        for name in ordered_names:
            root_to_tip[name] = tree.distance(tree.root, name)

        heights = list(root_to_tip.values())
        max_dist = max(heights)
        mean_height = float(np.mean(heights))
        D = np.array([max_dist - root_to_tip[name] for name in ordered_names])
        Tmax = max_dist
        return D, Tmax, mean_height

    def _logistic_starting_values(
        self, y: np.ndarray, X: np.ndarray, btol: float
    ) -> np.ndarray:
        """Standard logistic regression via IRLS as starting values."""
        n, p = X.shape
        beta = np.zeros(p)

        for _ in range(50):
            eta = X @ beta
            eta = np.clip(eta, -btol, btol)
            mu = 1.0 / (1.0 + np.exp(-eta))
            mu = np.clip(mu, 1e-6, 1 - 1e-6)
            W = np.diag(mu * (1 - mu))
            z = eta + (y - mu) / (mu * (1 - mu))
            try:
                beta_new = np.linalg.solve(X.T @ W @ X, X.T @ W @ z)
            except np.linalg.LinAlgError:
                break
            if np.max(np.abs(beta_new - beta)) < 1e-8:
                beta = beta_new
                break
            beta = beta_new

        return beta

    @staticmethod
    def _bernoulli_log_likelihood(
        y: np.ndarray, mu: np.ndarray,
    ) -> float:
        """Standard Bernoulli log-likelihood: sum(y*log(mu) + (1-y)*log(1-mu)).

        In phylolm's logistic_MPLE, the phylogenetic correction enters
        only through the Firth penalty, not through the likelihood itself.
        """
        mu_safe = np.clip(mu, 1e-10, 1 - 1e-10)
        return float(np.sum(
            y * np.log(mu_safe) + (1 - y) * np.log(1 - mu_safe)
        ))

    def _pruning_log_likelihood(
        self, tree, y: np.ndarray, mu: np.ndarray,
        alpha: float, ordered_names: List[str], dk: int,
    ) -> float:
        """Log-likelihood via 2-state CTMC pruning (Ives & Garland 2010).

        Implements the same algorithm as phylolm's logistreglikelihood.c:
        - Two-state Markov chain with rates q01=meanp*alpha, q10=meanq*alpha
        - Stationary distribution: (meanq, meanp) where meanp = mean(mu)
        - Tips have covariate-adjusted emission probabilities when dk > 1
        """
        name_to_idx = {name: i for i, name in enumerate(ordered_names)}

        meanp = float(np.mean(mu))
        meanq = 1.0 - meanp

        NEG_INF = float("-inf")

        llk0 = {}  # log P(data below node | state=0)
        llk1 = {}  # log P(data below node | state=1)

        def _logsumexp(a, b):
            if a == NEG_INF:
                return b
            if b == NEG_INF:
                return a
            mx = max(a, b)
            return mx + np.log(np.exp(a - mx) + np.exp(b - mx))

        def _postorder(node):
            nid = id(node)
            if not node.clades:
                # Tip node
                i = name_to_idx[node.name]
                yi = int(y[i])
                mui = float(mu[i])

                # Base: deterministic observation
                if yi == 1:
                    llk0[nid] = NEG_INF
                    llk1[nid] = 0.0
                else:
                    llk0[nid] = 0.0
                    llk1[nid] = NEG_INF

                # Covariate-adjusted emission probabilities (dk > 1)
                if dk > 1:
                    if mui < meanp:
                        et = mui / meanp
                        llk1[nid] = (
                            np.log(max(et, 1e-300)) if yi == 1
                            else np.log(max(1.0 - et, 1e-300))
                        )
                    else:
                        et = (1.0 - mui) / meanq
                        llk0[nid] = (
                            np.log(max(1.0 - et, 1e-300)) if yi == 1
                            else np.log(max(et, 1e-300))
                        )
            else:
                # Internal node
                llk0[nid] = 0.0
                llk1[nid] = 0.0

                for child in node.clades:
                    _postorder(child)
                    cid = id(child)
                    bl = child.branch_length or 0.0

                    e = np.exp(-bl * alpha)
                    pe = meanp * e
                    qe = meanq * e

                    # P(child=0|parent=0) = meanq + meanp*exp(-alpha*t)
                    # P(child=1|parent=0) = meanp - meanp*exp(-alpha*t)
                    # P(child=0|parent=1) = meanq - meanq*exp(-alpha*t)
                    # P(child=1|parent=1) = meanp + meanq*exp(-alpha*t)
                    c0 = _logsumexp(
                        np.log(max(meanq + pe, 1e-300)) + llk0[cid],
                        np.log(max(meanp - pe, 1e-300)) + llk1[cid],
                    )
                    c1 = _logsumexp(
                        np.log(max(meanq - qe, 1e-300)) + llk0[cid],
                        np.log(max(meanp + qe, 1e-300)) + llk1[cid],
                    )

                    llk0[nid] += c0
                    llk1[nid] += c1

        _postorder(tree.root)
        rid = id(tree.root)

        return _logsumexp(
            np.log(meanq) + llk0[rid], np.log(meanp) + llk1[rid]
        )

    def _three_point_compute(
        self, tree, ordered_names: List[str],
        alpha: float, mu: np.ndarray, D: np.ndarray, Tmax: float,
        X_mat: np.ndarray,
    ) -> Tuple[float, np.ndarray]:
        """Compute log|V| and X'V^{-1}X via postorder tree traversal.

        Matches R's three.point.compute: O(n) algorithm that computes
        V^{-1} products without forming V explicitly.

        Uses Ives & Garland 2010 branch length transformation:
        - distFromRoot[node] = exp(-2*alpha*(Tmax - h_node))
        - Internal edges: el = distFromRoot[desc] - distFromRoot[anc]
        - External edges: el = above + d2 (mu-dependent tip variance)
        - root_edge = exp(-2*alpha*Tmax)

        X_mat should already be scaled: X_mat = diag(mu*(1-mu)/dia) @ X

        Returns (logdet_V, X'V^{-1}X).
        """
        n = len(ordered_names)
        dk = X_mat.shape[1]
        name_to_idx = {name: i for i, name in enumerate(ordered_names)}

        meanp = float(np.mean(mu))
        meanq = 1.0 - meanp

        # Compute node heights and distFromRoot
        node_height = {}
        root = tree.root

        def _set_heights(node, h):
            node_height[id(node)] = h
            for child in node.clades:
                _set_heights(child, h + (child.branch_length or 0.0))
        _set_heights(root, 0.0)

        distFromRoot = {
            nid: np.exp(-2.0 * alpha * (Tmax - h))
            for nid, h in node_height.items()
        }
        root_edge = np.exp(-2.0 * alpha * Tmax)

        # Compute d2 for each tip
        d2 = np.zeros(n)
        for name in ordered_names:
            i = name_to_idx[name]
            mui = float(mu[i])
            Di = float(D[i])
            if mui < meanp:
                m = mui * np.sqrt(meanq / max(meanp, 1e-300))
            else:
                m = (1.0 - mui) * np.sqrt(meanp / max(meanq, 1e-300))
            m2 = max(m * m, 1e-300)
            d2[i] = np.exp(-2.0 * alpha * Di) * mui * (1.0 - mui) / m2

        # Per-node accumulators for the tree traversal
        logd = {}    # log|V_subtree|
        vec11 = {}   # 1'V_subtree^{-1}1
        X1 = {}      # X'V_subtree^{-1}1 (dk-vector)
        XX = {}      # X'V_subtree^{-1}X (dk x dk matrix)

        def _postorder(node):
            nid = id(node)

            if not node.clades:
                # Tip: will be initialized when parent processes this edge
                return

            # Initialize accumulators for this internal node
            logd[nid] = 0.0
            vec11[nid] = 0.0
            X1[nid] = np.zeros(dk)
            XX[nid] = np.zeros((dk, dk))

            # Process children
            for child in node.clades:
                _postorder(child)
                cid = id(child)
                is_tip = not child.clades

                if is_tip:
                    # External edge: initialize tip then accumulate
                    tip_idx = name_to_idx[child.name]
                    el = (
                        distFromRoot[cid] - distFromRoot[nid] + d2[tip_idx]
                    )
                    el = max(el, 1e-300)
                    xi = X_mat[tip_idx]

                    logd[cid] = np.log(el)
                    vec11[cid] = 1.0 / el
                    X1[cid] = xi / el
                    XX[cid] = np.outer(xi, xi) / el
                else:
                    # Internal edge: update subtree with edge contribution
                    el = distFromRoot[cid] - distFromRoot[nid]
                    el = max(el, 1e-300)

                    denom = 1.0 + el * vec11[cid]
                    ev = el / denom
                    ev2 = 1.0 / denom

                    logd[cid] += np.log(denom)
                    XX[cid] -= ev * np.outer(X1[cid], X1[cid])
                    X1[cid] *= ev2
                    vec11[cid] *= ev2

                # Accumulate child into parent
                logd[nid] += logd[cid]
                vec11[nid] += vec11[cid]
                X1[nid] += X1[cid]
                XX[nid] += XX[cid]

        _postorder(root)

        # Process root edge
        rid = id(root)
        denom = 1.0 + root_edge * vec11[rid]
        ev = root_edge / denom
        logd[rid] += np.log(denom)
        XX[rid] -= ev * np.outer(X1[rid], X1[rid])

        return logd[rid], XX[rid]

    def _compute_dia(
        self, ordered_names: List[str],
        alpha: float, mu: np.ndarray, D: np.ndarray,
    ) -> np.ndarray:
        """Compute diagonal scaling factors for Fisher info.

        dia[i] = m_i * exp(alpha * D_i) where m_i depends on mu_i
        relative to meanp.
        """
        n = len(ordered_names)
        meanp = float(np.mean(mu))
        meanq = 1.0 - meanp

        dia = np.zeros(n)
        for i in range(n):
            mui = float(mu[i])
            Di = float(D[i])
            if mui < meanp:
                m = mui * np.sqrt(meanq / max(meanp, 1e-300))
            else:
                m = (1.0 - mui) * np.sqrt(meanp / max(meanq, 1e-300))
            m2 = max(m * m, 1e-300)
            dia[i] = np.sqrt(m2) * np.exp(alpha * Di)

        return dia

    def _fit_logistic_mple(
        self, tree, y: np.ndarray, X: np.ndarray, ordered_names: List[str]
    ) -> Dict:
        """Fit logistic phylogenetic GLM via MPLE (Ives & Garland 2010).

        Minimizes -(LL + 0.5*log(det(I))) where:
        - LL = pruning log-likelihood for 2-state CTMC on the phylogeny
        - I = Fisher information via three_point_compute (tree-based V^{-1})
        """
        n, p = X.shape
        dk = p
        k = p - 1

        D, Tmax, mean_height = self._make_ultrametric(tree, ordered_names)

        # Starting values (matching R's phylolm)
        beta0 = self._logistic_starting_values(y, X, self.btol)
        eta0 = X @ beta0
        if np.any(np.abs(eta0) >= self.btol):
            beta0 = np.zeros(p)
            # Try log-odds as intermediate fallback
            n1 = np.sum(y == 1)
            n0 = np.sum(y == 0)
            if n1 > 0 and n0 > 0:
                beta0[0] = np.log(n1 / n0)
                if np.any(np.abs(X @ beta0) >= self.btol):
                    beta0[0] = 0.0

        # lL = log(1/alpha); starting alpha = 1/Tmax, so lL = log(Tmax)
        lL_init = np.log(Tmax)

        def objective(params):
            """Negative penalized log-likelihood."""
            beta = params[:p]
            lL = params[p]

            # Bound check (R returns 1e10 for out-of-bounds)
            if abs(lL - np.log(Tmax)) >= self.log_alpha_bound:
                return 1e10

            alpha = np.exp(-lL)

            eta = X @ beta
            if np.any(np.abs(eta) >= self.btol):
                return 1e10

            mu = 1.0 / (1.0 + np.exp(-eta))
            mu = np.clip(mu, 1e-10, 1 - 1e-10)

            # 2-state CTMC pruning log-likelihood
            ll = self._pruning_log_likelihood(
                tree, y, mu, alpha, ordered_names, dk,
            )
            if not np.isfinite(ll):
                return 1e10

            # Fisher info via three_point_compute (no dense V inversion)
            dia = self._compute_dia(ordered_names, alpha, mu, D)
            W = mu * (1 - mu)
            X_scaled = (W / dia)[:, np.newaxis] * X

            try:
                logdet_V, infoM = self._three_point_compute(
                    tree, ordered_names, alpha, mu, D, Tmax, X_scaled,
                )
            except Exception:
                return 1e10

            try:
                if dk == 1:
                    logdet_I = np.log(abs(infoM[0, 0]))
                else:
                    sign, logdet_I = np.linalg.slogdet(infoM)
                    if sign <= 0:
                        return 1e10
                penalty = 0.5 * logdet_I
            except (np.linalg.LinAlgError, ValueError):
                return 1e10

            return -(ll + penalty)

        x0 = np.concatenate([beta0, [lL_init]])

        # R uses unconstrained L-BFGS-B with factr=1e12 (loose tolerance)
        result = minimize(
            objective,
            x0,
            method="L-BFGS-B",
            options={"maxiter": 500, "ftol": 2.2e-4, "gtol": 1e-5},
        )

        beta_hat = result.x[:p]
        alpha_hat = np.exp(-result.x[p])

        # Final quantities
        eta = X @ beta_hat
        mu = 1.0 / (1.0 + np.exp(-eta))
        mu = np.clip(mu, 1e-10, 1 - 1e-10)

        # Reported LL is the unpenalized CTMC pruning LL
        ll = self._pruning_log_likelihood(
            tree, y, mu, alpha_hat, ordered_names, dk,
        )

        # Standard errors from Fisher information
        dia = self._compute_dia(ordered_names, alpha_hat, mu, D)
        W = mu * (1 - mu)
        X_scaled = (W / dia)[:, np.newaxis] * X

        try:
            _, infoM = self._three_point_compute(
                tree, ordered_names, alpha_hat, mu, D, Tmax, X_scaled,
            )
            I_inv = np.linalg.inv(infoM)
            se = np.sqrt(np.abs(np.diag(I_inv)))
        except np.linalg.LinAlgError:
            se = np.full(p, np.nan)

        z_stats = beta_hat / se
        p_values = 2.0 * norm_dist.sf(np.abs(z_stats))

        # AIC: -2*ll + 2*(p+1) where p coefficients + alpha
        aic = -2.0 * ll + 2.0 * (p + 1)

        coef_names = ["(Intercept)"] + list(self.predictors)
        formula = f"{self.response} ~ {' + '.join(self.predictors)}"
        fitted_dict = {ordered_names[i]: float(mu[i]) for i in range(n)}

        return self._format_result(
            family="binomial",
            method="logistic_MPLE",
            coef_names=coef_names,
            beta_hat=beta_hat,
            se=se,
            z_stats=z_stats,
            p_values=p_values,
            ll=ll,
            aic=aic,
            formula=formula,
            n=n,
            k=k,
            ordered_names=ordered_names,
            fitted=fitted_dict,
            alpha=alpha_hat,
            overdispersion=None,
        )

    # ----------------------------------------------------------------
    # Poisson GEE (Paradis & Claude 2002)
    # ----------------------------------------------------------------

    def _poisson_starting_values(
        self, y: np.ndarray, X: np.ndarray
    ) -> np.ndarray:
        """Standard Poisson GLM via IRLS as starting values."""
        n, p = X.shape
        beta = np.zeros(p)
        beta[0] = np.log(max(np.mean(y), 0.1))

        for _ in range(50):
            eta = X @ beta
            eta = np.clip(eta, -20, 20)
            mu = np.exp(eta)
            mu = np.clip(mu, 1e-6, 1e10)
            W = np.diag(mu)
            z = eta + (y - mu) / mu
            try:
                beta_new = np.linalg.solve(X.T @ W @ X, X.T @ W @ z)
            except np.linalg.LinAlgError:
                break
            if np.max(np.abs(beta_new - beta)) < 1e-8:
                beta = beta_new
                break
            beta = beta_new

        return beta

    def _fit_poisson_gee(
        self, tree, y: np.ndarray, X: np.ndarray, ordered_names: List[str]
    ) -> Dict:
        """Fit Poisson phylogenetic GLM via GEE with Fisher scoring.

        GEE for Poisson with log link and phylogenetic correlation R:
          D = diag(mu) X  (derivative of mu w.r.t. beta)
          V = diag(sqrt(mu)) R diag(sqrt(mu))
          Info I = D' V^{-1} D = X' diag(sqrt(mu)) R^{-1} diag(sqrt(mu)) X
          Score = D' V^{-1} (y-mu) = X' diag(sqrt(mu)) R^{-1} (y-mu)/sqrt(mu)
        """
        n, p = X.shape
        k = p - 1

        if getattr(self, "_precomputed_vcv", None) is not None:
            vcv = self._precomputed_vcv
        else:
            vcv = self._build_vcv_matrix(tree, ordered_names)
        d = np.sqrt(np.diag(vcv))
        R = vcv / np.outer(d, d)

        try:
            R_inv = np.linalg.inv(R)
        except np.linalg.LinAlgError:
            raise PhykitUserError(
                [
                    "Singular correlation matrix: cannot invert.",
                    "Check that the tree has valid branch lengths.",
                ],
                code=2,
            )

        beta = self._poisson_starting_values(y, X)

        max_iter = 10000
        tol = 1e-10
        converged = False

        for iteration in range(max_iter):
            eta = X @ beta
            eta = np.clip(eta, -20, 20)
            mu = np.exp(eta)
            mu = np.clip(mu, 1e-10, 1e10)

            sqrt_mu = np.sqrt(mu)
            inv_sqrt_mu = 1.0 / sqrt_mu

            WX = np.diag(sqrt_mu) @ X
            I_mat = WX.T @ R_inv @ WX
            score = WX.T @ R_inv @ (inv_sqrt_mu * (y - mu))

            try:
                delta = np.linalg.solve(I_mat, score)
            except np.linalg.LinAlgError:
                break

            beta = beta + delta

            if np.sum(np.abs(delta)) < tol:
                converged = True
                break

        if not converged:
            print(
                f"Warning: Poisson GEE did not converge after {max_iter} iterations.",
                file=sys.stderr,
            )

        # Final values
        eta = X @ beta
        eta = np.clip(eta, -20, 20)
        mu = np.exp(eta)

        # Overdispersion
        pearson_resid = (y - mu) / np.sqrt(mu)
        phi = float(np.sum(pearson_resid**2)) / (n - p)

        # Log-likelihood (Poisson)
        from scipy.special import gammaln
        ll = float(np.sum(
            y * np.log(np.clip(mu, 1e-300, None)) - mu - gammaln(y + 1)
        ))

        # Covariance: phi * I^{-1}
        sqrt_mu_final = np.sqrt(mu)
        WX_final = np.diag(sqrt_mu_final) @ X
        I_final = WX_final.T @ R_inv @ WX_final

        try:
            I_inv = np.linalg.inv(I_final)
            var_beta = phi * I_inv
            se = np.sqrt(np.abs(np.diag(var_beta)))
        except np.linalg.LinAlgError:
            se = np.full(p, np.nan)

        z_stats = beta / se
        p_values = 2.0 * norm_dist.sf(np.abs(z_stats))

        aic = -2.0 * ll + 2.0 * p

        coef_names = ["(Intercept)"] + list(self.predictors)
        formula = f"{self.response} ~ {' + '.join(self.predictors)}"
        fitted_dict = {ordered_names[i]: float(mu[i]) for i in range(n)}

        return self._format_result(
            family="poisson",
            method="poisson_GEE",
            coef_names=coef_names,
            beta_hat=beta,
            se=se,
            z_stats=z_stats,
            p_values=p_values,
            ll=ll,
            aic=aic,
            formula=formula,
            n=n,
            k=k,
            ordered_names=ordered_names,
            fitted=fitted_dict,
            alpha=None,
            overdispersion=phi,
        )

    # ----------------------------------------------------------------
    # Output formatting
    # ----------------------------------------------------------------

    def _print_text_output(self, result: Dict) -> None:
        family = result["family"]
        method = result["method"]

        if family == "binomial":
            print("Phylogenetic GLM (Logistic MPLE)")
        else:
            print("Phylogenetic GLM (Poisson GEE)")

        print(f"\nFormula: {result['formula']}")
        print(f"Family: {family}, Method: {method}")

        if result.get("alpha") is not None:
            print(f"\nEstimated alpha: {result['alpha']:.4f}")
        if result.get("overdispersion") is not None:
            print(f"\nOverdispersion (phi): {result['overdispersion']:.4f}")

        print("\nCoefficients:")
        print(
            f"{'':20s}{'Estimate':>12s}{'Std.Error':>12s}"
            f"{'z-value':>12s}{'p-value':>12s}"
        )

        coefficients = result["coefficients"]
        for name in coefficients:
            coef = coefficients[name]
            sig = self._signif_code(coef["p_value"])
            print(
                f"{name:20s}{coef['estimate']:12.4f}{coef['std_error']:12.4f}"
                f"{coef['z_value']:12.4f}{coef['p_value']:12.6f}    {sig}"
            )

        print("---")
        print("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1")

        print(
            f"\nLog-likelihood: {result['log_likelihood']:.4f}    "
            f"AIC: {result['aic']:.4f}"
        )
        print(f"Number of observations: {result['n_observations']}")

    def _format_result(
        self,
        *,
        family: str,
        method: str,
        coef_names: List[str],
        beta_hat: np.ndarray,
        se: np.ndarray,
        z_stats: np.ndarray,
        p_values: np.ndarray,
        ll: float,
        aic: float,
        formula: str,
        n: int,
        k: int,
        ordered_names: List[str],
        fitted: Dict,
        alpha: float,
        overdispersion: float,
    ) -> Dict:
        coefficients = {}
        for i, name in enumerate(coef_names):
            coefficients[name] = {
                "estimate": float(beta_hat[i]),
                "std_error": float(se[i]),
                "z_value": float(z_stats[i]),
                "p_value": float(p_values[i]),
            }

        result = {
            "family": family,
            "method": method,
            "formula": formula,
            "coefficients": coefficients,
            "log_likelihood": float(ll),
            "aic": float(aic),
            "n_observations": n,
            "n_predictors": k,
            "fitted_values": fitted,
        }

        if alpha is not None:
            result["alpha"] = float(alpha)
        if overdispersion is not None:
            result["overdispersion"] = float(overdispersion)

        return result
