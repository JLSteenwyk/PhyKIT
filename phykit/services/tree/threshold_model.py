from __future__ import annotations

import math
import sys

from .base import Tree
from ...errors import PhykitUserError


class _LazyNumpy:
    def __getattr__(self, name):
        import numpy as _np

        return getattr(_np, name)


np = _LazyNumpy()


_NDTR = None
_NDTRI = None
_CHO_FACTOR = None
_CHO_SOLVE = None
_INF = math.inf
_NEG_INF = -math.inf
_LOG_2PI = math.log(2.0 * math.pi)
_INV_SQRT2 = 1.0 / math.sqrt(2.0)
_TEXT_REPORT_TEMPLATE = (
    "Trait 1: %s (%s%s)\n"
    "Trait 2: %s (%s%s)\n"
    "MCMC: %d generations, sampled every %d, burn-in %d%%\n"
    "---\n"
    "Posterior correlation (r): %.4f (95%% HPD: %.3f, %.3f)\n"
    "Posterior sigma2_1: %.4f (95%% HPD: %.3f, %.3f)\n"
    "Posterior sigma2_2: %.4f (95%% HPD: %.3f, %.3f)\n"
    "Acceptance rates: r=%.3f, sigma2_1=%.3f, sigma2_2=%.3f, "
    "a1=%.3f, a2=%.3f"
)


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


def ndtr(*args, **kwargs):
    global _NDTR
    if _NDTR is None:
        from scipy.special import ndtr as _ndtr

        _NDTR = _ndtr

    return _NDTR(*args, **kwargs)


def ndtri(*args, **kwargs):
    global _NDTRI
    if _NDTRI is None:
        from scipy.special import ndtri as _ndtri

        _NDTRI = _ndtri

    return _NDTRI(*args, **kwargs)


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


def _normal_cdf(z):
    return 0.5 * math.erfc(-z * _INV_SQRT2)


def _truncnorm_rvs(lower_z, upper_z, loc, scale, rng):
    from scipy.stats import truncnorm

    return truncnorm.rvs(
        lower_z,
        upper_z,
        loc=loc,
        scale=scale,
        random_state=rng,
    )


class ThresholdModel(Tree):
    _FLOAT_TINY = sys.float_info.min
    _ONE_MINUS_EPS = 1.0 - sys.float_info.epsilon

    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"])
        self.trait_data_path = parsed["trait_data_path"]
        self.traits = parsed["traits"]
        self.types = parsed["types"]
        self.ngen = parsed["ngen"]
        self.sample = parsed["sample"]
        self.burnin = parsed["burnin"]
        self.seed = parsed["seed"]
        self.plot_output = parsed["plot_output"]
        self.json_output = parsed["json_output"]
        self.plot_config = parsed["plot_config"]

    def process_args(self, args) -> dict[str, str]:
        from ...helpers.plot_config import PlotConfig

        traits_str = getattr(args, "traits", None)
        types_str = getattr(args, "types", None)
        if not traits_str:
            raise PhykitUserError(
                [
                    "Two trait names must be specified with --traits.",
                    "Usage: --traits trait1,trait2",
                ],
                code=2,
            )
        traits = [t.strip() for t in traits_str.split(",")]
        if len(traits) != 2:
            raise PhykitUserError(
                [
                    f"Expected exactly 2 trait names, got {len(traits)}.",
                    "Usage: --traits trait1,trait2",
                ],
                code=2,
            )
        if not types_str:
            raise PhykitUserError(
                [
                    "Two trait types must be specified with --types.",
                    "Usage: --types discrete,continuous",
                ],
                code=2,
            )
        types = [t.strip() for t in types_str.split(",")]
        if len(types) != 2:
            raise PhykitUserError(
                [
                    f"Expected exactly 2 trait types, got {len(types)}.",
                    "Usage: --types discrete,continuous",
                ],
                code=2,
            )
        for t in types:
            if t not in ("discrete", "continuous"):
                raise PhykitUserError(
                    [
                        f"Invalid trait type '{t}'.",
                        "Each type must be 'discrete' or 'continuous'.",
                    ],
                    code=2,
                )
        return dict(
            tree_file_path=args.tree,
            trait_data_path=args.trait_data,
            traits=traits,
            types=types,
            ngen=getattr(args, "ngen", 100000),
            sample=getattr(args, "sample", 100),
            burnin=getattr(args, "burnin", 0.2),
            seed=getattr(args, "seed", None),
            plot_output=getattr(args, "plot_output", None),
            json_output=getattr(args, "json", False),
            plot_config=PlotConfig.from_args(args),
        )

    def run(self) -> None:
        tree = self.read_tree_file_unmodified()
        self.validate_tree(tree, min_tips=3, require_branch_lengths=True, context="threshold model")

        tree_tips = self.get_tip_names_from_tree(tree)
        trait1_dict, trait2_dict, ordered_names = self._parse_multi_trait_file(
            self.trait_data_path,
            self.traits[0],
            self.traits[1],
            self.types[0],
            self.types[1],
            tree_tips,
        )

        # Prune tree to shared taxa
        tips_to_prune = self._tips_to_prune_for_ordered_names(
            tree_tips,
            ordered_names,
        )
        if tips_to_prune:
            tree = self._fast_copy(tree)
            self._prune_tree_to_taxa(tree, ordered_names)

        C = self._build_vcv_matrix(tree, ordered_names)

        rng = np.random.default_rng(seed=self.seed)

        mcmc_result = self._run_mcmc(
            trait1_dict,
            trait2_dict,
            self.types[0],
            self.types[1],
            ordered_names,
            C,
            self.ngen,
            self.sample,
            self.burnin,
            rng,
        )

        summary = self._summarize_posterior(mcmc_result)

        if self.json_output:
            self._output_json(mcmc_result, summary)
            return

        self._output_text(summary)

        if self.plot_output:
            self._plot_trace(mcmc_result, self.plot_output)

    @staticmethod
    def _parse_multi_trait_file(path, trait1, trait2, type1, type2, tree_tips):
        try:
            with open(path) as f:
                header_line = f.readline()
                if not header_line:
                    raise PhykitUserError(
                        ["Trait file is empty."],
                        code=2,
                    )

                # Parse header
                header = header_line.strip().split("\t")
                if len(header) < 2:
                    raise PhykitUserError(
                        [
                            "Trait file header must have at least 2 tab-separated columns.",
                            f"Got: {header_line.strip()}",
                        ],
                        code=2,
                    )

                # Find column indices (first column is taxon name)
                col_names = header[1:]  # skip taxon column
                try:
                    idx1 = col_names.index(trait1) + 1
                except ValueError:
                    raise PhykitUserError(
                        [
                            f"Trait '{trait1}' not found in header.",
                            f"Available traits: {', '.join(col_names)}",
                        ],
                        code=2,
                    )
                try:
                    idx2 = col_names.index(trait2) + 1
                except ValueError:
                    raise PhykitUserError(
                        [
                            f"Trait '{trait2}' not found in header.",
                            f"Available traits: {', '.join(col_names)}",
                        ],
                        code=2,
                    )

                required_columns = len(header)
                trait1_dict = {}
                trait2_dict = {}
                trait1_discrete_values = set() if type1 == "discrete" else None
                trait2_discrete_values = set() if type2 == "discrete" else None
                parse_float = float
                for line_num, line in enumerate(f, 2):
                    line = line.strip()
                    if not line or line[0] == "#":
                        continue
                    parts = line.split("\t", required_columns)
                    if len(parts) < required_columns:
                        raise PhykitUserError(
                            [
                                f"Line {line_num} has {len(parts)} columns; "
                                f"expected {required_columns}.",
                            ],
                            code=2,
                        )
                    taxon = parts[0]
                    val1_str = parts[idx1]
                    val2_str = parts[idx2]
                    try:
                        val1 = parse_float(val1_str)
                        trait1_dict[taxon] = val1
                        if trait1_discrete_values is not None:
                            trait1_discrete_values.add(val1)
                    except ValueError:
                        raise PhykitUserError(
                            [
                                f"Non-numeric value '{val1_str}' for trait '{trait1}' "
                                f"taxon '{taxon}' on line {line_num}.",
                            ],
                            code=2,
                        )
                    try:
                        val2 = parse_float(val2_str)
                        trait2_dict[taxon] = val2
                        if trait2_discrete_values is not None:
                            trait2_discrete_values.add(val2)
                    except ValueError:
                        raise PhykitUserError(
                            [
                                f"Non-numeric value '{val2_str}' for trait '{trait2}' "
                                f"taxon '{taxon}' on line {line_num}.",
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

        # Validate discrete traits are binary
        for tname, vals in [
            (trait1, trait1_discrete_values),
            (trait2, trait2_discrete_values),
        ]:
            if vals is not None:
                if vals != {0.0, 1.0} and vals != {0.0} and vals != {1.0}:
                    raise PhykitUserError(
                        [
                            f"Discrete trait '{tname}' must have values 0 and 1.",
                            f"Found values: {sorted(vals)}",
                        ],
                        code=2,
                    )

        tree_tip_set = set(tree_tips)
        if (
            len(tree_tip_set) >= 3
            and len(tree_tip_set) == len(trait1_dict)
            and tree_tip_set == trait1_dict.keys()
            and tree_tip_set == trait2_dict.keys()
        ):
            return trait1_dict, trait2_dict, sorted(tree_tip_set)

        trait_taxa_set = set(trait1_dict)
        shared = sorted(tree_tip_set & trait_taxa_set)

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

        # Filter to shared taxa
        t1 = {t: trait1_dict[t] for t in shared}
        t2 = {t: trait2_dict[t] for t in shared}

        return t1, t2, shared

    @staticmethod
    def _build_vcv_matrix(tree, ordered_names: list[str]) -> np.ndarray:
        from .vcv_utils import build_vcv_matrix

        return build_vcv_matrix(tree, ordered_names)

    @staticmethod
    def _prune_tree_to_taxa(tree, taxa_to_keep) -> None:
        taxa_to_keep = set(taxa_to_keep)
        terminal_by_name = Tree._terminal_by_name_fast(tree)
        if terminal_by_name is None:
            for tip in list(tree.get_terminals()):
                if tip.name not in taxa_to_keep:
                    tree.prune(tip)
            return

        targets = [
            tip for name, tip in terminal_by_name.items() if name not in taxa_to_keep
        ]
        if not targets:
            return

        if (
            len(targets) < len(terminal_by_name)
            and all(hasattr(target, "clades") for target in targets)
            and Tree._prune_terminal_objects_batch_standard_tree(
                tree,
                {id(target) for target in targets},
            )
        ):
            return

        for target in targets:
            tree.prune(target)

    @staticmethod
    def _initialize_liabilities(trait_values, trait_type, ordered_names, rng):
        n = len(ordered_names)
        if trait_type == "continuous":
            return np.fromiter(
                (trait_values[name] for name in ordered_names),
                dtype=float,
                count=n,
            )

        states = np.fromiter(
            (int(trait_values[name]) for name in ordered_names),
            dtype=np.int8,
            count=n,
        )
        lower_z = np.where(states == 0, -np.inf, 0.0)
        upper_z = np.where(states == 0, 0.0, np.inf)
        return np.asarray(_truncnorm_rvs(lower_z, upper_z, 0, 1, rng), dtype=float)

    @staticmethod
    def _sample_truncated_normal(mu, sd, lower, upper, rng):
        """Sample one normal value truncated to [lower, upper] by inverse CDF."""
        if lower == _NEG_INF and upper == 0.0:
            return ThresholdModel._sample_truncated_normal_below_zero(mu, sd, rng)
        if lower == 0.0 and upper == _INF:
            return ThresholdModel._sample_truncated_normal_above_zero(mu, sd, rng)

        global _NDTRI
        if _NDTRI is None:
            from scipy.special import ndtri as _ndtri

            _NDTRI = _ndtri
        _ndtri = _NDTRI

        lower_z = (lower - mu) / sd
        upper_z = (upper - mu) / sd
        lower_cdf = 0.0 if lower_z == _NEG_INF else _normal_cdf(lower_z)
        upper_cdf = 1.0 if upper_z == _INF else _normal_cdf(upper_z)

        interval = upper_cdf - lower_cdf
        if interval <= 0 or not math.isfinite(interval):
            return _truncnorm_rvs(lower_z, upper_z, mu, sd, rng)

        u = lower_cdf + rng.random() * interval
        if u < ThresholdModel._FLOAT_TINY:
            u = ThresholdModel._FLOAT_TINY
        elif u > ThresholdModel._ONE_MINUS_EPS:
            u = ThresholdModel._ONE_MINUS_EPS
        return float(mu + sd * _ndtri(u))

    @staticmethod
    def _sample_truncated_normal_below_zero(mu, sd, rng):
        global _NDTRI
        if _NDTRI is None:
            from scipy.special import ndtri as _ndtri

            _NDTRI = _ndtri
        _ndtri = _NDTRI

        upper_z = -mu / sd
        upper_cdf = _normal_cdf(upper_z)
        if upper_cdf <= 0.0 or not math.isfinite(upper_cdf):
            return _truncnorm_rvs(_NEG_INF, upper_z, mu, sd, rng)

        u = rng.random() * upper_cdf
        if u < ThresholdModel._FLOAT_TINY:
            u = ThresholdModel._FLOAT_TINY
        elif u > ThresholdModel._ONE_MINUS_EPS:
            u = ThresholdModel._ONE_MINUS_EPS
        return float(mu + sd * _ndtri(u))

    @staticmethod
    def _sample_truncated_normal_above_zero(mu, sd, rng):
        global _NDTRI
        if _NDTRI is None:
            from scipy.special import ndtri as _ndtri

            _NDTRI = _ndtri
        _ndtri = _NDTRI

        lower_z = -mu / sd
        lower_cdf = _normal_cdf(lower_z)
        interval = 1.0 - lower_cdf
        if interval <= 0.0 or not math.isfinite(interval):
            return _truncnorm_rvs(lower_z, _INF, mu, sd, rng)

        u = lower_cdf + rng.random() * interval
        if u < ThresholdModel._FLOAT_TINY:
            u = ThresholdModel._FLOAT_TINY
        elif u > ThresholdModel._ONE_MINUS_EPS:
            u = ThresholdModel._ONE_MINUS_EPS
        return float(mu + sd * _ndtri(u))

    @staticmethod
    def _log_likelihood_bivariate_bm(
        x1, x2, C_inv, sigma2_1, sigma2_2, r, a1, a2, logdet_C, n
    ):
        stats = ThresholdModel._bivariate_quadratic_stats(x1, x2, C_inv)
        return ThresholdModel._log_likelihood_bivariate_bm_from_stats(
            stats, sigma2_1, sigma2_2, r, a1, a2, logdet_C, n
        )

    @staticmethod
    def _bivariate_quadratic_stats(x1, x2, C_inv):
        """Precompute C_inv quadratic terms independent of BM parameters."""
        C_inv_x1 = C_inv @ x1
        C_inv_x2 = C_inv @ x2
        one_C_one = float(C_inv.sum())
        return ThresholdModel._bivariate_quadratic_stats_from_products(
            x1, x2, C_inv_x1, C_inv_x2, one_C_one
        )

    @staticmethod
    def _bivariate_quadratic_stats_from_products(
        x1, x2, C_inv_x1, C_inv_x2, one_C_one
    ):
        return (
            float(x1 @ C_inv_x1),
            float(x1 @ C_inv_x2),
            float(x2 @ C_inv_x2),
            float(C_inv_x1.sum()),
            float(C_inv_x2.sum()),
            one_C_one,
        )

    @staticmethod
    def _log_likelihood_bivariate_bm_from_stats(
        stats, sigma2_1, sigma2_2, r, a1, a2, logdet_C, n
    ):
        # Guard against non-finite parameters
        if not (math.isfinite(sigma2_1) and math.isfinite(sigma2_2)):
            return -np.inf
        if sigma2_1 <= 0 or sigma2_2 <= 0:
            return -np.inf

        # Bivariate BM: V = Sigma kron C
        # Sigma = [[sigma2_1, r*sqrt(s1*s2)], [r*sqrt(s1*s2), sigma2_2]]
        # Use log-space arithmetic to avoid overflow
        log_s1s2 = math.log(sigma2_1) + math.log(sigma2_2)
        if not math.isfinite(log_s1s2):
            return -np.inf

        cov12 = r * math.exp(0.5 * log_s1s2)
        det_Sigma = sigma2_1 * sigma2_2 * (1.0 - r * r)
        if det_Sigma <= 0 or not math.isfinite(det_Sigma):
            return -np.inf

        # Inverse of Sigma
        S_inv_00 = sigma2_2 / det_Sigma
        S_inv_01 = -cov12 / det_Sigma
        S_inv_11 = sigma2_1 / det_Sigma

        if not (
            math.isfinite(S_inv_00)
            and math.isfinite(S_inv_01)
            and math.isfinite(S_inv_11)
        ):
            return -np.inf

        # Quadratic form: sum_ij S_inv[i,j] * r_i' C_inv r_j
        x1_C_x1, x1_C_x2, x2_C_x2, one_C_x1, one_C_x2, one_C_one = stats
        q1 = x1_C_x1 - 2.0 * a1 * one_C_x1 + a1 * a1 * one_C_one
        q12 = (
            x1_C_x2
            - a2 * one_C_x1
            - a1 * one_C_x2
            + a1 * a2 * one_C_one
        )
        q2 = x2_C_x2 - 2.0 * a2 * one_C_x2 + a2 * a2 * one_C_one

        if not (math.isfinite(q1) and math.isfinite(q12) and math.isfinite(q2)):
            return -np.inf

        quad = S_inv_00 * q1 + 2.0 * S_inv_01 * q12 + S_inv_11 * q2
        if not math.isfinite(quad):
            return -np.inf

        # Log determinant: log|V| = log|Sigma kron C| = n*log|Sigma| + 2*logdet_C
        log_det_Sigma = math.log(det_Sigma)
        log_det_V = n * log_det_Sigma + 2.0 * logdet_C

        ll = -0.5 * (2 * n * _LOG_2PI + log_det_V + quad)

        if not math.isfinite(ll):
            return -np.inf

        return float(ll)

    @staticmethod
    def _sample_liabilities_gibbs(
        liabilities, states, C, C_inv, sigma2, r, other_liabs,
        other_sigma2, a, other_a, rng, gibbs_context=None
    ):
        """Sample liabilities from their full conditional (truncated normal).

        Following phytools::threshBayes, the Gibbs step uses the univariate
        tree-conditional distribution for each tip's liability.  The bivariate
        structure is captured through the Metropolis updates of the correlation
        and variance parameters.

        For tip i given x_{-i}:
          var_cond  = sigma2 / C_inv[i,i]
          mu_cond   = a - (1/C_inv[i,i]) * sum_{j!=i} C_inv[i,j]*(x_j - a)
        Truncated to (-inf, 0) for state 0, or (0, inf) for state 1.
        """
        n = len(liabilities)
        new_liabilities = liabilities.copy()
        residuals = new_liabilities - a
        global _NDTR, _NDTRI
        if _NDTR is None:
            from scipy.special import ndtr as _ndtr_loaded

            _NDTR = _ndtr_loaded
        if _NDTRI is None:
            from scipy.special import ndtri as _ndtri_loaded

            _NDTRI = _ndtri_loaded
        _ndtr = _NDTR
        _ndtri = _NDTRI
        float_tiny = ThresholdModel._FLOAT_TINY
        one_minus_eps = ThresholdModel._ONE_MINUS_EPS
        if gibbs_context is None:
            c_inv_diag = C_inv.diagonal()
            sd_cond_by_tip = None
            valid_tip = None
        else:
            c_inv_diag = gibbs_context["c_inv_diag"]
            sd_cond_by_tip = gibbs_context["sd_cond"]
            valid_tip = gibbs_context["valid"]

        for i in range(n):
            if valid_tip is not None and not valid_tip[i]:
                continue
            c_ii_inv = c_inv_diag[i]
            if c_ii_inv <= 1e-15:
                continue

            # Conditional variance and mean from the tree MVN
            if sd_cond_by_tip is None:
                var_cond = sigma2 / c_ii_inv
                if var_cond <= 0 or not math.isfinite(var_cond):
                    continue
                sd_cond = math.sqrt(var_cond)
                if sd_cond <= 0 or not math.isfinite(sd_cond):
                    continue
            else:
                sd_cond = sd_cond_by_tip[i]

            row_sum = float(C_inv[i, :] @ residuals) - c_ii_inv * residuals[i]
            mu_cond = a - row_sum / c_ii_inv

            state = int(states[i])
            if state == 0:
                # Truncated to (-inf, 0)
                upper_z = -mu_cond / sd_cond
                upper_cdf = _ndtr(upper_z)
                if upper_cdf <= 0.0 or not math.isfinite(upper_cdf):
                    new_value = _truncnorm_rvs(
                        _NEG_INF, upper_z, mu_cond, sd_cond, rng
                    )
                else:
                    u = rng.random() * upper_cdf
                    if u < float_tiny:
                        u = float_tiny
                    elif u > one_minus_eps:
                        u = one_minus_eps
                    new_value = float(mu_cond + sd_cond * _ndtri(u))
            else:
                # Truncated to (0, inf)
                lower_z = -mu_cond / sd_cond
                lower_cdf = _ndtr(lower_z)
                interval = 1.0 - lower_cdf
                if interval <= 0.0 or not math.isfinite(interval):
                    new_value = _truncnorm_rvs(
                        lower_z, _INF, mu_cond, sd_cond, rng
                    )
                else:
                    u = lower_cdf + rng.random() * interval
                    if u < float_tiny:
                        u = float_tiny
                    elif u > one_minus_eps:
                        u = one_minus_eps
                    new_value = float(mu_cond + sd_cond * _ndtri(u))
            new_liabilities[i] = new_value
            residuals[i] = new_value - a

        return new_liabilities

    @staticmethod
    def _prepare_gibbs_context(C_inv, sigma2):
        c_inv_diag = C_inv.diagonal().copy()
        valid = (c_inv_diag > 1e-15) & np.isfinite(c_inv_diag)
        sd_cond = np.zeros_like(c_inv_diag, dtype=float)
        if sigma2 > 0 and np.isfinite(sigma2):
            sd_cond[valid] = np.sqrt(sigma2 / c_inv_diag[valid])
            valid &= np.isfinite(sd_cond) & (sd_cond > 0)
        else:
            valid[:] = False
        return {
            "c_inv_diag": c_inv_diag,
            "sd_cond": sd_cond,
            "valid": valid,
        }

    @staticmethod
    def _vcv_inverse_and_logdet(C):
        try:
            factor = cho_factor(C, lower=True, check_finite=False)
            C_inv = cho_solve(
                factor,
                np.eye(C.shape[0], dtype=C.dtype),
                check_finite=False,
            )
            logdet_C = 2.0 * float(np.log(np.diagonal(factor[0])).sum())
            return C_inv, logdet_C
        except (np.linalg.LinAlgError, FloatingPointError, ValueError):
            C_inv = np.linalg.inv(C)
            sign, logdet_C = np.linalg.slogdet(C)
            if sign <= 0:
                raise PhykitUserError(
                    ["VCV matrix is not positive definite."],
                    code=2,
                )
            return C_inv, logdet_C

    @staticmethod
    def _run_mcmc(
        trait1_dict, trait2_dict, type1, type2,
        ordered_names, C, ngen, sample, burnin_frac, rng
    ):
        n = len(ordered_names)
        C_inv, logdet_C = ThresholdModel._vcv_inverse_and_logdet(C)

        # Initialize liabilities
        x1 = ThresholdModel._initialize_liabilities(
            trait1_dict, type1, ordered_names, rng
        )
        x2 = ThresholdModel._initialize_liabilities(
            trait2_dict, type2, ordered_names, rng
        )

        # States for discrete traits (for Gibbs sampling)
        states1 = None
        states2 = None
        if type1 == "discrete":
            states1 = np.array([int(trait1_dict[name]) for name in ordered_names])
        if type2 == "discrete":
            states2 = np.array([int(trait2_dict[name]) for name in ordered_names])

        # For discrete traits, sigma2 is fixed at 1 and a is fixed at 0
        # for identifiability (threshold=0, sigma2=1 pins the liability scale).
        # Only continuous traits have their sigma2 and a estimated.
        mean_diag_C = float(np.trace(C) / n)
        if mean_diag_C <= 0:
            mean_diag_C = 1.0

        # Discrete: sigma2=1, a=0 (fixed); Continuous: estimate from data
        if type1 == "discrete":
            sigma2_1 = 1.0
            a1 = 0.0
        else:
            var_x1 = float(x1.var())
            var_x1 = var_x1 if var_x1 > 0 else 1.0
            sigma2_1 = var_x1 / mean_diag_C
            sigma2_1 = max(sigma2_1, 1e-10)
            a1 = float(x1.mean())

        if type2 == "discrete":
            sigma2_2 = 1.0
            a2 = 0.0
        else:
            var_x2 = float(x2.var())
            var_x2 = var_x2 if var_x2 > 0 else 1.0
            sigma2_2 = var_x2 / mean_diag_C
            sigma2_2 = max(sigma2_2, 1e-10)
            a2 = float(x2.mean())

        r = 0.0
        log_sigma2_1 = math.log(sigma2_1)
        log_sigma2_2 = math.log(sigma2_2)

        # Proposal variances (only used for estimated parameters)
        prop_var = {
            "log_sigma2_1": 0.2,
            "log_sigma2_2": 0.2,
            "r": 0.1,
            "a1": 0.5,
            "a2": 0.5,
        }
        prop_sd = {key: math.sqrt(value) for key, value in prop_var.items()}

        # Which parameters to update via MH
        update_sigma2_1 = (type1 == "continuous")
        update_sigma2_2 = (type2 == "continuous")
        update_a1 = (type1 == "continuous")
        update_a2 = (type2 == "continuous")
        gibbs_context1 = (
            ThresholdModel._prepare_gibbs_context(C_inv, sigma2_1)
            if type1 == "discrete"
            else None
        )
        gibbs_context2 = (
            ThresholdModel._prepare_gibbs_context(C_inv, sigma2_2)
            if type2 == "discrete"
            else None
        )

        # Acceptance counters
        accept = {k: 0 for k in prop_var}
        attempt = {k: 0 for k in prop_var}

        burnin = int(ngen * burnin_frac)
        n_samples = (ngen - burnin) // sample
        if n_samples < 1:
            n_samples = 1

        # Storage
        samples_r = []
        samples_sigma2_1 = []
        samples_sigma2_2 = []
        samples_a1 = []
        samples_a2 = []

        one_C_one = float(C_inv.sum())
        C_inv_x1 = C_inv @ x1
        C_inv_x2 = C_inv @ x2
        quad_stats = ThresholdModel._bivariate_quadratic_stats_from_products(
            x1, x2, C_inv_x1, C_inv_x2, one_C_one
        )
        current_ll = ThresholdModel._log_likelihood_bivariate_bm_from_stats(
            quad_stats, sigma2_1, sigma2_2, r, a1, a2, logdet_C, n
        )

        # Adaptive tuning interval
        tune_interval = max(100, sample)

        for gen in range(ngen):
            # Gibbs step: sample liabilities for discrete traits
            if type1 == "discrete":
                x1 = ThresholdModel._sample_liabilities_gibbs(
                    x1, states1, C, C_inv, sigma2_1, r,
                    x2, sigma2_2, a1, a2, rng, gibbs_context1
                )
                C_inv_x1 = C_inv @ x1
                quad_stats = ThresholdModel._bivariate_quadratic_stats_from_products(
                    x1, x2, C_inv_x1, C_inv_x2, one_C_one
                )
                current_ll = ThresholdModel._log_likelihood_bivariate_bm_from_stats(
                    quad_stats, sigma2_1, sigma2_2, r, a1, a2, logdet_C, n
                )

            if type2 == "discrete":
                x2 = ThresholdModel._sample_liabilities_gibbs(
                    x2, states2, C, C_inv, sigma2_2, r,
                    x1, sigma2_1, a2, a1, rng, gibbs_context2
                )
                C_inv_x2 = C_inv @ x2
                quad_stats = ThresholdModel._bivariate_quadratic_stats_from_products(
                    x1, x2, C_inv_x1, C_inv_x2, one_C_one
                )
                current_ll = ThresholdModel._log_likelihood_bivariate_bm_from_stats(
                    quad_stats, sigma2_1, sigma2_2, r, a1, a2, logdet_C, n
                )

            # Metropolis step: sigma2_1 (only for continuous traits)
            if update_sigma2_1:
                attempt["log_sigma2_1"] += 1
                log_s1_prop = log_sigma2_1 + rng.normal(0, prop_sd["log_sigma2_1"])
                s1_prop = math.exp(log_s1_prop)
                if math.isfinite(s1_prop) and s1_prop > 0:
                    prop_ll = ThresholdModel._log_likelihood_bivariate_bm_from_stats(
                        quad_stats, s1_prop, sigma2_2, r, a1, a2, logdet_C, n
                    )
                    log_alpha = prop_ll - current_ll + log_s1_prop - log_sigma2_1
                    if math.isfinite(log_alpha) and math.log(rng.random()) < log_alpha:
                        sigma2_1 = s1_prop
                        log_sigma2_1 = log_s1_prop
                        current_ll = prop_ll
                        accept["log_sigma2_1"] += 1

            # Metropolis step: sigma2_2 (only for continuous traits)
            if update_sigma2_2:
                attempt["log_sigma2_2"] += 1
                log_s2_prop = log_sigma2_2 + rng.normal(0, prop_sd["log_sigma2_2"])
                s2_prop = math.exp(log_s2_prop)
                if math.isfinite(s2_prop) and s2_prop > 0:
                    prop_ll = ThresholdModel._log_likelihood_bivariate_bm_from_stats(
                        quad_stats, sigma2_1, s2_prop, r, a1, a2, logdet_C, n
                    )
                    log_alpha = prop_ll - current_ll + log_s2_prop - log_sigma2_2
                    if math.isfinite(log_alpha) and math.log(rng.random()) < log_alpha:
                        sigma2_2 = s2_prop
                        log_sigma2_2 = log_s2_prop
                        current_ll = prop_ll
                        accept["log_sigma2_2"] += 1

            # Metropolis step: r (normal proposal with reflection on [-1,1])
            attempt["r"] += 1
            r_prop = r + rng.normal(0, prop_sd["r"])
            # Reflect into [-1, 1]
            while r_prop < -1 or r_prop > 1:
                if r_prop < -1:
                    r_prop = -2 - r_prop
                if r_prop > 1:
                    r_prop = 2 - r_prop
            prop_ll = ThresholdModel._log_likelihood_bivariate_bm_from_stats(
                quad_stats, sigma2_1, sigma2_2, r_prop, a1, a2, logdet_C, n
            )
            log_alpha = prop_ll - current_ll
            if math.isfinite(log_alpha) and math.log(rng.random()) < log_alpha:
                r = r_prop
                current_ll = prop_ll
                accept["r"] += 1

            # Metropolis step: a1 (only for continuous traits)
            if update_a1:
                attempt["a1"] += 1
                a1_prop = a1 + rng.normal(0, prop_sd["a1"])
                prop_ll = ThresholdModel._log_likelihood_bivariate_bm_from_stats(
                    quad_stats, sigma2_1, sigma2_2, r, a1_prop, a2, logdet_C, n
                )
                log_alpha = prop_ll - current_ll
                if math.isfinite(log_alpha) and math.log(rng.random()) < log_alpha:
                    a1 = a1_prop
                    current_ll = prop_ll
                    accept["a1"] += 1

            # Metropolis step: a2 (only for continuous traits)
            if update_a2:
                attempt["a2"] += 1
                a2_prop = a2 + rng.normal(0, prop_sd["a2"])
                prop_ll = ThresholdModel._log_likelihood_bivariate_bm_from_stats(
                    quad_stats, sigma2_1, sigma2_2, r, a1, a2_prop, logdet_C, n
                )
                log_alpha = prop_ll - current_ll
                if math.isfinite(log_alpha) and math.log(rng.random()) < log_alpha:
                    a2 = a2_prop
                    current_ll = prop_ll
                    accept["a2"] += 1

            # Adaptive tuning during burnin
            if gen < burnin and gen > 0 and gen % tune_interval == 0:
                for key in prop_var:
                    if attempt[key] > 0:
                        rate = accept[key] / attempt[key]
                        if rate < 0.15:
                            prop_var[key] *= 0.5
                        elif rate > 0.35:
                            prop_var[key] *= 1.5
                        prop_sd[key] = math.sqrt(prop_var[key])
                        accept[key] = 0
                        attempt[key] = 0

            # Store samples after burnin
            if gen >= burnin and (gen - burnin) % sample == 0:
                samples_r.append(r)
                samples_sigma2_1.append(sigma2_1)
                samples_sigma2_2.append(sigma2_2)
                samples_a1.append(a1)
                samples_a2.append(a2)

        # Compute final acceptance rates over full chain
        final_accept = {}
        for key in prop_var:
            total = accept[key]
            total_att = attempt[key]
            if total_att > 0:
                final_accept[key] = total / total_att
            else:
                final_accept[key] = 0.0

        return {
            "r": np.array(samples_r),
            "sigma2_1": np.array(samples_sigma2_1),
            "sigma2_2": np.array(samples_sigma2_2),
            "a1": np.array(samples_a1),
            "a2": np.array(samples_a2),
            "acceptance_rates": final_accept,
            "ngen": ngen,
            "sample": sample,
            "burnin": burnin,
            "trait1": None,  # filled in by caller if needed
            "trait2": None,
        }

    @staticmethod
    def _compute_hpd(samples, credible_mass=0.95):
        sorted_samples = np.sort(samples)
        return ThresholdModel._compute_hpd_from_sorted(
            sorted_samples, credible_mass
        )

    @staticmethod
    def _compute_hpd_from_sorted(sorted_samples, credible_mass=0.95):
        n = len(sorted_samples)
        interval_size = int(np.ceil(credible_mass * n))
        if interval_size >= n:
            return float(sorted_samples[0]), float(sorted_samples[-1])

        widths = sorted_samples[interval_size:] - sorted_samples[:n - interval_size]
        best = int(np.argmin(widths))
        return float(sorted_samples[best]), float(sorted_samples[best + interval_size])

    @staticmethod
    def _median_from_sorted(sorted_samples):
        n = len(sorted_samples)
        mid = n // 2
        if n % 2:
            return float(sorted_samples[mid])
        return float((sorted_samples[mid - 1] + sorted_samples[mid]) / 2.0)

    @staticmethod
    def _sample_mean(samples):
        if hasattr(samples, "mean") and len(samples) < 10_000:
            return float(samples.mean())
        return float(np.mean(samples))

    @staticmethod
    def _summarize_posterior(mcmc_result):
        summary = {}
        for param in ("r", "sigma2_1", "sigma2_2", "a1", "a2"):
            samples = mcmc_result[param]
            if len(samples) == 0:
                summary[param] = {
                    "mean": float("nan"),
                    "median": float("nan"),
                    "hpd_lower": float("nan"),
                    "hpd_upper": float("nan"),
                }
                continue
            sorted_samples = np.sort(samples)
            hpd_lo, hpd_hi = ThresholdModel._compute_hpd_from_sorted(
                sorted_samples
            )
            summary[param] = {
                "mean": ThresholdModel._sample_mean(samples),
                "median": ThresholdModel._median_from_sorted(sorted_samples),
                "hpd_lower": hpd_lo,
                "hpd_upper": hpd_hi,
            }
        summary["acceptance_rates"] = mcmc_result["acceptance_rates"]
        return summary

    def _output_text(self, summary):
        try:
            r_s = summary["r"]
            s1 = summary["sigma2_1"]
            s2 = summary["sigma2_2"]
            rates = summary["acceptance_rates"]
            trait1_type = self.types[0]
            trait2_type = self.types[1]
            trait1_states = (
                ", 2 states: 0, 1" if trait1_type == "discrete" else ""
            )
            trait2_states = (
                ", 2 states: 0, 1" if trait2_type == "discrete" else ""
            )
            print(
                _TEXT_REPORT_TEMPLATE % (
                    self.traits[0],
                    trait1_type,
                    trait1_states,
                    self.traits[1],
                    trait2_type,
                    trait2_states,
                    self.ngen,
                    self.sample,
                    int(self.burnin * 100),
                    r_s["mean"],
                    r_s["hpd_lower"],
                    r_s["hpd_upper"],
                    s1["mean"],
                    s1["hpd_lower"],
                    s1["hpd_upper"],
                    s2["mean"],
                    s2["hpd_lower"],
                    s2["hpd_upper"],
                    rates.get("r", 0),
                    rates.get("log_sigma2_1", 0),
                    rates.get("log_sigma2_2", 0),
                    rates.get("a1", 0),
                    rates.get("a2", 0),
                )
            )
        except BrokenPipeError:
            pass

    def _output_json(self, mcmc_result, summary):
        result = {
            "metadata": {
                "trait1": self.traits[0],
                "trait2": self.traits[1],
                "type1": self.types[0],
                "type2": self.types[1],
                "ngen": self.ngen,
                "sample": self.sample,
                "burnin_fraction": self.burnin,
            },
            "summary": {
                param: {
                    "mean": summary[param]["mean"],
                    "median": summary[param]["median"],
                    "hpd_95_lower": summary[param]["hpd_lower"],
                    "hpd_95_upper": summary[param]["hpd_upper"],
                }
                for param in ("r", "sigma2_1", "sigma2_2", "a1", "a2")
            },
            "acceptance_rates": summary["acceptance_rates"],
            "posterior_samples": {
                "r": mcmc_result["r"].tolist(),
                "sigma2_1": mcmc_result["sigma2_1"].tolist(),
                "sigma2_2": mcmc_result["sigma2_2"].tolist(),
                "a1": mcmc_result["a1"].tolist(),
                "a2": mcmc_result["a2"].tolist(),
            },
        }
        print_json(result)

    def _plot_trace(self, mcmc_result, output_path):
        """Generate a 3-row x 2-column figure.

        Left column: MCMC trace plots (sample index vs value).
        Right column: posterior density histograms with 95% HPD shading.

        Parameters plotted: r, sigma2_1, sigma2_2.  For discrete traits
        whose sigma2 is fixed at 1, the trace is flat — the histogram
        still renders but is trivial.
        """
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        params = [
            ("r", "Correlation (r)"),
            ("sigma2_1", r"$\sigma^2_1$"),
            ("sigma2_2", r"$\sigma^2_2$"),
        ]

        config = self.plot_config
        config.resolve(n_rows=None, n_cols=None)
        fig, axes = plt.subplots(3, 2, figsize=(config.fig_width, config.fig_height))

        for row, (param, label) in enumerate(params):
            samples = mcmc_result[param]
            mean_val = ThresholdModel._sample_mean(samples)
            hpd_lo, hpd_hi = ThresholdModel._compute_hpd(samples)

            # --- Left: trace plot ---
            ax_trace = axes[row, 0]
            ax_trace.plot(samples, linewidth=0.4, color="black", alpha=0.6)
            ax_trace.axhline(
                y=mean_val, color="red", linestyle="--",
                linewidth=1, alpha=0.8,
            )
            ax_trace.set_ylabel(label, fontsize=11)
            if row == 2:
                ax_trace.set_xlabel("Sample", fontsize=11)

            # --- Right: posterior density histogram ---
            ax_hist = axes[row, 1]
            unique_vals = np.unique(samples)
            if len(unique_vals) <= 1:
                # Fixed parameter (e.g. discrete sigma2 = 1)
                ax_hist.axvline(
                    x=unique_vals[0], color="black", linewidth=2,
                )
                ax_hist.set_xlim(unique_vals[0] - 0.5, unique_vals[0] + 0.5)
                ax_hist.text(
                    0.5, 0.5, "fixed", transform=ax_hist.transAxes,
                    ha="center", va="center", fontsize=11, color="gray",
                )
            else:
                n_bins = min(50, max(20, len(samples) // 20))
                ax_hist.hist(
                    samples, bins=n_bins, density=True,
                    color="steelblue", alpha=0.7, edgecolor="white",
                    linewidth=0.3,
                )
                # Shade 95% HPD region
                ax_hist.axvspan(
                    hpd_lo, hpd_hi, alpha=0.2, color="red",
                    label="95% HPD",
                )
                ax_hist.axvline(
                    x=mean_val, color="red", linestyle="--",
                    linewidth=1, alpha=0.8,
                )
            ax_hist.set_ylabel("Density", fontsize=10)
            if row == 0:
                ax_hist.legend(fontsize=9, loc="upper right")

        if config.show_title:
            axes[0, 0].set_title(config.title or "Trace", fontsize=config.title_fontsize)
            axes[0, 1].set_title(config.title or "Posterior", fontsize=config.title_fontsize)

        fig.tight_layout()
        fig.savefig(output_path, dpi=config.dpi, bbox_inches="tight")
        plt.close(fig)
