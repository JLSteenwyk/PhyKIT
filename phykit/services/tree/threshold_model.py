import sys
from typing import Dict, List, Tuple

import numpy as np
from scipy.stats import truncnorm

from .base import Tree
from ...helpers.json_output import print_json
from ...errors import PhykitUserError


class ThresholdModel(Tree):
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

    def process_args(self, args) -> Dict[str, str]:
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
        )

    def run(self) -> None:
        tree = self.read_tree_file()
        self._validate_tree(tree)

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
        shared_set = set(ordered_names)
        tips_to_prune = [t for t in tree_tips if t not in shared_set]
        if tips_to_prune:
            for tip_name in tips_to_prune:
                tips = [t for t in tree.get_terminals() if t.name == tip_name]
                if tips:
                    tree.prune(tips[0])

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

    def _validate_tree(self, tree) -> None:
        tips = list(tree.get_terminals())
        if len(tips) < 3:
            raise PhykitUserError(
                ["Tree must have at least 3 tips for threshold model."],
                code=2,
            )
        for clade in tree.find_clades():
            if clade.branch_length is None and clade != tree.root:
                raise PhykitUserError(
                    ["All branches in the tree must have lengths."],
                    code=2,
                )

    @staticmethod
    def _parse_multi_trait_file(path, trait1, trait2, type1, type2, tree_tips):
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

        if not lines:
            raise PhykitUserError(
                ["Trait file is empty."],
                code=2,
            )

        # Parse header
        header = lines[0].strip().split("\t")
        if len(header) < 2:
            raise PhykitUserError(
                [
                    "Trait file header must have at least 2 tab-separated columns.",
                    f"Got: {lines[0].strip()}",
                ],
                code=2,
            )

        # Find column indices (first column is taxon name)
        col_names = header[1:]  # skip taxon column
        try:
            idx1 = col_names.index(trait1)
        except ValueError:
            raise PhykitUserError(
                [
                    f"Trait '{trait1}' not found in header.",
                    f"Available traits: {', '.join(col_names)}",
                ],
                code=2,
            )
        try:
            idx2 = col_names.index(trait2)
        except ValueError:
            raise PhykitUserError(
                [
                    f"Trait '{trait2}' not found in header.",
                    f"Available traits: {', '.join(col_names)}",
                ],
                code=2,
            )

        trait1_dict = {}
        trait2_dict = {}
        for line_num, line in enumerate(lines[1:], 2):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < len(header):
                raise PhykitUserError(
                    [
                        f"Line {line_num} has {len(parts)} columns; "
                        f"expected {len(header)}.",
                    ],
                    code=2,
                )
            taxon = parts[0]
            val1_str = parts[idx1 + 1]  # +1 for taxon column
            val2_str = parts[idx2 + 1]
            try:
                trait1_dict[taxon] = float(val1_str)
            except ValueError:
                raise PhykitUserError(
                    [
                        f"Non-numeric value '{val1_str}' for trait '{trait1}' "
                        f"taxon '{taxon}' on line {line_num}.",
                    ],
                    code=2,
                )
            try:
                trait2_dict[taxon] = float(val2_str)
            except ValueError:
                raise PhykitUserError(
                    [
                        f"Non-numeric value '{val2_str}' for trait '{trait2}' "
                        f"taxon '{taxon}' on line {line_num}.",
                    ],
                    code=2,
                )

        # Validate discrete traits are binary
        for tname, ttype, tdict in [
            (trait1, type1, trait1_dict),
            (trait2, type2, trait2_dict),
        ]:
            if ttype == "discrete":
                vals = set(tdict.values())
                if vals != {0.0, 1.0} and vals != {0.0} and vals != {1.0}:
                    raise PhykitUserError(
                        [
                            f"Discrete trait '{tname}' must have values 0 and 1.",
                            f"Found values: {sorted(vals)}",
                        ],
                        code=2,
                    )

        tree_tip_set = set(tree_tips)
        trait_taxa_set = set(trait1_dict.keys())
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
    def _build_vcv_matrix(tree, ordered_names: List[str]) -> np.ndarray:
        n = len(ordered_names)
        vcv = np.zeros((n, n))

        root_to_tip = {}
        for name in ordered_names:
            root_to_tip[name] = tree.distance(tree.root, name)

        for i in range(n):
            for j in range(i, n):
                if i == j:
                    vcv[i, j] = root_to_tip[ordered_names[i]]
                else:
                    d_ij = tree.distance(ordered_names[i], ordered_names[j])
                    shared_path = (
                        root_to_tip[ordered_names[i]]
                        + root_to_tip[ordered_names[j]]
                        - d_ij
                    ) / 2.0
                    vcv[i, j] = shared_path
                    vcv[j, i] = shared_path

        return vcv

    @staticmethod
    def _initialize_liabilities(trait_values, trait_type, ordered_names, rng):
        n = len(ordered_names)
        liabilities = np.zeros(n)
        if trait_type == "continuous":
            for i, name in enumerate(ordered_names):
                liabilities[i] = trait_values[name]
        else:
            # discrete: state 0 -> negative, state 1 -> positive
            for i, name in enumerate(ordered_names):
                state = int(trait_values[name])
                if state == 0:
                    # Draw from truncated normal on (-inf, 0)
                    liabilities[i] = truncnorm.rvs(
                        -np.inf, 0, loc=0, scale=1, random_state=rng
                    )
                else:
                    # Draw from truncated normal on (0, inf)
                    liabilities[i] = truncnorm.rvs(
                        0, np.inf, loc=0, scale=1, random_state=rng
                    )
        return liabilities

    @staticmethod
    def _log_likelihood_bivariate_bm(
        x1, x2, C_inv, sigma2_1, sigma2_2, r, a1, a2, logdet_C, n
    ):
        # Guard against non-finite parameters
        if not (np.isfinite(sigma2_1) and np.isfinite(sigma2_2)):
            return -np.inf
        if sigma2_1 <= 0 or sigma2_2 <= 0:
            return -np.inf

        # Bivariate BM: V = Sigma kron C
        # Sigma = [[sigma2_1, r*sqrt(s1*s2)], [r*sqrt(s1*s2), sigma2_2]]
        # Use log-space arithmetic to avoid overflow
        log_s1s2 = np.log(sigma2_1) + np.log(sigma2_2)
        if not np.isfinite(log_s1s2):
            return -np.inf

        cov12 = r * np.exp(0.5 * log_s1s2)
        det_Sigma = sigma2_1 * sigma2_2 * (1.0 - r * r)
        if det_Sigma <= 0 or not np.isfinite(det_Sigma):
            return -np.inf

        # Inverse of Sigma
        S_inv_00 = sigma2_2 / det_Sigma
        S_inv_01 = -cov12 / det_Sigma
        S_inv_11 = sigma2_1 / det_Sigma

        if not (np.isfinite(S_inv_00) and np.isfinite(S_inv_01) and np.isfinite(S_inv_11)):
            return -np.inf

        # Residuals
        r1 = x1 - a1
        r2 = x2 - a2

        # Quadratic form: sum_ij S_inv[i,j] * r_i' C_inv r_j
        q1 = float(r1 @ C_inv @ r1)
        q12 = float(r1 @ C_inv @ r2)
        q2 = float(r2 @ C_inv @ r2)

        if not (np.isfinite(q1) and np.isfinite(q12) and np.isfinite(q2)):
            return -np.inf

        quad = S_inv_00 * q1 + 2.0 * S_inv_01 * q12 + S_inv_11 * q2
        if not np.isfinite(quad):
            return -np.inf

        # Log determinant: log|V| = log|Sigma kron C| = n*log|Sigma| + 2*logdet_C
        log_det_Sigma = np.log(det_Sigma)
        log_det_V = n * log_det_Sigma + 2.0 * logdet_C

        ll = -0.5 * (2 * n * np.log(2 * np.pi) + log_det_V + quad)

        if not np.isfinite(ll):
            return -np.inf

        return float(ll)

    @staticmethod
    def _sample_liabilities_gibbs(
        liabilities, states, C, C_inv, sigma2, r, other_liabs,
        other_sigma2, a, other_a, rng
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

        for i in range(n):
            c_ii_inv = C_inv[i, i]
            if c_ii_inv <= 1e-15:
                continue

            # Conditional variance and mean from the tree MVN
            var_cond = sigma2 / c_ii_inv
            if var_cond <= 0 or not np.isfinite(var_cond):
                continue

            residuals = new_liabilities - a
            row_sum = float(C_inv[i, :] @ residuals) - c_ii_inv * residuals[i]
            mu_cond = a - row_sum / c_ii_inv

            sd_cond = np.sqrt(var_cond)
            if sd_cond <= 0 or not np.isfinite(sd_cond):
                continue

            state = int(states[i])
            if state == 0:
                # Truncated to (-inf, 0)
                b_std = (0 - mu_cond) / sd_cond
                new_liabilities[i] = truncnorm.rvs(
                    -np.inf, b_std, loc=mu_cond, scale=sd_cond,
                    random_state=rng
                )
            else:
                # Truncated to (0, inf)
                a_std = (0 - mu_cond) / sd_cond
                new_liabilities[i] = truncnorm.rvs(
                    a_std, np.inf, loc=mu_cond, scale=sd_cond,
                    random_state=rng
                )

        return new_liabilities

    @staticmethod
    def _run_mcmc(
        trait1_dict, trait2_dict, type1, type2,
        ordered_names, C, ngen, sample, burnin_frac, rng
    ):
        n = len(ordered_names)
        C_inv = np.linalg.inv(C)
        sign, logdet_C = np.linalg.slogdet(C)
        if sign <= 0:
            raise PhykitUserError(
                ["VCV matrix is not positive definite."],
                code=2,
            )

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
        mean_diag_C = float(np.mean(np.diag(C)))
        if mean_diag_C <= 0:
            mean_diag_C = 1.0

        # Discrete: sigma2=1, a=0 (fixed); Continuous: estimate from data
        if type1 == "discrete":
            sigma2_1 = 1.0
            a1 = 0.0
        else:
            var_x1 = float(np.var(x1)) if float(np.var(x1)) > 0 else 1.0
            sigma2_1 = var_x1 / mean_diag_C
            sigma2_1 = max(sigma2_1, 1e-10)
            a1 = float(np.mean(x1))

        if type2 == "discrete":
            sigma2_2 = 1.0
            a2 = 0.0
        else:
            var_x2 = float(np.var(x2)) if float(np.var(x2)) > 0 else 1.0
            sigma2_2 = var_x2 / mean_diag_C
            sigma2_2 = max(sigma2_2, 1e-10)
            a2 = float(np.mean(x2))

        r = 0.0

        # Proposal variances (only used for estimated parameters)
        prop_var = {
            "log_sigma2_1": 0.2,
            "log_sigma2_2": 0.2,
            "r": 0.1,
            "a1": 0.5,
            "a2": 0.5,
        }

        # Which parameters to update via MH
        update_sigma2_1 = (type1 == "continuous")
        update_sigma2_2 = (type2 == "continuous")
        update_a1 = (type1 == "continuous")
        update_a2 = (type2 == "continuous")

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

        current_ll = ThresholdModel._log_likelihood_bivariate_bm(
            x1, x2, C_inv, sigma2_1, sigma2_2, r, a1, a2, logdet_C, n
        )

        # Adaptive tuning interval
        tune_interval = max(100, sample)

        for gen in range(ngen):
            # Gibbs step: sample liabilities for discrete traits
            if type1 == "discrete":
                x1 = ThresholdModel._sample_liabilities_gibbs(
                    x1, states1, C, C_inv, sigma2_1, r,
                    x2, sigma2_2, a1, a2, rng
                )
                current_ll = ThresholdModel._log_likelihood_bivariate_bm(
                    x1, x2, C_inv, sigma2_1, sigma2_2, r, a1, a2, logdet_C, n
                )

            if type2 == "discrete":
                x2 = ThresholdModel._sample_liabilities_gibbs(
                    x2, states2, C, C_inv, sigma2_2, r,
                    x1, sigma2_1, a2, a1, rng
                )
                current_ll = ThresholdModel._log_likelihood_bivariate_bm(
                    x1, x2, C_inv, sigma2_1, sigma2_2, r, a1, a2, logdet_C, n
                )

            # Metropolis step: sigma2_1 (only for continuous traits)
            if update_sigma2_1:
                attempt["log_sigma2_1"] += 1
                log_s1_prop = np.log(sigma2_1) + rng.normal(0, np.sqrt(prop_var["log_sigma2_1"]))
                s1_prop = np.exp(log_s1_prop)
                if np.isfinite(s1_prop) and s1_prop > 0:
                    prop_ll = ThresholdModel._log_likelihood_bivariate_bm(
                        x1, x2, C_inv, s1_prop, sigma2_2, r, a1, a2, logdet_C, n
                    )
                    log_alpha = prop_ll - current_ll + log_s1_prop - np.log(sigma2_1)
                    if np.isfinite(log_alpha) and np.log(rng.random()) < log_alpha:
                        sigma2_1 = s1_prop
                        current_ll = prop_ll
                        accept["log_sigma2_1"] += 1

            # Metropolis step: sigma2_2 (only for continuous traits)
            if update_sigma2_2:
                attempt["log_sigma2_2"] += 1
                log_s2_prop = np.log(sigma2_2) + rng.normal(0, np.sqrt(prop_var["log_sigma2_2"]))
                s2_prop = np.exp(log_s2_prop)
                if np.isfinite(s2_prop) and s2_prop > 0:
                    prop_ll = ThresholdModel._log_likelihood_bivariate_bm(
                        x1, x2, C_inv, sigma2_1, s2_prop, r, a1, a2, logdet_C, n
                    )
                    log_alpha = prop_ll - current_ll + log_s2_prop - np.log(sigma2_2)
                    if np.isfinite(log_alpha) and np.log(rng.random()) < log_alpha:
                        sigma2_2 = s2_prop
                        current_ll = prop_ll
                        accept["log_sigma2_2"] += 1

            # Metropolis step: r (normal proposal with reflection on [-1,1])
            attempt["r"] += 1
            r_prop = r + rng.normal(0, np.sqrt(prop_var["r"]))
            # Reflect into [-1, 1]
            while r_prop < -1 or r_prop > 1:
                if r_prop < -1:
                    r_prop = -2 - r_prop
                if r_prop > 1:
                    r_prop = 2 - r_prop
            prop_ll = ThresholdModel._log_likelihood_bivariate_bm(
                x1, x2, C_inv, sigma2_1, sigma2_2, r_prop, a1, a2, logdet_C, n
            )
            log_alpha = prop_ll - current_ll
            if np.isfinite(log_alpha) and np.log(rng.random()) < log_alpha:
                r = r_prop
                current_ll = prop_ll
                accept["r"] += 1

            # Metropolis step: a1 (only for continuous traits)
            if update_a1:
                attempt["a1"] += 1
                a1_prop = a1 + rng.normal(0, np.sqrt(prop_var["a1"]))
                prop_ll = ThresholdModel._log_likelihood_bivariate_bm(
                    x1, x2, C_inv, sigma2_1, sigma2_2, r, a1_prop, a2, logdet_C, n
                )
                log_alpha = prop_ll - current_ll
                if np.isfinite(log_alpha) and np.log(rng.random()) < log_alpha:
                    a1 = a1_prop
                    current_ll = prop_ll
                    accept["a1"] += 1

            # Metropolis step: a2 (only for continuous traits)
            if update_a2:
                attempt["a2"] += 1
                a2_prop = a2 + rng.normal(0, np.sqrt(prop_var["a2"]))
                prop_ll = ThresholdModel._log_likelihood_bivariate_bm(
                    x1, x2, C_inv, sigma2_1, sigma2_2, r, a1, a2_prop, logdet_C, n
                )
                log_alpha = prop_ll - current_ll
                if np.isfinite(log_alpha) and np.log(rng.random()) < log_alpha:
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
        n = len(sorted_samples)
        interval_size = int(np.ceil(credible_mass * n))
        if interval_size >= n:
            return float(sorted_samples[0]), float(sorted_samples[-1])

        widths = sorted_samples[interval_size:] - sorted_samples[:n - interval_size]
        best = int(np.argmin(widths))
        return float(sorted_samples[best]), float(sorted_samples[best + interval_size])

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
            hpd_lo, hpd_hi = ThresholdModel._compute_hpd(samples)
            summary[param] = {
                "mean": float(np.mean(samples)),
                "median": float(np.median(samples)),
                "hpd_lower": hpd_lo,
                "hpd_upper": hpd_hi,
            }
        summary["acceptance_rates"] = mcmc_result["acceptance_rates"]
        return summary

    def _output_text(self, summary):
        try:
            print(
                f"Trait 1: {self.traits[0]} ({self.types[0]}"
                + (", 2 states: 0, 1" if self.types[0] == "discrete" else "")
                + ")"
            )
            print(
                f"Trait 2: {self.traits[1]} ({self.types[1]}"
                + (", 2 states: 0, 1" if self.types[1] == "discrete" else "")
                + ")"
            )
            print(
                f"MCMC: {self.ngen} generations, sampled every "
                f"{self.sample}, burn-in {int(self.burnin * 100)}%"
            )
            print("---")
            r_s = summary["r"]
            print(
                f"Posterior correlation (r): {r_s['mean']:.4f} "
                f"(95% HPD: {r_s['hpd_lower']:.3f}, {r_s['hpd_upper']:.3f})"
            )
            s1 = summary["sigma2_1"]
            print(
                f"Posterior sigma2_1: {s1['mean']:.4f} "
                f"(95% HPD: {s1['hpd_lower']:.3f}, {s1['hpd_upper']:.3f})"
            )
            s2 = summary["sigma2_2"]
            print(
                f"Posterior sigma2_2: {s2['mean']:.4f} "
                f"(95% HPD: {s2['hpd_lower']:.3f}, {s2['hpd_upper']:.3f})"
            )
            rates = summary["acceptance_rates"]
            print(
                f"Acceptance rates: r={rates.get('r', 0):.3f}, "
                f"sigma2_1={rates.get('log_sigma2_1', 0):.3f}, "
                f"sigma2_2={rates.get('log_sigma2_2', 0):.3f}, "
                f"a1={rates.get('a1', 0):.3f}, "
                f"a2={rates.get('a2', 0):.3f}"
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

    @staticmethod
    def _plot_trace(mcmc_result, output_path):
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

        fig, axes = plt.subplots(3, 2, figsize=(12, 8))

        for row, (param, label) in enumerate(params):
            samples = mcmc_result[param]
            mean_val = float(np.mean(samples))
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

        axes[0, 0].set_title("Trace", fontsize=13)
        axes[0, 1].set_title("Posterior", fontsize=13)

        fig.tight_layout()
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close(fig)
