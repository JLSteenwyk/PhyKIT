"""
Phylogenetic ANOVA / MANOVA using RRPP.

Implements the Residual Randomization Permutation Procedure (RRPP)
of Adams & Collyer (2018) for testing group differences while
accounting for phylogenetic non-independence.

Auto-detects univariate (ANOVA) vs multivariate (MANOVA) based on
the number of response trait columns, with optional user override.
"""
import sys
from typing import Dict, List, Tuple

import numpy as np

from .base import Tree
from ...helpers.json_output import print_json
from ...helpers.plot_config import PlotConfig
from ...errors import PhykitUserError


class PhyloAnova(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"])
        self.trait_data_path = parsed["trait_data_path"]
        self.group_column = parsed["group_column"]
        self.method = parsed["method"]
        self.permutations = parsed["permutations"]
        self.pairwise = parsed["pairwise"]
        self.plot_output = parsed["plot_output"]
        self.plot_type = parsed["plot_type"]
        self.json_output = parsed["json_output"]
        self.seed = parsed["seed"]
        self.plot_config = parsed["plot_config"]

    def process_args(self, args) -> Dict:
        return dict(
            tree_file_path=args.tree,
            trait_data_path=args.traits,
            group_column=getattr(args, "group_column", None),
            method=getattr(args, "method", "auto"),
            permutations=getattr(args, "permutations", 1000),
            pairwise=getattr(args, "pairwise", False),
            plot_output=getattr(args, "plot_output", None),
            plot_type=getattr(args, "plot_type", "auto"),
            json_output=getattr(args, "json", False),
            seed=getattr(args, "seed", None),
            plot_config=PlotConfig.from_args(args),
        )

    def run(self) -> None:
        from .vcv_utils import build_vcv_matrix

        tree = self.read_tree_file()
        self._validate_tree(tree)

        tree_tips = self.get_tip_names_from_tree(tree)
        header, traits = self._parse_trait_file(self.trait_data_path, tree_tips)

        # Identify group column
        group_col = self.group_column or header[0]
        if group_col not in header:
            raise PhykitUserError(
                [
                    f"Group column '{group_col}' not found in trait file.",
                    f"Available columns: {', '.join(header)}",
                ],
                code=2,
            )
        group_idx = header.index(group_col)

        # Separate response columns (numeric) from group column
        response_names = [h for i, h in enumerate(header) if i != group_idx]
        if not response_names:
            raise PhykitUserError(
                ["No response trait columns found (only the group column)."],
                code=2,
            )

        # Determine method
        n_responses = len(response_names)
        if self.method == "auto":
            method = "manova" if n_responses > 1 else "anova"
        else:
            method = self.method
            if method == "anova" and n_responses > 1:
                raise PhykitUserError(
                    [
                        f"ANOVA requires 1 response trait but found {n_responses}.",
                        "Use --method manova for multiple response traits.",
                    ],
                    code=2,
                )

        # Build ordered data
        ordered_names = sorted(traits.keys())
        n = len(ordered_names)

        groups = []
        Y_list = []
        for name in ordered_names:
            row = traits[name]
            groups.append(row[group_idx])
            Y_list.append(
                [row[i] for i in range(len(header)) if i != group_idx]
            )

        Y = np.array(Y_list, dtype=float)
        unique_groups = sorted(set(groups))
        n_groups = len(unique_groups)

        if n_groups < 2:
            raise PhykitUserError(
                ["At least 2 groups are required for ANOVA/MANOVA."],
                code=2,
            )

        # Build design matrix (intercept + dummy variables)
        X_full = self._build_design_matrix(groups, unique_groups, n)
        X_reduced = np.ones((n, 1))

        # Build phylogenetic VCV
        vcv = build_vcv_matrix(tree, ordered_names)

        # Cholesky transform
        try:
            L = np.linalg.cholesky(vcv)
            L_inv = np.linalg.inv(L)
        except np.linalg.LinAlgError:
            raise PhykitUserError(
                [
                    "VCV matrix is not positive definite.",
                    "Check that the tree has valid branch lengths.",
                ],
                code=2,
            )

        # Transform data
        Y_star = L_inv @ Y
        X_full_star = L_inv @ X_full
        X_red_star = L_inv @ X_reduced

        # Fit models
        if method == "anova":
            result = self._run_anova(
                Y_star, X_full_star, X_red_star, n, n_groups,
            )
        else:
            result = self._run_manova(
                Y_star, X_full_star, X_red_star, n, n_groups,
            )

        # Pairwise comparisons
        pairwise_results = None
        if self.pairwise:
            pairwise_results = self._run_pairwise(
                Y_star, groups, unique_groups, X_red_star,
            )

        # Output
        self._print_results(
            method, result, unique_groups, response_names,
            pairwise_results, n, n_groups,
        )

        if self.json_output:
            self._print_json(
                method, result, unique_groups, response_names,
                pairwise_results, n, n_groups,
            )

        # Plot
        if self.plot_output:
            plot_type = self.plot_type
            if plot_type == "auto":
                plot_type = "boxplot" if method == "anova" else "phylomorphospace"
            if plot_type == "boxplot":
                self._plot_boxplot(
                    tree, ordered_names, Y, groups, unique_groups,
                    response_names, self.plot_output,
                )
            else:
                self._plot_phylomorphospace(
                    tree, ordered_names, Y, groups, unique_groups,
                    response_names, self.plot_output,
                )

    def _validate_tree(self, tree) -> None:
        for clade in tree.find_clades():
            if clade.branch_length is None and clade != tree.root:
                clade.branch_length = 1e-8

    def _parse_trait_file(
        self, path: str, tree_tips: List[str]
    ) -> Tuple[List[str], Dict[str, list]]:
        """Parse TSV with header. Group column kept as string; others as float."""
        try:
            with open(path) as f:
                lines = f.readlines()
        except FileNotFoundError:
            raise PhykitUserError(
                [f"{path} corresponds to no such file or directory."],
                code=2,
            )

        data_lines = [
            l.strip() for l in lines
            if l.strip() and not l.strip().startswith("#")
        ]
        if len(data_lines) < 2:
            raise PhykitUserError(
                ["Trait file must have a header row and at least one data row."],
                code=2,
            )

        header_parts = data_lines[0].split("\t")
        n_cols = len(header_parts)
        if n_cols < 3:
            raise PhykitUserError(
                ["Header must have at least 3 columns (taxon + group + trait)."],
                code=2,
            )
        header = header_parts[1:]

        # Detect which column is the group (categorical) column
        group_col_name = self.group_column or header[0]

        traits = {}
        for line_idx, line in enumerate(data_lines[1:], 2):
            parts = line.split("\t")
            if len(parts) != n_cols:
                raise PhykitUserError(
                    [f"Line {line_idx} has {len(parts)} columns; expected {n_cols}."],
                    code=2,
                )
            taxon = parts[0]
            values = []
            for i, val_str in enumerate(parts[1:]):
                col_name = header[i]
                if col_name == group_col_name:
                    values.append(val_str.strip())
                else:
                    try:
                        values.append(float(val_str))
                    except ValueError:
                        raise PhykitUserError(
                            [
                                f"Non-numeric trait value '{val_str}' for taxon "
                                f"'{taxon}' (trait '{col_name}') on line {line_idx}.",
                            ],
                            code=2,
                        )
            traits[taxon] = values

        tree_tip_set = set(tree_tips)
        trait_taxa_set = set(traits.keys())
        shared = tree_tip_set & trait_taxa_set

        if tree_tip_set - trait_taxa_set:
            print(
                f"Warning: {len(tree_tip_set - trait_taxa_set)} taxa in tree "
                f"but not in trait file",
                file=sys.stderr,
            )
        if trait_taxa_set - tree_tip_set:
            print(
                f"Warning: {len(trait_taxa_set - tree_tip_set)} taxa in trait "
                f"file but not in tree",
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
        return header, filtered

    @staticmethod
    def _build_design_matrix(
        groups: List[str], unique_groups: List[str], n: int
    ) -> np.ndarray:
        """Build design matrix: intercept + (g-1) dummy columns."""
        X = np.ones((n, len(unique_groups)))
        for i, g in enumerate(groups):
            col = unique_groups.index(g)
            X[i, col] = 1.0
        # Use treatment coding: intercept + dummies for groups 1..g-1
        X = np.zeros((n, len(unique_groups)))
        X[:, 0] = 1.0  # intercept
        for i, g in enumerate(groups):
            idx = unique_groups.index(g)
            if idx > 0:
                X[i, idx] = 1.0
        return X

    def _run_anova(
        self,
        Y_star: np.ndarray,
        X_full_star: np.ndarray,
        X_red_star: np.ndarray,
        n: int,
        n_groups: int,
    ) -> Dict:
        """Univariate phylogenetic ANOVA with RRPP."""
        y = Y_star.ravel() if Y_star.ndim > 1 else Y_star
        if Y_star.ndim > 1:
            y = Y_star[:, 0]

        # Fit models
        hat_red = X_red_star @ np.linalg.lstsq(X_red_star, y, rcond=None)[0]
        resid_red = y - hat_red

        hat_full = X_full_star @ np.linalg.lstsq(X_full_star, y, rcond=None)[0]
        resid_full = y - hat_full

        ss_total = float(resid_red @ resid_red)
        ss_resid = float(resid_full @ resid_full)
        ss_model = ss_total - ss_resid

        df_model = n_groups - 1
        df_resid = n - n_groups
        ms_model = ss_model / df_model if df_model > 0 else 0.0
        ms_resid = ss_resid / df_resid if df_resid > 0 else 0.0
        f_obs = ms_model / ms_resid if ms_resid > 0 else 0.0

        # RRPP permutations
        rng = np.random.RandomState(self.seed)
        f_perms = np.zeros(self.permutations)
        for p in range(self.permutations):
            perm_idx = rng.permutation(n)
            resid_perm = resid_red[perm_idx]
            y_perm = hat_red + resid_perm

            hat_full_p = X_full_star @ np.linalg.lstsq(
                X_full_star, y_perm, rcond=None
            )[0]
            resid_full_p = y_perm - hat_full_p

            ss_resid_p = float(resid_full_p @ resid_full_p)
            ss_model_p = ss_total - ss_resid_p
            ms_model_p = ss_model_p / df_model if df_model > 0 else 0.0
            ms_resid_p = ss_resid_p / df_resid if df_resid > 0 else 0.0
            f_perms[p] = ms_model_p / ms_resid_p if ms_resid_p > 0 else 0.0

        p_value = float(np.mean(f_perms >= f_obs))
        z_score = (
            (f_obs - np.mean(f_perms)) / np.std(f_perms)
            if np.std(f_perms) > 0
            else 0.0
        )

        return dict(
            ss_model=ss_model,
            ss_resid=ss_resid,
            ss_total=ss_total,
            df_model=df_model,
            df_resid=df_resid,
            df_total=n - 1,
            ms_model=ms_model,
            ms_resid=ms_resid,
            f_stat=f_obs,
            z_score=z_score,
            p_value=p_value,
            test_statistic_name="F",
        )

    def _run_manova(
        self,
        Y_star: np.ndarray,
        X_full_star: np.ndarray,
        X_red_star: np.ndarray,
        n: int,
        n_groups: int,
    ) -> Dict:
        """Multivariate phylogenetic MANOVA with RRPP using Pillai's trace."""
        # Fit models
        beta_red = np.linalg.lstsq(X_red_star, Y_star, rcond=None)[0]
        hat_red = X_red_star @ beta_red
        resid_red = Y_star - hat_red

        beta_full = np.linalg.lstsq(X_full_star, Y_star, rcond=None)[0]
        hat_full = X_full_star @ beta_full
        resid_full = Y_star - hat_full

        # SSCP matrices
        SS_total = resid_red.T @ resid_red
        SS_resid = resid_full.T @ resid_full
        SS_model = SS_total - SS_resid

        df_model = n_groups - 1
        df_resid = n - n_groups

        # Pillai's trace
        try:
            SS_resid_inv = np.linalg.inv(SS_resid)
            H_E_inv = SS_model @ SS_resid_inv
            eigenvalues = np.real(np.linalg.eigvals(H_E_inv))
            pillai = float(np.sum(eigenvalues / (1.0 + eigenvalues)))
        except np.linalg.LinAlgError:
            pillai = 0.0
            eigenvalues = np.array([0.0])

        # Approximate F for Pillai's trace
        p = Y_star.shape[1]  # number of response variables
        s = min(df_model, p)
        m = (abs(df_model - p) - 1) / 2.0
        nn = (df_resid - p - 1) / 2.0
        if s > 0 and (s + 2 * nn + 2) > 0:
            f_approx = (
                (pillai / s)
                * ((2 * nn + s + 1) / (2 * m + s + 1))
            )
        else:
            f_approx = 0.0

        # RRPP permutations
        rng = np.random.RandomState(self.seed)
        pillai_perms = np.zeros(self.permutations)
        for perm in range(self.permutations):
            perm_idx = rng.permutation(n)
            resid_perm = resid_red[perm_idx]
            Y_perm = hat_red + resid_perm

            beta_full_p = np.linalg.lstsq(X_full_star, Y_perm, rcond=None)[0]
            hat_full_p = X_full_star @ beta_full_p
            resid_full_p = Y_perm - hat_full_p

            SS_resid_p = resid_full_p.T @ resid_full_p
            SS_model_p = SS_total - SS_resid_p

            try:
                SS_resid_inv_p = np.linalg.inv(SS_resid_p)
                H_E_inv_p = SS_model_p @ SS_resid_inv_p
                eig_p = np.real(np.linalg.eigvals(H_E_inv_p))
                pillai_perms[perm] = float(
                    np.sum(eig_p / (1.0 + eig_p))
                )
            except np.linalg.LinAlgError:
                pillai_perms[perm] = 0.0

        p_value = float(np.mean(pillai_perms >= pillai))
        z_score = (
            (pillai - np.mean(pillai_perms)) / np.std(pillai_perms)
            if np.std(pillai_perms) > 0
            else 0.0
        )

        # Scalar SS for display (trace of SSCP)
        ss_model = float(np.trace(SS_model))
        ss_resid = float(np.trace(SS_resid))
        ss_total = float(np.trace(SS_total))

        return dict(
            ss_model=ss_model,
            ss_resid=ss_resid,
            ss_total=ss_total,
            df_model=df_model,
            df_resid=df_resid,
            df_total=n - 1,
            ms_model=ss_model / df_model if df_model > 0 else 0.0,
            ms_resid=ss_resid / df_resid if df_resid > 0 else 0.0,
            f_stat=f_approx,
            pillai_trace=pillai,
            z_score=z_score,
            p_value=p_value,
            test_statistic_name="Pillai's trace",
        )

    def _run_pairwise(
        self,
        Y_star: np.ndarray,
        groups: List[str],
        unique_groups: List[str],
        X_red_star: np.ndarray,
    ) -> List[Dict]:
        """Pairwise group comparisons via RRPP."""
        rng = np.random.RandomState(self.seed)
        results = []

        # Reduced model residuals
        beta_red = np.linalg.lstsq(X_red_star, Y_star, rcond=None)[0]
        hat_red = X_red_star @ beta_red
        resid_red = Y_star - hat_red

        groups_arr = np.array(groups)

        for i in range(len(unique_groups)):
            for j in range(i + 1, len(unique_groups)):
                g1, g2 = unique_groups[i], unique_groups[j]
                mask = (groups_arr == g1) | (groups_arr == g2)
                idx = np.where(mask)[0]

                Y_sub = Y_star[idx]
                groups_sub = groups_arr[idx]

                # Observed distance between group means
                mean1 = Y_sub[groups_sub == g1].mean(axis=0)
                mean2 = Y_sub[groups_sub == g2].mean(axis=0)
                d_obs = float(np.linalg.norm(mean1 - mean2))

                # Permute residuals for pairwise test
                resid_sub = resid_red[idx]
                hat_sub = hat_red[idx]
                d_perms = np.zeros(self.permutations)

                for p in range(self.permutations):
                    perm_idx = rng.permutation(len(idx))
                    Y_perm = hat_sub + resid_sub[perm_idx]
                    mean1_p = Y_perm[groups_sub == g1].mean(axis=0)
                    mean2_p = Y_perm[groups_sub == g2].mean(axis=0)
                    d_perms[p] = float(np.linalg.norm(mean1_p - mean2_p))

                p_val = float(np.mean(d_perms >= d_obs))
                z_val = (
                    (d_obs - np.mean(d_perms)) / np.std(d_perms)
                    if np.std(d_perms) > 0
                    else 0.0
                )

                results.append(dict(
                    group1=g1, group2=g2,
                    distance=d_obs, z_score=z_val, p_value=p_val,
                ))

        return results

    def _print_results(
        self, method, result, unique_groups, response_names,
        pairwise_results, n, n_groups,
    ) -> None:
        method_label = "Phylogenetic ANOVA" if method == "anova" else "Phylogenetic MANOVA"
        try:
            print(f"{method_label} (RRPP, {self.permutations} permutations)")
            if method == "anova":
                print(f"Response: {response_names[0]}")
            else:
                print(f"Response traits: {', '.join(response_names)}")
            print(f"Groups ({n_groups}): {', '.join(unique_groups)}")
            print(f"N taxa: {n}")
            print()
            print(f"{'Source':<15} {'Df':>5} {'SS':>12} {'MS':>12} "
                  f"{'F':>10} {'Z':>10} {'p-value':>10}")
            print("-" * 76)
            print(
                f"{'group':<15} {result['df_model']:>5d} "
                f"{result['ss_model']:>12.4f} {result['ms_model']:>12.4f} "
                f"{result['f_stat']:>10.4f} {result['z_score']:>10.4f} "
                f"{result['p_value']:>10.4f}"
            )
            print(
                f"{'residuals':<15} {result['df_resid']:>5d} "
                f"{result['ss_resid']:>12.4f} {result['ms_resid']:>12.4f}"
            )
            print(
                f"{'total':<15} {result['df_total']:>5d} "
                f"{result['ss_total']:>12.4f}"
            )

            if "pillai_trace" in result:
                print(f"\nPillai's trace: {result['pillai_trace']:.4f}")

            if pairwise_results:
                print("\nPairwise comparisons:")
                print(f"  {'Comparison':<30} {'d':>10} {'Z':>10} {'p-value':>10}")
                print("  " + "-" * 62)
                for pw in pairwise_results:
                    label = f"{pw['group1']} vs {pw['group2']}"
                    print(
                        f"  {label:<30} {pw['distance']:>10.4f} "
                        f"{pw['z_score']:>10.4f} {pw['p_value']:>10.4f}"
                    )
        except BrokenPipeError:
            return

    def _print_json(
        self, method, result, unique_groups, response_names,
        pairwise_results, n, n_groups,
    ) -> None:
        payload = {
            "method": method,
            "permutations": self.permutations,
            "n_taxa": n,
            "n_groups": n_groups,
            "groups": unique_groups,
            "response_traits": response_names,
            "anova_table": {
                "df_model": result["df_model"],
                "df_resid": result["df_resid"],
                "df_total": result["df_total"],
                "ss_model": round(result["ss_model"], 6),
                "ss_resid": round(result["ss_resid"], 6),
                "ss_total": round(result["ss_total"], 6),
                "ms_model": round(result["ms_model"], 6),
                "ms_resid": round(result["ms_resid"], 6),
                "f_stat": round(result["f_stat"], 6),
                "z_score": round(result["z_score"], 4),
                "p_value": round(result["p_value"], 4),
            },
        }
        if "pillai_trace" in result:
            payload["anova_table"]["pillai_trace"] = round(
                result["pillai_trace"], 6
            )
        if pairwise_results:
            payload["pairwise"] = [
                {
                    "group1": pw["group1"],
                    "group2": pw["group2"],
                    "distance": round(pw["distance"], 6),
                    "z_score": round(pw["z_score"], 4),
                    "p_value": round(pw["p_value"], 4),
                }
                for pw in pairwise_results
            ]
        print_json(payload, sort_keys=False)

    def _plot_boxplot(
        self, tree, ordered_names, Y, groups, unique_groups,
        response_names, output_path,
    ) -> None:
        """Violin + boxplot for univariate ANOVA."""
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        except ImportError:
            print("matplotlib is required for plotting.")
            return

        config = self.plot_config
        config.resolve(n_rows=len(unique_groups), n_cols=None)

        fig, ax = plt.subplots(figsize=(config.fig_width, config.fig_height))

        y = Y[:, 0] if Y.ndim > 1 else Y
        default_colors = ["#2b8cbe", "#d62728", "#969696", "#31a354",
                          "#756bb1", "#e6550d", "#636363", "#de2d26"]
        colors = config.merge_colors(default_colors)

        data_by_group = []
        group_labels = []
        group_colors = []
        for i, g in enumerate(unique_groups):
            mask = [groups[j] == g for j in range(len(groups))]
            vals = y[mask]
            data_by_group.append(vals)
            group_labels.append(g)
            group_colors.append(colors[i % len(colors)])

        positions = list(range(1, len(unique_groups) + 1))

        # Violin
        vp = ax.violinplot(
            data_by_group, positions=positions,
            showmeans=False, showmedians=False, showextrema=False,
        )
        for i, body in enumerate(vp["bodies"]):
            body.set_facecolor(group_colors[i])
            body.set_alpha(0.3)
            body.set_edgecolor(group_colors[i])

        # Boxplot
        bp = ax.boxplot(
            data_by_group, positions=positions,
            widths=0.15, patch_artist=True,
            showfliers=True, zorder=3,
        )
        for i, box in enumerate(bp["boxes"]):
            box.set_facecolor(group_colors[i])
            box.set_alpha(0.7)
        for element in ["whiskers", "caps", "medians"]:
            for line in bp[element]:
                line.set_color("black")

        ax.set_xticks(positions)
        ax.set_xticklabels(group_labels)
        label_fs = config.xlabel_fontsize or 10
        ax.tick_params(axis="x", labelsize=label_fs)

        ylabel = response_names[0] if response_names else "Trait value"
        ylabel_fs = config.ylabel_fontsize if config.ylabel_fontsize else 10
        ax.set_ylabel(ylabel, fontsize=ylabel_fs)

        if config.show_title:
            ax.set_title(
                config.title or "Phylogenetic ANOVA",
                fontsize=config.title_fontsize,
            )

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        fig.tight_layout()
        fig.savefig(output_path, dpi=config.dpi)
        plt.close(fig)
        print(f"Plot saved: {output_path}")

    def _plot_phylomorphospace(
        self, tree, ordered_names, Y, groups, unique_groups,
        response_names, output_path,
    ) -> None:
        """Phylomorphospace colored by group for multivariate MANOVA."""
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        except ImportError:
            print("matplotlib is required for plotting.")
            return

        config = self.plot_config
        config.resolve(n_rows=len(ordered_names), n_cols=None)

        fig, ax = plt.subplots(figsize=(config.fig_width, config.fig_height))

        default_colors = ["#2b8cbe", "#d62728", "#969696", "#31a354",
                          "#756bb1", "#e6550d", "#636363", "#de2d26"]
        colors = config.merge_colors(default_colors)

        # Use first 2 traits as axes (or PCA if more)
        if Y.shape[1] == 2:
            pc = Y
            xlabel = response_names[0]
            ylabel = response_names[1]
        else:
            # PCA
            Y_centered = Y - Y.mean(axis=0)
            U, S, Vt = np.linalg.svd(Y_centered, full_matrices=False)
            pc = Y_centered @ Vt[:2].T
            var_explained = S[:2] ** 2 / np.sum(S ** 2) * 100
            xlabel = f"PC1 ({var_explained[0]:.1f}%)"
            ylabel = f"PC2 ({var_explained[1]:.1f}%)"

        # Build parent map and reconstruct ancestors for phylo overlay
        parent_map = {}
        for clade in tree.find_clades(order="preorder"):
            for child in clade.clades:
                parent_map[id(child)] = clade

        # Map tip names to indices
        name_to_idx = {name: i for i, name in enumerate(ordered_names)}

        # Simple ancestral reconstruction (weighted mean of descendants)
        node_coords = {}
        for tip in tree.get_terminals():
            if tip.name in name_to_idx:
                idx = name_to_idx[tip.name]
                node_coords[id(tip)] = pc[idx]

        for clade in tree.find_clades(order="postorder"):
            if id(clade) in node_coords:
                continue
            child_coords = [
                node_coords[id(c)] for c in clade.clades
                if id(c) in node_coords
            ]
            if child_coords:
                node_coords[id(clade)] = np.mean(child_coords, axis=0)

        # Draw phylogeny branches
        for clade in tree.find_clades(order="preorder"):
            if clade == tree.root:
                continue
            pid = parent_map.get(id(clade))
            if pid is None or id(pid) not in node_coords or id(clade) not in node_coords:
                continue
            p_coord = node_coords[id(pid)]
            c_coord = node_coords[id(clade)]
            ax.plot(
                [p_coord[0], c_coord[0]],
                [p_coord[1], c_coord[1]],
                color="gray", lw=0.5, alpha=0.5, zorder=1,
            )

        # Plot points colored by group
        for i, g in enumerate(unique_groups):
            mask = np.array([groups[j] == g for j in range(len(groups))])
            color = colors[i % len(colors)]
            ax.scatter(
                pc[mask, 0], pc[mask, 1],
                c=color, label=g, s=50, edgecolors="black",
                linewidths=0.5, zorder=3,
            )

        ax.set_xlabel(xlabel, fontsize=config.axis_fontsize or 10)
        ax.set_ylabel(ylabel, fontsize=config.axis_fontsize or 10)
        ax.legend(fontsize=8, frameon=True)

        if config.show_title:
            ax.set_title(
                config.title or "Phylogenetic MANOVA — Phylomorphospace",
                fontsize=config.title_fontsize,
            )

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        fig.tight_layout()
        fig.savefig(output_path, dpi=config.dpi)
        plt.close(fig)
        print(f"Plot saved: {output_path}")
