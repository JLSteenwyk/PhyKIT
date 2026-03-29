"""
Phylogenetic imputation of missing trait values using conditional
multivariate normal distributions that capture both phylogenetic
relationships and between-trait correlations.
"""
import sys
from typing import Dict, List, Tuple

import numpy as np

from .base import Tree
from ...helpers.json_output import print_json
from ...errors import PhykitUserError


class PhyloImpute(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"])
        self.trait_data_path = parsed["trait_data_path"]
        self.output_path = parsed["output_path"]
        self.gene_trees_path = parsed["gene_trees_path"]
        self.json_output = parsed["json_output"]

    def run(self) -> None:
        from .vcv_utils import build_vcv_matrix, build_discordance_vcv, parse_gene_trees

        tree = self.read_tree_file()
        self.validate_tree(tree, min_tips=3, require_branch_lengths=True, context="phylogenetic imputation")

        tree_tips = self.get_tip_names_from_tree(tree)
        trait_names, trait_data, missing_info = self._parse_trait_file_with_na(
            self.trait_data_path, tree_tips
        )

        ordered_names = sorted(trait_data.keys())
        n = len(ordered_names)
        p = len(trait_names)

        # Build data matrix (n x p) with NaN for missing
        Y = np.full((n, p), np.nan)
        for i, name in enumerate(ordered_names):
            for j in range(p):
                Y[i, j] = trait_data[name][j]

        # Build VCV
        if self.gene_trees_path:
            gene_trees = parse_gene_trees(self.gene_trees_path)
            vcv, vcv_meta = build_discordance_vcv(tree, gene_trees, ordered_names)
            shared = vcv_meta["shared_taxa"]
            if set(shared) != set(ordered_names):
                # Subset to shared taxa
                keep_idx = [i for i, name in enumerate(ordered_names) if name in set(shared)]
                ordered_names = [ordered_names[i] for i in keep_idx]
                Y = Y[keep_idx, :]
                n = len(ordered_names)
            vcv_type = "discordance"
        else:
            vcv = build_vcv_matrix(tree, ordered_names)
            vcv_type = "BM"

        # Find complete cases (taxa with no missing values)
        complete_mask = ~np.any(np.isnan(Y), axis=1)
        n_complete = int(np.sum(complete_mask))

        if n_complete < 2:
            raise PhykitUserError(
                [
                    f"Only {n_complete} taxa have complete data for all traits.",
                    "At least 2 taxa with complete data are required for imputation.",
                ],
                code=2,
            )

        # Step 1: Estimate phylogenetic mean and trait covariance from complete cases
        complete_idx = np.where(complete_mask)[0]
        Y_complete = Y[complete_idx, :]
        C_complete = vcv[np.ix_(complete_idx, complete_idx)]

        try:
            C_complete_inv = np.linalg.inv(C_complete)
        except np.linalg.LinAlgError:
            C_complete_inv = np.linalg.pinv(C_complete)

        ones = np.ones(n_complete)
        denom = ones @ C_complete_inv @ ones

        # GLS phylogenetic mean for each trait
        a_hat = np.array(
            [(ones @ C_complete_inv @ Y_complete[:, j]) / denom for j in range(p)]
        )

        # GLS residuals
        E = Y_complete - a_hat

        # Trait covariance from phylogenetic residuals
        if n_complete > 1:
            Sigma_trait = (E.T @ C_complete_inv @ E) / (n_complete - 1)
        else:
            Sigma_trait = np.eye(p)

        # Step 2: Impute each missing value
        imputed_results = []
        Y_imputed = Y.copy()

        for i in range(n):
            missing_traits = np.where(np.isnan(Y[i, :]))[0]
            if len(missing_traits) == 0:
                continue

            observed_traits = np.where(~np.isnan(Y[i, :]))[0]

            imp_values, imp_ses = self._impute_taxon(
                Y_imputed, vcv, i, list(observed_traits), list(missing_traits),
                a_hat, Sigma_trait
            )

            for j in missing_traits:
                val = imp_values[j]
                se = imp_ses[j]
                ci_lower = val - 1.96 * se
                ci_upper = val + 1.96 * se

                Y_imputed[i, j] = val

                imputed_results.append({
                    "taxon": ordered_names[i],
                    "trait": trait_names[j],
                    "value": round(val, 6),
                    "se": round(se, 6),
                    "ci_lower": round(ci_lower, 6),
                    "ci_upper": round(ci_upper, 6),
                })

        # Write output TSV
        self._write_output_tsv(
            self.output_path, trait_names, ordered_names, Y_imputed
        )

        # Print summary
        if self.json_output:
            self._print_json(
                n, p, trait_names, vcv_type, imputed_results, self.output_path
            )
        else:
            self._print_text(
                n, p, trait_names, vcv_type, imputed_results,
                self.trait_data_path, self.output_path,
            )

    def process_args(self, args) -> Dict:
        return dict(
            tree_file_path=args.tree,
            trait_data_path=args.trait_data,
            output_path=args.output,
            gene_trees_path=getattr(args, "gene_trees", None),
            json_output=getattr(args, "json", False),
        )

    def _parse_trait_file_with_na(
        self, path: str, tree_tips: List[str]
    ) -> Tuple[List[str], Dict[str, List[float]], List[dict]]:
        """Parse multi-trait TSV allowing NA/? for missing values."""
        missing_markers = {"NA", "na", "Na", "?", ""}

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
        missing_info = []
        for line_idx, line in enumerate(data_lines[1:], 2):
            parts = line.split("\t")
            if len(parts) != n_cols:
                raise PhykitUserError(
                    [
                        f"Line {line_idx} has {len(parts)} columns; expected {n_cols}.",
                    ],
                    code=2,
                )
            taxon = parts[0]
            values = []
            for i, val_str in enumerate(parts[1:]):
                if val_str.strip() in missing_markers:
                    values.append(float("nan"))
                    missing_info.append({
                        "taxon": taxon,
                        "trait": trait_names[i],
                        "line": line_idx,
                    })
                else:
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
        filtered_missing = [m for m in missing_info if m["taxon"] in shared]
        return trait_names, filtered, filtered_missing

    @staticmethod
    def _impute_taxon(
        Y: np.ndarray,
        C_phylo: np.ndarray,
        taxon_idx: int,
        observed_trait_indices: List[int],
        missing_trait_indices: List[int],
        a_hat: np.ndarray,
        Sigma_trait: np.ndarray,
    ) -> Tuple[Dict[int, float], Dict[int, float]]:
        """Impute missing traits for one taxon using phylogeny + trait correlations."""
        n = C_phylo.shape[0]

        other_idx = [k for k in range(n) if k != taxon_idx]

        C_io = C_phylo[taxon_idx, other_idx]
        C_oo = C_phylo[np.ix_(other_idx, other_idx)]

        imputed = {}
        imputed_se = {}

        for j in missing_trait_indices:
            # Get observed values of trait j from other taxa
            obs_other = [
                (k_idx, Y[other_idx[k_idx], j])
                for k_idx in range(len(other_idx))
                if not np.isnan(Y[other_idx[k_idx], j])
            ]

            if not obs_other:
                # No other taxa have this trait -- use phylogenetic mean
                imputed[j] = float(a_hat[j])
                imputed_se[j] = float(np.sqrt(max(Sigma_trait[j, j], 0.0)))
                continue

            obs_k_indices = [x[0] for x in obs_other]
            obs_values = np.array([x[1] for x in obs_other])

            # Phylogenetic prediction
            c_io_sub = C_io[obs_k_indices]
            C_oo_sub = C_oo[np.ix_(obs_k_indices, obs_k_indices)]
            try:
                C_oo_sub_inv = np.linalg.inv(C_oo_sub)
            except np.linalg.LinAlgError:
                C_oo_sub_inv = np.linalg.pinv(C_oo_sub)

            weights = c_io_sub @ C_oo_sub_inv
            y_hat = float(a_hat[j] + weights @ (obs_values - a_hat[j]))

            # Conditional variance
            c_ii = C_phylo[taxon_idx, taxon_idx]
            cond_var = (c_ii - c_io_sub @ C_oo_sub_inv @ c_io_sub) * Sigma_trait[j, j]
            se = float(np.sqrt(max(cond_var, 0.0)))

            # If taxon has observed traits that correlate with trait j,
            # refine the prediction using trait correlations
            if len(observed_trait_indices) > 0:
                obs_traits_values = np.array(
                    [Y[taxon_idx, k] for k in observed_trait_indices]
                )
                obs_means = a_hat[observed_trait_indices]

                Sigma_jo = Sigma_trait[j, observed_trait_indices]
                Sigma_oo_traits = Sigma_trait[
                    np.ix_(observed_trait_indices, observed_trait_indices)
                ]
                try:
                    Sigma_oo_inv = np.linalg.inv(Sigma_oo_traits)
                    # Conditional expectation adjustment
                    trait_adjustment = float(
                        Sigma_jo @ Sigma_oo_inv @ (obs_traits_values - obs_means)
                    )
                    y_hat += trait_adjustment

                    # Refined variance
                    cond_var_trait = Sigma_trait[j, j] - Sigma_jo @ Sigma_oo_inv @ Sigma_jo
                    se = float(
                        np.sqrt(max(cond_var * cond_var_trait / Sigma_trait[j, j], 0.0))
                    )
                except np.linalg.LinAlgError:
                    pass  # keep phylogeny-only prediction

            imputed[j] = float(y_hat)
            imputed_se[j] = float(se)

        return imputed, imputed_se

    @staticmethod
    def _write_output_tsv(
        path: str,
        trait_names: List[str],
        ordered_names: List[str],
        Y: np.ndarray,
    ) -> None:
        """Write imputed data as TSV (drop-in replacement for the input)."""
        with open(path, "w") as f:
            header = "taxon\t" + "\t".join(trait_names) + "\n"
            f.write(header)
            for i, name in enumerate(ordered_names):
                vals = "\t".join(f"{Y[i, j]:.6f}" for j in range(Y.shape[1]))
                f.write(f"{name}\t{vals}\n")

    def _print_text(
        self, n, p, trait_names, vcv_type, imputed_results, data_path, output_path
    ) -> None:
        print("Phylogenetic Imputation")
        print(f"Data: {data_path}")
        print(f"Taxa: {n}")
        print(f"Traits: {p} ({', '.join(trait_names)})")
        print(f"Missing values: {len(imputed_results)}")
        vcv_label = "standard" if vcv_type == "BM" else "discordance-aware"
        print(f"VCV: {vcv_type} ({vcv_label})")

        if imputed_results:
            # Determine column widths
            max_taxon = max(len(r["taxon"]) for r in imputed_results)
            max_trait = max(len(r["trait"]) for r in imputed_results)
            max_taxon = max(max_taxon, len("Taxon"))
            max_trait = max(max_trait, len("Trait"))

            print("\nImputed values:")
            header = (
                f"  {'Taxon':<{max_taxon}}  "
                f"{'Trait':<{max_trait}}  "
                f"{'Imputed':>10}  "
                f"{'SE':>10}  "
                f"{'95% CI':>20}"
            )
            print(header)

            for r in imputed_results:
                ci_str = f"[{r['ci_lower']:.4f}, {r['ci_upper']:.4f}]"
                row = (
                    f"  {r['taxon']:<{max_taxon}}  "
                    f"  {r['trait']:<{max_trait}}  "
                    f"{r['value']:>10.4f}  "
                    f"{r['se']:>10.4f}  "
                    f"{ci_str:>20}"
                )
                print(row)

        print(f"\nOutput: {output_path}")

    def _print_json(self, n, p, trait_names, vcv_type, imputed_results, output_path):
        payload = {
            "n_taxa": n,
            "n_traits": p,
            "trait_names": trait_names,
            "n_missing": len(imputed_results),
            "vcv_type": vcv_type,
            "imputed": imputed_results,
            "output_file": output_path,
        }
        print_json(payload)
