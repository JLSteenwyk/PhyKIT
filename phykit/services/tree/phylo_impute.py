"""
Phylogenetic imputation of missing trait values using conditional
multivariate normal distributions that capture both phylogenetic
relationships and between-trait correlations.
"""
from __future__ import annotations

import sys

from .base import Tree
from ...errors import PhykitUserError


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


class _LazyNumpy:
    def __getattr__(self, name):
        import numpy as _np

        return getattr(_np, name)


np = _LazyNumpy()


def cho_factor(*args, **kwargs):
    from scipy.linalg import cho_factor as _cho_factor

    return _cho_factor(*args, **kwargs)


def cho_solve(*args, **kwargs):
    from scipy.linalg import cho_solve as _cho_solve

    return _cho_solve(*args, **kwargs)


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

        tree = self.read_tree_file_unmodified()
        self.validate_tree(tree, min_tips=3, require_branch_lengths=True, context="phylogenetic imputation")

        tree_tips = self.get_tip_names_from_tree(tree)
        trait_names, trait_data, missing_info = self._parse_trait_file_with_na(
            self.trait_data_path, tree_tips
        )

        ordered_names = sorted(trait_data.keys())
        n = len(ordered_names)
        p = len(trait_names)

        # Build data matrix (n x p) with NaN for missing
        Y = self._build_data_matrix(ordered_names, trait_data)

        # Build VCV
        if self.gene_trees_path:
            gene_trees = parse_gene_trees(self.gene_trees_path)
            vcv, vcv_meta = build_discordance_vcv(tree, gene_trees, ordered_names)
            shared = vcv_meta["shared_taxa"]
            ordered_names, Y = self._subset_to_shared_taxa(ordered_names, Y, shared)
            n = len(ordered_names)
            vcv_type = "discordance"
        else:
            vcv = build_vcv_matrix(tree, ordered_names)
            vcv_type = "BM"

        # Find complete cases (taxa with no missing values)
        missing_mask = np.isnan(Y)
        complete_mask = ~np.any(missing_mask, axis=1)
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

        a_hat, Sigma_trait = self._estimate_complete_case_stats(
            Y_complete, C_complete
        )

        # Step 2: Impute each missing value
        imputed_results = []
        Y_imputed = Y.copy()

        for i in np.flatnonzero(~complete_mask):
            row_missing = missing_mask[i]
            missing_traits = np.flatnonzero(row_missing)
            observed_traits = np.flatnonzero(~row_missing)

            imp_values, imp_ses = self._impute_taxon(
                Y_imputed, vcv, i, observed_traits, missing_traits,
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

    def process_args(self, args) -> dict:
        return dict(
            tree_file_path=args.tree,
            trait_data_path=args.trait_data,
            output_path=args.output,
            gene_trees_path=getattr(args, "gene_trees", None),
            json_output=getattr(args, "json", False),
        )

    @staticmethod
    def _build_data_matrix(ordered_names: list[str], trait_data: dict) -> np.ndarray:
        """Build the trait matrix in the exact requested taxon order."""
        return np.array([trait_data[name] for name in ordered_names], dtype=float)

    def _parse_trait_file_with_na(
        self, path: str, tree_tips: list[str]
    ) -> tuple[list[str], dict[str, list[float]], list[dict]]:
        """Parse multi-trait TSV allowing NA/? for missing values."""
        missing_markers = {"NA", "na", "Na", "?", ""}
        tree_tip_set = set(tree_tips)

        try:
            with open(path) as f:
                header_parts = None
                for line in f:
                    stripped = line.strip()
                    if not stripped or stripped.startswith("#"):
                        continue
                    header_parts = stripped.split("\t")
                    break

                if header_parts is None:
                    raise PhykitUserError(
                        [
                            "Multi-trait file must have a header row and at least one data row.",
                        ],
                        code=2,
                    )

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
                trait_taxa_set = set()
                missing_info = []
                saw_data = False
                data_line_idx = 2
                for line in f:
                    line = line.strip()
                    if not line or line.startswith("#"):
                        continue
                    saw_data = True
                    parts = line.split("\t")
                    if len(parts) != n_cols:
                        raise PhykitUserError(
                            [
                                f"Line {data_line_idx} has {len(parts)} columns; expected {n_cols}.",
                            ],
                            code=2,
                        )
                    taxon = parts[0]
                    trait_taxa_set.add(taxon)
                    keep_taxon = taxon in tree_tip_set
                    value_parts = parts[1:]
                    if not keep_taxon:
                        for i, val_str in enumerate(value_parts):
                            if val_str.strip() in missing_markers:
                                continue
                            try:
                                float(val_str)
                            except ValueError:
                                raise PhykitUserError(
                                    [
                                        f"Non-numeric trait value '{val_str}' for taxon '{taxon}' "
                                        f"(trait '{trait_names[i]}') on line {data_line_idx}.",
                                    ],
                                    code=2,
                                )
                        data_line_idx += 1
                        continue
                    try:
                        values = [float(val_str) for val_str in value_parts]
                    except ValueError:
                        values = []
                        for i, val_str in enumerate(value_parts):
                            if val_str.strip() in missing_markers:
                                values.append(float("nan"))
                                missing_info.append({
                                    "taxon": taxon,
                                    "trait": trait_names[i],
                                    "line": data_line_idx,
                                })
                                continue
                            try:
                                value = float(val_str)
                            except ValueError:
                                raise PhykitUserError(
                                    [
                                        f"Non-numeric trait value '{val_str}' for taxon '{taxon}' "
                                        f"(trait '{trait_names[i]}') on line {data_line_idx}.",
                                    ],
                                    code=2,
                                )
                            values.append(value)
                    traits[taxon] = values
                    data_line_idx += 1

                if not saw_data:
                    raise PhykitUserError(
                        [
                            "Multi-trait file must have a header row and at least one data row.",
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

        if (
            len(tree_tip_set) >= 3
            and len(tree_tip_set) == len(trait_taxa_set)
            and tree_tip_set == trait_taxa_set
        ):
            return trait_names, traits, missing_info

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
        return trait_names, filtered, missing_info

    @staticmethod
    def _subset_to_shared_taxa(ordered_names, Y, shared):
        shared_set = set(shared)
        if shared_set == set(ordered_names):
            return ordered_names, Y

        keep_idx = [
            index for index, name in enumerate(ordered_names) if name in shared_set
        ]
        return [ordered_names[index] for index in keep_idx], Y[keep_idx, :]

    @staticmethod
    def _estimate_complete_case_stats(
        Y_complete: np.ndarray,
        C_complete: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray]:
        try:
            return PhyloImpute._estimate_complete_case_stats_cholesky(
                Y_complete, C_complete
            )
        except (np.linalg.LinAlgError, FloatingPointError, ValueError):
            return PhyloImpute._estimate_complete_case_stats_inverse(
                Y_complete, C_complete
            )

    @staticmethod
    def _estimate_complete_case_stats_cholesky(
        Y_complete: np.ndarray,
        C_complete: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray]:
        n_complete, p = Y_complete.shape
        if n_complete <= 1:
            return np.mean(Y_complete, axis=0), np.eye(p)

        factor = cho_factor(C_complete, lower=True, check_finite=False)
        ones = np.ones(n_complete)
        solve_rhs = np.empty(
            (n_complete, p + 1),
            dtype=np.result_type(Y_complete, C_complete),
        )
        solve_rhs[:, 0] = 1.0
        solve_rhs[:, 1:] = Y_complete
        solved = cho_solve(factor, solve_rhs, check_finite=False)
        C_inv_ones = solved[:, 0]
        C_inv_Y = solved[:, 1:]

        denom = ones @ C_inv_ones
        a_hat = (ones @ C_inv_Y) / denom

        E = Y_complete - a_hat
        C_inv_E = C_inv_Y - C_inv_ones[:, None] * a_hat
        Sigma_trait = (E.T @ C_inv_E) / (n_complete - 1)
        return a_hat, Sigma_trait

    @staticmethod
    def _estimate_complete_case_stats_inverse(
        Y_complete: np.ndarray,
        C_complete: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray]:
        n_complete, p = Y_complete.shape

        try:
            C_complete_inv = np.linalg.inv(C_complete)
        except np.linalg.LinAlgError:
            C_complete_inv = np.linalg.pinv(C_complete)

        ones = np.ones(n_complete)
        C_inv_ones = C_complete_inv @ ones
        C_inv_Y = C_complete_inv @ Y_complete
        denom = ones @ C_inv_ones
        a_hat = (ones @ C_inv_Y) / denom

        E = Y_complete - a_hat

        if n_complete > 1:
            C_inv_E = C_inv_Y - C_inv_ones[:, None] * a_hat
            Sigma_trait = (E.T @ C_inv_E) / (n_complete - 1)
        else:
            Sigma_trait = np.eye(p)

        return a_hat, Sigma_trait

    @staticmethod
    def _impute_taxon(
        Y: np.ndarray,
        C_phylo: np.ndarray,
        taxon_idx: int,
        observed_trait_indices: list[int],
        missing_trait_indices: list[int],
        a_hat: np.ndarray,
        Sigma_trait: np.ndarray,
    ) -> tuple[dict[int, float], dict[int, float]]:
        """Impute missing traits for one taxon using phylogeny + trait correlations."""
        try:
            return PhyloImpute._impute_taxon_cholesky(
                Y,
                C_phylo,
                taxon_idx,
                observed_trait_indices,
                missing_trait_indices,
                a_hat,
                Sigma_trait,
            )
        except (np.linalg.LinAlgError, FloatingPointError, ValueError):
            return PhyloImpute._impute_taxon_inverse(
                Y,
                C_phylo,
                taxon_idx,
                observed_trait_indices,
                missing_trait_indices,
                a_hat,
                Sigma_trait,
            )

    @staticmethod
    def _impute_taxon_cholesky(
        Y: np.ndarray,
        C_phylo: np.ndarray,
        taxon_idx: int,
        observed_trait_indices: list[int],
        missing_trait_indices: list[int],
        a_hat: np.ndarray,
        Sigma_trait: np.ndarray,
    ) -> tuple[dict[int, float], dict[int, float]]:
        n = C_phylo.shape[0]

        other_idx_array = PhyloImpute._other_taxon_indices(n, taxon_idx)

        C_io = C_phylo[taxon_idx, other_idx_array]
        C_oo = C_phylo[np.ix_(other_idx_array, other_idx_array)]

        imputed = {}
        imputed_se = {}
        trait_context = None
        trait_context_ready = False

        for j in missing_trait_indices:
            # Get observed values of trait j from other taxa
            obs_mask = ~np.isnan(Y[other_idx_array, j])

            if not np.any(obs_mask):
                # No other taxa have this trait -- use phylogenetic mean
                imputed[j] = float(a_hat[j])
                imputed_se[j] = float(np.sqrt(max(Sigma_trait[j, j], 0.0)))
                continue

            obs_k_indices = np.flatnonzero(obs_mask)
            obs_values = Y[other_idx_array[obs_k_indices], j]

            # Phylogenetic prediction
            c_io_sub = C_io[obs_k_indices]
            C_oo_sub = C_oo[np.ix_(obs_k_indices, obs_k_indices)]
            phylo_factor = cho_factor(C_oo_sub, lower=True, check_finite=False)

            phylo_rhs = np.empty(
                (obs_values.size, 2),
                dtype=np.result_type(obs_values, c_io_sub, C_oo_sub),
            )
            phylo_rhs[:, 0] = obs_values - a_hat[j]
            phylo_rhs[:, 1] = c_io_sub
            phylo_solved = cho_solve(phylo_factor, phylo_rhs, check_finite=False)
            phylo_delta = phylo_solved[:, 0]
            weights_for_var = phylo_solved[:, 1]
            y_hat = float(a_hat[j] + c_io_sub @ phylo_delta)

            # Conditional variance
            c_ii = C_phylo[taxon_idx, taxon_idx]
            cond_var = (c_ii - c_io_sub @ weights_for_var) * Sigma_trait[j, j]
            se = float(np.sqrt(max(cond_var, 0.0)))

            # If taxon has observed traits that correlate with trait j,
            # refine the prediction using trait correlations
            if len(observed_trait_indices) > 0:
                if not trait_context_ready:
                    obs_traits_values = Y[taxon_idx, observed_trait_indices]
                    obs_means = a_hat[observed_trait_indices]
                    Sigma_oo_traits = Sigma_trait[
                        np.ix_(observed_trait_indices, observed_trait_indices)
                    ]
                    trait_factor = cho_factor(
                        Sigma_oo_traits, lower=True, check_finite=False
                    )
                    trait_delta = cho_solve(
                        trait_factor,
                        obs_traits_values - obs_means,
                        check_finite=False,
                    )
                    trait_context = (trait_factor, trait_delta)
                    trait_context_ready = True

                trait_factor, trait_delta = trait_context
                Sigma_jo = Sigma_trait[j, observed_trait_indices]
                trait_adjustment = float(Sigma_jo @ trait_delta)
                y_hat += trait_adjustment

                sigma_jo_solved = cho_solve(
                    trait_factor,
                    Sigma_jo,
                    check_finite=False,
                )
                cond_var_trait = Sigma_trait[j, j] - Sigma_jo @ sigma_jo_solved
                se = float(
                    np.sqrt(max(cond_var * cond_var_trait / Sigma_trait[j, j], 0.0))
                )

            imputed[j] = float(y_hat)
            imputed_se[j] = float(se)

        return imputed, imputed_se

    @staticmethod
    def _impute_taxon_inverse(
        Y: np.ndarray,
        C_phylo: np.ndarray,
        taxon_idx: int,
        observed_trait_indices: list[int],
        missing_trait_indices: list[int],
        a_hat: np.ndarray,
        Sigma_trait: np.ndarray,
    ) -> tuple[dict[int, float], dict[int, float]]:
        n = C_phylo.shape[0]

        other_idx_array = PhyloImpute._other_taxon_indices(n, taxon_idx)

        C_io = C_phylo[taxon_idx, other_idx_array]
        C_oo = C_phylo[np.ix_(other_idx_array, other_idx_array)]

        imputed = {}
        imputed_se = {}
        trait_context = None
        trait_context_ready = False

        for j in missing_trait_indices:
            obs_mask = ~np.isnan(Y[other_idx_array, j])

            if not np.any(obs_mask):
                imputed[j] = float(a_hat[j])
                imputed_se[j] = float(np.sqrt(max(Sigma_trait[j, j], 0.0)))
                continue

            obs_k_indices = np.flatnonzero(obs_mask)
            obs_values = Y[other_idx_array[obs_k_indices], j]

            c_io_sub = C_io[obs_k_indices]
            C_oo_sub = C_oo[np.ix_(obs_k_indices, obs_k_indices)]
            try:
                C_oo_sub_inv = np.linalg.inv(C_oo_sub)
            except np.linalg.LinAlgError:
                C_oo_sub_inv = np.linalg.pinv(C_oo_sub)

            weights = c_io_sub @ C_oo_sub_inv
            y_hat = float(a_hat[j] + weights @ (obs_values - a_hat[j]))

            c_ii = C_phylo[taxon_idx, taxon_idx]
            cond_var = (c_ii - c_io_sub @ C_oo_sub_inv @ c_io_sub) * Sigma_trait[j, j]
            se = float(np.sqrt(max(cond_var, 0.0)))

            if len(observed_trait_indices) > 0:
                if not trait_context_ready:
                    obs_traits_values = Y[taxon_idx, observed_trait_indices]
                    obs_means = a_hat[observed_trait_indices]
                    Sigma_oo_traits = Sigma_trait[
                        np.ix_(observed_trait_indices, observed_trait_indices)
                    ]
                    try:
                        Sigma_oo_inv = np.linalg.inv(Sigma_oo_traits)
                        trait_delta = Sigma_oo_inv @ (obs_traits_values - obs_means)
                        trait_context = (Sigma_oo_inv, trait_delta)
                    except np.linalg.LinAlgError:
                        trait_context = None
                    trait_context_ready = True

                if trait_context is not None:
                    Sigma_oo_inv, trait_delta = trait_context
                    Sigma_jo = Sigma_trait[j, observed_trait_indices]
                    y_hat += float(Sigma_jo @ trait_delta)

                    cond_var_trait = (
                        Sigma_trait[j, j] - Sigma_jo @ Sigma_oo_inv @ Sigma_jo
                    )
                    se = float(
                        np.sqrt(max(cond_var * cond_var_trait / Sigma_trait[j, j], 0.0))
                    )

            imputed[j] = float(y_hat)
            imputed_se[j] = float(se)

        return imputed, imputed_se

    @staticmethod
    def _other_taxon_indices(n_taxa: int, taxon_idx: int) -> np.ndarray:
        if taxon_idx == 0:
            return np.arange(1, n_taxa, dtype=np.intp)
        if taxon_idx == n_taxa - 1:
            return np.arange(0, n_taxa - 1, dtype=np.intp)
        return np.concatenate(
            (
                np.arange(taxon_idx, dtype=np.intp),
                np.arange(taxon_idx + 1, n_taxa, dtype=np.intp),
            )
        )

    @staticmethod
    def _write_output_tsv(
        path: str,
        trait_names: list[str],
        ordered_names: list[str],
        Y: np.ndarray,
    ) -> None:
        """Write imputed data as TSV (drop-in replacement for the input)."""
        n_traits = len(trait_names)
        fmt = ["%s"] + ["%.6f"] * n_traits
        chunk_size = 1000
        with open(path, "w") as f:
            header = "taxon\t" + "\t".join(trait_names) + "\n"
            f.write(header)
            for start in range(0, len(ordered_names), chunk_size):
                stop = min(start + chunk_size, len(ordered_names))
                block = np.empty((stop - start, n_traits + 1), dtype=object)
                block[:, 0] = ordered_names[start:stop]
                block[:, 1:] = Y[start:stop]
                np.savetxt(f, block, delimiter="\t", fmt=fmt)

    def _print_text(
        self, n, p, trait_names, vcv_type, imputed_results, data_path, output_path
    ) -> None:
        vcv_label = "standard" if vcv_type == "BM" else "discordance-aware"
        lines = [
            "Phylogenetic Imputation",
            f"Data: {data_path}",
            f"Taxa: {n}",
            f"Traits: {p} ({', '.join(trait_names)})",
            f"Missing values: {len(imputed_results)}",
            f"VCV: {vcv_type} ({vcv_label})",
        ]

        if imputed_results:
            # Determine column widths
            max_taxon = len("Taxon")
            max_trait = len("Trait")
            for r in imputed_results:
                taxon_len = len(r["taxon"])
                if taxon_len > max_taxon:
                    max_taxon = taxon_len
                trait_len = len(r["trait"])
                if trait_len > max_trait:
                    max_trait = trait_len

            lines.append("\nImputed values:")
            header = (
                f"  {'Taxon':<{max_taxon}}  "
                f"{'Trait':<{max_trait}}  "
                f"{'Imputed':>10}  "
                f"{'SE':>10}  "
                f"{'95% CI':>20}"
            )
            lines.append(header)

            for r in imputed_results:
                ci_str = f"[{r['ci_lower']:.4f}, {r['ci_upper']:.4f}]"
                row = (
                    f"  {r['taxon']:<{max_taxon}}  "
                    f"  {r['trait']:<{max_trait}}  "
                    f"{r['value']:>10.4f}  "
                    f"{r['se']:>10.4f}  "
                    f"{ci_str:>20}"
                )
                lines.append(row)

        lines.append(f"\nOutput: {output_path}")
        print("\n".join(lines))

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
