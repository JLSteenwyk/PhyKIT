import copy
import math
import sys
from typing import Dict, List, Tuple

import numpy as np

from .base import Tree
from ...helpers.json_output import print_json
from ...errors import PhykitUserError


class AncestralReconstruction(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"])
        self.trait_data_path = parsed["trait_data_path"]
        self.trait_column = parsed["trait_column"]
        self.method = parsed["method"]
        self.ci = parsed["ci"]
        self.plot_output = parsed["plot_output"]
        self.json_output = parsed["json_output"]

    def run(self) -> None:
        tree = self.read_tree_file()
        self._validate_tree(tree)

        tree_tips = self.get_tip_names_from_tree(tree)

        if self.trait_column is not None:
            trait_values = self._parse_multi_trait_data(
                self.trait_data_path, tree_tips, self.trait_column
            )
            trait_name = self.trait_column
        else:
            trait_values = self._parse_single_trait_data(
                self.trait_data_path, tree_tips
            )
            trait_name = "trait"

        ordered_names = sorted(trait_values.keys())
        n = len(ordered_names)
        x = np.array([trait_values[name] for name in ordered_names])

        # Prune tree to shared taxa
        tree_copy = copy.deepcopy(tree)
        tip_names_in_tree = [t.name for t in tree_copy.get_terminals()]
        tips_to_prune = [t for t in tip_names_in_tree if t not in trait_values]
        if tips_to_prune:
            tree_copy = self.prune_tree_using_taxa_list(tree_copy, tips_to_prune)

        # Label internal nodes
        node_labels = self._label_internal_nodes(tree_copy)

        if self.method == "ml":
            (
                node_estimates, node_cis, sigma2, log_likelihood
            ) = self._anc_ml(tree_copy, x, ordered_names, node_labels)
        else:
            (
                node_estimates, node_cis, sigma2, log_likelihood
            ) = self._fast_anc(tree_copy, x, ordered_names, node_labels)

        # Build result
        result = self._format_result(
            method=self.method,
            trait_name=trait_name,
            n_tips=n,
            sigma2=sigma2,
            log_likelihood=log_likelihood,
            node_estimates=node_estimates,
            node_cis=node_cis,
            node_labels=node_labels,
            tree=tree_copy,
            trait_values=trait_values,
        )

        if self.plot_output:
            self._plot_contmap(
                tree_copy, node_estimates, node_labels,
                trait_values, trait_name, self.plot_output
            )
            result["plot_output"] = self.plot_output

        if self.json_output:
            print_json(result)
        else:
            self._print_text_output(
                method=self.method,
                trait_name=trait_name,
                n_tips=n,
                sigma2=sigma2,
                log_likelihood=log_likelihood,
                node_estimates=node_estimates,
                node_cis=node_cis,
                node_labels=node_labels,
                tree=tree_copy,
            )

    def process_args(self, args) -> Dict:
        return dict(
            tree_file_path=args.tree,
            trait_data_path=args.trait_data,
            trait_column=getattr(args, "trait", None),
            method=getattr(args, "method", "fast"),
            ci=getattr(args, "ci", False),
            plot_output=getattr(args, "plot", None),
            json_output=getattr(args, "json", False),
        )

    def _validate_tree(self, tree) -> None:
        tips = list(tree.get_terminals())
        if len(tips) < 3:
            raise PhykitUserError(
                ["Tree must have at least 3 tips for ancestral reconstruction."],
                code=2,
            )
        for clade in tree.find_clades():
            if clade.branch_length is None and clade != tree.root:
                raise PhykitUserError(
                    ["All branches in the tree must have lengths."],
                    code=2,
                )

    def _parse_single_trait_data(
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

    def _parse_multi_trait_data(
        self, path: str, tree_tips: List[str], trait_column: str
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
        if len(header_parts) < 2:
            raise PhykitUserError(
                [
                    "Header must have at least 2 columns (taxon + at least 1 trait).",
                ],
                code=2,
            )
        trait_names = header_parts[1:]

        if trait_column not in trait_names:
            raise PhykitUserError(
                [
                    f"Column '{trait_column}' not found in trait file.",
                    f"Available columns: {', '.join(trait_names)}",
                ],
                code=2,
            )

        col_idx = trait_names.index(trait_column)

        traits = {}
        for line_idx, line in enumerate(data_lines[1:], 2):
            parts = line.split("\t")
            if len(parts) != len(header_parts):
                raise PhykitUserError(
                    [
                        f"Line {line_idx} has {len(parts)} columns; expected {len(header_parts)}.",
                    ],
                    code=2,
                )
            taxon = parts[0]
            val_str = parts[1 + col_idx]
            try:
                traits[taxon] = float(val_str)
            except ValueError:
                raise PhykitUserError(
                    [
                        f"Non-numeric trait value '{val_str}' for taxon '{taxon}' "
                        f"(trait '{trait_column}') on line {line_idx}.",
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

    def _label_internal_nodes(self, tree) -> Dict:
        """Assign N1, N2, ... labels to unnamed internal nodes.

        Returns dict mapping id(clade) -> label for all internal nodes.
        """
        labels = {}
        counter = 1
        for clade in tree.find_clades(order="preorder"):
            if not clade.is_terminal():
                if clade.name:
                    labels[id(clade)] = clade.name
                else:
                    labels[id(clade)] = f"N{counter}"
                    counter += 1
        return labels

    def _get_descendant_tips(self, tree, node) -> List[str]:
        tips = []
        for tip in node.get_terminals():
            tips.append(tip.name)
        return sorted(tips)

    def _build_parent_map(self, tree) -> Dict:
        parent_map = {}
        for clade in tree.find_clades(order="preorder"):
            for child in clade.clades:
                parent_map[id(child)] = clade
        return parent_map

    def _fast_anc(
        self, tree, x: np.ndarray, ordered_names: List[str],
        node_labels: Dict,
    ) -> Tuple[Dict, Dict, float, float]:
        """Two-pass ML ancestral reconstruction (Felsenstein's algorithm).

        Pass 1 (postorder): compute subtree estimates and variances.
        Pass 2 (preorder): incorporate information from the rest of the
        tree so that every internal node gets the full-tree ML estimate.
        """
        name_to_idx = {name: i for i, name in enumerate(ordered_names)}
        parent_map = self._build_parent_map(tree)

        # --- Pass 1: postorder (subtree estimates) ---
        est_down = {}
        var_down = {}

        for clade in tree.find_clades(order="postorder"):
            if clade.is_terminal():
                name = clade.name
                if name in name_to_idx:
                    est_down[id(clade)] = x[name_to_idx[name]]
                    var_down[id(clade)] = 0.0
            else:
                children = clade.clades
                precisions = []
                weighted_vals = []
                for child in children:
                    child_id = id(child)
                    if child_id not in est_down:
                        continue
                    v_i = child.branch_length if child.branch_length else 0.0
                    child_var = var_down[child_id]
                    denom = v_i + child_var
                    if denom == 0:
                        denom = 1e-10
                    prec = 1.0 / denom
                    precisions.append(prec)
                    weighted_vals.append(prec * est_down[child_id])

                if precisions:
                    total_prec = sum(precisions)
                    est_down[id(clade)] = sum(weighted_vals) / total_prec
                    var_down[id(clade)] = 1.0 / total_prec

        # --- Pass 2: preorder (incorporate full-tree information) ---
        # For each node, compute the "message from above": an estimate
        # and variance representing all information NOT in the subtree
        # below that node.
        est_up = {}
        var_up = {}
        est_final = {}
        var_final = {}

        root = tree.root
        # Root has no information from above
        est_final[id(root)] = est_down[id(root)]
        var_final[id(root)] = var_down[id(root)]

        for clade in tree.find_clades(order="preorder"):
            if clade == root:
                continue
            if id(clade) not in parent_map:
                continue
            parent = parent_map[id(clade)]
            bl = clade.branch_length if clade.branch_length else 0.0

            # Gather precision from parent's "above" message
            prec_above = 0.0
            est_above_weighted = 0.0
            if id(parent) in est_up:
                v_up_parent = var_up[id(parent)]
                if v_up_parent > 0:
                    p = 1.0 / v_up_parent
                    prec_above += p
                    est_above_weighted += p * est_up[id(parent)]
                elif v_up_parent == 0 and id(parent) != id(root):
                    p = 1.0 / 1e-10
                    prec_above += p
                    est_above_weighted += p * est_up[id(parent)]

            # Gather precision from siblings
            for sibling in parent.clades:
                if id(sibling) == id(clade):
                    continue
                if id(sibling) not in est_down:
                    continue
                sib_bl = sibling.branch_length if sibling.branch_length else 0.0
                sib_var = var_down[id(sibling)] + sib_bl
                if sib_var == 0:
                    sib_var = 1e-10
                p = 1.0 / sib_var
                prec_above += p
                est_above_weighted += p * est_down[id(sibling)]

            # The "above" message for this node
            if prec_above > 0:
                est_up[id(clade)] = est_above_weighted / prec_above
                var_up[id(clade)] = bl + 1.0 / prec_above
            else:
                # No information from above (shouldn't happen in a valid tree)
                est_up[id(clade)] = est_down.get(id(clade), 0.0)
                var_up[id(clade)] = bl + 1e10

            # Combine down and up for final estimate
            if id(clade) in est_down:
                vd = var_down[id(clade)]
                vu = var_up[id(clade)]
                if vd == 0:
                    vd = 1e-10
                if vu == 0:
                    vu = 1e-10
                prec_d = 1.0 / vd
                prec_u = 1.0 / vu
                total = prec_d + prec_u
                est_final[id(clade)] = (
                    prec_d * est_down[id(clade)] + prec_u * est_up[id(clade)]
                ) / total
                var_final[id(clade)] = 1.0 / total

        # Compute sigma^2 from phylogenetic independent contrasts
        sigma2 = self._compute_sigma2_from_contrasts(
            tree, x, est_down, var_down, ordered_names
        )

        # Compute log-likelihood
        log_likelihood = self._compute_log_likelihood(
            tree, x, ordered_names, sigma2
        )

        # Build output dicts keyed by node label
        node_estimates = {}
        node_cis = {}
        for clade in tree.find_clades(order="preorder"):
            if not clade.is_terminal() and id(clade) in node_labels:
                label = node_labels[id(clade)]
                if id(clade) in est_final:
                    est = est_final[id(clade)]
                    node_estimates[label] = est
                    if self.ci and id(clade) in var_final:
                        var = var_final[id(clade)]
                        se = math.sqrt(sigma2 * var) if sigma2 * var > 0 else 0.0
                        node_cis[label] = (est - 1.96 * se, est + 1.96 * se)

        return node_estimates, node_cis, sigma2, log_likelihood

    def _compute_sigma2_from_contrasts(
        self, tree, x: np.ndarray, node_estimates: Dict,
        node_variances: Dict, ordered_names: List[str],
    ) -> float:
        """Compute BM rate (sigma^2) from phylogenetic independent contrasts.

        For nodes with >2 children (polytomies), iteratively compute
        pairwise contrasts and update the combined estimate/variance,
        producing one contrast per additional child (n_children - 1
        contrasts per node). This ensures the total number of contrasts
        equals n_tips - 1 for a fully resolved tree.
        """
        contrasts_sq = []
        for clade in tree.find_clades(order="postorder"):
            if clade.is_terminal():
                continue
            children = clade.clades
            valid_children = [
                c for c in children if id(c) in node_estimates
            ]
            if len(valid_children) < 2:
                continue

            # Process first two children
            c1, c2 = valid_children[0], valid_children[1]
            v1 = (c1.branch_length or 0.0) + node_variances[id(c1)]
            v2 = (c2.branch_length or 0.0) + node_variances[id(c2)]
            denom = v1 + v2
            if denom > 0:
                contrast = (
                    node_estimates[id(c1)] - node_estimates[id(c2)]
                ) / math.sqrt(denom)
                contrasts_sq.append(contrast ** 2)

            # Combined estimate and variance after first two children
            if denom > 0:
                combined_est = (
                    v2 * node_estimates[id(c1)] + v1 * node_estimates[id(c2)]
                ) / denom
                combined_var = (v1 * v2) / denom
            else:
                combined_est = node_estimates[id(c1)]
                combined_var = 0.0

            # Process additional children (for polytomies)
            for k in range(2, len(valid_children)):
                ck = valid_children[k]
                vk = (ck.branch_length or 0.0) + node_variances[id(ck)]
                vc = combined_var
                denom_k = vk + vc
                if denom_k > 0:
                    contrast_k = (
                        node_estimates[id(ck)] - combined_est
                    ) / math.sqrt(denom_k)
                    contrasts_sq.append(contrast_k ** 2)
                    combined_est = (
                        vc * node_estimates[id(ck)] + vk * combined_est
                    ) / denom_k
                    combined_var = (vk * vc) / denom_k

        if contrasts_sq:
            return sum(contrasts_sq) / len(contrasts_sq)
        return 0.0

    def _compute_log_likelihood(
        self, tree, x: np.ndarray, ordered_names: List[str], sigma2: float,
    ) -> float:
        """Compute log-likelihood under BM model."""
        n = len(ordered_names)
        if sigma2 <= 0 or n == 0:
            return float("-inf")

        vcv = self._build_vcv_matrix(tree, ordered_names)
        C = sigma2 * vcv

        try:
            sign, logdet = np.linalg.slogdet(C)
            if sign <= 0:
                return float("-inf")
            C_inv = np.linalg.inv(C)
        except np.linalg.LinAlgError:
            return float("-inf")

        ones = np.ones(n)
        a_hat = float(ones @ C_inv @ x) / float(ones @ C_inv @ ones)
        residuals = x - a_hat

        ll = -0.5 * (
            n * np.log(2 * np.pi) + logdet + float(residuals @ C_inv @ residuals)
        )
        return float(ll)

    def _build_vcv_matrix(
        self, tree, ordered_names: List[str]
    ) -> np.ndarray:
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

    def _anc_ml(
        self, tree, x: np.ndarray, ordered_names: List[str],
        node_labels: Dict,
    ) -> Tuple[Dict, Dict, float, float]:
        """VCV-based ML ancestral reconstruction with exact conditional CIs."""
        n = len(ordered_names)
        vcv = self._build_vcv_matrix(tree, ordered_names)

        try:
            C_inv = np.linalg.inv(vcv)
        except np.linalg.LinAlgError:
            raise PhykitUserError(
                [
                    "Singular VCV matrix: cannot invert.",
                    "Check that the tree has valid branch lengths.",
                ],
                code=2,
            )

        ones = np.ones(n)
        # ML root estimate
        a_hat = float(ones @ C_inv @ x) / float(ones @ C_inv @ ones)
        # ML rate (sigma^2)
        residuals = x - a_hat
        sigma2_ml = float(residuals @ C_inv @ residuals) / n
        # REML rate for CIs (matches R's fastAnc)
        sigma2_reml = float(residuals @ C_inv @ residuals) / (n - 1)

        # Log-likelihood (uses ML sigma^2)
        sign, logdet = np.linalg.slogdet(vcv)
        ll = -0.5 * (
            n * np.log(2 * np.pi) + n * np.log(sigma2_ml) + logdet + n
        )

        # Variance of the root estimate
        root_var = 1.0 / float(ones @ C_inv @ ones)

        # Compute cross-covariance and conditional estimates for each internal node
        node_estimates = {}
        node_cis = {}

        for clade in tree.find_clades(order="preorder"):
            if clade.is_terminal():
                continue
            if id(clade) not in node_labels:
                continue

            label = node_labels[id(clade)]
            d_i = tree.distance(tree.root, clade)

            # Cross-covariance: c_i[j] = distance from root to MRCA(node_i, tip_j)
            c_i = np.zeros(n)
            for j, tip_name in enumerate(ordered_names):
                tip_clade = None
                for t in tree.get_terminals():
                    if t.name == tip_name:
                        tip_clade = t
                        break
                mrca = tree.common_ancestor(clade, tip_clade)
                c_i[j] = tree.distance(tree.root, mrca)

            # Conditional mean: E[x_i | x] = a_hat + c_i' C^{-1} (x - a_hat)
            cond_mean = a_hat + float(c_i @ C_inv @ residuals)
            node_estimates[label] = cond_mean

            if self.ci:
                if clade == tree.root:
                    # Root variance from GLS estimation
                    cond_var = sigma2_reml * root_var
                else:
                    # Conditional variance using REML sigma^2
                    cond_var = sigma2_reml * (d_i - float(c_i @ C_inv @ c_i))
                if cond_var < 0:
                    cond_var = 0.0
                se = math.sqrt(cond_var)
                node_cis[label] = (cond_mean - 1.96 * se, cond_mean + 1.96 * se)

        return node_estimates, node_cis, sigma2_reml, float(ll)

    def _format_result(
        self, *, method, trait_name, n_tips, sigma2, log_likelihood,
        node_estimates, node_cis, node_labels, tree, trait_values,
    ) -> Dict:
        ancestral_estimates = {}
        root_id = id(tree.root)

        for clade in tree.find_clades(order="preorder"):
            if clade.is_terminal():
                continue
            if id(clade) not in node_labels:
                continue
            label = node_labels[id(clade)]
            if label not in node_estimates:
                continue

            entry = {
                "estimate": float(node_estimates[label]),
                "descendants": self._get_descendant_tips(tree, clade),
                "is_root": id(clade) == root_id,
            }
            if label in node_cis:
                entry["ci_lower"] = float(node_cis[label][0])
                entry["ci_upper"] = float(node_cis[label][1])
            ancestral_estimates[label] = entry

        result = {
            "method": method,
            "trait": trait_name,
            "n_tips": n_tips,
            "log_likelihood": float(log_likelihood),
            "sigma2": float(sigma2),
            "ancestral_estimates": ancestral_estimates,
            "tip_values": {k: float(v) for k, v in sorted(trait_values.items())},
        }
        return result

    def _print_text_output(
        self, *, method, trait_name, n_tips, sigma2, log_likelihood,
        node_estimates, node_cis, node_labels, tree,
    ) -> None:
        print("Ancestral State Reconstruction")

        method_desc = (
            "fast (Felsenstein's contrasts)" if method == "fast"
            else "ml (maximum likelihood, VCV-based)"
        )
        print(f"\nMethod: {method_desc}")
        print(f"Trait: {trait_name}")
        print(f"Number of tips: {n_tips}")
        print(f"\nLog-likelihood: {log_likelihood:.4f}")
        print(f"Sigma-squared (BM rate): {sigma2:.6f}")

        print("\nAncestral estimates:")
        root_id = id(tree.root)

        if node_cis:
            print(
                f"  {'Node':<12s}{'Descendants':>12s}{'Estimate':>12s}{'95% CI':>24s}"
            )
            for clade in tree.find_clades(order="preorder"):
                if clade.is_terminal():
                    continue
                if id(clade) not in node_labels:
                    continue
                label = node_labels[id(clade)]
                if label not in node_estimates:
                    continue
                n_desc = len(self._get_descendant_tips(tree, clade))
                est = node_estimates[label]
                root_tag = " (root)" if id(clade) == root_id else ""
                if label in node_cis:
                    ci_lo, ci_hi = node_cis[label]
                    ci_str = f"[{ci_lo:.4f}, {ci_hi:.4f}]"
                else:
                    ci_str = ""
                print(
                    f"  {label + root_tag:<12s}{n_desc:>12d}{est:>12.4f}{ci_str:>24s}"
                )
        else:
            print(
                f"  {'Node':<12s}{'Descendants':>12s}{'Estimate':>12s}"
            )
            for clade in tree.find_clades(order="preorder"):
                if clade.is_terminal():
                    continue
                if id(clade) not in node_labels:
                    continue
                label = node_labels[id(clade)]
                if label not in node_estimates:
                    continue
                n_desc = len(self._get_descendant_tips(tree, clade))
                est = node_estimates[label]
                root_tag = " (root)" if id(clade) == root_id else ""
                print(
                    f"  {label + root_tag:<12s}{n_desc:>12d}{est:>12.4f}"
                )

    def _plot_contmap(
        self, tree, node_estimates, node_labels, trait_values,
        trait_name, output_path,
    ) -> None:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            from matplotlib.colors import Normalize
        except ImportError:
            print(
                "matplotlib is required for contMap plotting. "
                "Install matplotlib and retry."
            )
            raise SystemExit(2)

        parent_map = self._build_parent_map(tree)

        # Build estimates dict keyed by id(clade) for all nodes
        all_estimates = {}
        for clade in tree.find_clades(order="preorder"):
            if clade.is_terminal():
                if clade.name in trait_values:
                    all_estimates[id(clade)] = trait_values[clade.name]
            else:
                if id(clade) in node_labels:
                    label = node_labels[id(clade)]
                    if label in node_estimates:
                        all_estimates[id(clade)] = node_estimates[label]

        tips = list(tree.get_terminals())
        node_x = {}
        node_y = {}

        # Assign tip y-positions (evenly spaced)
        for i, tip in enumerate(tips):
            node_y[id(tip)] = i

        # Assign x-positions via preorder traversal
        root = tree.root
        for clade in tree.find_clades(order="preorder"):
            if clade == root:
                node_x[id(clade)] = 0.0
            else:
                if id(clade) in parent_map:
                    parent = parent_map[id(clade)]
                    t = clade.branch_length if clade.branch_length else 0.0
                    node_x[id(clade)] = node_x[id(parent)] + t

        # Assign internal y-positions (mean of children)
        for clade in tree.find_clades(order="postorder"):
            if not clade.is_terminal() and id(clade) not in node_y:
                child_ys = [
                    node_y[id(c)] for c in clade.clades if id(c) in node_y
                ]
                if child_ys:
                    node_y[id(clade)] = np.mean(child_ys)
                else:
                    node_y[id(clade)] = 0.0

        # Normalize trait values for colormap
        all_vals = list(all_estimates.values())
        if not all_vals:
            return
        vmin, vmax = min(all_vals), max(all_vals)
        if vmin == vmax:
            vmax = vmin + 1.0
        norm = Normalize(vmin=vmin, vmax=vmax)
        cmap = plt.get_cmap("coolwarm")

        fig, ax = plt.subplots(figsize=(10, max(4, len(tips) * 0.4)))

        n_seg = 50

        for clade in tree.find_clades(order="preorder"):
            if clade == root:
                continue
            if id(clade) not in parent_map:
                continue

            parent = parent_map[id(clade)]
            if id(parent) not in node_x or id(clade) not in node_x:
                continue
            if id(parent) not in all_estimates or id(clade) not in all_estimates:
                continue

            x0 = node_x[id(parent)]
            x1 = node_x[id(clade)]
            y0 = node_y[id(parent)]
            y1 = node_y[id(clade)]
            val0 = all_estimates[id(parent)]
            val1 = all_estimates[id(clade)]

            # Horizontal segment (at child's y): colored gradient
            for s in range(n_seg):
                frac_s = s / n_seg
                frac_e = (s + 1) / n_seg
                xs = x0 + frac_s * (x1 - x0)
                xe = x0 + frac_e * (x1 - x0)
                v = val0 + ((frac_s + frac_e) / 2) * (val1 - val0)
                ax.plot(
                    [xs, xe], [y1, y1],
                    color=cmap(norm(v)), lw=3, solid_capstyle="butt",
                )

            # Vertical connector (at parent's x): parent's trait color
            ax.plot(
                [x0, x0], [y0, y1],
                color=cmap(norm(val0)), lw=3, solid_capstyle="butt",
            )

        # Tip labels
        max_x = max(node_x.values()) if node_x else 0
        offset = max_x * 0.02
        for tip in tips:
            ax.text(
                node_x[id(tip)] + offset, node_y[id(tip)],
                tip.name, va="center", fontsize=9,
            )

        # Colorbar
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax, pad=0.15)
        cbar.set_label(trait_name)

        ax.set_xlabel("Branch length")
        ax.set_yticks([])
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)
        ax.set_title("Continuous Trait Map (contMap)")
        fig.tight_layout()
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved contMap plot: {output_path}")
