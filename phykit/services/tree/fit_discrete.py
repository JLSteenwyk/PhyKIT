"""
Compare models of discrete trait evolution on a phylogeny.

Fits ER (Equal Rates), SYM (Symmetric), and ARD (All Rates Different)
Mk models and compares them using AIC/BIC. Analogous to R's
geiger::fitDiscrete().
"""
from __future__ import annotations

import math

from .base import Tree
from ...errors import PhykitUserError


VALID_DISCRETE_MODELS = frozenset(["ER", "SYM", "ARD"])
_PREPARE_FELSENSTEIN_CONTEXT = None
_COUNT_PARAMS = None
_FIT_Q_MATRIX = None
_PARSE_DISCRETE_TRAITS = None


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


def _prepare_felsenstein_context(*args, **kwargs):
    global _PREPARE_FELSENSTEIN_CONTEXT

    if _PREPARE_FELSENSTEIN_CONTEXT is None:
        from ...helpers.discrete_models import (
            _prepare_felsenstein_context as _prepare_context,
        )

        _PREPARE_FELSENSTEIN_CONTEXT = _prepare_context

    return _PREPARE_FELSENSTEIN_CONTEXT(*args, **kwargs)


def count_params(*args, **kwargs):
    global _COUNT_PARAMS

    if _COUNT_PARAMS is None:
        from ...helpers.discrete_models import count_params as _count_params

        _COUNT_PARAMS = _count_params

    return _COUNT_PARAMS(*args, **kwargs)


def fit_q_matrix(*args, **kwargs):
    global _FIT_Q_MATRIX

    if _FIT_Q_MATRIX is None:
        from ...helpers.discrete_models import fit_q_matrix as _fit_q_matrix

        _FIT_Q_MATRIX = _fit_q_matrix

    return _FIT_Q_MATRIX(*args, **kwargs)


def parse_discrete_traits(*args, **kwargs):
    global _PARSE_DISCRETE_TRAITS

    if _PARSE_DISCRETE_TRAITS is None:
        from ...helpers.discrete_models import parse_discrete_traits as _parse_traits

        _PARSE_DISCRETE_TRAITS = _parse_traits

    return _PARSE_DISCRETE_TRAITS(*args, **kwargs)


class FitDiscrete(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"])
        self.trait_data_path = parsed["trait_data_path"]
        self.trait_column = parsed["trait_column"]
        self.selected_models = parsed["models"]
        self.json_output = parsed["json_output"]

    def run(self) -> None:
        tree = self.read_tree_file_unmodified()
        copied_tree = False
        if self._needs_default_branch_lengths(tree):
            tree = self._fast_copy(tree)
            copied_tree = True
        self.validate_tree(tree, min_tips=3, assign_default_branch_length=1e-8, context="model fitting")

        tree_tips = self.get_tip_names_from_tree(tree)
        tip_states = parse_discrete_traits(
            self.trait_data_path, tree_tips, trait_column=self.trait_column
        )

        tips_to_prune = self._tips_to_prune_for_states(tree_tips, tip_states)
        if tips_to_prune:
            if not copied_tree:
                tree = self._fast_copy(tree)
            tree = self.prune_tree_using_taxa_list(tree, tips_to_prune)

        states = sorted(set(tip_states.values()))
        n_obs = len(tip_states)
        pruning_context = _prepare_felsenstein_context(tree, tip_states, states)

        results = []
        for model_name in self.selected_models:
            result = self._fit_model(
                tree, tip_states, states, model_name, n_obs,
                pruning_context=pruning_context,
            )
            results.append(result)

        results = self._compute_model_comparison(results)

        if self.json_output:
            self._print_json(results, n_obs, states)
        else:
            self._print_text(results, n_obs, states)

    @staticmethod
    def _tips_to_prune_for_states(tree_tips, tip_states):
        if len(tree_tips) == len(tip_states):
            for tip_name, state_name in zip(tree_tips, tip_states):
                if tip_name != state_name:
                    break
            else:
                return []

        shared_taxa = set(tip_states)
        return [tip for tip in tree_tips if tip not in shared_taxa]

    @staticmethod
    def _needs_default_branch_lengths(tree) -> bool:
        try:
            root = tree.root
            root.clades
        except AttributeError:
            return True

        stack = [root]
        pop = stack.pop
        extend = stack.extend
        try:
            while stack:
                clade = pop()
                children = clade.clades
                if clade is not root and clade.branch_length is None:
                    return True
                if children:
                    extend(children)
        except AttributeError:
            return True

        return False

    def process_args(self, args) -> dict:
        models_str = getattr(args, "models", None)
        if models_str is None:
            models = ["ER", "SYM", "ARD"]
        else:
            models = [m.strip().upper() for m in models_str.split(",")]
            for m in models:
                if m not in VALID_DISCRETE_MODELS:
                    raise PhykitUserError(
                        [
                            f"Unknown model '{m}'.",
                            f"Valid models: {', '.join(sorted(VALID_DISCRETE_MODELS))}",
                        ],
                        code=2,
                    )

        return dict(
            tree_file_path=args.tree,
            trait_data_path=args.trait_data,
            trait_column=args.trait,
            models=models,
            json_output=getattr(args, "json", False),
        )

    def _fit_model(
        self, tree, tip_states, states, model_name, n_obs, pruning_context=None
    ):
        k = len(states)
        Q, lnL = fit_q_matrix(
            tree, tip_states, states, model_name,
            pruning_context=pruning_context,
        )
        n_params = count_params(k, model_name)
        aic = -2.0 * lnL + 2.0 * n_params
        bic = -2.0 * lnL + n_params * math.log(n_obs)
        return {
            "model": model_name,
            "lnL": lnL,
            "aic": aic,
            "bic": bic,
            "n_params": n_params,
            "n_states": k,
            "q_matrix": Q.tolist(),
        }

    def _compute_model_comparison(self, results):
        if len(results) == 1:
            result = results[0]
            result["delta_aic"] = 0.0
            result["aic_weight"] = 1.0
            result["delta_bic"] = 0.0
            return results

        # Delta-AIC and AIC weights
        min_aic = min(r["aic"] for r in results)
        for r in results:
            r["delta_aic"] = r["aic"] - min_aic

        raw_weights = [math.exp(-0.5 * r["delta_aic"]) for r in results]
        total_weight = sum(raw_weights)
        for r, w in zip(results, raw_weights):
            r["aic_weight"] = w / total_weight if total_weight > 0 else 0.0

        # Delta-BIC
        min_bic = min(r["bic"] for r in results)
        for r in results:
            r["delta_bic"] = r["bic"] - min_bic

        # Sort by AIC (best first)
        results.sort(key=lambda r: r["aic"])
        return results

    def _print_text(self, results, n_obs, states):
        header = f"{'Model':<8}{'lnL':>12}{'AIC':>12}{'dAIC':>10}{'AIC_w':>10}{'BIC':>12}{'dBIC':>10}{'n_params':>10}"
        lines = [
            f"Number of taxa: {n_obs}",
            f"Number of states: {len(states)}",
            f"States: {', '.join(states)}",
            "",
            header,
            "-" * len(header),
        ]
        for r in results:
            lines.append(
                f"{r['model']:<8}"
                f"{r['lnL']:>12.4f}"
                f"{r['aic']:>12.4f}"
                f"{r['delta_aic']:>10.4f}"
                f"{r['aic_weight']:>10.4f}"
                f"{r['bic']:>12.4f}"
                f"{r['delta_bic']:>10.4f}"
                f"{r['n_params']:>10d}"
            )
        print("\n".join(lines))

    def _print_json(self, results, n_obs, states):
        # Remove q_matrix numpy arrays (already converted to lists in _fit_model)
        payload = {
            "n_taxa": n_obs,
            "n_states": len(states),
            "states": states,
            "models": [
                {
                    "model": r["model"],
                    "lnL": round(r["lnL"], 4),
                    "aic": round(r["aic"], 4),
                    "delta_aic": round(r["delta_aic"], 4),
                    "aic_weight": round(r["aic_weight"], 4),
                    "bic": round(r["bic"], 4),
                    "delta_bic": round(r["delta_bic"], 4),
                    "n_params": r["n_params"],
                    "q_matrix": [
                        [round(v, 6) for v in row] for row in r["q_matrix"]
                    ],
                }
                for r in results
            ],
        }
        print_json(payload)
