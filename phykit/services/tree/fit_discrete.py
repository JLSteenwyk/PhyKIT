"""
Compare models of discrete trait evolution on a phylogeny.

Fits ER (Equal Rates), SYM (Symmetric), and ARD (All Rates Different)
Mk models and compares them using AIC/BIC. Analogous to R's
geiger::fitDiscrete().
"""
import math
from typing import Dict, List

import numpy as np

from .base import Tree
from ...helpers.json_output import print_json
from ...helpers.discrete_models import (
    VALID_DISCRETE_MODELS,
    count_params,
    fit_q_matrix,
    parse_discrete_traits,
)
from ...errors import PhykitUserError


class FitDiscrete(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"])
        self.trait_data_path = parsed["trait_data_path"]
        self.trait_column = parsed["trait_column"]
        self.selected_models = parsed["models"]
        self.json_output = parsed["json_output"]

    def run(self) -> None:
        tree = self.read_tree_file()
        self._validate_tree(tree)

        tree_tips = self.get_tip_names_from_tree(tree)
        tip_states = parse_discrete_traits(
            self.trait_data_path, tree_tips, trait_column=self.trait_column
        )

        # Prune tree to shared taxa
        shared_taxa = set(tip_states.keys())
        tips_to_prune = [t for t in tree_tips if t not in shared_taxa]
        if tips_to_prune:
            tree = self.prune_tree_using_taxa_list(tree, tips_to_prune)

        states = sorted(set(tip_states.values()))
        n_obs = len(tip_states)

        results = []
        for model_name in self.selected_models:
            result = self._fit_model(tree, tip_states, states, model_name, n_obs)
            results.append(result)

        results = self._compute_model_comparison(results)

        if self.json_output:
            self._print_json(results, n_obs, states)
        else:
            self._print_text(results, n_obs, states)

    def process_args(self, args) -> Dict:
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

    def _validate_tree(self, tree) -> None:
        tips = list(tree.get_terminals())
        if len(tips) < 3:
            raise PhykitUserError(
                ["Tree must have at least 3 tips."], code=2
            )
        for clade in tree.find_clades():
            if clade.branch_length is None and clade != tree.root:
                clade.branch_length = 1e-8

    def _fit_model(self, tree, tip_states, states, model_name, n_obs):
        k = len(states)
        Q, lnL = fit_q_matrix(tree, tip_states, states, model_name)
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
        print(f"Number of taxa: {n_obs}")
        print(f"Number of states: {len(states)}")
        print(f"States: {', '.join(states)}")
        print()

        header = f"{'Model':<8}{'lnL':>12}{'AIC':>12}{'dAIC':>10}{'AIC_w':>10}{'BIC':>12}{'dBIC':>10}{'n_params':>10}"
        print(header)
        print("-" * len(header))
        for r in results:
            print(
                f"{r['model']:<8}"
                f"{r['lnL']:>12.4f}"
                f"{r['aic']:>12.4f}"
                f"{r['delta_aic']:>10.4f}"
                f"{r['aic_weight']:>10.4f}"
                f"{r['bic']:>12.4f}"
                f"{r['delta_bic']:>10.4f}"
                f"{r['n_params']:>10d}"
            )

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
