"""
Cross-validate PhyKIT threshold_model against R/phytools::threshBayes.

Requires: R with phytools and ape packages installed.
Run:  python -m pytest tests/validation/test_threshold_model_vs_phytools.py -v
"""
import subprocess
import shutil
import tempfile
from argparse import Namespace
from pathlib import Path

import numpy as np
import pytest

from phykit.services.tree.threshold_model import ThresholdModel


TREE_SIMPLE = "tests/sample_files/tree_simple.tre"
TRAITS_FILE = "tests/sample_files/tree_simple_threshold_traits.tsv"


R_SCRIPT = r"""
suppressPackageStartupMessages({library(ape); library(phytools)})

args <- commandArgs(trailingOnly = TRUE)
tree_file <- args[1]
trait_file <- args[2]
seed_val <- as.integer(args[3])

tree <- read.tree(tree_file)
data <- read.delim(trait_file, row.names=1)
X <- data[tree$tip.label, , drop=FALSE]

set.seed(seed_val)
result <- threshBayes(tree, X, types=c("discrete","continuous"),
                      ngen=100000, control=list(sample=100))

# Discard first 20% as burnin
n_total <- nrow(result$par)
burnin <- ceiling(0.2 * n_total)
post <- result$par[(burnin+1):n_total, ]

cat(sprintf("r_mean=%.6f\n", mean(post[,"r"])))
cat(sprintf("r_q025=%.6f\n", quantile(post[,"r"], 0.025)))
cat(sprintf("r_q975=%.6f\n", quantile(post[,"r"], 0.975)))
"""


def _has_r_phytools():
    if not shutil.which("Rscript"):
        return False
    try:
        result = subprocess.run(
            ["Rscript", "-e", "library(phytools)"],
            capture_output=True, timeout=10,
        )
        return result.returncode == 0
    except Exception:
        return False


def _run_r_threshbayes(tree_file, trait_file, seed=42):
    """Run R threshBayes and return posterior summary."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".R", delete=False) as f:
        f.write(R_SCRIPT)
        r_script_path = f.name

    try:
        result = subprocess.run(
            ["Rscript", r_script_path, tree_file, trait_file, str(seed)],
            capture_output=True, text=True, timeout=600,
        )
        assert result.returncode == 0, f"R failed: {result.stderr}"
    finally:
        Path(r_script_path).unlink()

    summary = {}
    for line in result.stdout.strip().splitlines():
        key, val = line.split("=")
        summary[key] = float(val)
    return summary


def _run_phykit_threshold(tree_file, trait_file, seed=42):
    """Run PhyKIT threshold_model and return posterior summary."""
    args = Namespace(
        tree=tree_file,
        trait_data=trait_file,
        traits="habitat,body_mass",
        types="discrete,continuous",
        ngen=100000,
        sample=100,
        burnin=0.2,
        seed=seed,
        plot_output=None,
        json=False,
    )
    svc = ThresholdModel(args)
    svc.run = lambda: None  # Don't run output

    from Bio import Phylo
    tree = Phylo.read(tree_file, "newick")
    tree_tips = [t.name for t in tree.get_terminals()]
    t1, t2, ordered = ThresholdModel._parse_multi_trait_file(
        trait_file, "habitat", "body_mass",
        "discrete", "continuous", tree_tips,
    )
    shared_set = set(ordered)
    tips_to_prune = [t for t in tree_tips if t not in shared_set]
    if tips_to_prune:
        for tip_name in tips_to_prune:
            tips = [t for t in tree.get_terminals() if t.name == tip_name]
            if tips:
                tree.prune(tips[0])

    C = ThresholdModel._build_vcv_matrix(tree, ordered)
    rng = np.random.default_rng(seed)

    result = ThresholdModel._run_mcmc(
        t1, t2, "discrete", "continuous",
        ordered, C, 100000, 100, 0.2, rng,
    )

    r_samples = result["r"]
    return {
        "r_mean": float(np.mean(r_samples)),
        "r_q025": float(np.quantile(r_samples, 0.025)),
        "r_q975": float(np.quantile(r_samples, 0.975)),
    }


@pytest.mark.skipif(not _has_r_phytools(), reason="R with phytools not available")
class TestThresholdModelVsPhytools:
    def test_correlation_direction_agrees(self):
        """Both implementations should agree on the sign of r."""
        r_summary = _run_r_threshbayes(TREE_SIMPLE, TRAITS_FILE)
        pk_summary = _run_phykit_threshold(TREE_SIMPLE, TRAITS_FILE)

        # Both should have the same sign for r_mean
        r_sign = np.sign(r_summary["r_mean"])
        pk_sign = np.sign(pk_summary["r_mean"])
        assert r_sign == pk_sign, (
            f"Sign mismatch: R r_mean={r_summary['r_mean']:.4f}, "
            f"PK r_mean={pk_summary['r_mean']:.4f}"
        )

    def test_intervals_overlap(self):
        """95% intervals from both implementations should overlap."""
        r_summary = _run_r_threshbayes(TREE_SIMPLE, TRAITS_FILE)
        pk_summary = _run_phykit_threshold(TREE_SIMPLE, TRAITS_FILE)

        # Check that intervals overlap
        overlap = (
            max(r_summary["r_q025"], pk_summary["r_q025"])
            < min(r_summary["r_q975"], pk_summary["r_q975"])
        )
        assert overlap, (
            f"No overlap: R=[{r_summary['r_q025']:.4f}, {r_summary['r_q975']:.4f}], "
            f"PK=[{pk_summary['r_q025']:.4f}, {pk_summary['r_q975']:.4f}]"
        )

    def test_means_within_tolerance(self):
        """Posterior means should be within 0.3 of each other.

        MCMC samplers won't produce identical results, but with enough
        iterations the posterior means should converge to similar values.
        The tolerance is generous because the dataset is small (8 taxa).
        """
        r_summary = _run_r_threshbayes(TREE_SIMPLE, TRAITS_FILE)
        pk_summary = _run_phykit_threshold(TREE_SIMPLE, TRAITS_FILE)

        diff = abs(r_summary["r_mean"] - pk_summary["r_mean"])
        assert diff < 0.3, (
            f"Means too different: R={r_summary['r_mean']:.4f}, "
            f"PK={pk_summary['r_mean']:.4f}, diff={diff:.4f}"
        )
