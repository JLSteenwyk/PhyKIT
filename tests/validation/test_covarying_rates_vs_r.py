from argparse import Namespace
from pathlib import Path
import shutil
import subprocess

import pytest

from phykit.services.tree.covarying_evolutionary_rates import (
    CovaryingEvolutionaryRates,
    _pearsonr,
)


ROOT = Path(__file__).resolve().parents[2]
SAMPLES = ROOT / "tests" / "sample_files"
R_VALIDATOR = Path(__file__).with_name("validate_covarying_rates.R")


@pytest.mark.validation
def test_covarying_rates_match_r_ape_and_stats():
    if shutil.which("Rscript") is None:
        pytest.skip("Rscript is not installed")
    dependency = subprocess.run(
        ["Rscript", "-e", "quit(status=!requireNamespace('ape', quietly=TRUE))"],
        check=False,
    )
    if dependency.returncode:
        pytest.skip("R package 'ape' is not installed")

    tree_zero = SAMPLES / "tree_simple.tre"
    tree_one = SAMPLES / "tree_simple_1.tre"
    reference = SAMPLES / "tree_simple_2.tre"
    completed = subprocess.run(
        [
            "Rscript",
            str(R_VALIDATOR),
            str(tree_zero),
            str(tree_one),
            str(reference),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    r_correlation, r_p_value, r_branch_count = completed.stdout.strip().split("\t")

    service = CovaryingEvolutionaryRates(
        Namespace(
            tree_zero=str(tree_zero),
            tree_one=str(tree_one),
            reference=str(reference),
            verbose=False,
        )
    )
    trees = (
        service.read_tree_file_unmodified(),
        service.read_tree1_file_unmodified(),
        service.read_reference_tree_file_unmodified(),
    )
    rates_zero, rates_one, _ = service.correct_branch_lengths(*trees)
    outliers = service.get_indices_of_outlier_branch_lengths(rates_zero, [])
    outliers = service.get_indices_of_outlier_branch_lengths(rates_one, outliers)
    rates_zero = service.remove_outliers_based_on_indices(rates_zero, outliers)
    rates_one = service.remove_outliers_based_on_indices(rates_one, outliers)
    correlation, p_value = _pearsonr(rates_zero, rates_one)

    assert len(rates_zero) == int(r_branch_count) == 13
    assert correlation == pytest.approx(float(r_correlation), abs=1e-12)
    assert p_value == pytest.approx(float(r_p_value), abs=1e-12)
