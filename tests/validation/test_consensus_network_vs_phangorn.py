"""
Cross-validate PhyKIT consensus_network against R/phangorn as.splits().

Requires: R with phangorn package installed.
Run:  python -m pytest tests/validation/test_consensus_network_vs_phangorn.py -v
"""
import json
import subprocess
import shutil
import tempfile
from argparse import Namespace
from pathlib import Path

import pytest

from phykit.services.tree.consensus_network import ConsensusNetwork


R_SCRIPT = r"""
suppressPackageStartupMessages({library(ape); library(phangorn)})

args <- commandArgs(trailingOnly = TRUE)
tree_file <- args[1]

trees <- read.tree(tree_file)
sp <- as.splits(trees)
labs <- attr(sp, "labels")
wts <- attr(sp, "weights")
confs <- attr(sp, "confidences")
n_taxa <- length(labs)
n_trees <- length(trees)

results <- list()
idx <- 1
for (i in seq_along(sp)) {
  taxa <- sort(labs[sp[[i]]])
  n_in <- length(taxa)
  if (n_in <= 1 || n_in >= (n_taxa - 1)) next
  comp <- sort(setdiff(labs, taxa))
  if (length(taxa) > length(comp)) { display <- comp }
  else if (length(taxa) < length(comp)) { display <- taxa }
  else { display <- if (paste(taxa,collapse=",") < paste(comp,collapse=",")) taxa else comp }
  cat(sprintf("%s\t%d\t%.6f\n", paste(display, collapse=","), wts[i], confs[i]))
}
"""


def _has_r_phangorn():
    if not shutil.which("Rscript"):
        return False
    try:
        result = subprocess.run(
            ["Rscript", "-e", "library(phangorn)"],
            capture_output=True, timeout=10
        )
        return result.returncode == 0
    except Exception:
        return False


def _run_r_splits(tree_file_path):
    """Run R script and return dict of {frozenset: (count, freq)}."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".R", delete=False) as f:
        f.write(R_SCRIPT)
        r_script_path = f.name

    try:
        result = subprocess.run(
            ["Rscript", r_script_path, tree_file_path],
            capture_output=True, text=True, timeout=30
        )
        assert result.returncode == 0, f"R failed: {result.stderr}"
    finally:
        Path(r_script_path).unlink()

    splits = {}
    for line in result.stdout.strip().splitlines():
        parts = line.split("\t")
        taxa = frozenset(parts[0].split(","))
        count = int(parts[1])
        freq = float(parts[2])
        splits[taxa] = (count, freq)
    return splits


def _run_phykit_splits(tree_file_path, threshold=0.0):
    """Run PhyKIT and return dict of {frozenset: (count, freq)}."""
    captured = {}

    import phykit.services.tree.consensus_network as mod
    original_print_json = mod.print_json

    def capture_json(payload):
        captured["payload"] = payload

    mod.print_json = capture_json
    try:
        svc = ConsensusNetwork(Namespace(
            trees=tree_file_path,
            threshold=threshold,
            missing_taxa="error",
            plot_output=None,
            json=True,
        ))
        svc.run()
    finally:
        mod.print_json = original_print_json

    splits = {}
    for s in captured["payload"]["splits"]:
        taxa = frozenset(s["split"])
        splits[taxa] = (s["count"], s["frequency"])
    return splits


@pytest.mark.skipif(not _has_r_phangorn(), reason="R with phangorn not available")
class TestConsensusNetworkVsPhangorn:
    def test_sample_gene_trees(self):
        tree_file = "tests/sample_files/gene_trees_for_network.nwk"
        r_splits = _run_r_splits(tree_file)
        pk_splits = _run_phykit_splits(tree_file)

        assert set(r_splits.keys()) == set(pk_splits.keys()), (
            f"Split sets differ.\nR only: {set(r_splits) - set(pk_splits)}\n"
            f"PK only: {set(pk_splits) - set(r_splits)}"
        )
        for split in r_splits:
            r_count, r_freq = r_splits[split]
            pk_count, pk_freq = pk_splits[split]
            assert r_count == pk_count, (
                f"Count mismatch for {sorted(split)}: R={r_count}, PK={pk_count}"
            )
            assert abs(r_freq - pk_freq) < 0.001, (
                f"Freq mismatch for {sorted(split)}: R={r_freq}, PK={pk_freq}"
            )

    def test_simple_four_taxon(self, tmp_path):
        tree_file = tmp_path / "trees.nwk"
        tree_file.write_text(
            "((A,B),(C,D));\n((A,B),(C,D));\n((A,C),(B,D));\n"
        )
        r_splits = _run_r_splits(str(tree_file))
        pk_splits = _run_phykit_splits(str(tree_file))

        assert set(r_splits.keys()) == set(pk_splits.keys())
        for split in r_splits:
            assert r_splits[split][0] == pk_splits[split][0]

    def test_identical_trees(self, tmp_path):
        tree_file = tmp_path / "trees.nwk"
        tree_file.write_text("((A,B),(C,D));\n" * 5)
        r_splits = _run_r_splits(str(tree_file))
        pk_splits = _run_phykit_splits(str(tree_file))

        assert set(r_splits.keys()) == set(pk_splits.keys())
        for split in r_splits:
            assert r_splits[split][0] == pk_splits[split][0]
