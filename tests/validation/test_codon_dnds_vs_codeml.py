"""Cross-validate PhyKIT codon dN/dS ML estimates against PAML codeml."""

from argparse import Namespace
from pathlib import Path
import re
import shutil
import subprocess

import pytest

from phykit.services.alignment.codon_dnds import CodonDnDs


pytestmark = pytest.mark.validation

REFERENCE_ALIGNMENT = (
    Path(__file__).parents[1] / "sample_files" / "codon_dnds_reference.fa"
)


def _read_fasta_records(path):
    records = []
    name = None
    sequence = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                records.append((name, "".join(sequence)))
            name = line[1:]
            sequence = []
        else:
            sequence.append(line.strip())
    if name is not None:
        records.append((name, "".join(sequence)))
    return records


def _run_codeml_pairwise(tmp_path):
    executable = shutil.which("codeml")
    if executable is None:
        pytest.skip("PAML codeml is not available")

    records = _read_fasta_records(REFERENCE_ALIGNMENT)
    alignment_length = len(records[0][1])
    (tmp_path / "pair.phy").write_text(
        f"{len(records)} {alignment_length}\n"
        + "".join(f"{name:<10} {sequence}\n" for name, sequence in records)
    )
    (tmp_path / "pair.tree").write_text("(rna1:0.1,rna2:0.1);\n")
    (tmp_path / "codeml.ctl").write_text(
        """seqfile = pair.phy
treefile = pair.tree
outfile = results.txt
noisy = 0
verbose = 0
runmode = -2
seqtype = 1
CodonFreq = 2
model = 0
NSsites = 0
icode = 0
fix_kappa = 0
kappa = 2
fix_omega = 0
omega = 0.4
cleandata = 1
"""
    )
    completed = subprocess.run(
        [executable, "codeml.ctl"],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr
    output = (tmp_path / "results.txt").read_text()
    match = re.search(
        r"dN/dS=\s*([0-9.eE+-]+)\s+dN\s*=\s*([0-9.eE+-]+)"
        r"\s+dS\s*=\s*([0-9.eE+-]+)",
        output,
    )
    assert match is not None, output
    omega, d_n, d_s = (float(value) for value in match.groups())
    return {"dN": d_n, "dS": d_s, "omega": omega}


def _run_phykit_ml():
    service = CodonDnDs(
        Namespace(
            alignment=str(REFERENCE_ALIGNMENT),
            method="ML",
            genetic_code=1,
            kappa=1.0,
            codon_frequency="F3x4",
            stop_policy="error",
            reference=None,
            verbose=False,
            json=False,
        )
    )
    alignment, _, is_protein = service.get_alignment_and_format()
    table = service._get_codon_table()
    records = service._validate_alignment(alignment, is_protein, table)
    pair = next(iter(service._select_pairs(records)))
    return service._estimate_pair(*pair, table)


def test_ml_f3x4_agrees_with_paml_codeml(tmp_path):
    codeml = _run_codeml_pairwise(tmp_path)
    phykit = _run_phykit_ml()

    assert phykit["status"] == "ok"
    assert phykit["dN"] == pytest.approx(codeml["dN"], abs=0.003)
    assert phykit["dS"] == pytest.approx(codeml["dS"], abs=0.003)
    assert phykit["omega"] == pytest.approx(codeml["omega"], abs=0.005)
