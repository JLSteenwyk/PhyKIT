import json
import subprocess
import sys
from argparse import Namespace
from pathlib import Path
from types import SimpleNamespace

import pytest

import phykit.services.alignment.codon_dnds as codon_dnds_module
from phykit.errors import PhykitUserError
from phykit.services.alignment.codon_dnds import CodonDnDs


SAMPLE_DIR = Path(__file__).parents[3] / "sample_files"
REFERENCE_ALIGNMENT = SAMPLE_DIR / "codon_dnds_reference.fa"


def _make_args(**overrides):
    defaults = dict(
        alignment=str(REFERENCE_ALIGNMENT),
        method="NG86",
        genetic_code=1,
        kappa=1.0,
        codon_frequency="F3x4",
        stop_policy="error",
        reference=None,
        verbose=False,
        json=False,
    )
    defaults.update(overrides)
    return Namespace(**defaults)


def _write_alignment(tmp_path, records):
    path = tmp_path / "codons.fa"
    path.write_text("".join(f">{name}\n{sequence}\n" for name, sequence in records))
    return path


def _single_result(service):
    alignment, _, is_protein = service.get_alignment_and_format()
    table = service._get_codon_table()
    records = service._validate_alignment(alignment, is_protein, table)
    pairs = list(service._select_pairs(records))
    assert len(pairs) == 1
    return service._estimate_pair(*pairs[0], table)


def test_module_import_defers_biopython_codon_api():
    code = """
import sys
import phykit.services.alignment.codon_dnds as module
assert module._CALCULATE_DN_DS is None
assert "Bio.Align.analysis" not in sys.modules
assert "Bio.Data.CodonTable" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


@pytest.mark.parametrize(
    ("method", "expected_d_n", "expected_d_s"),
    [
        ("NG86", 0.0677151350, 0.2011979899),
        ("LWL85", 0.0687284919, 0.2075511181),
        ("YN00", 0.0814688234, 0.1277062093),
    ],
)
def test_reference_estimates_match_biopython_documentation(
    method, expected_d_n, expected_d_s
):
    service = CodonDnDs(_make_args(method=method))

    result = _single_result(service)

    assert result["dN"] == pytest.approx(expected_d_n, rel=1e-7)
    assert result["dS"] == pytest.approx(expected_d_s, rel=1e-7)
    assert result["omega"] == pytest.approx(expected_d_n / expected_d_s, rel=1e-7)
    assert result["codons_used"] == 23
    assert result["codons_skipped"] == 0
    assert result["status"] == "ok"


def test_run_prints_pairwise_tsv(capsys):
    CodonDnDs(_make_args(verbose=True)).run()

    lines = capsys.readouterr().out.splitlines()
    assert lines[0] == (
        "taxon_a\ttaxon_b\tdN\tdS\tomega\tcodons_used\t"
        "codons_skipped\tstatus"
    )
    fields = lines[1].split("\t")
    assert fields[:2] == ["rna1", "rna2"]
    assert float(fields[2]) == pytest.approx(0.0677151350)
    assert float(fields[3]) == pytest.approx(0.2011979899)
    assert fields[-3:] == ["23", "0", "ok"]


def test_run_prints_structured_json(capsys):
    CodonDnDs(_make_args(verbose=True, json=True)).run()

    payload = json.loads(capsys.readouterr().out)
    assert payload["method"] == "NG86"
    assert payload["genetic_code"] == 1
    assert payload["summary"]["pairs_total"] == 1
    assert payload["summary"]["pairs_with_omega"] == 1
    assert payload["pairs"][0]["status"] == "ok"


def test_reference_limits_pairs_and_preserves_reference_first(tmp_path):
    path = _write_alignment(
        tmp_path,
        [
            ("a", "ATGGCTGAA" * 4),
            ("b", "ATGGCCGAA" * 4),
            ("c", "ATGGCTGAC" * 4),
        ],
    )
    service = CodonDnDs(_make_args(alignment=str(path), reference="b"))
    alignment, _, is_protein = service.get_alignment_and_format()
    table = service._get_codon_table()
    records = service._validate_alignment(alignment, is_protein, table)

    pairs = list(service._select_pairs(records))

    assert [(left[0], right[0]) for left, right in pairs] == [("b", "a"), ("b", "c")]


def test_pairwise_deletion_excludes_ambiguous_and_terminal_stop_codons(tmp_path):
    path = _write_alignment(
        tmp_path,
        [
            ("a", "ATGGCTNNNTAA"),
            ("b", "ATGGCNGACTAG"),
        ],
    )
    service = CodonDnDs(_make_args(alignment=str(path)))

    result = _single_result(service)

    assert result["codons_used"] == 1
    assert result["codons_skipped"] == 3


def test_identical_sequences_report_zero_synonymous_distance(tmp_path):
    sequence = "ATGGCTGAA" * 10
    path = _write_alignment(tmp_path, [("a", sequence), ("b", sequence)])
    service = CodonDnDs(_make_args(alignment=str(path)))

    result = _single_result(service)

    assert result["dN"] == pytest.approx(0.0)
    assert result["dS"] == pytest.approx(0.0)
    assert result["omega"] is None
    assert result["status"] == "zero_synonymous_distance"


def test_alignment_length_must_preserve_frame(tmp_path):
    path = _write_alignment(tmp_path, [("a", "ATGG"), ("b", "ATGG")])
    service = CodonDnDs(_make_args(alignment=str(path)))

    with pytest.raises(PhykitUserError) as error:
        service.run()
    assert "not divisible by 3" in "\n".join(error.value.messages)


def test_alignment_records_must_have_equal_length():
    class UnequalAlignment(list):
        @staticmethod
        def get_alignment_length():
            return 6

    alignment = UnequalAlignment(
        [
            SimpleNamespace(id="a", seq="ATGGCT"),
            SimpleNamespace(id="b", seq="ATG"),
        ]
    )
    service = CodonDnDs(_make_args())

    with pytest.raises(PhykitUserError) as error:
        service._validate_alignment(alignment, False, service._get_codon_table())
    assert "must have equal length" in "\n".join(error.value.messages)


def test_partial_codon_gap_is_rejected(tmp_path):
    path = _write_alignment(tmp_path, [("a", "ATGA--"), ("b", "ATGGCT")])
    service = CodonDnDs(_make_args(alignment=str(path)))

    with pytest.raises(PhykitUserError) as error:
        service.run()
    assert "frame-disrupting gap" in "\n".join(error.value.messages)


def test_internal_stop_is_rejected_by_default(tmp_path):
    path = _write_alignment(tmp_path, [("a", "ATGTAGGCT"), ("b", "ATGGCTGCT")])
    service = CodonDnDs(_make_args(alignment=str(path)))

    with pytest.raises(PhykitUserError) as error:
        service.run()
    assert "internal stop codon" in "\n".join(error.value.messages)


def test_internal_stop_can_be_skipped_explicitly(tmp_path):
    path = _write_alignment(
        tmp_path,
        [("a", "ATGTAGGCTGAA"), ("b", "ATGGCTGCTGAC")],
    )
    service = CodonDnDs(
        _make_args(alignment=str(path), stop_policy="skip")
    )

    result = _single_result(service)

    assert result["codons_used"] == 3
    assert result["codons_skipped"] == 1


def test_duplicate_taxon_ids_are_rejected(tmp_path):
    path = _write_alignment(tmp_path, [("a", "ATGGCT"), ("a", "ATGGCC")])
    service = CodonDnDs(_make_args(alignment=str(path)))

    with pytest.raises(PhykitUserError) as error:
        service.run()
    assert "duplicate taxon ID" in "\n".join(error.value.messages)


def test_unknown_genetic_code_is_rejected():
    service = CodonDnDs(_make_args(genetic_code=999))

    with pytest.raises(PhykitUserError) as error:
        service.run()
    assert "Unknown NCBI genetic code" in "\n".join(error.value.messages)


def test_missing_reference_is_rejected():
    service = CodonDnDs(_make_args(reference="absent"))

    with pytest.raises(PhykitUserError) as error:
        service.run()
    assert "is not present" in "\n".join(error.value.messages)


@pytest.mark.parametrize("kappa", [0.0, -1.0, float("inf")])
def test_kappa_must_be_positive_and_finite(kappa):
    with pytest.raises(PhykitUserError) as error:
        CodonDnDs(_make_args(kappa=kappa))
    assert "--kappa" in "\n".join(error.value.messages)


def test_ml_codon_frequency_is_forwarded(monkeypatch):
    captured = {}

    def fake_calculate(alignment, **kwargs):
        captured.update(kwargs)
        return 0.2, 0.4

    monkeypatch.setattr(codon_dnds_module, "_CALCULATE_DN_DS", fake_calculate)
    monkeypatch.setattr(
        codon_dnds_module, "_BIO_ALIGNMENT", lambda sequences: sequences
    )
    monkeypatch.setattr(codon_dnds_module, "_BIO_SEQ", lambda sequence: sequence)
    service = CodonDnDs(_make_args(method="ML", codon_frequency="F61"))

    result = _single_result(service)

    assert captured["method"] == "ML"
    assert captured["cfreq"] == "F61"
    assert "k" not in captured
    assert result["omega"] == pytest.approx(0.5)
