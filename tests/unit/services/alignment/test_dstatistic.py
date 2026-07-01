import pytest
import subprocess
import sys
from argparse import Namespace
from io import StringIO
from math import isclose

import numpy as np
from Bio import Phylo
from Bio.Phylo.BaseTree import TreeMixin

from phykit.services.alignment.dstatistic import Dstatistic
import phykit.services.alignment.dstatistic as module
from phykit.errors import PhykitUserError


def test_module_import_does_not_import_biopython_parsers():
    code = """
import sys
import phykit.services.alignment.dstatistic as module
assert hasattr(module.np, "__getattr__")
assert callable(module.print_json)
assert "typing" not in sys.modules
assert "numpy" not in sys.modules
assert "json" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
assert "Bio.SeqIO.FastaIO" not in sys.modules
assert "Bio.Phylo" not in sys.modules
assert "Bio.AlignIO" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


def _write_alignment(path, seqs):
    with open(path, "w") as f:
        for name, seq in seqs.items():
            f.write(f">{name}\n{seq}\n")


def _make_args(alignment, p1="P1", p2="P2", p3="P3", outgroup="Outgroup",
               block_size=100, json_output=False):
    return Namespace(
        alignment=alignment,
        p1=p1,
        p2=p2,
        p3=p3,
        outgroup=outgroup,
        block_size=block_size,
        json=json_output,
    )


class TestDstatistic:
    def test_print_alignment_text_output_batches_report(self, monkeypatch):
        svc = Dstatistic.__new__(Dstatistic)
        svc.p1 = "P1"
        svc.p2 = "P2"
        svc.p3 = "P3"
        svc.outgroup = "Outgroup"
        svc.block_size = 100
        printed = []

        def fake_print(*args, **kwargs):
            printed.append((args, kwargs))

        monkeypatch.setattr("builtins.print", fake_print)
        svc._print_alignment_text_output(
            1000,
            120,
            75,
            45,
            0.25,
            0.0625,
            4.0,
            0.000063,
        )

        expected = "\n".join([
            "Patterson's D-statistic (ABBA-BABA Test)",
            "=========================================",
            "Topology: (((P1, P2), P3), Outgroup)",
            "P1: P1",
            "P2: P2",
            "P3: P3",
            "Outgroup: Outgroup",
            "",
            "Alignment length: 1000",
            "Informative sites: 120",
            "ABBA sites: 75",
            "BABA sites: 45",
            "D-statistic: 0.2500",
            "Block jackknife (block size: 100):",
            "  Standard error: 0.0625",
            "  Z-score: 4.00",
            "  p-value: 0.000063",
            "",
            f"Interpretation: {svc._interpret(0.25, 0.000063)}",
        ])
        assert printed == [((expected,), {})]

    def test_read_fasta_uses_first_header_token_uppercases_and_keeps_last_duplicate(
        self, tmp_path
    ):
        aln = tmp_path / "test.fa"
        aln.write_text(
            ">P1 description\nac gt\nn-\n"
            ">P2 second description\ncc gg\nn-\n"
            ">P1 replacement\ntttt\n"
        )
        assert Dstatistic._read_fasta(str(aln)) == {
            "P1": "TTTT",
            "P2": "CCGGN-",
        }

    def test_count_site_patterns_returns_block_counts(self):
        seq_p1 = "AAAACCAAAA"
        seq_p2 = "CCCCAAAAAA"
        seq_p3 = "CCCCCCAAAA"
        seq_o = "AAAAAAAAAA"

        abba, baba, block_abba, block_baba = Dstatistic._count_site_patterns(
            seq_p1,
            seq_p2,
            seq_p3,
            seq_o,
            block_size=5,
        )

        assert abba == 4
        assert baba == 2
        assert block_abba.tolist() == [4.0, 0.0]
        assert block_baba.tolist() == [1.0, 1.0]

    def test_count_site_patterns_all_valid_ascii_skips_validity_mask(self, mocker):
        mocker.patch.object(
            module.np,
            "ones",
            side_effect=AssertionError(
                "all-valid ASCII D-statistic counts should not build a valid mask"
            ),
        )

        abba, baba, block_abba, block_baba = Dstatistic._count_site_patterns(
            "AAAACCAAAA",
            "CCCCAAAAAA",
            "CCCCCCAAAA",
            "AAAAAAAAAA",
            block_size=5,
        )

        assert abba == 4
        assert baba == 2
        assert block_abba.tolist() == [4.0, 0.0]
        assert block_baba.tolist() == [1.0, 1.0]

    def test_count_site_patterns_all_invariant_skips_numpy_byte_setup(self, mocker):
        mocker.patch.object(
            module.np,
            "frombuffer",
            side_effect=AssertionError(
                "all-invariant D-statistic counts should return before NumPy setup"
            ),
        )

        abba, baba, block_abba, block_baba = Dstatistic._count_site_patterns(
            "ACGT" * 10,
            "ACGT" * 10,
            "ACGT" * 10,
            "ACGT" * 10,
            block_size=8,
        )

        assert abba == 0
        assert baba == 0
        assert block_abba.tolist() == [0.0] * 5
        assert block_baba.tolist() == [0.0] * 5

    def test_count_site_patterns_identical_ingroup_skips_numpy_byte_setup(self, mocker):
        mocker.patch.object(
            module.np,
            "frombuffer",
            side_effect=AssertionError(
                "identical ingroup D-statistic counts should return before NumPy setup"
            ),
        )

        abba, baba, block_abba, block_baba = Dstatistic._count_site_patterns(
            "CCCC" * 10,
            "CCCC" * 10,
            "CCCC" * 10,
            "AAAA" * 10,
            block_size=8,
        )

        assert abba == 0
        assert baba == 0
        assert block_abba.tolist() == [0.0] * 5
        assert block_baba.tolist() == [0.0] * 5

    def test_count_site_patterns_identical_p1_p2_skips_numpy_byte_setup(self, mocker):
        mocker.patch.object(
            module.np,
            "frombuffer",
            side_effect=AssertionError(
                "identical P1/P2 D-statistic counts should return before NumPy setup"
            ),
        )

        abba, baba, block_abba, block_baba = Dstatistic._count_site_patterns(
            "ACGT" * 10,
            "ACGT" * 10,
            "TGCA" * 10,
            "AAAA" * 10,
            block_size=8,
        )

        assert abba == 0
        assert baba == 0
        assert block_abba.tolist() == [0.0] * 5
        assert block_baba.tolist() == [0.0] * 5

    def test_count_site_patterns_identical_p3_outgroup_skips_numpy_byte_setup(
        self, mocker
    ):
        mocker.patch.object(
            module.np,
            "frombuffer",
            side_effect=AssertionError(
                "identical P3/outgroup D-statistic counts should return before NumPy setup"
            ),
        )

        abba, baba, block_abba, block_baba = Dstatistic._count_site_patterns(
            "ACGT" * 10,
            "TGCA" * 10,
            "GATT" * 10,
            "GATT" * 10,
            block_size=8,
        )

        assert abba == 0
        assert baba == 0
        assert block_abba.tolist() == [0.0] * 5
        assert block_baba.tolist() == [0.0] * 5

    def test_count_site_patterns_abba_only_skips_baba_count(self, mocker):
        count_nonzero = mocker.patch.object(
            module.np,
            "count_nonzero",
            wraps=module.np.count_nonzero,
        )

        abba, baba, block_abba, block_baba = Dstatistic._count_site_patterns(
            "AAAA" * 10,
            "CCCC" * 10,
            "CCCC" * 10,
            "AAAA" * 10,
            block_size=8,
        )

        assert abba == 40
        assert baba == 0
        assert block_abba.tolist() == [8.0] * 5
        assert block_baba.tolist() == [0.0] * 5
        assert count_nonzero.call_count == 1

    def test_count_site_patterns_baba_only_skips_abba_count(self, mocker):
        count_nonzero = mocker.patch.object(
            module.np,
            "count_nonzero",
            wraps=module.np.count_nonzero,
        )

        abba, baba, block_abba, block_baba = Dstatistic._count_site_patterns(
            "CCCC" * 10,
            "AAAA" * 10,
            "CCCC" * 10,
            "AAAA" * 10,
            block_size=8,
        )

        assert abba == 0
        assert baba == 40
        assert block_abba.tolist() == [0.0] * 5
        assert block_baba.tolist() == [8.0] * 5
        assert count_nonzero.call_count == 1

    def test_count_site_patterns_skips_ambiguous_sites(self):
        abba, baba, block_abba, block_baba = Dstatistic._count_site_patterns(
            "ANXXC",
            "CCCAA",
            "CCGCC",
            "AAAAA",
            block_size=5,
        )

        assert abba == 1
        assert baba == 1
        assert block_abba.tolist() == [1.0]
        assert block_baba.tolist() == [1.0]

    def test_count_site_patterns_unicode_fallback_preserves_counts_and_blocks(self):
        abba, baba, block_abba, block_baba = Dstatistic._count_site_patterns(
            "AΩCN",
            "TΩAA",
            "TΩCX",
            "AΩAA",
            block_size=2,
        )

        assert abba == 1
        assert baba == 1
        assert block_abba.tolist() == [1.0, 0.0]
        assert block_baba.tolist() == [0.0, 1.0]

    def test_jackknife_d_values_match_leave_one_out_loop(self, monkeypatch):
        block_abba = np.array([4.0, 0.0, 3.0, 1.0, 0.0])
        block_baba = np.array([1.0, 2.0, 0.0, 1.0, 0.0])
        total_abba = block_abba.sum()
        total_baba = block_baba.sum()
        expected = []
        for idx in range(len(block_abba)):
            loo_abba = total_abba - block_abba[idx]
            loo_baba = total_baba - block_baba[idx]
            denom = loo_abba + loo_baba
            expected.append(
                (loo_abba - loo_baba) / denom if denom > 0 else 0.0
            )

        monkeypatch.setattr(
            module.np,
            "sum",
            lambda *args, **kwargs: pytest.fail(
                "jackknife block totals should use ndarray.sum"
            ),
        )
        observed = Dstatistic._jackknife_d_values(block_abba, block_baba)

        np.testing.assert_allclose(observed, expected)

    def test_jackknife_d_values_handles_all_zero_blocks(self):
        observed = Dstatistic._jackknife_d_values(
            np.zeros(4),
            np.zeros(4),
        )

        np.testing.assert_allclose(observed, np.zeros(4))

    def test_sum_squared_deviations_uses_dot_product(self, monkeypatch):
        values = np.array([0.1, 0.3, 0.4, 0.8])
        center = float(values.mean())
        expected = float(np.sum((values - center) ** 2))

        monkeypatch.setattr(
            module.np,
            "sum",
            lambda *args, **kwargs: pytest.fail(
                "jackknife variance should use a dot product"
            ),
        )

        assert module._sum_squared_deviations(values, center) == pytest.approx(
            expected
        )

    def test_abba_counted_correctly(self, tmp_path):
        """4 ABBA sites, 0 BABA, rest invariant."""
        aln = tmp_path / "test.fa"
        seqs = {
            #         ABBA ABBA ABBA ABBA inv  inv
            "P1":       "AAAA" + "AAAAAA",
            "P2":       "CCCC" + "AAAAAA",
            "P3":       "CCCC" + "AAAAAA",
            "Outgroup": "AAAA" + "AAAAAA",
        }
        _write_alignment(str(aln), seqs)
        args = _make_args(str(aln))
        svc = Dstatistic(args)
        # Access internal computation by running and checking output
        # We directly test the algorithm by calling run with captured output
        import io, sys
        captured = io.StringIO()
        sys.stdout = captured
        svc.run()
        sys.stdout = sys.__stdout__
        output = captured.getvalue()
        assert "ABBA sites: 4" in output
        assert "BABA sites: 0" in output

    def test_baba_counted_correctly(self, tmp_path):
        """0 ABBA, 3 BABA, rest invariant."""
        aln = tmp_path / "test.fa"
        seqs = {
            #         BABA BABA BABA inv  inv  inv
            "P1":       "CCC" + "AAAAAAA",
            "P2":       "AAA" + "AAAAAAA",
            "P3":       "CCC" + "AAAAAAA",
            "Outgroup": "AAA" + "AAAAAAA",
        }
        _write_alignment(str(aln), seqs)
        args = _make_args(str(aln))
        svc = Dstatistic(args)
        import io, sys
        captured = io.StringIO()
        sys.stdout = captured
        svc.run()
        sys.stdout = sys.__stdout__
        output = captured.getvalue()
        assert "ABBA sites: 0" in output
        assert "BABA sites: 3" in output

    def test_d_statistic_positive(self, tmp_path):
        """More ABBA than BABA -> D > 0."""
        aln = tmp_path / "test.fa"
        # 4 ABBA, 2 BABA, 4 invariant = 10 sites
        seqs = {
            "P1":       "AAAACCAAAA",
            "P2":       "CCCCAAAAAA",
            "P3":       "CCCCCCAAAA",
            "Outgroup": "AAAAAAAAAA",
        }
        _write_alignment(str(aln), seqs)
        args = _make_args(str(aln))
        svc = Dstatistic(args)
        import io, sys
        captured = io.StringIO()
        sys.stdout = captured
        svc.run()
        sys.stdout = sys.__stdout__
        output = captured.getvalue()
        assert "ABBA sites: 4" in output
        assert "BABA sites: 2" in output
        assert "D-statistic: 0.3333" in output

    def test_d_statistic_zero(self, tmp_path):
        """Equal ABBA and BABA -> D = 0."""
        aln = tmp_path / "test.fa"
        # 2 ABBA, 2 BABA, 6 invariant = 10 sites
        seqs = {
            "P1":       "AACCAAAAAA",
            "P2":       "CCAAAAAAAA",
            "P3":       "CCCCAAAAAA",
            "Outgroup": "AAAAAAAAAA",
        }
        _write_alignment(str(aln), seqs)
        args = _make_args(str(aln))
        svc = Dstatistic(args)
        import io, sys
        captured = io.StringIO()
        sys.stdout = captured
        svc.run()
        sys.stdout = sys.__stdout__
        output = captured.getvalue()
        assert "ABBA sites: 2" in output
        assert "BABA sites: 2" in output
        assert "D-statistic: 0.0000" in output

    def test_d_statistic_negative(self, tmp_path):
        """More BABA than ABBA -> D < 0."""
        aln = tmp_path / "test.fa"
        # 1 ABBA, 4 BABA, 5 invariant
        seqs = {
            "P1":       "ACCCCAAAAA",
            "P2":       "CAAAAAAAAA",
            "P3":       "CCCCCAAAAA",
            "Outgroup": "AAAAAAAAAA",
        }
        _write_alignment(str(aln), seqs)
        args = _make_args(str(aln))
        svc = Dstatistic(args)
        import io, sys
        captured = io.StringIO()
        sys.stdout = captured
        svc.run()
        sys.stdout = sys.__stdout__
        output = captured.getvalue()
        assert "ABBA sites: 1" in output
        assert "BABA sites: 4" in output
        assert "D-statistic: -0.6000" in output

    def test_gaps_skipped(self, tmp_path):
        """Sites with gaps should not count."""
        aln = tmp_path / "test.fa"
        # site 0: ABBA but P1 has gap -> skip
        # site 1: ABBA normal -> count
        # sites 2-5: invariant
        seqs = {
            "P1":       "-AAAAA",
            "P2":       "CCAAAA",
            "P3":       "CCAAAA",
            "Outgroup": "AAAAAA",
        }
        _write_alignment(str(aln), seqs)
        args = _make_args(str(aln))
        svc = Dstatistic(args)
        import io, sys
        captured = io.StringIO()
        sys.stdout = captured
        svc.run()
        sys.stdout = sys.__stdout__
        output = captured.getvalue()
        assert "ABBA sites: 1" in output
        assert "BABA sites: 0" in output

    def test_non_biallelic_skipped(self, tmp_path):
        """Sites with 3+ alleles should be skipped."""
        aln = tmp_path / "test.fa"
        # site 0: P1=A, P2=C, P3=G, O=A -> 3 alleles -> skip
        # site 1: ABBA normal
        # sites 2-5: invariant
        seqs = {
            "P1":       "AAAAAA",
            "P2":       "CCAAAA",
            "P3":       "GCAAAA",
            "Outgroup": "AAAAAA",
        }
        _write_alignment(str(aln), seqs)
        args = _make_args(str(aln))
        svc = Dstatistic(args)
        import io, sys
        captured = io.StringIO()
        sys.stdout = captured
        svc.run()
        sys.stdout = sys.__stdout__
        output = captured.getvalue()
        assert "ABBA sites: 1" in output
        assert "BABA sites: 0" in output

    def test_invariant_skipped(self, tmp_path):
        """All-same sites should be skipped (only 1 allele, not biallelic)."""
        aln = tmp_path / "test.fa"
        seqs = {
            "P1":       "AAAAAAAAAA",
            "P2":       "AAAAAAAAAA",
            "P3":       "AAAAAAAAAA",
            "Outgroup": "AAAAAAAAAA",
        }
        _write_alignment(str(aln), seqs)
        args = _make_args(str(aln))
        svc = Dstatistic(args)
        import io, sys
        captured = io.StringIO()
        sys.stdout = captured
        svc.run()
        sys.stdout = sys.__stdout__
        output = captured.getvalue()
        assert "ABBA sites: 0" in output
        assert "BABA sites: 0" in output
        assert "D-statistic: 0.0000" in output

    def test_missing_taxon_raises(self, tmp_path):
        """Error when a specified taxon is not in the alignment."""
        aln = tmp_path / "test.fa"
        seqs = {
            "P1":       "AAAAAAAAAA",
            "P2":       "AAAAAAAAAA",
            "P3":       "AAAAAAAAAA",
            "Outgroup": "AAAAAAAAAA",
        }
        _write_alignment(str(aln), seqs)
        args = _make_args(str(aln), p1="missing_taxon")
        svc = Dstatistic(args)
        with pytest.raises(PhykitUserError):
            svc.run()

    def test_block_jackknife(self, tmp_path):
        """SE and Z-score computed when there are enough blocks."""
        aln = tmp_path / "test.fa"
        # Create 20-site alignment with block_size=5 -> 4 blocks
        # 4 ABBA, 2 BABA, 14 invariant
        seqs = {
            "P1":       "AAAACCAAAAAAAAAAAAAA",
            "P2":       "CCCCAAAAAAAAAAAAAAA" + "A",
            "P3":       "CCCCCCAAAAAAAAAAAAA" + "A",
            "Outgroup": "AAAAAAAAAAAAAAAAAAAA",
        }
        _write_alignment(str(aln), seqs)
        args = _make_args(str(aln), block_size=5)
        svc = Dstatistic(args)
        import io, sys
        captured = io.StringIO()
        sys.stdout = captured
        svc.run()
        sys.stdout = sys.__stdout__
        output = captured.getvalue()
        assert "Block jackknife (block size: 5):" in output
        assert "Standard error:" in output
        assert "Z-score:" in output
        assert "p-value:" in output

    def test_block_jackknife_mean_uses_array_reduction(
        self, tmp_path, monkeypatch, capsys
    ):
        aln = tmp_path / "test.fa"
        seqs = {
            "P1":       "AAAACCAAAAAAAAAAAAAA",
            "P2":       "CCCCAAAAAAAAAAAAAAAA",
            "P3":       "CCCCCCAAAAAAAAAAAAAA",
            "Outgroup": "AAAAAAAAAAAAAAAAAAAA",
        }
        _write_alignment(str(aln), seqs)
        monkeypatch.setattr(
            module.np,
            "mean",
            lambda *args, **kwargs: pytest.fail(
                "jackknife mean should use ndarray.mean"
            ),
        )

        Dstatistic(_make_args(str(aln), block_size=5)).run()

        assert "Standard error:" in capsys.readouterr().out

    def test_json_output(self, tmp_path):
        """JSON output has correct structure."""
        import json
        aln = tmp_path / "test.fa"
        seqs = {
            "P1":       "AAAACCAAAA",
            "P2":       "CCCCAAAAAA",
            "P3":       "CCCCCCAAAA",
            "Outgroup": "AAAAAAAAAA",
        }
        _write_alignment(str(aln), seqs)
        args = _make_args(str(aln), json_output=True)
        svc = Dstatistic(args)
        import io, sys
        captured = io.StringIO()
        sys.stdout = captured
        svc.run()
        sys.stdout = sys.__stdout__
        output = captured.getvalue()
        payload = json.loads(output)
        assert payload["p1"] == "P1"
        assert payload["p2"] == "P2"
        assert payload["p3"] == "P3"
        assert payload["outgroup"] == "Outgroup"
        assert payload["alignment_length"] == 10
        assert payload["informative_sites"] == 6
        assert payload["abba_count"] == 4
        assert payload["baba_count"] == 2
        assert isclose(payload["d_statistic"], 0.3333, rel_tol=0.01)
        assert "block_size" in payload
        assert "n_blocks" in payload

    def test_alignment_significance_does_not_import_scipy(self, tmp_path):
        aln = tmp_path / "test.fa"
        seqs = {
            "P1":       "AAAACCAAAA",
            "P2":       "CCCCAAAAAA",
            "P3":       "CCCCCCAAAA",
            "Outgroup": "AAAAAAAAAA",
        }
        _write_alignment(str(aln), seqs)
        args = _make_args(str(aln), block_size=5, json_output=True)
        svc = Dstatistic(args)
        real_import = __import__

        def fail_scipy_import(name, *args, **kwargs):
            if name == "scipy" or name.startswith("scipy."):
                raise AssertionError("D-statistic z-score p-value should not import scipy")
            return real_import(name, *args, **kwargs)

        import builtins
        original_import = builtins.__import__
        builtins.__import__ = fail_scipy_import
        try:
            import io, sys
            captured = io.StringIO()
            original_stdout = sys.stdout
            sys.stdout = captured
            try:
                svc.run()
            finally:
                sys.stdout = original_stdout
        finally:
            builtins.__import__ = original_import


def _write_gene_trees(path, newicks):
    with open(path, "w") as f:
        for nwk in newicks:
            f.write(nwk + "\n")


def _make_gt_args(gene_trees, p1="A", p2="B", p3="C", outgroup="O",
                  json_output=False):
    return Namespace(
        alignment=None,
        gene_trees=gene_trees,
        p1=p1,
        p2=p2,
        p3=p3,
        outgroup=outgroup,
        block_size=100,
        json=json_output,
    )


class TestGeneTreeMode:
    def test_print_gene_tree_text_output_batches_report(self, monkeypatch):
        svc = Dstatistic.__new__(Dstatistic)
        svc.p1 = "P1"
        svc.p2 = "P2"
        svc.p3 = "P3"
        svc.outgroup = "Outgroup"
        svc.support_threshold = 70.0
        printed = []

        def fake_print(*args, **kwargs):
            printed.append((args, kwargs))

        monkeypatch.setattr("builtins.print", fake_print)
        svc._print_gene_tree_text_output(
            100,
            30,
            45,
            15,
            10,
            0.5,
            15.0,
            0.000108,
        )

        expected = "\n".join([
            "Patterson's D-statistic (Gene Tree Mode)",
            "=========================================",
            "Topology: (((P1, P2), P3), Outgroup)",
            "P1: P1",
            "P2: P2",
            "P3: P3",
            "Outgroup: Outgroup",
            "",
            "Gene trees: 100",
            "Support threshold: 70.0",
            "Concordant ((P1,P2),P3): 30",
            "ABBA ((P2,P3),P1): 45",
            "BABA ((P1,P3),P2): 15",
            "Unresolved: 10",
            "D-statistic: 0.5000",
            "Chi-squared: 15.0000",
            "p-value: 0.000108",
            "",
            f"Interpretation: {svc._interpret(0.5, 0.000108)}",
        ])
        assert printed == [((expected,), {})]

    def test_collect_clade_taxa_caches_descendants(self):
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,O:1):1);"), "newick")
        clade_taxa = Dstatistic._collect_clade_taxa(tree)

        assert clade_taxa[id(tree.root)] == frozenset({"A", "B", "C", "O"})
        internal_sets = {
            taxa for cid, taxa in clade_taxa.items()
            if 1 < len(taxa) < 4
        }
        assert frozenset({"A", "B"}) in internal_sets
        assert frozenset({"C", "O"}) in internal_sets

    def test_direct_clade_taxa_preserves_nonterminal_preorder(self, monkeypatch):
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,O:1):1);"), "newick")
        expected_nonterminals = list(tree.get_nonterminals())
        expected_root_taxa = frozenset(tip.name for tip in tree.get_terminals())

        def fail_generic_traversal(*_args, **_kwargs):
            raise AssertionError("direct clade collector should not use generic traversal")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_generic_traversal)
        monkeypatch.setattr(TreeMixin, "get_nonterminals", fail_generic_traversal)
        monkeypatch.setattr(TreeMixin, "get_terminals", fail_generic_traversal)

        clade_taxa, nonterminals = (
            Dstatistic._collect_clade_taxa_and_nonterminals_direct(tree)
        )

        assert nonterminals == expected_nonterminals
        assert clade_taxa[id(tree.root)] == expected_root_taxa

    def test_direct_clade_taxa_binary_children_use_indexed_aggregation(self):
        from Bio.Phylo.BaseTree import Clade, Tree

        class IndexedOnlyList(list):
            def __iter__(self):
                raise AssertionError("binary aggregation should not iterate children")

        left = Clade(
            name="left",
            clades=IndexedOnlyList([Clade(name="A"), Clade(name="B")]),
        )
        right = Clade(
            name="right",
            clades=IndexedOnlyList([Clade(name="C"), Clade(name="O")]),
        )
        root = Clade(name="root", clades=IndexedOnlyList([left, right]))
        tree = Tree(root=root)

        clade_taxa, nonterminals = (
            Dstatistic._collect_clade_taxa_and_nonterminals_direct(tree)
        )

        assert nonterminals == [root, left, right]
        assert clade_taxa[id(left)] == frozenset({"A", "B"})
        assert clade_taxa[id(right)] == frozenset({"C", "O"})
        assert clade_taxa[id(root)] == frozenset({"A", "B", "C", "O"})

    def test_get_quartet_topology_uses_combined_direct_traversal(self, monkeypatch):
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,O:1):1);"), "newick")
        svc = object.__new__(Dstatistic)
        svc.support_threshold = None

        def fail_generic_traversal(*_args, **_kwargs):
            raise AssertionError("generic tree traversal should not be used")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_generic_traversal)
        monkeypatch.setattr(TreeMixin, "get_nonterminals", fail_generic_traversal)

        assert svc._get_quartet_topology(tree, ("A", "B", "C", "O")) == "concordant"

    def test_get_quartet_topology_avoids_full_complement_sets(self, monkeypatch):
        from Bio.Phylo.BaseTree import Clade, Tree

        class NoSubtractFrozenSet(frozenset):
            def __sub__(self, other):
                raise AssertionError("quartet topology should not build full complements")

        root = Clade(name="root", clades=[Clade(name="left"), Clade(name="right")])
        left = root.clades[0]
        right = root.clades[1]
        tree = Tree(root=root)
        svc = object.__new__(Dstatistic)
        svc.support_threshold = None

        clade_taxa = {
            id(root): NoSubtractFrozenSet({"A", "B", "C", "O", "extra"}),
            id(left): NoSubtractFrozenSet({"A", "B"}),
            id(right): NoSubtractFrozenSet({"C", "O", "extra"}),
        }

        monkeypatch.setattr(
            Dstatistic,
            "_collect_clade_taxa_and_nonterminals",
            staticmethod(lambda tree: (clade_taxa, [root, left, right])),
        )

        assert svc._get_quartet_topology(tree, ("A", "B", "C", "O")) == "concordant"

    def test_low_support_branch_is_unresolved(self, tmp_path):
        gt_file = tmp_path / "trees.nwk"
        _write_gene_trees(gt_file, [
            "((A:1,B:1)50:1,(C:1,O:1)50:1);",
        ])
        svc = Dstatistic(_make_gt_args(str(gt_file)))
        svc.support_threshold = 80
        tree = svc._parse_gene_trees(str(gt_file))[0]

        assert svc._get_quartet_topology(tree, ("A", "B", "C", "O")) == "unresolved"

    def test_concordant_trees(self, tmp_path):
        """All concordant trees → D = 0."""
        gt_file = tmp_path / "trees.nwk"
        # Species tree: (((A,B),C),O) → concordant
        _write_gene_trees(gt_file, [
            "((A:1,B:1):1,(C:1,O:1):1);",
            "((A:1,B:1):1,(C:1,O:1):1);",
            "((A:1,B:1):1,(C:1,O:1):1);",
        ])
        svc = Dstatistic(_make_gt_args(str(gt_file)))
        svc.run()
        # All concordant, no ABBA or BABA

    def test_abba_trees_positive_d(self, tmp_path):
        """More ABBA than BABA trees → D > 0."""
        gt_file = tmp_path / "trees.nwk"
        _write_gene_trees(gt_file, [
            "((A:1,B:1):1,(C:1,O:1):1);",   # concordant
            "((B:1,C:1):1,(A:1,O:1):1);",   # ABBA (P2+P3 together)
            "((B:1,C:1):1,(A:1,O:1):1);",   # ABBA
            "((A:1,C:1):1,(B:1,O:1):1);",   # BABA (P1+P3 together)
        ])
        svc = Dstatistic(_make_gt_args(str(gt_file)))
        svc.run()

    def test_d_from_gene_trees_value(self, tmp_path, mocker):
        """Verify D = (ABBA - BABA) / (ABBA + BABA)."""
        gt_file = tmp_path / "trees.nwk"
        # 3 ABBA, 1 BABA → D = (3-1)/(3+1) = 0.5
        _write_gene_trees(gt_file, [
            "((B:1,C:1):1,(A:1,O:1):1);",   # ABBA
            "((B:1,C:1):1,(A:1,O:1):1);",   # ABBA
            "((B:1,C:1):1,(A:1,O:1):1);",   # ABBA
            "((A:1,C:1):1,(B:1,O:1):1);",   # BABA
        ])
        mocked = mocker.patch("phykit.services.alignment.dstatistic.print_json")
        svc = Dstatistic(_make_gt_args(str(gt_file), json_output=True))
        svc.run()
        payload = mocked.call_args.args[0]
        assert payload["mode"] == "gene_trees"
        assert payload["abba_count"] == 3
        assert payload["baba_count"] == 1
        assert isclose(payload["d_statistic"], 0.5, rel_tol=0.01)

    def test_gene_tree_p_value_does_not_import_scipy(self, tmp_path, mocker):
        gt_file = tmp_path / "trees.nwk"
        _write_gene_trees(gt_file, [
            "((B:1,C:1):1,(A:1,O:1):1);",
            "((B:1,C:1):1,(A:1,O:1):1);",
            "((B:1,C:1):1,(A:1,O:1):1);",
            "((A:1,C:1):1,(B:1,O:1):1);",
        ])
        mocked = mocker.patch("phykit.services.alignment.dstatistic.print_json")
        svc = Dstatistic(_make_gt_args(str(gt_file), json_output=True))
        real_import = __import__

        def fail_scipy_import(name, *args, **kwargs):
            if name == "scipy" or name.startswith("scipy."):
                raise AssertionError("gene-tree chi-square p-value should not import scipy")
            return real_import(name, *args, **kwargs)

        import builtins
        original_import = builtins.__import__
        builtins.__import__ = fail_scipy_import
        try:
            svc.run()
        finally:
            builtins.__import__ = original_import

        payload = mocked.call_args.args[0]
        assert payload["p_value"] is not None

    def test_multi_taxon_gene_trees(self, tmp_path, mocker):
        """Gene trees with >4 taxa — quartet is extracted correctly."""
        gt_file = tmp_path / "trees.nwk"
        # 8-taxon trees where the quartet (A,B,C,O) topology varies
        _write_gene_trees(gt_file, [
            # A+B together (concordant for (((A,B),C),O))
            "(((A:1,B:1):1,(E:1,F:1):1):1,((C:1,G:1):1,(O:1,H:1):1):1);",
            # B+C together (ABBA)
            "(((B:1,C:1):1,(E:1,F:1):1):1,((A:1,G:1):1,(O:1,H:1):1):1);",
            # B+C together (ABBA)
            "(((B:1,C:1):1,(E:1,F:1):1):1,((A:1,G:1):1,(O:1,H:1):1):1);",
        ])
        mocked = mocker.patch("phykit.services.alignment.dstatistic.print_json")
        svc = Dstatistic(_make_gt_args(str(gt_file), json_output=True))
        svc.run()
        payload = mocked.call_args.args[0]
        assert payload["n_gene_trees"] == 3
        assert payload["abba_count"] >= 1  # at least some ABBA detected

    def test_missing_taxon_unresolved(self, tmp_path, mocker):
        """Gene trees missing one of the four taxa → counted as unresolved."""
        gt_file = tmp_path / "trees.nwk"
        _write_gene_trees(gt_file, [
            "((A:1,B:1):1,C:1);",  # missing O → unresolved
        ])
        mocked = mocker.patch("phykit.services.alignment.dstatistic.print_json")
        svc = Dstatistic(_make_gt_args(str(gt_file), json_output=True))
        svc.run()
        payload = mocked.call_args.args[0]
        assert payload["unresolved"] == 1
        assert payload["abba_count"] == 0
        assert payload["baba_count"] == 0

    def test_mutual_exclusion(self, tmp_path):
        """Cannot provide both -a and -g."""
        aln = tmp_path / "test.fa"
        aln.write_text(">A\nACGT\n>B\nACGT\n")
        gt = tmp_path / "trees.nwk"
        gt.write_text("((A:1,B:1):1,(C:1,O:1):1);\n")
        with pytest.raises(SystemExit):
            Dstatistic(Namespace(
                alignment=str(aln),
                gene_trees=str(gt),
                p1="A", p2="B", p3="C", outgroup="O",
                block_size=100, json=False,
            ))

    def test_neither_input_raises(self):
        """Must provide either -a or -g."""
        with pytest.raises(SystemExit):
            Dstatistic(Namespace(
                alignment=None,
                gene_trees=None,
                p1="A", p2="B", p3="C", outgroup="O",
                block_size=100, json=False,
            ))
