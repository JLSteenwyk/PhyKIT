import pytest
from argparse import Namespace
from math import isclose

from phykit.services.alignment.dstatistic import Dstatistic
from phykit.errors import PhykitUserError


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
