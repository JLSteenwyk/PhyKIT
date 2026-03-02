import pytest
from phykit.services.tree.relative_rate_test import RelativeRateTest


class TestTajimaTest:
    def test_equal_sequences(self):
        """Identical ingroup sequences: m1=m2=0, chi2=0, p=1.0."""
        result = RelativeRateTest._tajima_test(
            seq_x="ACGTACGT",
            seq_y="ACGTACGT",
            seq_out="ACGTACGT",
        )
        assert result["m1"] == 0
        assert result["m2"] == 0
        assert result["chi2"] == 0.0
        assert result["p_value"] == 1.0

    def test_one_lineage_faster(self):
        """One lineage has unique changes, the other doesn't."""
        # x matches outgroup, y differs at 2 sites
        result = RelativeRateTest._tajima_test(
            seq_x="ACGTACGT",
            seq_y="ACTTACTT",
            seq_out="ACGTACGT",
        )
        assert result["m1"] == 0
        assert result["m2"] == 2
        assert result["chi2"] == 2.0
        assert abs(result["p_value"] - 0.1573) < 0.001

    def test_symmetric_changes(self):
        """Equal unique changes on both lineages: chi2=0, p=1."""
        result = RelativeRateTest._tajima_test(
            seq_x="TCGTACGT",
            seq_y="ACGTACTT",
            seq_out="ACGTACGT",
        )
        assert result["m1"] == 1
        assert result["m2"] == 1
        assert result["chi2"] == 0.0
        assert result["p_value"] == 1.0

    def test_highly_unequal_rates(self):
        """Many changes on one lineage: should reject clock."""
        # seq_y differs from outgroup at 6 sites (positions 0,1,2,4,5,6)
        # positions 3,7 have T in all three sequences (no difference)
        result = RelativeRateTest._tajima_test(
            seq_x="ACGTACGTACGT",
            seq_y="TTTTTTTTACGT",
            seq_out="ACGTACGTACGT",
        )
        assert result["m1"] == 0
        assert result["m2"] == 6
        assert result["chi2"] == 6.0
        assert result["p_value"] < 0.02

    def test_gaps_are_skipped(self):
        """Sites with gaps in any sequence should be excluded."""
        result = RelativeRateTest._tajima_test(
            seq_x="A-GTACGT",
            seq_y="ACTTACGT",
            seq_out="ACGTACGT",
        )
        # Position 1 has gap in x, skip it. Position 2: y=T, out=G, x=G -> m2.
        # Only non-gap informative site where y differs: pos 2 (y=T, out=G, x=G)
        assert result["m1"] == 0
        assert result["m2"] == 1

    def test_shared_changes_not_counted(self):
        """Sites where BOTH ingroup differ from outgroup are uninformative."""
        result = RelativeRateTest._tajima_test(
            seq_x="TCGTACGT",
            seq_y="TCGTACGT",
            seq_out="ACGTACGT",
        )
        # Position 0: both x and y differ from outgroup -> shared, skip
        assert result["m1"] == 0
        assert result["m2"] == 0
        assert result["chi2"] == 0.0


class TestIdentifyOutgroup:
    def test_simple_rooted_tree(self):
        """Outgroup should be the earliest-diverging taxon."""
        from io import StringIO
        from Bio import Phylo

        tree = Phylo.read(StringIO("(((A:0.1,B:0.1):0.05,C:0.15):0.1,O:0.2);"), "newick")
        outgroup = RelativeRateTest._identify_outgroup(tree)
        assert outgroup == "O"

    def test_two_taxa_at_root(self):
        """With two clades at root, outgroup is the smaller one."""
        from io import StringIO
        from Bio import Phylo

        tree = Phylo.read(StringIO("((A:0.1,B:0.1):0.1,O:0.2);"), "newick")
        outgroup = RelativeRateTest._identify_outgroup(tree)
        assert outgroup == "O"


class TestMultipleTestingCorrection:
    def test_bonferroni(self):
        p_values = [0.01, 0.04, 0.06]
        corrected = RelativeRateTest._bonferroni(p_values)
        assert corrected[0] == 0.03  # 0.01 * 3
        assert corrected[1] == 0.12  # 0.04 * 3
        assert corrected[2] == 0.18  # 0.06 * 3

    def test_bonferroni_capped_at_one(self):
        p_values = [0.5, 0.8]
        corrected = RelativeRateTest._bonferroni(p_values)
        assert corrected[0] == 1.0
        assert corrected[1] == 1.0

    def test_fdr(self):
        p_values = [0.01, 0.04, 0.06]
        corrected = RelativeRateTest._fdr(p_values)
        # BH-FDR: rank 1: 0.01*3/1=0.03, rank 2: 0.04*3/2=0.06, rank 3: 0.06*3/3=0.06
        assert abs(corrected[0] - 0.03) < 1e-10
        assert abs(corrected[1] - 0.06) < 1e-10
        assert abs(corrected[2] - 0.06) < 1e-10


class TestMatchesRReference:
    """Validate against R's pegas::rr.test() reference values.

    R code used to generate reference values:
        library(pegas); library(ape)
        aln <- read.FASTA("rrt_test.fa")
        rr.test(aln[["A"]], aln[["B"]], aln[["O"]])
    """

    def test_a_vs_b_outgroup_o(self):
        """R: Chi=4.0, Pval=0.0455002639."""
        result = RelativeRateTest._tajima_test(
            seq_x="ACGTACGTACGTACGTACGT",
            seq_y="ACGTACGTACTAACTAACGT",
            seq_out="ACGTACGTACGTACGTACGT",
        )
        assert result["m1"] == 0
        assert result["m2"] == 4
        assert result["chi2"] == 4.0
        assert abs(result["p_value"] - 0.0455002639) < 1e-8

    def test_a_vs_c_outgroup_o(self):
        """R: Chi=2.0, Pval=0.1572992071."""
        result = RelativeRateTest._tajima_test(
            seq_x="ACGTACGTACGTACGTACGT",
            seq_y="ACTAACGTACGTACGTACGT",
            seq_out="ACGTACGTACGTACGTACGT",
        )
        assert result["m1"] == 0
        assert result["m2"] == 2
        assert result["chi2"] == 2.0
        assert abs(result["p_value"] - 0.1572992071) < 1e-8

    def test_b_vs_c_outgroup_o(self):
        """R: Chi=0.6667, Pval=0.4142161782."""
        result = RelativeRateTest._tajima_test(
            seq_x="ACGTACGTACTAACTAACGT",
            seq_y="ACTAACGTACGTACGTACGT",
            seq_out="ACGTACGTACGTACGTACGT",
        )
        assert result["m1"] == 4
        assert result["m2"] == 2
        assert abs(result["chi2"] - 0.6666666667) < 1e-8
        assert abs(result["p_value"] - 0.4142161782) < 1e-8

    def test_one_lineage_faster_matches_r(self):
        """R: Chi=2.0, Pval=0.1572992071."""
        result = RelativeRateTest._tajima_test(
            seq_x="ACGTACGT",
            seq_y="ACTTACTT",
            seq_out="ACGTACGT",
        )
        assert result["m1"] == 0
        assert result["m2"] == 2
        assert result["chi2"] == 2.0
        assert abs(result["p_value"] - 0.1572992071) < 1e-8

    def test_equal_sequences_edge_case(self):
        """R returns NaN for m1+m2=0; we return chi2=0, p=1.0."""
        result = RelativeRateTest._tajima_test(
            seq_x="ACGTACGT",
            seq_y="ACGTACGT",
            seq_out="ACGTACGT",
        )
        assert result["chi2"] == 0.0
        assert result["p_value"] == 1.0
