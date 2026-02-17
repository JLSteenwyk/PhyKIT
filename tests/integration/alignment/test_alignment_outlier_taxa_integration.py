import json
import sys
from pathlib import Path

import pytest
from mock import patch

from phykit.phykit import Phykit

here = Path(__file__)


@pytest.mark.integration
class TestAlignmentOutlierTaxa:
    @patch("builtins.print")
    def test_alignment_outlier_taxa(self, mocked_print):
        testargs = [
            "phykit",
            "alignment_outlier_taxa",
            f"{here.parent.parent.parent}/sample_files/alignment_outlier_taxa.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        printed = [call.args[0] for call in mocked_print.call_args_list]
        assert (
            printed[0]
            == "features_evaluated\tgap_rate,occupancy,composition_distance,long_branch_proxy,rcvt,entropy_burden"
        )
        assert printed[1].startswith("thresholds\t")
        assert any(line.startswith("taxon_d\t") for line in printed)
        assert any(line.startswith("taxon_e\t") for line in printed)
        assert any("composition_distance" in line for line in printed if line.startswith("taxon_d\t"))
        assert any("long_branch_proxy" in line for line in printed if line.startswith("taxon_d\t"))
        assert any("gap_rate" in line for line in printed if line.startswith("taxon_e\t"))

    @patch("builtins.print")
    def test_alignment_outlier_taxa_alias(self, mocked_print):
        testargs = [
            "phykit",
            "aot",
            f"{here.parent.parent.parent}/sample_files/alignment_outlier_taxa.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        printed = [call.args[0] for call in mocked_print.call_args_list]
        assert any(line.startswith("taxon_d\t") for line in printed)

    @patch("builtins.print")
    def test_alignment_outlier_taxa_json(self, mocked_print):
        testargs = [
            "phykit",
            "alignment_outlier_taxa",
            f"{here.parent.parent.parent}/sample_files/alignment_outlier_taxa.fa",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["features_evaluated"] == [
            "gap_rate",
            "occupancy",
            "composition_distance",
            "long_branch_proxy",
            "rcvt",
            "entropy_burden",
        ]
        assert payload["rows"][0] == payload["taxa"][0]

        flagged_taxa = {row["taxon"] for row in payload["outliers"]}
        assert "taxon_d" in flagged_taxa
        assert "taxon_e" in flagged_taxa

        outlier_d = [row for row in payload["outliers"] if row["taxon"] == "taxon_d"][0]
        reason_features_d = {reason["feature"] for reason in outlier_d["reasons"]}
        assert "composition_distance" in reason_features_d
        assert "long_branch_proxy" in reason_features_d

        outlier_e = [row for row in payload["outliers"] if row["taxon"] == "taxon_e"][0]
        reason_features_e = {reason["feature"] for reason in outlier_e["reasons"]}
        assert "gap_rate" in reason_features_e

    @patch("builtins.print")
    def test_alignment_outlier_taxa_bad_path(self, mocked_print):
        testargs = [
            "phykit",
            "alignment_outlier_taxa",
            f"{here.parent.parent.parent}/sample_files/alignment_outlier_taxa.fasta.nope",
        ]
        with patch.object(sys, "argv", testargs):
            with pytest.raises(SystemExit) as pytest_wrapped_e:
                Phykit()

        assert pytest_wrapped_e.type is SystemExit
        assert pytest_wrapped_e.value.code == 2
