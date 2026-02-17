from argparse import Namespace
from pathlib import Path

from Bio import AlignIO

from phykit.services.alignment.alignment_outlier_taxa import AlignmentOutlierTaxa


here = Path(__file__)


class TestAlignmentOutlierTaxa:
    def test_init(self):
        args = Namespace(
            alignment="test.fa",
            gap_z=3.0,
            composition_z=3.0,
            distance_z=3.0,
            rcvt_z=3.0,
            occupancy_z=3.0,
            entropy_z=3.0,
            json=False,
        )
        service = AlignmentOutlierTaxa(args)
        assert service.alignment_file_path == "test.fa"
        assert service.gap_z == 3.0
        assert service.composition_z == 3.0
        assert service.distance_z == 3.0
        assert service.rcvt_z == 3.0
        assert service.occupancy_z == 3.0
        assert service.entropy_z == 3.0

    def test_calculate_outliers_flags_expected_taxa(self):
        args = Namespace(
            alignment="unused.fa",
            gap_z=3.0,
            composition_z=3.0,
            distance_z=3.0,
            rcvt_z=3.0,
            occupancy_z=3.0,
            entropy_z=3.0,
            json=False,
        )
        service = AlignmentOutlierTaxa(args)
        alignment = AlignIO.read(
            f"{here.parent.parent.parent.parent}/sample_files/alignment_outlier_taxa.fa",
            "fasta",
        )

        result = service.calculate_outliers(alignment, is_protein=False)

        flagged = {row["taxon"] for row in result["outliers"]}
        assert "taxon_d" in flagged
        assert "taxon_e" in flagged

        outlier_d = [row for row in result["outliers"] if row["taxon"] == "taxon_d"][0]
        features_d = {reason["feature"] for reason in outlier_d["reasons"]}
        assert "composition_distance" in features_d
        assert "long_branch_proxy" in features_d

        outlier_e = [row for row in result["outliers"] if row["taxon"] == "taxon_e"][0]
        features_e = {reason["feature"] for reason in outlier_e["reasons"]}
        assert "gap_rate" in features_e
