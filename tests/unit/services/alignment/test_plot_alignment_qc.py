from argparse import Namespace

from phykit.services.alignment.plot_alignment_qc import PlotAlignmentQC


class TestPlotAlignmentQC:
    def test_init(self):
        args = Namespace(
            alignment="x.fa",
            output="out.png",
            width=14.0,
            height=10.0,
            dpi=300,
            gap_z=3.0,
            composition_z=3.0,
            distance_z=3.0,
            rcvt_z=3.0,
            occupancy_z=3.0,
            entropy_z=3.0,
            json=False,
        )
        svc = PlotAlignmentQC(args)
        assert svc.alignment_file_path == "x.fa"
        assert svc.output == "out.png"
        assert svc.width == 14.0
        assert svc.height == 10.0
        assert svc.dpi == 300
        assert svc.gap_z == 3.0
        assert svc.composition_z == 3.0
        assert svc.distance_z == 3.0
        assert svc.rcvt_z == 3.0
        assert svc.occupancy_z == 3.0
        assert svc.entropy_z == 3.0
