from argparse import Namespace
import builtins

import pytest

from phykit.services.alignment.plot_alignment_qc import PlotAlignmentQC
import phykit.services.alignment.plot_alignment_qc as plot_alignment_qc_module


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

    def test_process_args_defaults_json_false(self):
        args = Namespace(
            alignment="x.fa",
            output="out.png",
            width=10.0,
            height=8.0,
            dpi=200,
            gap_z=2.0,
            composition_z=2.0,
            distance_z=2.0,
            rcvt_z=2.0,
            occupancy_z=2.0,
            entropy_z=2.0,
        )
        parsed = PlotAlignmentQC(args).process_args(args)
        assert parsed["json_output"] is False

    def test_get_outlier_result_forwards_thresholds(self, monkeypatch):
        captured = {}

        class FakeOutlierService:
            def __init__(self, ns):
                captured["ns"] = ns

            def calculate_outliers(self, alignment, is_protein):
                captured["alignment"] = alignment
                captured["is_protein"] = is_protein
                return {"rows": [], "outliers": [], "thresholds": {}, "features": []}

        monkeypatch.setattr(plot_alignment_qc_module, "AlignmentOutlierTaxa", FakeOutlierService)
        svc = PlotAlignmentQC(
            Namespace(
                alignment="my.fa",
                output="out.png",
                width=14.0,
                height=10.0,
                dpi=300,
                gap_z=1.1,
                composition_z=1.2,
                distance_z=1.3,
                rcvt_z=1.4,
                occupancy_z=1.5,
                entropy_z=1.6,
                json=False,
            )
        )

        result = svc._get_outlier_result("aln", True)

        assert result["rows"] == []
        assert captured["alignment"] == "aln"
        assert captured["is_protein"] is True
        assert captured["ns"].alignment == "my.fa"
        assert captured["ns"].gap_z == 1.1
        assert captured["ns"].composition_z == 1.2
        assert captured["ns"].distance_z == 1.3
        assert captured["ns"].rcvt_z == 1.4
        assert captured["ns"].occupancy_z == 1.5
        assert captured["ns"].entropy_z == 1.6

    def test_run_importerror_exits(self, monkeypatch, capsys):
        original_import = builtins.__import__

        def fake_import(name, globals=None, locals=None, fromlist=(), level=0):
            if name.startswith("matplotlib"):
                raise ImportError("no matplotlib")
            return original_import(name, globals, locals, fromlist, level)

        monkeypatch.setattr(builtins, "__import__", fake_import)
        svc = PlotAlignmentQC(
            Namespace(
                alignment="x.fa",
                output="out.png",
                width=10.0,
                height=8.0,
                dpi=100,
                gap_z=2.0,
                composition_z=2.0,
                distance_z=2.0,
                rcvt_z=2.0,
                occupancy_z=2.0,
                entropy_z=2.0,
                json=False,
            )
        )

        with pytest.raises(SystemExit) as exc:
            svc.run()
        assert exc.value.code == 2
        out, _ = capsys.readouterr()
        assert "matplotlib is required for plot_alignment_qc" in out

    def test_run_outputs_json(self, tmp_path, monkeypatch):
        pytest.importorskip("matplotlib")
        output_path = tmp_path / "qc.png"
        captured = {}

        class FakeAlignment:
            def __len__(self):
                return 3

            def get_alignment_length(self):
                return 12

        svc = PlotAlignmentQC(
            Namespace(
                alignment="x.fa",
                output=str(output_path),
                width=8.0,
                height=6.0,
                dpi=100,
                gap_z=3.0,
                composition_z=3.0,
                distance_z=3.0,
                rcvt_z=3.0,
                occupancy_z=3.0,
                entropy_z=3.0,
                json=True,
            )
        )
        monkeypatch.setattr(
            svc,
            "get_alignment_and_format",
            lambda: (FakeAlignment(), None, True),
        )
        monkeypatch.setattr(
            svc,
            "_get_outlier_result",
            lambda _alignment, _is_protein: {
                "rows": [
                    {
                        "taxon": "t1",
                        "occupancy": 1.0,
                        "gap_rate": 0.0,
                        "composition_distance": 0.1,
                        "long_branch_proxy": 0.2,
                        "rcvt": 0.01,
                        "entropy_burden": 0.02,
                    },
                    {
                        "taxon": "t2",
                        "occupancy": 0.7,
                        "gap_rate": 0.3,
                        "composition_distance": 0.6,
                        "long_branch_proxy": None,
                        "rcvt": 0.1,
                        "entropy_burden": 0.2,
                    },
                    {
                        "taxon": "t3",
                        "occupancy": 0.8,
                        "gap_rate": 0.2,
                        "composition_distance": 0.5,
                        "long_branch_proxy": 0.4,
                        "rcvt": 0.05,
                        "entropy_burden": 0.1,
                    },
                ],
                "outliers": [{"taxon": "t2"}],
                "thresholds": {
                    "occupancy": 0.75,
                    "gap_rate": 0.25,
                    "composition_distance": 0.4,
                    "long_branch_proxy": 0.3,
                },
                "features": [
                    "gap_rate",
                    "occupancy",
                    "composition_distance",
                    "long_branch_proxy",
                    "rcvt",
                    "entropy_burden",
                ],
            },
        )
        monkeypatch.setattr(plot_alignment_qc_module, "print_json", lambda payload: captured.setdefault("payload", payload))

        svc.run()

        assert output_path.exists()
        assert captured["payload"]["output_file"] == str(output_path)
        assert captured["payload"]["n_taxa"] == 3
        assert captured["payload"]["alignment_length"] == 12
        assert captured["payload"]["outlier_count"] == 1

    def test_run_outputs_terminal_message(self, tmp_path, monkeypatch, capsys):
        pytest.importorskip("matplotlib")
        output_path = tmp_path / "qc.png"

        class FakeAlignment:
            def __len__(self):
                return 1

            def get_alignment_length(self):
                return 4

        svc = PlotAlignmentQC(
            Namespace(
                alignment="x.fa",
                output=str(output_path),
                width=7.0,
                height=5.0,
                dpi=80,
                gap_z=2.0,
                composition_z=2.0,
                distance_z=2.0,
                rcvt_z=2.0,
                occupancy_z=2.0,
                entropy_z=2.0,
                json=False,
            )
        )
        monkeypatch.setattr(svc, "get_alignment_and_format", lambda: (FakeAlignment(), None, True))
        monkeypatch.setattr(
            svc,
            "_get_outlier_result",
            lambda _alignment, _is_protein: {
                "rows": [
                    {
                        "taxon": "t1",
                        "occupancy": 1.0,
                        "gap_rate": 0.0,
                        "composition_distance": 0.1,
                        "long_branch_proxy": 0.2,
                        "rcvt": 0.01,
                        "entropy_burden": 0.02,
                    }
                ],
                "outliers": [],
                "thresholds": {
                    "occupancy": None,
                    "gap_rate": None,
                    "composition_distance": None,
                    "long_branch_proxy": None,
                },
                "features": [
                    "gap_rate",
                    "occupancy",
                    "composition_distance",
                    "long_branch_proxy",
                    "rcvt",
                    "entropy_burden",
                ],
            },
        )

        svc.run()

        out, _ = capsys.readouterr()
        assert output_path.exists()
        assert "Saved alignment QC plot:" in out
