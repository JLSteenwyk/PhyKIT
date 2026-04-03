from argparse import Namespace
import json

import pytest

from phykit.services.tree.chronogram import Chronogram


TREE = "tests/sample_files/ultrametric_tree.tre"


def _make_args(**overrides):
    defaults = dict(
        tree=TREE, root_age=70.0, plot_output=None, timescale="auto",
        node_ages=False, json=False,
        fig_width=None, fig_height=None, dpi=150, no_title=False,
        title=None, legend_position=None, ylabel_fontsize=None,
        xlabel_fontsize=None, title_fontsize=None, axis_fontsize=None,
        colors=None, ladderize=False, cladogram=False, circular=False,
        color_file=None,
    )
    defaults.update(overrides)
    return Namespace(**defaults)


class TestChronogram:
    def test_rectangular_creates_file(self, tmp_path):
        out = str(tmp_path / "chrono.png")
        svc = Chronogram(_make_args(plot_output=out))
        svc.run()
        assert (tmp_path / "chrono.png").exists()

    def test_circular_creates_file(self, tmp_path):
        out = str(tmp_path / "chrono_circ.png")
        svc = Chronogram(_make_args(plot_output=out, circular=True))
        svc.run()
        assert (tmp_path / "chrono_circ.png").exists()

    def test_node_ages_flag(self, tmp_path):
        out = str(tmp_path / "chrono_ages.png")
        svc = Chronogram(_make_args(plot_output=out, node_ages=True))
        svc.run()
        assert (tmp_path / "chrono_ages.png").exists()

    def test_timescale_epoch(self, tmp_path):
        out = str(tmp_path / "chrono_epoch.png")
        svc = Chronogram(_make_args(plot_output=out, timescale="epoch"))
        svc.run()
        assert (tmp_path / "chrono_epoch.png").exists()

    def test_timescale_period(self, tmp_path):
        out = str(tmp_path / "chrono_period.png")
        svc = Chronogram(_make_args(
            plot_output=out, timescale="period", root_age=250.0,
        ))
        svc.run()
        assert (tmp_path / "chrono_period.png").exists()

    def test_timescale_era(self, tmp_path):
        out = str(tmp_path / "chrono_era.png")
        svc = Chronogram(_make_args(
            plot_output=out, timescale="era", root_age=500.0,
        ))
        svc.run()
        assert (tmp_path / "chrono_era.png").exists()

    def test_json_output(self, tmp_path, capsys):
        out = str(tmp_path / "chrono.png")
        svc = Chronogram(_make_args(plot_output=out, json=True))
        svc.run()
        captured = capsys.readouterr()
        payload = json.loads(captured.out)
        assert payload["root_age"] == 70.0
        assert payload["n_tips"] == 8
        assert len(payload["node_ages"]) > 0
        # Check root node age is ~70
        root_entry = [n for n in payload["node_ages"] if n["n_descendants"] == 8]
        assert len(root_entry) == 1
        assert root_entry[0]["age_ma"] == pytest.approx(70.0, abs=0.1)

    def test_ladderize(self, tmp_path):
        out = str(tmp_path / "chrono_ladder.png")
        svc = Chronogram(_make_args(plot_output=out, ladderize=True))
        svc.run()
        assert (tmp_path / "chrono_ladder.png").exists()

    def test_non_ultrametric_tree(self, tmp_path):
        """Non-ultrametric tree should still work (scaled by root age)."""
        out = str(tmp_path / "chrono_nonultra.png")
        svc = Chronogram(_make_args(
            tree="tests/sample_files/tree_simple.tre",
            plot_output=out, root_age=100.0,
        ))
        svc.run()
        assert (tmp_path / "chrono_nonultra.png").exists()

    def test_beast_hpd_auto_detected(self, tmp_path):
        """BEAST-style HPD annotations should be auto-detected and plotted."""
        out = str(tmp_path / "chrono_beast.png")
        svc = Chronogram(_make_args(
            tree="tests/sample_files/beast_annotated_tree.tre",
            plot_output=out, root_age=70.0,
        ))
        svc.run()
        assert (tmp_path / "chrono_beast.png").exists()

    def test_hpd_in_json(self, tmp_path, capsys):
        """JSON should include HPD intervals when present."""
        out = str(tmp_path / "chrono_beast.png")
        svc = Chronogram(_make_args(
            tree="tests/sample_files/beast_annotated_tree.tre",
            plot_output=out, root_age=70.0, json=True,
        ))
        svc.run()
        captured = capsys.readouterr()
        payload = json.loads(captured.out)
        # At least some nodes should have HPD
        hpd_nodes = [n for n in payload["node_ages"] if "hpd_lower" in n]
        assert len(hpd_nodes) > 0
        for n in hpd_nodes:
            assert n["hpd_lower"] < n["hpd_upper"]

    def test_beast_hpd_circular(self, tmp_path):
        """BEAST HPD should work in circular mode too."""
        out = str(tmp_path / "chrono_beast_circ.png")
        svc = Chronogram(_make_args(
            tree="tests/sample_files/beast_annotated_tree.tre",
            plot_output=out, root_age=70.0, circular=True,
        ))
        svc.run()
        assert (tmp_path / "chrono_beast_circ.png").exists()
