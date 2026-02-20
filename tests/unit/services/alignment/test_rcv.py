import pytest
from argparse import Namespace
from math import isclose

from phykit.services.alignment.rcv import RelativeCompositionVariability
import phykit.services.alignment.rcv as rcv_module


@pytest.fixture
def args():
    kwargs = dict(alignment="/some/path/to/file.fa")
    return Namespace(**kwargs)


class TestRelativeCompositionVariability(object):
    def test_init_sets_alignment_file_path(self, args):
        rcv = RelativeCompositionVariability(args)
        assert rcv.alignment_file_path == args.alignment
        assert rcv.output_file_path is None

    def test_relative_composition_variability(self, mocker, alignment_simple, args):
        mocker.patch(
            "phykit.services.alignment.rcv.RelativeCompositionVariability.get_alignment_and_format",
            return_value=(alignment_simple, "fa", True),
        )
        rcv = RelativeCompositionVariability(args)
        relative_composition_variability = rcv.calculate_rcv()
        assert isinstance(relative_composition_variability, float)
        assert isclose(relative_composition_variability, 0.292, rel_tol=0.001)

    def test_process_args_defaults_json_false(self):
        parsed = RelativeCompositionVariability(Namespace(alignment="x.fa")).process_args(
            Namespace(alignment="x.fa")
        )
        assert parsed["json_output"] is False

    def test_run_prints_rcv(self, mocker, capsys):
        rcv = RelativeCompositionVariability(Namespace(alignment="x.fa", json=False))
        mocker.patch.object(RelativeCompositionVariability, "calculate_rcv", return_value=0.123456)
        rcv.run()
        out, _ = capsys.readouterr()
        assert out.strip() == "0.1235"

    def test_run_json_output(self, mocker):
        rcv = RelativeCompositionVariability(Namespace(alignment="x.fa", json=True))
        mocker.patch.object(RelativeCompositionVariability, "calculate_rcv", return_value=0.56789)
        mocked_json = mocker.patch.object(rcv_module, "print_json")
        rcv.run()
        mocked_json.assert_called_once_with({"rcv": 0.5679})
