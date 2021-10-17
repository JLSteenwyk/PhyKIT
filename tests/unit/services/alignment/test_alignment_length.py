import pytest
from argparse import Namespace
from mock import patch, call

from phykit.services.alignment.alignment_length import AlignmentLength
from phykit.services.alignment.base import Alignment


@pytest.fixture
def args():
    kwargs = dict(alignment="/some/path/to/file.fa")
    return Namespace(**kwargs)


class TestAlignmentLength(object):
    def test_init_sets_alignment_file_path(self, args):
        aln = AlignmentLength(args)
        assert aln.alignment_file_path == args.alignment
        assert aln.output_file_path is None

    def test_alignment_length_is_printed(self, mocker, args):
        expected_length = "6"
        aln = mocker.MagicMock(
            get_alignment_length=mocker.MagicMock(return_value=expected_length)
        )
        mocked_print = mocker.patch("builtins.print")
        mocked_get_alignment_and_format = mocker.patch("phykit.services.alignment.alignment_length.AlignmentLength.get_alignment_and_format", return_value=(aln, ''))
        aln_len = AlignmentLength(args)
        res = aln_len.run()

        assert mocked_get_alignment_and_format.called
        assert mocked_print.mock_calls == [
            call(expected_length)
        ]
        
