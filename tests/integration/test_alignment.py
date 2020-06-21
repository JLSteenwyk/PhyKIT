import pytest
import sys
from mock import patch, call
from pathlib import Path

from phykit.phykit import Phykit

here = Path(__file__)

@pytest.mark.integration
class TestAlignment(object):

	@patch('builtins.print')
	def test_alignment_length(self, mocked_print):
		expected_length = '6'
		testargs = ["clipkit", "alignment_length", f"{here.parent.parent}/sample_files/simple.fa"]
		with patch.object(sys, 'argv', testargs):
			Phykit()
		assert mocked_print.mock_calls == [call(expected_length)]
