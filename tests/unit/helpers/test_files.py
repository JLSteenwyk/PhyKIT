import argparse

import pytest

from phykit.helpers.files import (get_alignment_and_format, read_single_column_file_to_list)

class TestFileErrorHandling(object):
    def test_get_alignment_and_format_error_handling(self):
        file_path = "not_real_file"
        with pytest.raises(SystemExit) as excinfo:
            get_alignment_and_format(file_path)
        assert excinfo.type == SystemExit

    def test_get_read_single_column_file_to_list_error_handling(self):
        file_path = "not_real_file"
        with pytest.raises(SystemExit) as excinfo:
            read_single_column_file_to_list(file_path)
        assert excinfo.type == SystemExit

        
        

