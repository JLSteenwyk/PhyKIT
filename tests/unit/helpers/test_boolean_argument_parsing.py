import argparse

import pytest

from phykit.helpers.boolean_argument_parsing import str2bool

class TestBooleanHandling(object):
    def test_str2bool_true_boolean(self):
        v = True
        v = str2bool(v)
        assert v == True

    def test_str2bool_false_boolean(self):
        v = False
        v = str2bool(v)
        assert v == False

    def test_str2bool_true_str_boolean(self):
        v = 'true'
        v = str2bool(v)
        assert v == True

    def test_str2bool_false_str_boolean(self):
        v = 'false'
        v = str2bool(v)
        assert v == False

    def test_str2bool_true_str_boolean(self):
        v = 't'
        v = str2bool(v)
        assert v == True

    def test_str2bool_false_str_boolean(self):
        v = 'f'
        v = str2bool(v)
        assert v == False

    def test_str2bool_1_str_boolean(self):
        v = '1'
        v = str2bool(v)
        assert v == True

    def test_str2bool_0_str_boolean(self):
        v = '0'
        v = str2bool(v)
        assert v == False

    def test_str2bool_argument_type_error(self):
        v = 'Not_valid'
        with pytest.raises(Exception) as excinfo:
            v = str2bool(v)
        assert "Boolean value expected." in str(excinfo.value)
        

