import argparse
import subprocess
import sys

import pytest

from phykit.helpers.boolean_argument_parsing import str2bool


def test_import_and_valid_parse_do_not_import_argparse():
    code = """
import sys
from phykit.helpers.boolean_argument_parsing import str2bool
assert str2bool("true") is True
assert str2bool("false") is False
assert "argparse" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


class LowerCountingStr(str):
    def __new__(cls, value):
        obj = str.__new__(cls, value)
        obj.lower_calls = 0
        return obj

    def lower(self):
        self.lower_calls += 1
        return super().lower()


def test_str2bool_normalizes_string_once():
    v = LowerCountingStr("false")
    assert str2bool(v) is False
    assert v.lower_calls == 1


class TestBooleanHandling(object):
    def test_str2bool_true_boolean(self):
        v = True
        v = str2bool(v)
        assert v is True

    def test_str2bool_false_boolean(self):
        v = False
        v = str2bool(v)
        assert v is False

    def test_str2bool_true_str_boolean(self):
        v = 'true'
        v = str2bool(v)
        assert v is True

    def test_str2bool_false_str_boolean(self):
        v = 'false'
        v = str2bool(v)
        assert v is False

    def test_str2bool_t_str_boolean(self):
        v = 't'
        v = str2bool(v)
        assert v is True

    def test_str2bool_f_str_boolean(self):
        v = 'f'
        v = str2bool(v)
        assert v is False

    def test_str2bool_1_str_boolean(self):
        v = '1'
        v = str2bool(v)
        assert v is True

    def test_str2bool_0_str_boolean(self):
        v = '0'
        v = str2bool(v)
        assert v is False

    def test_str2bool_argument_type_error(self):
        v = 'Not_valid'
        with pytest.raises(argparse.ArgumentTypeError) as excinfo:
            v = str2bool(v)
        assert "Boolean value expected." in str(excinfo.value)
