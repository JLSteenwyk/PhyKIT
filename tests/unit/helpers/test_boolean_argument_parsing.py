import argparse
import subprocess
import sys

import pytest

from phykit.helpers.boolean_argument_parsing import str2bool


YES_NO_CASES = [
    ("yes", True),
    ("YES", True),
    ("Yes", True),
    ("y", True),
    ("Y", True),
    ("no", False),
    ("NO", False),
    ("No", False),
    ("n", False),
    ("N", False),
]


def test_import_and_valid_parse_do_not_import_argparse():
    code = """
import sys
from phykit.helpers.boolean_argument_parsing import str2bool
for value in ("true", "t", "1", "yes", "YES", "Yes", "y", "Y", True):
    assert str2bool(value) is True
for value in ("false", "f", "0", "no", "NO", "No", "n", "N", False):
    assert str2bool(value) is False
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


@pytest.mark.parametrize(
    "value, expected",
    [("true", True), ("false", False)] + YES_NO_CASES,
)
def test_str2bool_normalizes_string_once(value, expected):
    v = LowerCountingStr(value)
    assert str2bool(v) is expected
    assert v.lower_calls == 1


@pytest.mark.parametrize("value, expected", YES_NO_CASES)
def test_str2bool_yes_no_spellings(value, expected):
    assert str2bool(value) is expected


@pytest.mark.parametrize("value, expected", YES_NO_CASES)
def test_argparse_yes_no_spellings(value, expected):
    parser = argparse.ArgumentParser()
    parser.add_argument("--flag", type=str2bool)
    args = parser.parse_args(["--flag", value])
    assert args.flag is expected


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

    @pytest.mark.parametrize("value", ["Not_valid", "maybe"])
    def test_str2bool_argument_type_error(self, value):
        v = LowerCountingStr(value)
        with pytest.raises(argparse.ArgumentTypeError) as excinfo:
            str2bool(v)
        assert str(excinfo.value) == "Boolean value expected."
        assert v.lower_calls == 1
