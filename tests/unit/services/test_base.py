import pytest

from phykit.services.base import process_args

class TestBase(object):
    def test_base(self):
        with pytest.raises(Exception) as excinfo:
            process_args()
        assert "Boolean value expected." in str(excinfo.value)