import pytest
from argparse import Namespace

from phykit.services.base import BaseService

class TestBaseService(object):
    def test_must_implement_process_args(self):
        with pytest.raises(NotImplementedError):
            BaseService().process_args(None)
    
    def test_must_implement_run(self):
        with pytest.raises(NotImplementedError):
            BaseService().run()