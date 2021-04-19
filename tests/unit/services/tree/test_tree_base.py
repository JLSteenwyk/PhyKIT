import pytest
from argparse import Namespace

from phykit.services.tree.base import Tree

class TestBaseTreeService(object):
    def test_file_not_found_read_tree(self):
        with pytest.raises(TypeError):
            Tree().read_tree_file()
    
    def test_file_not_found_read_tree1(self):
        with pytest.raises(TypeError):
            Tree().read_tree1_file()
    
    def test_file_not_found_read_ref_tree(self):
        with pytest.raises(TypeError):
            Tree().read_reference_tree_file()

    def test_no_common_tips(self):
        with pytest.raises(SystemExit):
            Tree().shared_tips(['a'], ['b'])