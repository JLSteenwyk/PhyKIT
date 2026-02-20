import pytest

import phykit.services.alignment as alignment_services
import phykit.services.tree as tree_services


class TestLazyExports:
    def test_alignment_lazy_getattr_success(self):
        cls = alignment_services.AlignmentLength
        assert cls.__name__ == "AlignmentLength"
        assert "AlignmentLength" in alignment_services.__all__

    def test_alignment_lazy_getattr_invalid(self):
        with pytest.raises(AttributeError):
            _ = alignment_services.NotARealAlignmentService

    def test_tree_lazy_getattr_success(self):
        cls = tree_services.DVMC
        assert cls.__name__ == "DVMC"
        assert "DVMC" in tree_services.__all__

    def test_tree_lazy_getattr_invalid(self):
        with pytest.raises(AttributeError):
            _ = tree_services.NotARealTreeService
