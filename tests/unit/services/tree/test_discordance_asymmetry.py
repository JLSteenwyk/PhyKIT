import os
import sys
import tempfile

import pytest
from mock import patch
from argparse import Namespace

from phykit.errors import PhykitUserError

TREE_SIMPLE = "tests/sample_files/tree_simple.tre"
GENE_TREES = "tests/sample_files/gene_trees_simple.nwk"


class TestProcessArgs:
    def test_default_args(self):
        from phykit.services.tree.discordance_asymmetry import DiscordanceAsymmetry

        args = Namespace(
            tree=TREE_SIMPLE,
            gene_trees=GENE_TREES,
            verbose=False,
            json=False,
            plot_output=None,
        )
        svc = DiscordanceAsymmetry(args)
        assert svc.gene_trees_path == GENE_TREES
        assert svc.verbose is False
        assert svc.json_output is False
        assert svc.plot_output is None

    def test_verbose_arg(self):
        from phykit.services.tree.discordance_asymmetry import DiscordanceAsymmetry

        args = Namespace(
            tree=TREE_SIMPLE,
            gene_trees=GENE_TREES,
            verbose=True,
            json=False,
            plot_output=None,
        )
        svc = DiscordanceAsymmetry(args)
        assert svc.verbose is True

    def test_json_arg(self):
        from phykit.services.tree.discordance_asymmetry import DiscordanceAsymmetry

        args = Namespace(
            tree=TREE_SIMPLE,
            gene_trees=GENE_TREES,
            verbose=False,
            json=True,
            plot_output=None,
        )
        svc = DiscordanceAsymmetry(args)
        assert svc.json_output is True

    def test_plot_output_arg(self):
        from phykit.services.tree.discordance_asymmetry import DiscordanceAsymmetry

        args = Namespace(
            tree=TREE_SIMPLE,
            gene_trees=GENE_TREES,
            verbose=False,
            json=False,
            plot_output="/tmp/test.png",
        )
        svc = DiscordanceAsymmetry(args)
        assert svc.plot_output == "/tmp/test.png"
