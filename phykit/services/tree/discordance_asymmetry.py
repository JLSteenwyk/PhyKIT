from io import StringIO
from pathlib import Path
from typing import Dict, List, Optional

from Bio import Phylo

from .base import Tree
from ...helpers.json_output import print_json
from ...errors import PhykitUserError


class DiscordanceAsymmetry(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"])
        self.gene_trees_path = parsed["gene_trees_path"]
        self.verbose = parsed["verbose"]
        self.json_output = parsed["json_output"]
        self.plot_output = parsed["plot_output"]

    def process_args(self, args) -> Dict:
        return dict(
            tree_file_path=args.tree,
            gene_trees_path=args.gene_trees,
            verbose=getattr(args, "verbose", False),
            json_output=getattr(args, "json", False),
            plot_output=getattr(args, "plot_output", None),
        )

    def run(self) -> None:
        pass
