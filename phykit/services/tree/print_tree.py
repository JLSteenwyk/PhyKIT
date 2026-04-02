from typing import Dict
import pickle

from Bio import Phylo

from .base import Tree
from ...helpers.json_output import print_json


class PrintTree(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(
            tree_file_path=parsed["tree_file_path"],
            remove=parsed["remove"],
        )
        self.json_output = parsed["json_output"]

    def run(self):
        tree = self.read_tree_file()

        if self.remove:
            # Make a deep copy to avoid modifying the cached tree
            tree = pickle.loads(pickle.dumps(tree, protocol=pickle.HIGHEST_PROTOCOL))
            for node in tree.get_terminals() + tree.get_nonterminals():
                node.branch_length = None

        if self.json_output:
            print_json(
                dict(
                    remove_branch_lengths=self.remove,
                    tree_newick=tree.format("newick").strip(),
                )
            )
            return

        try:
            Phylo.draw_ascii(tree)
        except BrokenPipeError:
            pass

    def process_args(self, args) -> Dict[str, str]:
        return dict(
            tree_file_path=args.tree,
            remove=args.remove,
            json_output=getattr(args, "json", False),
        )
