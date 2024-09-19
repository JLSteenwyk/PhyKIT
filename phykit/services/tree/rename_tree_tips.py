import sys
from typing import Dict

from Bio.Phylo import Newick

from .base import Tree


class RenameTreeTips(Tree):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        tree = self.read_tree_file()

        idmap = self.read_id_map()

        tree = self.replace_tip_names(tree, idmap)

        self.write_tree_file(tree, self.output_file_path)

    def process_args(self, args) -> Dict[str, str]:
        tree_file_path = args.tree

        output_file_path = \
            f"{args.output}" if args.output else f"{tree_file_path}.renamed"

        return dict(
            tree_file_path=tree_file_path,
            idmap=args.idmap,
            output_file_path=output_file_path,
        )

    def read_id_map(self) -> Dict[str, str]:
        idmap = dict()
        try:
            with open(self.idmap) as identifiers:
                for line in identifiers:
                    (key, val) = line.split()
                    idmap[key] = val
        except FileNotFoundError:
            try:
                print(f"{self.idmap} corresponds to no such file.")
                print("Please check file name and pathing")
                sys.exit()
            except BrokenPipeError:
                pass

        return idmap

    def replace_tip_names(
        self,
        tree: Tree,
        idmap: Dict[str, str]
    ) -> Newick.Tree:
        for term in tree.get_terminals():
            name = term.name
            if name in idmap:
                term.name = idmap[name]

        return tree
