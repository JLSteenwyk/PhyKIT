import pickle
import sys
from typing import Dict

from Bio.Phylo import Newick

from .base import Tree
from ...helpers.json_output import print_json


class RenameTreeTips(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(
            tree_file_path=parsed["tree_file_path"],
            idmap=parsed["idmap"],
            output_file_path=parsed["output_file_path"],
        )
        self.json_output = parsed["json_output"]

    def run(self):
        tree = self.read_tree_file()
        # Make a deep copy to avoid modifying the cached tree
        tree_copy = pickle.loads(pickle.dumps(tree, protocol=pickle.HIGHEST_PROTOCOL))

        idmap = self.read_id_map()

        tree_copy, renamed_count = self.replace_tip_names(tree_copy, idmap)

        self.write_tree_file(tree_copy, self.output_file_path)

        if self.json_output:
            print_json(
                dict(
                    input_tree=self.tree_file_path,
                    idmap=self.idmap,
                    output_file=self.output_file_path,
                    renamed_tips=renamed_count,
                )
            )

    def process_args(self, args) -> Dict[str, str]:
        tree_file_path = args.tree

        output_file_path = \
            f"{args.output}" if args.output else f"{tree_file_path}.renamed"

        return dict(
            tree_file_path=tree_file_path,
            idmap=args.idmap,
            output_file_path=output_file_path,
            json_output=getattr(args, "json", False),
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
                sys.exit(2)
            except BrokenPipeError:
                pass

        return idmap

    def replace_tip_names(
        self,
        tree: Tree,
        idmap: Dict[str, str]
    ) -> tuple[Newick.Tree, int]:
        renamed_count = 0
        for term in tree.get_terminals():
            name = term.name
            if name in idmap:
                term.name = idmap[name]
                renamed_count += 1

        return tree, renamed_count
