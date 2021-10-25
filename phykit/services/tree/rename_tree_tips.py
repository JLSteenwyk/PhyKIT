import logging
import sys

from .base import Tree

class RenameTreeTips(Tree):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        tree = self.read_tree_file()

        # save idmap to a dictionary
        idmap = self.read_id_map(self.idmap)
        
        tree = self.replace_tip_names(tree, idmap)
        
        self.write_tree_file(tree, self.output_file_path)
    
    def process_args(self, args):
        tree_file_path = args.tree

        if args.output is None:
            output_file_path = f"{tree_file_path}.renamed"
        else:
            output_file_path = f"{args.output}"

        return dict(
            tree_file_path=tree_file_path, 
            idmap=args.idmap,
            output_file_path=output_file_path
            )

    def read_id_map(self, idmap: str) -> dict:
        """
        read two column file to dictionary
        """
        idmap={}
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

    def replace_tip_names(self, tree: Tree, idmap: dict) -> Tree:
        """
        if a tip name is in the idmap, replace it
        """
        # for tips with a name as a key in the idmap,
        # replace that tip name with the value in the
        # idmap
        for term in tree.get_terminals():
            if term.name in idmap:
                term.name = idmap[term.name]
        
        return tree

