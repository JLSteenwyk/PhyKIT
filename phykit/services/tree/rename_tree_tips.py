import logging

from Bio import Phylo

from .base import Tree

class RenameTreeTips(Tree):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        tree = self.read_tree_file()

        # save idmap to a dictionary
        idmap={}
        with open(self.idmap) as identifiers:
            for line in identifiers:
                (key, val) = line.split()
                idmap[key] = val
        
        # for tips with a name as a key in the idmap,
        # replace that tip name with the value in the
        # idmap
        for term in tree.get_terminals():
            if term.name in idmap:
                term.name = idmap[term.name]
        
        self.write_tree_file(tree, self.output_file_path)
    
    def process_args(self, args):
        tree_file_path = args.tree
        return dict(
            tree_file_path=tree_file_path, 
            idmap=args.idmap,
            output_file_path=f"{tree_file_path}.renamed.tre"
            )