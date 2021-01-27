import statistics as stat
from typing import Tuple

from .base import Tree

class SpuriousSequence(Tree):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        tree = self.read_tree_file()
        factor = self.factor
        name_and_branch_len, threshold, median = self.identify_spurious_sequence(tree, factor)
        
        counter = 0
        for name, length in name_and_branch_len.items():
            if length >= threshold:
                try:
                    print(f"{name}\t{round(length, 4)}\t{round(threshold, 4)}\t{round(median, 4)}")
                except BrokenPipeError:
                    pass
                counter += 1
        
        # if no terminal branch is longer than the one specified
        # inform the user and print "None"
        if counter == 0:
            print("None")


    def process_args(self, args):
        return dict(tree_file_path=args.tree, factor=args.factor)

    def identify_spurious_sequence(self, tree, factor):
        branch_lengths, name_and_branch_len = self.get_branch_lengths_and_their_names(tree)

        median = stat.median(branch_lengths)

        if factor is None:
            factor = 20

        threshold = median*factor

        return name_and_branch_len, threshold, median
    
    def get_branch_lengths_and_their_names(
        self,
        tree
    ) -> Tuple[list, list]:
        # initialize a list to hold branch lengths and a
        # dictionary with terminal names and the branch
        # lengths that lead up to them.
        branch_lengths = []
        name_and_branch_len = {}
        
        # collect terminal branch lengths
        for terminal in tree.get_terminals():
            branch_lengths.append(terminal.branch_length)
            name_and_branch_len[terminal.name] = terminal.branch_length
        
        # collect internal branch lengths
        for internal in tree.get_nonterminals():
            branch_lengths.append(terminal.branch_length)

        return branch_lengths, name_and_branch_len      

