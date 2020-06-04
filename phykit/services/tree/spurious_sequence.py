import statistics as stat

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
                print(f"{name}\t{length}\t{threshold}\t{median}")
                counter += 1
        
        # if no terminal branch is longer than the one specified
        # inform the user and print "None"
        if counter == 0:
            print("None")

    def process_args(self, args):
        return dict(tree_file_path=args.tree, factor=args.factor)

    def identify_spurious_sequence(self, tree, factor):
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

        median = stat.median(branch_lengths)

        if factor is None:
            factor = 20

        threshold = median*factor

        return(name_and_branch_len, threshold, median)
            
            
