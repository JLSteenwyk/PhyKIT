import sys

from .base import Tree

from ...helpers.stats_summary import calculate_summary_statistics_from_arr, print_summary_statistics

class TerminalBranchStats(Tree):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        tree = self.read_tree_file()
        terminal_branch_lengths, stats = self.calculate_terminal_branch_stats(tree)

        if self.verbose:
            try:
                for terminal_branch_length in terminal_branch_lengths:
                    print(round(terminal_branch_length, 4))
            except BrokenPipeError:
                pass
        else:
            print_summary_statistics(stats)

    def process_args(self, args):
        return dict(tree_file_path=args.tree, verbose=args.verbose)

    def get_terminal_branch_lengths(self, tree) -> list:
        """
        loop through tree and get all terminal branch lengths
        """
        terminal_branch_lengths = []
        for terminal_branch in tree.get_terminals():
            if terminal_branch.branch_length != None:
                terminal_branch_lengths.append(terminal_branch.branch_length)

        return terminal_branch_lengths

    def check_tree_has_branch_lengths(self, terminal_branch_lengths:list) -> None:
        """
        if tree has no branch lengths, exit
        """
        if len(terminal_branch_lengths) == 0:
            print("Calculating terminal branch statistics requires a phylogeny with branch lengths.")
            sys.exit()


    def calculate_terminal_branch_stats(self, tree):
        # save terminal branch lengths to terminal_branch_lengths
        terminal_branch_lengths = self.get_terminal_branch_lengths(tree)
        
        # If the phylogeny had no branch lengths, inform user and quit
        self.check_tree_has_branch_lengths(terminal_branch_lengths)
        
        # calculate summary stats
        stats = calculate_summary_statistics_from_arr(terminal_branch_lengths)

        return terminal_branch_lengths, stats