from .base import Tree
from ...helpers.stats_summary import (
    calculate_summary_statistics_from_arr,
    print_summary_statistics,
)


class BipartitionSupportStats(Tree):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        tree = self.read_tree_file()
        bs_vals, term_names = self.get_bipartition_support_vals(tree)

        if self.verbose:
            try:
                for bs_val, names in zip(bs_vals, term_names):
                    print(bs_val, ";".join(names))
            except BrokenPipeError:
                pass
        else:
            stats = calculate_summary_statistics_from_arr(bs_vals)
            print_summary_statistics(stats)

    def process_args(self, args):
        return dict(tree_file_path=args.tree, verbose=args.verbose)

    def get_bipartition_support_vals(self, tree):
        # initialize list to hold bootstrap values
        bs_vals = []
        term_names = []

        # populate bs_vals with bootstrap values
        for nonterminal in tree.get_nonterminals():
            # only include if a bootstrap value is present
            if nonterminal.confidence != None:
                bs_vals.append(nonterminal.confidence)
                temp = []
                for term in nonterminal.get_terminals():
                    temp.append(term.name)
                term_names.append(temp)
        return bs_vals, term_names
