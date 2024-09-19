from multiprocessing import Pool
from typing import Dict, List, Tuple, Union

from Bio.Phylo import Newick

from .base import Tree
from ...helpers.stats_summary import (
    calculate_summary_statistics_from_arr,
    print_summary_statistics,
)


class BipartitionSupportStats(Tree):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self) -> None:
        tree = self.read_tree_file()
        bs_vals, term_names = self.get_bipartition_support_vals(tree)

        if self.verbose:
            try:
                for i in range(len(bs_vals)):
                    print(bs_vals[i], ";".join(term_names[i]))
            except BrokenPipeError:
                pass
        else:
            stats = calculate_summary_statistics_from_arr(bs_vals)
            print_summary_statistics(stats)

    def process_args(self, args) -> Dict[str, str]:
        return dict(
            tree_file_path=args.tree,
            verbose=args.verbose,
            cpu=args.cpu
        )

    def get_bipartition_support_vals(
        self,
        tree: Newick.Tree,
    ) -> Tuple[List[float], List[List[str]]]:
        nonterminals = tree.get_nonterminals()

        cpu = self.set_cpu()

        with Pool(cpu) as pool:
            results = pool.map(
                self.extract_support_and_terminals, nonterminals
            )

        bs_vals = [result[0] for result in results if result[0] is not None]
        term_names = [result[1] for result in results if result[0] is not None]

        return bs_vals, term_names

    def extract_support_and_terminals(
        self,
        nonterminal
    ) -> Tuple[
        Union[
                float, List[str], None
        ]
    ]:
        print(nonterminal, type(nonterminal))
        if nonterminal.confidence is None:
            return None, None
        term_names = [term.name for term in nonterminal.get_terminals()]
        return nonterminal.confidence, term_names
