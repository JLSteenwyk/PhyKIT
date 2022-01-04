import sys

from Bio.Phylo.BaseTree import TreeMixin

from .base import Tree


class TipToTipNodeDistance(Tree):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        tree_zero = self.read_tree_file()

        # determine is leaves exists
        self.check_leaves(tree_zero, self.tip_1, self.tip_2)

        # print distance
        print(len(TreeMixin.trace(tree_zero, self.tip_1, self.tip_2)))

    def check_leaves(self, tree_zero, tip_1, tip_2):
        leaf1 = TreeMixin.find_any(tree_zero, tip_1)
        if (bool(leaf1)) == False:
            print(tip_1, "not on tree\nExiting...")
            sys.exit()
        leaf2 = TreeMixin.find_any(tree_zero, tip_2)
        if (bool(leaf2)) == False:
            print(tip_2, "not on tree\nExiting...")
            sys.exit()

    def process_args(self, args):
        return dict(
            tree_file_path=args.tree_zero,
            tip_1=args.tip_1,
            tip_2=args.tip_2
        )
