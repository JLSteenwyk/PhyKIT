from .base import Tree


class TotalTreeLength(Tree):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        tree = self.read_tree_file()
        total_tree_length = self.calculate_total_tree_length(tree)
        print(round(total_tree_length, 4))


    def process_args(self, args):
        return dict(tree_file_path=args.tree)

    def calculate_total_tree_length(self, tree):
        total_len = float(0.0)
        total_len = tree.total_branch_length()

        if isinstance(total_len, (int, float)):
            return total_len
