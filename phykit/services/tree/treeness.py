from .base import Tree


class Treeness(Tree):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        tree = self.read_tree_file()
        treeness = self.calculate_treeness(tree)
        print(round(treeness, 4))


    def process_args(self, args):
        return dict(tree_file_path=args.tree)
