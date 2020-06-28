from .base import Tree


class InternodeLabeler(Tree):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        tree = self.read_tree_file()
        tree_with_labels = self.add_labels_to_tree(tree)
        self.write_tree_file(tree_with_labels, self.output_file_path)

    def process_args(self, args):
        tree_file_path = args.tree

        if args.output is None:
            output_file_path = f"{tree_file_path}.internode_labels.tre"
        else:
            output_file_path = f"{args.output}"

        return dict(
            tree_file_path=tree_file_path,
            output_file_path=output_file_path,
        )

    def add_labels_to_tree(self, tree):
        label = 1
        for node in tree.get_nonterminals():
            node.confidence = label
            label += 1
        return tree
