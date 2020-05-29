from phykit.services.tree.base import Tree


class Treeness(Tree):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        tree = self.read_tree_file()
        treeness = self.calculate_treeness(tree)
        if treeness:
            print(f"Treeness score: {treeness}")

    def process_args(self, args):
        return dict(tree_file_path=args.tree)

    def calculate_treeness(self, tree):
        inter_len = float(0.0)
        # determine internal branch lengths
        for interal in tree.get_nonterminals():
            # only include if a branch length value is present
            if interal.branch_length != None:
                inter_len += interal.branch_length
        # determine total branch length
        total_len = tree.total_branch_length()

        try:
            treeness = float(inter_len / total_len)
            return treeness
        except ZeroDivisionError:
            print("Invalid tree. Tree should contain branch lengths")
            return None

