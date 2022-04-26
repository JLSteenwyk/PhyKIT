from .base import Tree

from ...helpers.files import read_single_column_file_to_list

class PruneTree(Tree):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        tree = self.read_tree_file()
        
        taxa = read_single_column_file_to_list(self.list_of_taxa)

        if self.keep:
            tips_in_tree = []
            for term in tree.get_terminals():
                tips_in_tree.append(term.name)
            taxa = [x for x in tips_in_tree if x not in taxa]
            
        tree = self.prune_tree_using_taxa_list(tree, taxa)

        self.write_tree_file(tree, self.output_file_path)
    
    def process_args(self, args):
        tree_file_path = args.tree
        if args.output is None:
            output_file_path = f"{tree_file_path}.pruned"
        else:
            output_file_path = f"{args.output}"

        if args.keep is None:
            keep=True
        else:
            keep=args.keep

        return dict(
            tree_file_path=tree_file_path,
            list_of_taxa=args.list_of_taxa,
            output_file_path=output_file_path,
            keep=keep
        )