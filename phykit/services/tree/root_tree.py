from .base import Tree


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


def read_single_column_file_to_list(*args, **kwargs):
    from ...helpers.files import (
        read_single_column_file_to_list as _read_single_column_file_to_list,
    )

    return _read_single_column_file_to_list(*args, **kwargs)


class _LazyBaseTreeTree:
    def root_with_outgroup(self, *args, **kwargs):
        from Bio.Phylo.BaseTree import Tree as _Tree

        return _Tree.root_with_outgroup(*args, **kwargs)


class _LazyBaseTree:
    Tree = _LazyBaseTreeTree()


class _LazyPhylo:
    BaseTree = _LazyBaseTree()


Phylo = _LazyPhylo()


class RootTree(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(
            tree_file_path=parsed["tree_file_path"],
            outgroup_taxa_file_path=parsed["outgroup_taxa_file_path"],
            output_file_path=parsed["output_file_path"],
        )
        self.json_output = parsed["json_output"]

    def run(self):
        tree = self.read_tree_file()

        outgroup = \
            read_single_column_file_to_list(self.outgroup_taxa_file_path)

        Phylo.BaseTree.Tree.root_with_outgroup(tree, outgroup)

        self.write_tree_file(tree, self.output_file_path)

        if self.json_output:
            print_json(
                dict(
                    input_tree=self.tree_file_path,
                    outgroup_taxa_file=self.outgroup_taxa_file_path,
                    outgroup_taxa=outgroup,
                    output_file=self.output_file_path,
                )
            )

    def process_args(self, args):
        tree_file_path = args.tree

        output_file_path = \
            args.output if args.output else f"{tree_file_path}.rooted"

        return dict(
            tree_file_path=tree_file_path,
            outgroup_taxa_file_path=args.root,
            output_file_path=output_file_path,
            json_output=getattr(args, "json", False),
        )
