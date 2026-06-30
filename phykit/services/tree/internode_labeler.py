from __future__ import annotations

from .base import Tree


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


class InternodeLabeler(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(
            tree_file_path=parsed["tree_file_path"],
            output_file_path=parsed["output_file_path"],
        )
        self.json_output = parsed["json_output"]

    def run(self):
        tree = self.read_tree_file()
        labeled_count = self.add_labels_to_tree(tree)
        self.write_tree_file(tree, self.output_file_path)

        if self.json_output:
            print_json(
                dict(
                    input_tree=self.tree_file_path,
                    labeled_internodes=labeled_count,
                    output_file=self.output_file_path,
                )
            )

    def process_args(self, args) -> dict[str, str]:
        output_file_path = args.output or f"{args.tree}.internode_labels.tre"

        return dict(
            tree_file_path=args.tree,
            output_file_path=output_file_path,
            json_output=getattr(args, "json", False),
        )

    def add_labels_to_tree(
        self,
        tree: Newick.Tree
    ) -> int:
        count = self._add_labels_to_standard_tree(tree)
        if count is not None:
            return count

        labeled_count = 0
        for label, node in enumerate(tree.get_nonterminals(), start=1):
            node.confidence = label
            labeled_count += 1
        return labeled_count

    @staticmethod
    def _add_labels_to_standard_tree(tree: Newick.Tree):
        try:
            root = tree.root
            root.clades
        except AttributeError:
            return None

        labeled_count = 0
        stack = [root]
        pop = stack.pop
        append = stack.append
        try:
            while stack:
                node = pop()
                children = node.clades
                if children:
                    labeled_count += 1
                    node.confidence = labeled_count
                    child_count = len(children)
                    if child_count == 2:
                        append(children[1])
                        append(children[0])
                    else:
                        for idx in range(child_count - 1, -1, -1):
                            append(children[idx])
        except AttributeError:
            return None
        return labeled_count
