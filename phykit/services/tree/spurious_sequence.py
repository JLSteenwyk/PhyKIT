from __future__ import annotations

from .base import Tree


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


class _LazyNumpy:
    _module = None

    def __getattr__(self, name):
        module = self._module
        if module is None:
            import numpy as _np

            module = _np
            self._module = module

        value = getattr(module, name)
        setattr(self, name, value)
        return value


np = _LazyNumpy()

_MEDIAN_NUMPY_THRESHOLD = 1024


class SpuriousSequence(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"], factor=parsed["factor"])
        self.json_output = parsed["json_output"]

    def run(self) -> None:
        tree = self.read_tree_file_unmodified()
        name_and_branch_len, threshold, median = \
            self.identify_spurious_sequence(
                tree, self.factor
            )

        rounded_threshold = round(threshold, 4)
        rounded_median = round(median, 4)
        has_spurious_sequence = self._has_spurious_sequence(
            name_and_branch_len,
            threshold,
        )

        if self.json_output:
            if not has_spurious_sequence:
                print_json(dict(rows=[], spurious_sequences=[]))
                return

            spurious_rows = [
                {
                    "taxon": name,
                    "branch_length": round(length, 4),
                    "threshold": rounded_threshold,
                    "median": rounded_median,
                }
                for name, length in name_and_branch_len.items()
                if length >= threshold
            ]
            print_json(dict(rows=spurious_rows, spurious_sequences=spurious_rows))
            return

        if not has_spurious_sequence:
            print("None")
            return

        try:
            print(
                "\n".join(
                    f"{name}\t{round(length, 4)}\t{rounded_threshold}\t{rounded_median}"
                    for name, length in name_and_branch_len.items()
                    if length >= threshold
                )
            )
        except BrokenPipeError:
            pass

    def process_args(self, args) -> dict[str, str]:
        return dict(
            tree_file_path=args.tree,
            factor=args.factor or 20,
            json_output=getattr(args, "json", False),
        )

    def identify_spurious_sequence(
        self,
        tree: Newick.Tree,
        factor: float,
    ) -> tuple[
        dict[str, float],
        float,
        float,
    ]:
        branch_lengths, name_and_branch_len = \
            self.get_branch_lengths_and_their_names(tree)

        median = self._median_branch_length(branch_lengths)

        threshold = median * factor

        return name_and_branch_len, threshold, median

    @staticmethod
    def _has_spurious_sequence(name_and_branch_len, threshold) -> bool:
        if not name_and_branch_len:
            return False

        longest_branch_length = max(name_and_branch_len.values())
        if longest_branch_length >= threshold:
            return True

        if longest_branch_length != longest_branch_length:
            for length in name_and_branch_len.values():
                if length >= threshold:
                    return True

        return False

    @staticmethod
    def _median_branch_length(branch_lengths: list[float]) -> float:
        count = len(branch_lengths)
        if count == 0:
            raise ValueError("no median for empty data")

        middle = count // 2
        if count >= _MEDIAN_NUMPY_THRESHOLD:
            values = np.asarray(branch_lengths, dtype=float)
            if count % 2:
                values.partition(middle)
                return float(values[middle])
            values.partition((middle - 1, middle))
            return float((values[middle - 1] + values[middle]) / 2)

        values = sorted(branch_lengths)
        if count % 2:
            return values[middle]
        return (values[middle - 1] + values[middle]) / 2

    def get_branch_lengths_and_their_names(
        self,
        tree: Newick.Tree,
    ) -> tuple[
        list[float],
        dict[str, float],
    ]:
        direct_result = self._get_branch_lengths_and_names_standard_tree(tree)
        if direct_result is not None:
            return direct_result

        branch_lengths = []
        name_and_branch_len = {}

        # collect terminal branch lengths only for spurious sequence detection
        # (internal branches are not considered for spurious sequence detection)
        terminals = self._iter_terminal_clades(tree)
        for terminal in terminals:
            if terminal.branch_length is not None:
                branch_lengths.append(terminal.branch_length)
                name_and_branch_len[terminal.name] = terminal.branch_length

        return branch_lengths, name_and_branch_len

    @staticmethod
    def _get_branch_lengths_and_names_standard_tree(tree):
        try:
            root = tree.root
            root.clades
        except AttributeError:
            return None

        branch_lengths = []
        name_and_branch_len = {}
        stack = [root]
        pop = stack.pop
        append = stack.append
        append_length = branch_lengths.append
        try:
            while stack:
                clade = pop()
                children = clade.clades
                child_count = len(children)
                if child_count == 2:
                    append(children[1])
                    append(children[0])
                elif child_count:
                    for idx in range(child_count - 1, -1, -1):
                        append(children[idx])
                else:
                    branch_length = clade.branch_length
                    if branch_length is not None:
                        append_length(branch_length)
                        name_and_branch_len[clade.name] = branch_length
        except AttributeError:
            return None

        return branch_lengths, name_and_branch_len

    @staticmethod
    def _iter_terminal_clades(tree):
        try:
            root = tree.root
            root.clades
        except AttributeError:
            return tree.get_terminals()

        terminals = []
        stack = [root]
        pop = stack.pop
        append = stack.append
        append_terminal = terminals.append
        while stack:
            clade = pop()
            children = clade.clades
            child_count = len(children)
            if child_count == 2:
                append(children[1])
                append(children[0])
            elif child_count:
                for idx in range(child_count - 1, -1, -1):
                    append(children[idx])
            else:
                append_terminal(clade)
        return terminals
