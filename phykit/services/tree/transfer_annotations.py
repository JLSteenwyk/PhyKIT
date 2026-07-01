"""
Transfer annotations between phylogenies.

Copies internal node annotations (e.g., wASTRAL q1/q2/q3, pp1, f1)
from an annotated source tree onto a target tree with optimized branch
lengths. Nodes are matched by their bipartition (set of descendant taxa).

Typical use case: transfer wASTRAL support annotations onto a
branch-length-optimized topology from RAxML-NG, IQ-TREE, or any
other tool for use with quartet_pie or other visualization commands.
"""
from __future__ import annotations

from io import StringIO

from .base import Tree
from ...errors import PhykitUserError


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


class _LazyPhylo:
    def read(self, *args, **kwargs):
        from Bio import Phylo as _Phylo

        return _Phylo.read(*args, **kwargs)

    def write(self, *args, **kwargs):
        from Bio import Phylo as _Phylo

        return _Phylo.write(*args, **kwargs)


Phylo = _LazyPhylo()


class TransferAnnotations(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["target_path"])
        self.source_path = parsed["source_path"]
        self.output_path = parsed["output_path"]
        self.json_output = parsed["json_output"]

    def process_args(self, args) -> dict:
        return dict(
            target_path=args.target,
            source_path=args.source,
            output_path=getattr(args, "output", None),
            json_output=getattr(args, "json", False),
        )

    def run(self) -> None:
        # Read both trees
        target = self.read_tree_file()
        source = self._read_tree_with_error(
            self.source_path,
            "source_path",
            copy_tree=False,
        )

        # Validate taxa match
        target_tips = set(self.get_tip_names_from_tree(target))
        source_tips = set(self.get_tip_names_from_tree(source))

        if target_tips != source_tips:
            only_target = target_tips - source_tips
            only_source = source_tips - target_tips
            msgs = ["Taxa in source and target trees do not match."]
            if only_target:
                msgs.append(
                    f"Only in target: {', '.join(sorted(only_target)[:5])}"
                    + (f" ... ({len(only_target)} total)" if len(only_target) > 5 else "")
                )
            if only_source:
                msgs.append(
                    f"Only in source: {', '.join(sorted(only_source)[:5])}"
                    + (f" ... ({len(only_source)} total)" if len(only_source) > 5 else "")
                )
            raise PhykitUserError(msgs, code=2)

        # Build bipartition -> annotation map from source tree
        all_taxa = frozenset(source_tips)
        source_annotations = self._extract_annotations(source, all_taxa)

        # Match and transfer annotations to target tree
        n_transferred, n_unmatched = self._transfer(
            target, source_annotations, all_taxa
        )

        # Output
        output_path = self.output_path or f"{self.tree_file_path}.annotated"
        self._write_annotated_tree(target, output_path)

        if self.json_output:
            self._print_json(
                n_transferred, n_unmatched, output_path,
                len(source_annotations), target_tips,
            )
        else:
            try:
                self._print_text_output(n_transferred, n_unmatched, output_path)
            except BrokenPipeError:
                return

    def _print_text_output(
        self,
        n_transferred: int,
        n_unmatched: int,
        output_path: str,
    ) -> None:
        print(
            "\n".join(
                [
                    f"Annotations transferred: {n_transferred}",
                    f"Unmatched nodes: {n_unmatched}",
                    f"Output: {output_path}",
                ]
            )
        )

    @staticmethod
    def _collect_clade_taxa(tree) -> dict[int, frozenset]:
        clade_taxa: dict[int, frozenset] = {}
        postorder = TransferAnnotations._postorder_clades_direct(tree)
        if postorder is None:
            postorder = tree.find_clades(order="postorder")

        id_ = id
        get = clade_taxa.__getitem__
        for clade in postorder:
            children = clade.clades
            child_count = len(children)
            if child_count == 0:
                clade_taxa[id_(clade)] = frozenset((clade.name,))
            elif child_count == 1:
                clade_taxa[id_(clade)] = get(id_(children[0]))
            elif child_count == 2:
                clade_taxa[id_(clade)] = (
                    get(id_(children[0]))
                    | get(id_(children[1]))
                )
            else:
                clade_taxa[id_(clade)] = frozenset().union(
                    *(get(id_(child)) for child in children)
                )
        return clade_taxa

    @staticmethod
    def _canonical_bipartition(tips: frozenset, all_taxa: frozenset) -> frozenset:
        # Use the smaller side as the canonical key
        if len(tips) * 2 <= len(all_taxa):
            return tips
        return all_taxa - tips

    @staticmethod
    def _postorder_clades_direct(tree):
        try:
            root = tree.root
            root.clades
        except AttributeError:
            return None

        clades = []
        stack = [root]
        try:
            while stack:
                clade = stack.pop()
                clades.append(clade)
                stack.extend(clade.clades)
        except AttributeError:
            return None
        clades.reverse()
        return clades

    @staticmethod
    def _combine_child_taxa(children, clade_taxa: dict[int, frozenset]) -> frozenset:
        if len(children) == 2:
            return clade_taxa[id(children[0])] | clade_taxa[id(children[1])]

        taxa = set()
        for child in children:
            taxa.update(clade_taxa.get(id(child), ()))
        return frozenset(taxa)

    @staticmethod
    def _get_bipartition(clade, all_taxa: frozenset, clade_taxa=None) -> frozenset:
        """Get canonical bipartition for a clade (smaller side)."""
        if clade_taxa is None:
            tips = frozenset(t.name for t in clade.get_terminals())
        else:
            tips = clade_taxa.get(id(clade), frozenset())
        return TransferAnnotations._canonical_bipartition(tips, all_taxa)

    def _extract_annotations(
        self, tree, all_taxa: frozenset
    ) -> dict[frozenset, str]:
        """Extract bipartition -> annotation mapping from source tree."""
        annotations = {}
        clade_taxa: dict[int, frozenset] = {}
        root = tree.root
        postorder = self._postorder_clades_direct(tree)
        if postorder is None:
            postorder = tree.find_clades(order="postorder")

        for clade in postorder:
            children = clade.clades
            if not children:
                clade_taxa[id(clade)] = frozenset({clade.name})
                continue

            tips = self._combine_child_taxa(children, clade_taxa)
            clade_taxa[id(clade)] = tips

            if clade is root:
                continue
            annotation = clade.comment or clade.name or ""
            if not annotation:
                continue
            # wASTRAL --support 3 wraps annotations as '[key=val;...]'
            # which BioPython parses as a quoted node name including the
            # brackets.  Strip them so BioPython's comment writer can
            # re-add proper brackets on output.
            if annotation.startswith("[") and annotation.endswith("]"):
                annotation = annotation[1:-1]
            bp = self._canonical_bipartition(tips, all_taxa)
            if bp:
                annotations[bp] = annotation
        return annotations

    def _transfer(
        self, target, source_annotations: dict[frozenset, str],
        all_taxa: frozenset,
    ) -> tuple[int, int]:
        """Transfer annotations from source to target by bipartition matching."""
        transferred = 0
        unmatched = 0
        clade_taxa: dict[int, frozenset] = {}
        root = target.root
        postorder = self._postorder_clades_direct(target)
        if postorder is None:
            postorder = target.find_clades(order="postorder")

        for clade in postorder:
            children = clade.clades
            if not children:
                clade_taxa[id(clade)] = frozenset({clade.name})
                continue

            tips = self._combine_child_taxa(children, clade_taxa)
            clade_taxa[id(clade)] = tips

            if clade is root:
                continue
            bp = self._canonical_bipartition(tips, all_taxa)
            if bp in source_annotations:
                clade.comment = source_annotations[bp]
                transferred += 1
            else:
                unmatched += 1

        return transferred, unmatched

    @staticmethod
    def _write_annotated_tree(tree, output_path: str) -> None:
        """Write tree preserving comment annotations in brackets."""
        buffer = StringIO()
        Phylo.write(tree, buffer, "newick")
        content = buffer.getvalue()

        # BioPython writes comments as [&comment] but wASTRAL uses
        # [comment]. Fix the format before writing the output once.
        # Remove the & prefix that BioPython adds
        content = content.replace("[&", "[")
        # BioPython may also backslash-escape brackets inside comments
        content = content.replace("\\[", "[").replace("\\]", "]")
        with open(output_path, "w") as f:
            f.write(content)

    def _print_json(
        self, n_transferred, n_unmatched, output_path,
        n_source_annotations, target_tips,
    ):
        payload = {
            "annotations_transferred": n_transferred,
            "unmatched_nodes": n_unmatched,
            "source_annotations_available": n_source_annotations,
            "n_taxa": len(target_tips),
            "output": output_path,
        }
        print_json(payload, sort_keys=False)
