"""
Transfer annotations between phylogenies.

Copies internal node annotations (e.g., wASTRAL q1/q2/q3, pp1, f1)
from an annotated source tree onto a target tree with optimized branch
lengths. Nodes are matched by their bipartition (set of descendant taxa).

Typical use case: transfer wASTRAL support annotations onto a RAxML-NG
branch-length-optimized topology for use with quartet_pie or other
visualization commands.
"""
from io import StringIO
from typing import Dict, List, Set, Tuple

from Bio import Phylo

from .base import Tree
from ...helpers.json_output import print_json
from ...errors import PhykitUserError


class TransferAnnotations(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["target_path"])
        self.source_path = parsed["source_path"]
        self.output_path = parsed["output_path"]
        self.json_output = parsed["json_output"]

    def process_args(self, args) -> Dict:
        return dict(
            target_path=args.target,
            source_path=args.source,
            output_path=getattr(args, "output", None),
            json_output=getattr(args, "json", False),
        )

    def run(self) -> None:
        # Read both trees
        target = self.read_tree_file()
        try:
            source = Phylo.read(self.source_path, "newick")
        except FileNotFoundError:
            raise PhykitUserError(
                [
                    f"{self.source_path} corresponds to no such file or directory.",
                    "Please check filename and pathing",
                ],
                code=2,
            )

        # Validate taxa match
        target_tips = {t.name for t in target.get_terminals()}
        source_tips = {t.name for t in source.get_terminals()}

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
                print(f"Annotations transferred: {n_transferred}")
                print(f"Unmatched nodes: {n_unmatched}")
                print(f"Output: {output_path}")
            except BrokenPipeError:
                return

    @staticmethod
    def _get_bipartition(clade, all_taxa: frozenset) -> frozenset:
        """Get canonical bipartition for a clade (smaller side)."""
        tips = frozenset(t.name for t in clade.get_terminals())
        complement = all_taxa - tips
        # Use the smaller side as the canonical key
        return tips if len(tips) <= len(complement) else complement

    def _extract_annotations(
        self, tree, all_taxa: frozenset
    ) -> Dict[frozenset, str]:
        """Extract bipartition -> annotation mapping from source tree."""
        annotations = {}
        for clade in tree.find_clades(order="preorder"):
            if clade.is_terminal() or clade == tree.root:
                continue
            annotation = clade.comment or clade.name or ""
            if not annotation:
                continue
            bp = self._get_bipartition(clade, all_taxa)
            if bp:
                annotations[bp] = annotation
        return annotations

    def _transfer(
        self, target, source_annotations: Dict[frozenset, str],
        all_taxa: frozenset,
    ) -> Tuple[int, int]:
        """Transfer annotations from source to target by bipartition matching."""
        transferred = 0
        unmatched = 0

        for clade in target.find_clades(order="preorder"):
            if clade.is_terminal() or clade == target.root:
                continue
            bp = self._get_bipartition(clade, all_taxa)
            if bp in source_annotations:
                clade.comment = source_annotations[bp]
                transferred += 1
            else:
                unmatched += 1

        return transferred, unmatched

    @staticmethod
    def _write_annotated_tree(tree, output_path: str) -> None:
        """Write tree preserving comment annotations in brackets."""
        Phylo.write(tree, output_path, "newick")
        # BioPython writes comments as [&comment] but wASTRAL uses
        # [comment]. Fix the format.
        with open(output_path) as f:
            content = f.read()
        # Remove the & prefix that BioPython adds
        content = content.replace("[&", "[")
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
