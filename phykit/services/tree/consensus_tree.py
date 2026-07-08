from __future__ import annotations

from io import StringIO
import os
from pathlib import Path

from .base import Tree
from ...errors import PhykitUserError


_path_exists = os.path.exists
_path_isabs = os.path.isabs


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


class _LazyPhylo:
    def read(self, *args, **kwargs):
        from Bio import Phylo as _Phylo

        self.read = _Phylo.read
        return self.read(*args, **kwargs)


class _LazyConsensus:
    def strict_consensus(self, *args, **kwargs):
        from Bio.Phylo import Consensus as _Consensus

        self.strict_consensus = _Consensus.strict_consensus
        return self.strict_consensus(*args, **kwargs)

    def majority_consensus(self, *args, **kwargs):
        from Bio.Phylo import Consensus as _Consensus

        self.majority_consensus = _Consensus.majority_consensus
        return self.majority_consensus(*args, **kwargs)


Phylo = _LazyPhylo()
Consensus = _LazyConsensus()


def _all_tip_sets_identical(tip_sets) -> bool:
    iterator = iter(tip_sets)
    first_tip_set = next(iterator, None)
    if first_tip_set is None:
        return True
    for tip_set in iterator:
        if tip_set != first_tip_set:
            return False
    return True


def _shared_tip_set(tip_sets) -> set[str]:
    iterator = iter(tip_sets)
    shared = set(next(iterator, ()))
    for tip_set in iterator:
        shared.intersection_update(tip_set)
        if not shared:
            break
    return shared


class ConsensusTree(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(trees=parsed["trees"])
        self.method = parsed["method"]
        self.missing_taxa = parsed["missing_taxa"]
        self.json_output = parsed["json_output"]

    def process_args(self, args) -> dict[str, str]:
        return dict(
            trees=args.trees,
            method=args.method,
            missing_taxa=args.missing_taxa,
            json_output=getattr(args, "json", False),
        )

    def _parse_trees_from_source(self, trees_path: str):
        source = Path(trees_path)
        try:
            with source.open() as handle:
                cleaned = [
                    stripped
                    for line in handle
                    if (stripped := line.strip())
                    and stripped[0] != "#"
                ]
        except FileNotFoundError:
            raise PhykitUserError(
                [
                    f"{trees_path} corresponds to no such file or directory.",
                    "Please check filename and pathing",
                ],
                code=2,
            )

        if not cleaned:
            raise PhykitUserError(["No trees found in input."], code=2)

        # Multi-newick file: each non-empty line is a Newick tree.
        if all(line.startswith("(") for line in cleaned):
            trees = []
            try:
                for line in cleaned:
                    trees.append(Phylo.read(StringIO(line), "newick"))
            except Exception as exc:  # pragma: no cover - exercised via CLI validation paths
                raise PhykitUserError([f"Failed to parse Newick trees: {exc}"], code=2)
            return trees

        # Otherwise, treat as a file containing one tree path per line.
        trees = []
        parent_str = str(source.parent)
        parent_prefix = "" if parent_str == "." else parent_str + os.sep
        for line in cleaned:
            tree_path = line if _path_isabs(line) else parent_prefix + line
            if not _path_exists(tree_path):
                raise PhykitUserError(
                    [
                        f"{tree_path} corresponds to no such file or directory.",
                        "Please check filename and pathing",
                    ],
                    code=2,
                )
            trees.append(Phylo.read(tree_path, "newick"))

        return trees

    @staticmethod
    def _tips(tree) -> set[str]:
        tip_names = Tree.calculate_terminal_names_fast(tree)
        if tip_names is not None:
            return set(tip_names)
        return {tip.name for tip in tree.get_terminals()}

    @staticmethod
    def _prune_to_taxa(tree, taxa: set[str]):
        terminals = Tree.calculate_terminal_clades_fast(tree)
        if terminals is None:
            terminals = tree.get_terminals()
        remove = [tip for tip in terminals if tip.name not in taxa]
        if len(remove) < len(terminals):
            target_ids = {id(tip) for tip in remove}
            if Tree._prune_terminal_objects_batch_standard_tree(tree, target_ids):
                return tree
        for tip in remove:
            tree.prune(tip)
        return tree

    def _normalize_taxa(self, trees: list):
        tip_sets = [self._tips(tree) for tree in trees]
        first_tip_set = tip_sets[0]
        if _all_tip_sets_identical(tip_sets):
            return trees, False, len(first_tip_set)

        if self.missing_taxa == "error":
            raise PhykitUserError(
                [
                    "Input trees do not share an identical taxon set.",
                    "Use --missing-taxa shared to prune all trees to their shared taxa before inferring consensus.",
                ],
                code=2,
            )

        shared_taxa = _shared_tip_set(tip_sets)
        if len(shared_taxa) < 3:
            raise PhykitUserError(
                [
                    "Unable to infer consensus after pruning to shared taxa.",
                    "At least 3 shared taxa are required.",
                ],
                code=2,
            )

        pruned = [self._prune_to_taxa(tree, shared_taxa) for tree in trees]
        return pruned, True, len(shared_taxa)

    def _infer_consensus(self, trees: list):
        if self.method == "strict":
            return Consensus.strict_consensus(trees)
        return Consensus.majority_consensus(trees, cutoff=0.5)

    def run(self):
        trees = self._parse_trees_from_source(self.trees)
        if len(trees) == 0:
            raise PhykitUserError(["No trees were parsed from input."], code=2)

        trees, pruned, taxa_count = self._normalize_taxa(trees)
        consensus = self._infer_consensus(trees)
        newick = consensus.format("newick").strip()

        if self.json_output:
            print_json(
                dict(
                    method=self.method,
                    missing_taxa=self.missing_taxa,
                    pruned_to_shared_taxa=pruned,
                    input_tree_count=len(trees),
                    taxa_in_consensus=taxa_count,
                    tree_newick=newick,
                )
            )
            return

        print(newick)
