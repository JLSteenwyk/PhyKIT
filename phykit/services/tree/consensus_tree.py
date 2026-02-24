from io import StringIO
from pathlib import Path
from typing import Dict, List, Set

from Bio import Phylo
from Bio.Phylo import Consensus

from .base import Tree
from ...errors import PhykitUserError
from ...helpers.json_output import print_json


class ConsensusTree(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(trees=parsed["trees"])
        self.method = parsed["method"]
        self.missing_taxa = parsed["missing_taxa"]
        self.json_output = parsed["json_output"]

    def process_args(self, args) -> Dict[str, str]:
        return dict(
            trees=args.trees,
            method=args.method,
            missing_taxa=args.missing_taxa,
            json_output=getattr(args, "json", False),
        )

    def _parse_trees_from_source(self, trees_path: str):
        source = Path(trees_path)
        try:
            lines = source.read_text().splitlines()
        except FileNotFoundError:
            raise PhykitUserError(
                [
                    f"{trees_path} corresponds to no such file or directory.",
                    "Please check filename and pathing",
                ],
                code=2,
            )

        cleaned = [line.strip() for line in lines if line.strip() and not line.strip().startswith("#")]
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
        for line in cleaned:
            tree_path = Path(line)
            if not tree_path.is_absolute():
                tree_path = source.parent / tree_path
            if not tree_path.exists():
                raise PhykitUserError(
                    [
                        f"{tree_path} corresponds to no such file or directory.",
                        "Please check filename and pathing",
                    ],
                    code=2,
                )
            trees.append(Phylo.read(str(tree_path), "newick"))

        return trees

    @staticmethod
    def _tips(tree) -> Set[str]:
        return {tip.name for tip in tree.get_terminals()}

    @staticmethod
    def _prune_to_taxa(tree, taxa: Set[str]):
        remove = [tip.name for tip in tree.get_terminals() if tip.name not in taxa]
        for tip in remove:
            tree.prune(tip)
        return tree

    def _normalize_taxa(self, trees: List):
        tip_sets = [self._tips(tree) for tree in trees]
        shared_taxa = set.intersection(*tip_sets)

        identical = all(tip_set == tip_sets[0] for tip_set in tip_sets[1:])
        if identical:
            return trees, False, len(tip_sets[0])

        if self.missing_taxa == "error":
            raise PhykitUserError(
                [
                    "Input trees do not share an identical taxon set.",
                    "Use --missing-taxa shared to prune all trees to their shared taxa before inferring consensus.",
                ],
                code=2,
            )

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

    def _infer_consensus(self, trees: List):
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
