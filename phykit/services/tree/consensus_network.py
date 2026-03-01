import math
from collections import Counter
from io import StringIO
from pathlib import Path
from typing import Dict, List, Set, Tuple

from Bio import Phylo
from Bio.Phylo import Consensus

from .base import Tree
from ...errors import PhykitUserError
from ...helpers.json_output import print_json


class ConsensusNetwork(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(trees=parsed["trees"])
        self.threshold = parsed["threshold"]
        self.missing_taxa = parsed["missing_taxa"]
        self.plot_output = parsed["plot_output"]
        self.json_output = parsed["json_output"]

    def process_args(self, args) -> Dict[str, str]:
        return dict(
            trees=args.trees,
            threshold=args.threshold,
            missing_taxa=args.missing_taxa,
            plot_output=getattr(args, "plot_output", None),
            json_output=getattr(args, "json", False),
        )

    # ------------------------------------------------------------------
    # Tree parsing (copied from ConsensusTree, self-contained per convention)
    # ------------------------------------------------------------------

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

        if all(line.startswith("(") for line in cleaned):
            trees = []
            try:
                for line in cleaned:
                    trees.append(Phylo.read(StringIO(line), "newick"))
            except Exception as exc:
                raise PhykitUserError([f"Failed to parse Newick trees: {exc}"], code=2)
            return trees

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
            return trees, False, tip_sets[0]

        if self.missing_taxa == "error":
            raise PhykitUserError(
                [
                    "Input trees do not share an identical taxon set.",
                    "Use --missing-taxa shared to prune all trees to their shared taxa.",
                ],
                code=2,
            )

        if len(shared_taxa) < 3:
            raise PhykitUserError(
                [
                    "Unable to compute network after pruning to shared taxa.",
                    "At least 3 shared taxa are required.",
                ],
                code=2,
            )

        pruned = [self._prune_to_taxa(tree, shared_taxa) for tree in trees]
        return pruned, True, shared_taxa

    # ------------------------------------------------------------------
    # Split extraction
    # ------------------------------------------------------------------

    @staticmethod
    def _canonical_split(taxa_side: frozenset, all_taxa: frozenset) -> frozenset:
        complement = all_taxa - taxa_side
        if len(taxa_side) < len(complement):
            return taxa_side
        if len(taxa_side) > len(complement):
            return complement
        # Equal size: return the lexicographically smaller side
        if sorted(taxa_side) < sorted(complement):
            return taxa_side
        return complement

    @staticmethod
    def _extract_splits_from_tree(tree, all_taxa: frozenset) -> Set[frozenset]:
        splits = set()
        for clade in tree.get_nonterminals():
            tips = frozenset(tip.name for tip in clade.get_terminals())
            # Skip trivial splits (single taxon or all-but-one)
            if len(tips) <= 1 or len(tips) >= len(all_taxa) - 1:
                continue
            # Skip root-level clade that contains all taxa
            if tips == all_taxa:
                continue
            canonical = ConsensusNetwork._canonical_split(tips, all_taxa)
            splits.add(canonical)
        return splits

    @staticmethod
    def _count_splits(trees: List, all_taxa: frozenset) -> Counter:
        counter = Counter()
        for tree in trees:
            tree_splits = ConsensusNetwork._extract_splits_from_tree(tree, all_taxa)
            for split in tree_splits:
                counter[split] += 1
        return counter

    @staticmethod
    def _filter_splits(split_counts: Counter, n_trees: int, threshold: float) -> List[Tuple[frozenset, int, float]]:
        results = []
        for split, count in split_counts.items():
            freq = count / n_trees
            if freq >= threshold:
                results.append((split, count, freq))
        results.sort(key=lambda x: (-x[2], sorted(x[0])))
        return results

    # ------------------------------------------------------------------
    # Circular ordering for visualization
    # ------------------------------------------------------------------

    @staticmethod
    def _compute_circular_ordering(trees: List, all_taxa: frozenset) -> List[str]:
        try:
            consensus = Consensus.majority_consensus(trees, cutoff=0.5)
            ordering = [tip.name for tip in consensus.get_terminals()]
            if set(ordering) == set(all_taxa):
                return ordering
        except Exception:
            pass
        return sorted(all_taxa)

    # ------------------------------------------------------------------
    # Visualization
    # ------------------------------------------------------------------

    def _draw_network(self, ordering: List[str], filtered_splits, all_taxa: frozenset, output_path: str):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np

        n = len(ordering)
        angles = {taxon: 2 * math.pi * i / n for i, taxon in enumerate(ordering)}

        fig, ax = plt.subplots(1, 1, figsize=(8, 8))
        ax.set_aspect("equal")

        # Draw taxa labels on circle
        radius = 1.0
        label_radius = 1.15
        for taxon in ordering:
            angle = angles[taxon]
            x = radius * math.cos(angle)
            y = radius * math.sin(angle)
            ax.plot(x, y, "o", color="black", markersize=6, zorder=3)

            lx = label_radius * math.cos(angle)
            ly = label_radius * math.sin(angle)
            # Rotate label to face outward
            deg = math.degrees(angle)
            ha = "left" if -90 < deg < 90 or deg > 270 else "right"
            rotation = deg if -90 < deg < 90 else deg - 180
            if deg > 270:
                rotation = deg - 360
            ax.text(lx, ly, taxon, ha=ha, va="center", fontsize=10, rotation=rotation)

        # Draw circle outline
        theta = np.linspace(0, 2 * np.pi, 200)
        ax.plot(radius * np.cos(theta), radius * np.sin(theta), "-", color="lightgray", linewidth=0.5, zorder=1)

        # Draw split chords
        if filtered_splits:
            max_freq = max(fs[2] for fs in filtered_splits)
        else:
            max_freq = 1.0

        position_index = {taxon: i for i, taxon in enumerate(ordering)}

        for split, count, freq in filtered_splits:
            complement = all_taxa - split
            # Find boundary points: positions where a split-member and
            # complement-member are adjacent in the circular ordering.
            boundary_midpoints = []
            for i in range(n):
                curr = ordering[i]
                nxt = ordering[(i + 1) % n]
                if (curr in split and nxt in complement) or (curr in complement and nxt in split):
                    mid_angle = (angles[curr] + angles[nxt]) / 2
                    # Handle wrap-around
                    if abs(angles[curr] - angles[nxt]) > math.pi:
                        mid_angle += math.pi
                    boundary_midpoints.append(mid_angle)

            # Draw chords between boundary midpoints
            alpha = 0.3 + 0.7 * (freq / max_freq)
            linewidth = 1.0 + 3.0 * (freq / max_freq)
            for i in range(len(boundary_midpoints)):
                for j in range(i + 1, len(boundary_midpoints)):
                    a1, a2 = boundary_midpoints[i], boundary_midpoints[j]
                    x1 = radius * math.cos(a1)
                    y1 = radius * math.sin(a1)
                    x2 = radius * math.cos(a2)
                    y2 = radius * math.sin(a2)
                    ax.plot([x1, x2], [y1, y2], "-", color="steelblue", alpha=alpha, linewidth=linewidth, zorder=2)

        ax.set_xlim(-1.5, 1.5)
        ax.set_ylim(-1.5, 1.5)
        ax.axis("off")
        ax.set_title("Consensus Splits Network", fontsize=14, pad=20)

        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close()

    # ------------------------------------------------------------------
    # Output
    # ------------------------------------------------------------------

    def _format_split(self, split: frozenset) -> str:
        return "{" + ", ".join(sorted(split)) + "}"

    def run(self):
        trees = self._parse_trees_from_source(self.trees)
        if len(trees) == 0:
            raise PhykitUserError(["No trees were parsed from input."], code=2)

        trees, pruned, all_taxa_set = self._normalize_taxa(trees)
        all_taxa = frozenset(all_taxa_set)
        n_trees = len(trees)

        split_counts = self._count_splits(trees, all_taxa)
        filtered = self._filter_splits(split_counts, n_trees, self.threshold)

        if self.json_output:
            splits_list = [
                {
                    "split": sorted(split),
                    "count": count,
                    "frequency": round(freq, 4),
                }
                for split, count, freq in filtered
            ]
            print_json(
                dict(
                    input_tree_count=n_trees,
                    taxa=sorted(all_taxa),
                    taxa_count=len(all_taxa),
                    threshold=self.threshold,
                    pruned_to_shared_taxa=pruned,
                    total_unique_splits=len(split_counts),
                    filtered_split_count=len(filtered),
                    splits=splits_list,
                )
            )
        else:
            print(f"Number of input trees: {n_trees}")
            print(f"Number of taxa: {len(all_taxa)}")
            print(f"Threshold: {self.threshold}")
            if pruned:
                print("Pruned to shared taxa: yes")
            print(f"Total unique splits: {len(split_counts)}")
            print(f"Splits above threshold: {len(filtered)}")
            print("---")
            for split, count, freq in filtered:
                print(f"{self._format_split(split)}\t{count}/{n_trees}\t{freq:.4f}")

        if self.plot_output:
            ordering = self._compute_circular_ordering(trees, all_taxa)
            self._draw_network(ordering, filtered, all_taxa, self.plot_output)
