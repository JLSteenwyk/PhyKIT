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
    # Visualization: Planar Splits Graph
    # ------------------------------------------------------------------

    @staticmethod
    def _is_circular_split(split, ordering):
        """Check if a split has exactly 2 boundary gaps in the circular ordering."""
        n = len(ordering)
        gaps = 0
        for i in range(n):
            curr = ordering[i]
            nxt = ordering[(i + 1) % n]
            if (curr in split) != (nxt in split):
                gaps += 1
        return gaps == 2

    @staticmethod
    def _compute_split_directions(ordering, circular_splits):
        """Compute 2D direction vectors for each circular split."""
        n = len(ordering)
        angles = {taxon: 2 * math.pi * i / n for i, taxon in enumerate(ordering)}
        directions = {}
        for split, count, freq in circular_splits:
            gap_positions = []
            for i in range(n):
                curr = ordering[i]
                nxt = ordering[(i + 1) % n]
                if (curr in split) != (nxt in split):
                    gap_positions.append(i)
            if len(gap_positions) != 2:
                continue
            # Gap midpoint angles on the unit circle
            g1 = math.pi * (2 * gap_positions[0] + 1) / n
            g2 = math.pi * (2 * gap_positions[1] + 1) / n
            # Chord between gap midpoints
            cx = math.cos(g2) - math.cos(g1)
            cy = math.sin(g2) - math.sin(g1)
            # Perpendicular direction (rotated 90 degrees)
            dx = -cy
            dy = cx
            length = math.sqrt(dx * dx + dy * dy)
            if length > 1e-10:
                dx /= length
                dy /= length
            else:
                dx, dy = 1.0, 0.0
            # Orient toward the canonical (positive / smaller) side
            cx_split = sum(math.cos(angles[t]) for t in split) / len(split)
            cy_split = sum(math.sin(angles[t]) for t in split) / len(split)
            if dx * cx_split + dy * cy_split < 0:
                dx = -dx
                dy = -dy
            directions[split] = (dx, dy)
        return directions

    @staticmethod
    def _build_splits_graph(circular_splits, all_taxa):
        """Build the Buneman splits graph: valid sign vectors and edges.

        Returns (taxon_signs, valid_nodes, edges, splits_list, weights).
        """
        splits_list = [s[0] for s in circular_splits]
        weights = [s[2] for s in circular_splits]
        n_splits = len(splits_list)
        if n_splits == 0:
            return {}, set(), [], [], []
        # Each taxon's sign vector: +1 if in canonical side, -1 otherwise
        taxon_signs = {}
        for taxon in all_taxa:
            signs = tuple(1 if taxon in sp else -1 for sp in splits_list)
            taxon_signs[taxon] = signs
        # Pairwise forbidden sign combinations
        forbidden = {}
        for i in range(n_splits):
            for j in range(i + 1, n_splits):
                si_pos = splits_list[i]
                si_neg = all_taxa - si_pos
                sj_pos = splits_list[j]
                sj_neg = all_taxa - sj_pos
                fb = set()
                if not (si_pos & sj_pos):
                    fb.add((1, 1))
                if not (si_pos & sj_neg):
                    fb.add((1, -1))
                if not (si_neg & sj_pos):
                    fb.add((-1, 1))
                if not (si_neg & sj_neg):
                    fb.add((-1, -1))
                if fb:
                    forbidden[(i, j)] = fb
        # Enumerate valid sign vectors
        valid_nodes = set()
        for combo in range(2 ** n_splits):
            signs = tuple(1 if (combo >> j) & 1 else -1 for j in range(n_splits))
            valid = True
            for (i, j), fb_pairs in forbidden.items():
                if (signs[i], signs[j]) in fb_pairs:
                    valid = False
                    break
            if valid:
                valid_nodes.add(signs)
        # Ensure taxon sign vectors are included
        for signs in taxon_signs.values():
            valid_nodes.add(signs)
        # Edges: connect nodes differing in exactly one split
        edges = []
        valid_list = list(valid_nodes)
        for i in range(len(valid_list)):
            for j in range(i + 1, len(valid_list)):
                s1 = valid_list[i]
                s2 = valid_list[j]
                diff = [k for k in range(n_splits) if s1[k] != s2[k]]
                if len(diff) == 1:
                    edges.append((s1, s2, diff[0]))
        return taxon_signs, valid_nodes, edges, splits_list, weights

    def _draw_network(self, ordering: List[str], filtered_splits, all_taxa: frozenset, output_path: str):
        """Draw a planar splits graph using matplotlib."""
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        n = len(ordering)
        angles = {taxon: 2 * math.pi * i / n for i, taxon in enumerate(ordering)}

        # Filter to circular splits only
        circular_splits = [
            (split, count, freq)
            for split, count, freq in filtered_splits
            if self._is_circular_split(split, ordering)
        ]

        # Compute direction vectors
        directions = self._compute_split_directions(ordering, circular_splits)

        # Keep only splits that have computed directions
        circular_splits = [
            (split, count, freq)
            for split, count, freq in circular_splits
            if split in directions
        ]

        # Build splits graph
        taxon_signs, valid_nodes, edges, splits_list, weights = self._build_splits_graph(
            circular_splits, all_taxa
        )

        fig, ax = plt.subplots(1, 1, figsize=(8, 8))
        ax.set_aspect("equal")

        if not splits_list:
            # No splits: place taxa evenly on a circle
            for taxon in ordering:
                angle = angles[taxon]
                x = math.cos(angle)
                y = math.sin(angle)
                deg = math.degrees(angle)
                ha = "left" if -90 < deg < 90 or deg > 270 else "right"
                ax.text(x, y, taxon, ha=ha, va="center", fontsize=10)
        else:
            # Compute node positions: pos = sum(sign_i * weight_i * dir_i) / 2
            node_positions = {}
            for node in valid_nodes:
                x, y = 0.0, 0.0
                for i, sp in enumerate(splits_list):
                    dx, dy = directions[sp]
                    x += node[i] * weights[i] * dx / 2
                    y += node[i] * weights[i] * dy / 2
                node_positions[node] = (x, y)

            # Pendant edge length proportional to graph extent
            all_x = [p[0] for p in node_positions.values()]
            all_y = [p[1] for p in node_positions.values()]
            extent = max(
                max(all_x) - min(all_x) if all_x else 0,
                max(all_y) - min(all_y) if all_y else 0,
            )
            pendant_len = max(0.15, extent * 0.3)

            # Draw internal edges
            for s1, s2, split_idx in edges:
                x1, y1 = node_positions[s1]
                x2, y2 = node_positions[s2]
                ax.plot([x1, x2], [y1, y2], "-", color="black", linewidth=1.5, zorder=2)

            # Draw pendant edges and taxon labels
            for taxon in ordering:
                angle = angles[taxon]
                if taxon in taxon_signs and taxon_signs[taxon] in node_positions:
                    nx, ny = node_positions[taxon_signs[taxon]]
                else:
                    nx, ny = 0.0, 0.0
                tx = nx + pendant_len * math.cos(angle)
                ty = ny + pendant_len * math.sin(angle)
                ax.plot([nx, tx], [ny, ty], "-", color="black", linewidth=1.5, zorder=2)
                lx = tx + 0.03 * math.cos(angle)
                ly = ty + 0.03 * math.sin(angle)
                deg = math.degrees(angle)
                ha = "left" if -90 < deg < 90 or deg > 270 else "right"
                ax.text(lx, ly, taxon, ha=ha, va="center", fontsize=10, zorder=4)

            # Frequency labels on internal edges (one per split)
            labeled_splits = set()
            for s1, s2, split_idx in edges:
                if split_idx not in labeled_splits:
                    x1, y1 = node_positions[s1]
                    x2, y2 = node_positions[s2]
                    mx = (x1 + x2) / 2
                    my = (y1 + y2) / 2
                    freq = weights[split_idx]
                    ax.text(
                        mx, my, f"{freq:.2f}",
                        fontsize=7, ha="center", va="center",
                        color="gray", zorder=3,
                    )
                    labeled_splits.add(split_idx)

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
