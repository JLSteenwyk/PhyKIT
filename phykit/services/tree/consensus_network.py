from __future__ import annotations

import math
from collections import Counter
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


def _position_extent(positions):
    iterator = iter(positions)
    try:
        min_x, min_y = next(iterator)
    except StopIteration:
        return 0
    max_x = min_x
    max_y = min_y
    for x, y in iterator:
        if x < min_x:
            min_x = x
        elif x > max_x:
            max_x = x
        if y < min_y:
            min_y = y
        elif y > max_y:
            max_y = y
    x_extent = max_x - min_x
    y_extent = max_y - min_y
    return x_extent if x_extent >= y_extent else y_extent


def _all_tip_sets_identical(tip_sets) -> bool:
    iterator = iter(tip_sets)
    first_tip_set = next(iterator, None)
    if first_tip_set is None:
        return True
    for tip_set in iterator:
        if tip_set != first_tip_set:
            return False
    return True


class _LazyPhylo:
    def read(self, *args, **kwargs):
        from Bio import Phylo as _Phylo

        return _Phylo.read(*args, **kwargs)


class _LazyConsensus:
    def majority_consensus(self, *args, **kwargs):
        from Bio.Phylo import Consensus as _Consensus

        return _Consensus.majority_consensus(*args, **kwargs)


Phylo = _LazyPhylo()
Consensus = _LazyConsensus()


class ConsensusNetwork(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(trees=parsed["trees"])
        self.threshold = parsed["threshold"]
        self.missing_taxa = parsed["missing_taxa"]
        self.max_splits = parsed["max_splits"]
        self.histogram = parsed["histogram"]
        self.plot_output = parsed["plot_output"]
        self.json_output = parsed["json_output"]
        self.plot_config = parsed["plot_config"]

    def process_args(self, args) -> dict[str, str]:
        from ...helpers.plot_config import PlotConfig

        return dict(
            trees=args.trees,
            threshold=args.threshold,
            missing_taxa=args.missing_taxa,
            max_splits=getattr(args, "max_splits", 30),
            histogram=getattr(args, "histogram", None),
            plot_output=getattr(args, "plot_output", None),
            json_output=getattr(args, "json", False),
            plot_config=PlotConfig.from_args(args),
        )

    # ------------------------------------------------------------------
    # Tree parsing (copied from ConsensusTree, self-contained per convention)
    # ------------------------------------------------------------------

    def _parse_trees_from_source(self, trees_path: str):
        source = Path(trees_path)
        try:
            with source.open() as handle:
                cleaned = [
                    stripped
                    for line in handle
                    if (stripped := line.strip())
                    and not stripped.startswith("#")
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

        if all(line.startswith("(") for line in cleaned):
            trees = []
            try:
                for line in cleaned:
                    trees.append(Phylo.read(StringIO(line), "newick"))
            except Exception as exc:
                raise PhykitUserError([f"Failed to parse Newick trees: {exc}"], code=2)
            return trees

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
            return trees, False, first_tip_set

        if self.missing_taxa == "error":
            raise PhykitUserError(
                [
                    "Input trees do not share an identical taxon set.",
                    "Use --missing-taxa allow or --missing-taxa shared.",
                ],
                code=2,
            )

        if self.missing_taxa == "allow":
            # Use the union of all taxa; each tree contributes splits
            # using its own taxon set. Split frequencies are normalized
            # by how many trees could contain each split.
            union_taxa = set.union(*tip_sets)
            if len(union_taxa) < 3:
                raise PhykitUserError(
                    ["Fewer than 3 taxa found across all trees."], code=2
                )
            return trees, False, union_taxa

        # shared mode
        shared_taxa = set.intersection(*tip_sets)
        if len(shared_taxa) < 3:
            raise PhykitUserError(
                [
                    "Unable to compute network after pruning to shared taxa.",
                    "At least 3 shared taxa are required.",
                    "Consider using --missing-taxa allow instead.",
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
        taxa_side_len = len(taxa_side)
        all_taxa_len = len(all_taxa)
        if taxa_side_len * 2 < all_taxa_len:
            return taxa_side
        complement = all_taxa - taxa_side
        if taxa_side_len * 2 > all_taxa_len:
            return complement
        # Equal size: return the lexicographically smaller side
        if not taxa_side:
            return taxa_side
        if min(taxa_side) < min(complement):
            return taxa_side
        return complement

    @staticmethod
    def _mask_to_split(mask: int, names: tuple[str, ...]) -> frozenset:
        taxa = []
        while mask:
            bit = mask & -mask
            taxa.append(names[bit.bit_length() - 1])
            mask ^= bit
        return frozenset(taxa)

    @staticmethod
    def _canonical_split_mask(
        mask: int,
        all_mask: int,
        n_taxa: int,
        names: tuple[str, ...],
    ) -> int:
        complement = all_mask ^ mask
        side_size = mask.bit_count()
        complement_size = n_taxa - side_size
        if side_size < complement_size:
            return mask
        if side_size > complement_size:
            return complement

        # ``names`` is sorted and masks are disjoint, so the first differing
        # taxon is the side with the lower set bit.
        if (mask & -mask) < (complement & -complement):
            return mask
        return complement

    @staticmethod
    def _extract_split_masks_from_tree(
        tree,
        taxon_to_bit: dict[str, int],
        names: tuple[str, ...],
    ):
        try:
            root = tree.root
            root.clades
        except AttributeError:
            return None

        n_taxa = len(names)
        all_mask = (1 << n_taxa) - 1
        split_masks = set()
        masks_by_id = {}
        clades = []
        stack = [root]

        try:
            pop = stack.pop
            extend = stack.extend
            append_clade = clades.append
            while stack:
                clade = pop()
                append_clade(clade)
                children = clade.clades
                if children:
                    extend(children)

            for clade in reversed(clades):
                children = clade.clades
                if not children:
                    try:
                        masks_by_id[id(clade)] = 1 << taxon_to_bit[clade.name]
                    except KeyError:
                        return None
                    continue

                n_children = len(children)
                if n_children == 2:
                    mask = (
                        masks_by_id[id(children[0])]
                        | masks_by_id[id(children[1])]
                    )
                else:
                    mask = 0
                    for child in children:
                        mask |= masks_by_id[id(child)]
                masks_by_id[id(clade)] = mask

                # Skip polytomous nodes (>2 children = unresolved branching),
                # but allow trifurcating roots (standard unrooted Newick).
                if n_children > 2:
                    is_root = clade == root
                    if not (is_root and n_children == 3):
                        continue

                # Skip trivial splits (single taxon or all-but-one), plus the
                # root-level all-taxa clade.
                split_size = mask.bit_count()
                if split_size <= 1 or split_size >= n_taxa - 1:
                    continue
                if mask == all_mask:
                    continue

                split_masks.add(
                    ConsensusNetwork._canonical_split_mask(
                        mask,
                        all_mask,
                        n_taxa,
                        names,
                    )
                )
        except AttributeError:
            return None

        return split_masks

    @staticmethod
    def _extract_splits_from_tree_legacy(
        tree,
        all_taxa: frozenset,
    ) -> set[frozenset]:
        splits = set()
        clade_taxa = {}

        for clade in tree.find_clades(order="postorder"):
            if clade.is_terminal():
                tips = frozenset([clade.name])
                clade_taxa[id(clade)] = tips
                continue

            tips = frozenset().union(
                *(clade_taxa[id(child)] for child in clade.clades)
            )
            clade_taxa[id(clade)] = tips

            # Skip polytomous nodes (>2 children = unresolved branching),
            # but allow trifurcating roots (standard unrooted Newick).
            n_children = len(clade.clades)
            if n_children > 2:
                is_root = clade == tree.root
                if not (is_root and n_children == 3):
                    continue
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
    def _extract_splits_from_tree(tree, all_taxa: frozenset) -> set[frozenset]:
        names = tuple(sorted(all_taxa))
        split_masks = ConsensusNetwork._extract_split_masks_from_tree(
            tree,
            {name: index for index, name in enumerate(names)},
            names,
        )
        if split_masks is not None:
            return {
                ConsensusNetwork._mask_to_split(mask, names)
                for mask in split_masks
            }

        return ConsensusNetwork._extract_splits_from_tree_legacy(tree, all_taxa)

    @staticmethod
    def _count_splits(trees: list, all_taxa: frozenset,
                      allow_mode: bool = False) -> tuple[Counter, Counter]:
        """Count splits across trees.

        Returns (split_counts, split_possible) where split_possible[s]
        is the number of trees that contain ALL taxa in split s (and
        its complement). In allow mode, each tree uses its own taxon
        set; in shared mode, all trees use all_taxa.
        """
        counter = Counter()
        possible = Counter()

        if allow_mode:
            # Precompute taxon sets for all trees
            tree_taxa_list = [
                frozenset(ConsensusNetwork._tips(tree))
                for tree in trees
            ]

            # Extract splits from each tree using its own taxon set
            for tree, tree_taxa in zip(trees, tree_taxa_list):
                tree_splits = ConsensusNetwork._extract_splits_from_tree(
                    tree, tree_taxa
                )
                counter.update(tree_splits)

            # Actually, we should count how many trees COULD have
            # produced the split but didn't. A more accurate approach:
            # possible = number of trees containing all taxa in the
            # split's smaller side. But since splits are defined
            # relative to each tree's own taxon set, the split IS
            # the canonical smaller side from that tree. Different
            # trees may have different "all_taxa" so the same
            # bipartition in two trees means different things.
            #
            # The simplest correct normalization for incomplete
            # taxon sampling: frequency = count / n_trees.
            # This is what most software does.
            possible = Counter(dict.fromkeys(counter, len(trees)))
        else:
            names = tuple(sorted(all_taxa))
            taxon_to_bit = {name: index for index, name in enumerate(names)}
            mask_counter = Counter()
            use_mask_counter = True
            for tree in trees:
                split_masks = ConsensusNetwork._extract_split_masks_from_tree(
                    tree,
                    taxon_to_bit,
                    names,
                )
                if split_masks is None:
                    use_mask_counter = False
                    break
                mask_counter.update(split_masks)

            if use_mask_counter:
                counter.update(
                    {
                        ConsensusNetwork._mask_to_split(mask, names): count
                        for mask, count in mask_counter.items()
                    }
                )
            else:
                for tree in trees:
                    tree_splits = ConsensusNetwork._extract_splits_from_tree(
                        tree, all_taxa
                    )
                    counter.update(tree_splits)

            possible = Counter(dict.fromkeys(counter, len(trees)))

        return counter, possible

    @staticmethod
    def _filter_splits(
        split_counts: Counter, n_trees: int, threshold: float,
        split_possible: Counter = None,
    ) -> list[tuple[frozenset, int, float]]:
        results = []
        for split, count in split_counts.items():
            if split_possible and split in split_possible:
                denom = split_possible[split]
            else:
                denom = n_trees
            freq = count / denom if denom > 0 else 0.0
            if freq >= threshold:
                results.append((split, count, freq))
        results.sort(key=lambda x: (-x[2], sorted(x[0])))
        return results

    # ------------------------------------------------------------------
    # Circular ordering for visualization
    # ------------------------------------------------------------------

    @staticmethod
    def _compute_circular_ordering(trees: list, all_taxa: frozenset) -> list[str]:
        try:
            consensus = Consensus.majority_consensus(trees, cutoff=0.5)
            ordering = Tree.calculate_terminal_names_fast(consensus)
            if ordering is None:
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
    def _circular_gap_positions(split, ordering):
        """Return the two circular boundary gap positions, or None."""
        n = len(ordering)
        if n == 0:
            return None
        gap_positions = []
        first_in_split = ordering[0] in split
        previous_in_split = first_in_split
        for i in range(1, n):
            current_in_split = ordering[i] in split
            if previous_in_split != current_in_split:
                gap_positions.append(i - 1)
                if len(gap_positions) > 2:
                    return None
            previous_in_split = current_in_split
        if previous_in_split != first_in_split:
            gap_positions.append(n - 1)
        if len(gap_positions) == 2:
            return (gap_positions[0], gap_positions[1])
        return None

    @staticmethod
    def _is_circular_split(split, ordering):
        """Check if a split has exactly 2 boundary gaps in the circular ordering."""
        return ConsensusNetwork._circular_gap_positions(split, ordering) is not None

    @staticmethod
    def _compute_split_directions(ordering, circular_splits,
                                  gap_positions_by_split=None):
        """Compute 2D direction vectors for each circular split."""
        n = len(ordering)
        angles = {taxon: 2 * math.pi * i / n for i, taxon in enumerate(ordering)}
        directions = {}
        for split, count, freq in circular_splits:
            if gap_positions_by_split is None:
                gap_positions = ConsensusNetwork._circular_gap_positions(
                    split,
                    ordering,
                )
            else:
                gap_positions = gap_positions_by_split.get(split)
            if gap_positions is None:
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
        bit_values = [1 << idx for idx in range(n_splits)]
        # Each taxon's sign vector: +1 if in canonical side, -1 otherwise
        taxon_signs = {}
        taxon_masks = []
        for taxon in all_taxa:
            sign_mask = 0
            signs = []
            for split_idx, sp in enumerate(splits_list):
                if taxon in sp:
                    signs.append(1)
                    sign_mask |= bit_values[split_idx]
                else:
                    signs.append(-1)
            taxon_signs[taxon] = tuple(signs)
            taxon_masks.append(sign_mask)
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
        valid_masks = set()
        constraints_by_split = [[] for _ in range(n_splits)]
        for (i, j), fb_pairs in forbidden.items():
            constraints_by_split[j].append((i, fb_pairs))

        signs = [-1] * n_splits

        def add_valid_signs(split_idx, sign_mask):
            if split_idx == n_splits:
                valid_masks.add(sign_mask)
                return
            bit_value = bit_values[split_idx]
            for sign, next_mask in (
                (-1, sign_mask),
                (1, sign_mask | bit_value),
            ):
                valid = True
                for prev_idx, fb_pairs in constraints_by_split[split_idx]:
                    if (signs[prev_idx], sign) in fb_pairs:
                        valid = False
                        break
                if valid:
                    signs[split_idx] = sign
                    add_valid_signs(split_idx + 1, next_mask)

        add_valid_signs(0, 0)
        # Ensure taxon sign vectors are included
        valid_masks.update(taxon_masks)

        signs_cache = {}

        def mask_to_signs(mask):
            try:
                return signs_cache[mask]
            except KeyError:
                sign_tuple = tuple(1 if mask & bit else -1 for bit in bit_values)
                signs_cache[mask] = sign_tuple
                return sign_tuple

        valid_nodes = {mask_to_signs(mask) for mask in valid_masks}
        edges = []
        for mask in valid_masks:
            node = mask_to_signs(mask)
            for split_idx, bit_value in enumerate(bit_values):
                flipped = mask ^ bit_value
                if flipped in valid_masks and mask < flipped:
                    edges.append((node, mask_to_signs(flipped), split_idx))
        return taxon_signs, valid_nodes, edges, splits_list, weights

    def _draw_network(self, ordering: list[str], filtered_splits, all_taxa: frozenset, output_path: str):
        """Draw a planar splits graph using matplotlib."""
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.collections import LineCollection

        n = len(ordering)
        angles = {taxon: 2 * math.pi * i / n for i, taxon in enumerate(ordering)}

        # Filter to circular splits only
        circular_splits = []
        gap_positions_by_split = {}
        for split, count, freq in filtered_splits:
            gap_positions = self._circular_gap_positions(split, ordering)
            if gap_positions is not None:
                circular_splits.append((split, count, freq))
                gap_positions_by_split[split] = gap_positions

        # Compute direction vectors
        directions = self._compute_split_directions(
            ordering,
            circular_splits,
            gap_positions_by_split,
        )

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

        config = self.plot_config
        n = len(ordering)
        # Force square figure for network graphs
        if config.fig_width is None and config.fig_height is None:
            size = max(10, min(30, 8 + n * 0.05))
            config.fig_width = size
            config.fig_height = size
        config.resolve(n_rows=n, n_cols=None)
        fig, ax = plt.subplots(1, 1, figsize=(config.fig_width, config.fig_height))
        ax.set_aspect("equal")

        # Determine label fontsize: auto-suppress for large trees
        if config.ylabel_fontsize is not None:
            label_fontsize = config.ylabel_fontsize
        elif n > 100:
            label_fontsize = 0  # auto-hide for very large trees
        elif n > 50:
            label_fontsize = max(3, 8 - (n - 50) * 0.1)
        else:
            label_fontsize = 10

        if not splits_list:
            # No splits: place taxa evenly on a circle
            point_x = []
            point_y = []
            for taxon in ordering:
                angle = angles[taxon]
                x = math.cos(angle)
                y = math.sin(angle)
                if label_fontsize > 0:
                    deg = math.degrees(angle)
                    ha = "left" if -90 < deg < 90 or deg > 270 else "right"
                    rotation = deg if -90 < deg < 90 or deg > 270 else deg + 180
                    ax.text(x, y, taxon, ha=ha, va="center",
                            fontsize=label_fontsize, rotation=rotation,
                            rotation_mode="anchor")
                else:
                    point_x.append(x)
                    point_y.append(y)
            if point_x:
                ax.scatter(point_x, point_y, color="black", s=4)
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
            extent = _position_extent(node_positions.values())
            pendant_len = max(0.15, extent * 0.3)

            # Draw internal edges
            internal_segments = []
            for s1, s2, split_idx in edges:
                x1, y1 = node_positions[s1]
                x2, y2 = node_positions[s2]
                internal_segments.append(((x1, y1), (x2, y2)))
            if internal_segments:
                ax.add_collection(
                    LineCollection(
                        internal_segments,
                        colors="black",
                        linewidths=1.5,
                        zorder=2,
                    ),
                    autolim=True,
                )

            # Find taxa on split boundaries (participating in displayed splits)
            boundary_taxa = set()
            for split, count, freq in circular_splits:
                for i in range(n):
                    curr = ordering[i]
                    nxt = ordering[(i + 1) % n]
                    if (curr in split) != (nxt in split):
                        boundary_taxa.add(curr)
                        boundary_taxa.add(nxt)

            # Draw pendant edges and taxon labels
            pendant_segments = []
            pendant_widths = []
            pendant_colors = []
            for taxon in ordering:
                angle = angles[taxon]
                if taxon in taxon_signs and taxon_signs[taxon] in node_positions:
                    nx, ny = node_positions[taxon_signs[taxon]]
                else:
                    nx, ny = 0.0, 0.0

                # Only draw full pendant edges for boundary taxa;
                # non-boundary taxa get a short stub or are skipped
                if taxon in boundary_taxa or n <= 50:
                    plen = pendant_len
                else:
                    plen = pendant_len * 0.3  # short stub

                tx = nx + plen * math.cos(angle)
                ty = ny + plen * math.sin(angle)
                is_boundary = taxon in boundary_taxa
                pendant_segments.append(((nx, ny), (tx, ty)))
                pendant_widths.append(1.5 if is_boundary else 0.8)
                pendant_colors.append((0.0, 0.0, 0.0, 1.0 if is_boundary else 0.3))

                if label_fontsize > 0 and (taxon in boundary_taxa or n <= 50):
                    lx = tx + 0.03 * math.cos(angle)
                    ly = ty + 0.03 * math.sin(angle)
                    deg = math.degrees(angle)
                    ha = "left" if -90 < deg < 90 or deg > 270 else "right"
                    rotation = deg if -90 < deg < 90 or deg > 270 else deg + 180
                    ax.text(lx, ly, taxon, ha=ha, va="center",
                            fontsize=label_fontsize, rotation=rotation,
                            rotation_mode="anchor", zorder=4)
            if pendant_segments:
                ax.add_collection(
                    LineCollection(
                        pendant_segments,
                        colors=pendant_colors,
                        linewidths=pendant_widths,
                        zorder=2,
                    ),
                    autolim=True,
                )
            ax.autoscale_view()

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
        if config.show_title:
            ax.set_title(config.title or "Consensus Splits Network", fontsize=config.title_fontsize or 14, pad=20)

        plt.tight_layout()
        plt.savefig(output_path, dpi=config.dpi, bbox_inches="tight")
        plt.close()

    # ------------------------------------------------------------------
    # Output
    # ------------------------------------------------------------------

    def _draw_histogram(self, split_counts, n_trees, output_path):
        """Draw a histogram of split frequencies."""
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        config = self.plot_config
        config.resolve(n_rows=10, n_cols=None)

        frequencies = [count / n_trees for count in split_counts.values()]

        fig, ax = plt.subplots(figsize=(config.fig_width or 10, config.fig_height or 6))
        ax.hist(frequencies, bins=50, color="#377eb8", edgecolor="black", linewidth=0.5)
        ax.set_xlabel("Split frequency", fontsize=config.axis_fontsize or 12)
        ax.set_ylabel("Number of splits", fontsize=config.axis_fontsize or 12)
        ax.axvline(x=self.threshold, color="red", linestyle="--", lw=1,
                   label=f"Threshold ({self.threshold})")
        ax.legend(fontsize=9)

        if config.show_title:
            ax.set_title(
                config.title or "Split Frequency Distribution",
                fontsize=config.title_fontsize or 14,
            )

        fig.tight_layout()
        fig.savefig(output_path, dpi=config.dpi, bbox_inches="tight")
        plt.close(fig)

    def _format_split(self, split: frozenset) -> str:
        return "{" + ", ".join(sorted(split)) + "}"

    def run(self):
        trees = self._parse_trees_from_source(self.trees)
        if len(trees) == 0:
            raise PhykitUserError(["No trees were parsed from input."], code=2)

        trees, pruned, all_taxa_set = self._normalize_taxa(trees)
        all_taxa = frozenset(all_taxa_set)
        n_trees = len(trees)

        allow_mode = (self.missing_taxa == "allow")
        split_counts, split_possible = self._count_splits(
            trees, all_taxa, allow_mode=allow_mode
        )
        filtered = self._filter_splits(
            split_counts, n_trees, self.threshold,
            split_possible=split_possible,
        )

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
            lines = [
                f"Number of input trees: {n_trees}",
                f"Number of taxa: {len(all_taxa)}",
                f"Threshold: {self.threshold}",
            ]
            if pruned:
                lines.append("Pruned to shared taxa: yes")
            lines.extend(
                [
                    f"Total unique splits: {len(split_counts)}",
                    f"Splits above threshold: {len(filtered)}",
                    "---",
                ]
            )
            lines.extend(
                f"{self._format_split(split)}\t{count}/{n_trees}\t{freq:.4f}"
                for split, count, freq in filtered
            )
            print("\n".join(lines))

        if self.histogram:
            self._draw_histogram(split_counts, n_trees, self.histogram)
            if not self.json_output:
                print(f"Histogram saved: {self.histogram}")

        if self.plot_output:
            import sys
            n_taxa = len(all_taxa)

            if n_taxa > 100:
                print(
                    f"Warning: {n_taxa} taxa — network graph may not be "
                    f"informative at this scale. Consider using "
                    f"--histogram for a split frequency distribution instead.",
                    file=sys.stderr,
                )

            ordering = self._compute_circular_ordering(trees, all_taxa)

            # Cap splits for graph visualization to avoid exponential blowup
            plot_filtered = filtered
            if len(filtered) > self.max_splits:
                print(
                    f"Warning: {len(filtered)} splits above threshold; "
                    f"using top {self.max_splits} for network graph "
                    f"(use --max-splits to adjust).",
                    file=sys.stderr,
                )
                # filtered is already sorted by frequency (descending)
                plot_filtered = filtered[:self.max_splits]

            self._draw_network(ordering, plot_filtered, all_taxa, self.plot_output)
