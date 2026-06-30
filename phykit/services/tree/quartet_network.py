import itertools
import math
from collections import Counter
from io import StringIO
import os
from pathlib import Path

from .base import Tree
from ...errors import PhykitUserError


_path_exists = os.path.exists
_path_isabs = os.path.isabs


class _LazyPhylo:
    def read(self, *args, **kwargs):
        from Bio import Phylo as _Phylo

        return _Phylo.read(*args, **kwargs)


Phylo = _LazyPhylo()


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


def _chi2_sf(x, df):
    if df == 1:
        return math.erfc(math.sqrt(x / 2.0))
    if df == 2:
        return math.exp(-x / 2.0)

    from scipy.stats import chi2
    return float(chi2.sf(x, df))


class QuartetNetwork(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(trees=parsed["trees"])
        self.alpha = parsed["alpha"]
        self.beta = parsed["beta"]
        self.missing_taxa = parsed["missing_taxa"]
        self.plot_output = parsed["plot_output"]
        self.json_output = parsed["json_output"]
        self.plot_config = parsed["plot_config"]

    def process_args(self, args) -> dict[str, str]:
        from ...helpers.plot_config import PlotConfig

        return dict(
            trees=args.trees,
            alpha=args.alpha,
            beta=getattr(args, "beta", 0.95),
            missing_taxa=args.missing_taxa,
            plot_output=getattr(args, "plot_output", None),
            json_output=getattr(args, "json", False),
            plot_config=PlotConfig.from_args(args),
        )

    # ------------------------------------------------------------------
    # Tree parsing (copied per codebase convention)
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
        names = Tree.calculate_terminal_names_fast(tree)
        if names is not None:
            return set(names)
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
                    "Use --missing-taxa shared to prune all trees to their shared taxa.",
                ],
                code=2,
            )

        shared_taxa = set.intersection(*tip_sets)
        if len(shared_taxa) < 4:
            raise PhykitUserError(
                [
                    "Unable to compute quartet network after pruning to shared taxa.",
                    "At least 4 shared taxa are required.",
                ],
                code=2,
            )

        pruned = [self._prune_to_taxa(tree, shared_taxa) for tree in trees]
        return pruned, True, shared_taxa

    # ------------------------------------------------------------------
    # Bipartition extraction
    # ------------------------------------------------------------------

    @staticmethod
    def _extract_bipartitions(tree, all_taxa: frozenset) -> list[tuple[frozenset, frozenset]]:
        direct_result = QuartetNetwork._extract_bipartitions_direct(tree, all_taxa)
        if direct_result is not None:
            return direct_result

        bipartitions = []
        clade_taxa = {}

        for clade in tree.find_clades(order="postorder"):
            if clade.is_terminal():
                clade_taxa[id(clade)] = frozenset({clade.name})
            else:
                taxa = frozenset()
                for child in clade.clades:
                    taxa = taxa | clade_taxa.get(id(child), frozenset())
                clade_taxa[id(clade)] = taxa

        for clade in tree.get_nonterminals():
            # Skip polytomous nodes (>2 children = unresolved branching),
            # but allow trifurcating roots (standard unrooted Newick).
            n_children = len(clade.clades)
            if n_children > 2:
                is_root = (clade == tree.root)
                if not (is_root and n_children == 3):
                    continue
            tips = clade_taxa.get(id(clade), frozenset())
            if len(tips) <= 1 or len(tips) >= len(all_taxa) - 1:
                continue
            if tips == all_taxa:
                continue
            complement = all_taxa - tips
            bipartitions.append((tips, complement))
        return bipartitions

    @staticmethod
    def _extract_bipartitions_direct(
        tree,
        all_taxa: frozenset,
    ) -> list[tuple[frozenset, frozenset]] | None:
        try:
            root = tree.root
            root.clades
        except AttributeError:
            return None

        preorder = []
        stack = [root]
        pop = stack.pop
        append = preorder.append
        extend = stack.extend
        try:
            while stack:
                clade = pop()
                append(clade)
                children = clade.clades
                if children:
                    extend(reversed(children))
        except AttributeError:
            return None

        empty = frozenset()
        clade_taxa = {}
        try:
            for clade in reversed(preorder):
                children = clade.clades
                if not children:
                    clade_taxa[id(clade)] = frozenset({clade.name})
                elif len(children) == 1:
                    clade_taxa[id(clade)] = clade_taxa.get(id(children[0]), empty)
                elif len(children) == 2:
                    clade_taxa[id(clade)] = (
                        clade_taxa.get(id(children[0]), empty)
                        | clade_taxa.get(id(children[1]), empty)
                    )
                else:
                    clade_taxa[id(clade)] = empty.union(
                        *(clade_taxa.get(id(child), empty) for child in children)
                    )

            n_taxa = len(all_taxa)
            bipartitions = []
            for clade in preorder:
                children = clade.clades
                if not children:
                    continue

                n_children = len(children)
                if n_children > 2:
                    is_root = clade is root
                    if not (is_root and n_children == 3):
                        continue

                tips = clade_taxa.get(id(clade), empty)
                if len(tips) <= 1 or len(tips) >= n_taxa - 1:
                    continue
                if tips == all_taxa:
                    continue
                bipartitions.append((tips, all_taxa - tips))
        except AttributeError:
            return None

        return bipartitions

    # ------------------------------------------------------------------
    # Quartet topology determination
    # ------------------------------------------------------------------

    @staticmethod
    def _determine_quartet_topology(
        quartet: tuple[str, str, str, str],
        bipartitions: list[tuple[frozenset, frozenset]],
    ) -> int | None:
        """Determine which of the 3 unrooted topologies a gene tree displays
        for the given quartet.

        Returns:
            0 if ab|cd, 1 if ac|bd, 2 if ad|bc, or None if unresolved (star).
        """
        a, b, c, d = quartet
        quartet_set = {a, b, c, d}
        for side_a, side_b in bipartitions:
            in_a = quartet_set & side_a
            in_b = quartet_set & side_b
            if len(in_a) == 2 and len(in_b) == 2:
                pair = frozenset(in_a)
                if pair == frozenset({a, b}) or pair == frozenset({c, d}):
                    return 0  # ab|cd
                if pair == frozenset({a, c}) or pair == frozenset({b, d}):
                    return 1  # ac|bd
                if pair == frozenset({a, d}) or pair == frozenset({b, c}):
                    return 2  # ad|bc
        return None

    @staticmethod
    def _postorder_clades_direct(tree):
        try:
            root = tree.root
            root.clades
        except AttributeError:
            return None

        clades = []
        stack = [root]
        append = clades.append
        pop = stack.pop
        extend = stack.extend
        try:
            while stack:
                clade = pop()
                append(clade)
                children = clade.clades
                if children:
                    extend(children)
        except AttributeError:
            return None
        clades.reverse()
        return clades

    @staticmethod
    def _extract_bipartition_masks(tree, all_taxa: frozenset, taxa_index: dict[str, int]) -> list[int]:
        split_masks = []
        clade_masks = {}
        n_taxa = len(all_taxa)
        root = tree.root
        postorder = QuartetNetwork._postorder_clades_direct(tree)
        if postorder is None:
            postorder = tree.find_clades(order="postorder")

        for clade in postorder:
            children = clade.clades
            if not children:
                idx = taxa_index.get(clade.name)
                clade_masks[id(clade)] = 0 if idx is None else 1 << idx
                continue

            n_children = len(children)
            if n_children == 2:
                mask = (
                    clade_masks.get(id(children[0]), 0)
                    | clade_masks.get(id(children[1]), 0)
                )
            elif n_children == 1:
                mask = clade_masks.get(id(children[0]), 0)
            else:
                mask = 0
                for child in children:
                    mask |= clade_masks.get(id(child), 0)
            clade_masks[id(clade)] = mask

            if n_children > 2:
                is_root = (clade is root)
                if not (is_root and n_children == 3):
                    continue

            n_tips = mask.bit_count()
            if n_tips <= 1 or n_tips >= n_taxa - 1:
                continue
            split_masks.append(mask)

        return split_masks

    @staticmethod
    def _determine_quartet_topology_from_masks(
        quartet_bits: tuple[int, int, int, int],
        split_masks: list[int],
    ) -> int | None:
        bit_a, bit_b, bit_c, bit_d = quartet_bits
        quartet_mask = bit_a | bit_b | bit_c | bit_d
        ab = bit_a | bit_b
        ac = bit_a | bit_c
        ad = bit_a | bit_d
        cd = bit_c | bit_d
        bd = bit_b | bit_d
        bc = bit_b | bit_c

        for split_mask in split_masks:
            pair = quartet_mask & split_mask
            if pair.bit_count() != 2:
                continue
            if pair == ab or pair == cd:
                return 0
            if pair == ac or pair == bd:
                return 1
            if pair == ad or pair == bc:
                return 2
        return None

    @staticmethod
    def _determine_quartet_topology_from_precomputed_masks(
        quartet_data: tuple[int, int, int, int, int, int, int],
        split_masks: list[int],
    ) -> int | None:
        quartet_mask, ab, ac, ad, cd, bd, bc = quartet_data
        for split_mask in split_masks:
            pair = quartet_mask & split_mask
            if pair.bit_count() != 2:
                continue
            if pair == ab or pair == cd:
                return 0
            if pair == ac or pair == bd:
                return 1
            if pair == ad or pair == bc:
                return 2
        return None

    # ------------------------------------------------------------------
    # Compute quartet concordance factors
    # ------------------------------------------------------------------

    @staticmethod
    def _compute_quartet_cfs(
        trees: list,
        all_taxa: frozenset,
    ) -> dict[tuple[str, str, str, str], list[int]]:
        taxa_sorted = sorted(all_taxa)
        taxa_index = {taxon: idx for idx, taxon in enumerate(taxa_sorted)}
        all_split_masks = []
        for tree in trees:
            split_masks = QuartetNetwork._extract_bipartition_masks(
                tree,
                all_taxa,
                taxa_index,
            )
            all_split_masks.append(split_masks)

        quartet_data = []
        for quartet_indices in itertools.combinations(range(len(taxa_sorted)), 4):
            bit_a, bit_b, bit_c, bit_d = (1 << i for i in quartet_indices)
            quartet = tuple(taxa_sorted[i] for i in quartet_indices)
            quartet_data.append(
                (
                    quartet,
                    (
                        bit_a | bit_b | bit_c | bit_d,
                        bit_a | bit_b,
                        bit_a | bit_c,
                        bit_a | bit_d,
                        bit_c | bit_d,
                        bit_b | bit_d,
                        bit_b | bit_c,
                    ),
                )
            )

        results = {}
        topology_from_masks = (
            QuartetNetwork._determine_quartet_topology_from_precomputed_masks
        )
        for quartet, precomputed_masks in quartet_data:
            counts = [0, 0, 0]
            for split_masks in all_split_masks:
                topo = topology_from_masks(
                    precomputed_masks,
                    split_masks,
                )
                if topo is not None:
                    counts[topo] += 1
            results[quartet] = counts
        return results

    # ------------------------------------------------------------------
    # Statistical tests matching MSCquartets/NANUQ
    # ------------------------------------------------------------------

    @staticmethod
    def _star_test(counts: list[int]) -> float:
        """Star tree test (Pearson chi-squared, H0: all three topologies
        equally likely). Matches MSCquartets quartetStarTest.

        Returns p-value.
        """
        total = sum(counts)
        if total == 0:
            return 1.0
        expected = total / 3.0
        stat = sum((obs - expected) ** 2 / expected for obs in counts)
        return _chi2_sf(stat, df=2)

    @staticmethod
    def _tree_test(counts: list[int]) -> float:
        """T3 tree model test (G-test / likelihood ratio).

        Under H0 (any resolved quartet tree under MSC), the dominant
        topology count is the MLE and the two minor topologies should
        be equal.  Expected values: [max_count, (n-max_count)/2,
        (n-max_count)/2].  G-statistic with conservative chi-sq(df=1)
        p-value, matching MSCquartets quartetTreeTest(method="conservative").

        Returns p-value.
        """
        total = sum(counts)
        if total == 0:
            return 1.0

        sorted_counts = sorted(counts, reverse=True)
        major = sorted_counts[0]
        remainder = total - major

        if remainder == 0:
            # All counts in one topology — perfect tree, p=1
            return 1.0

        expected = [major, remainder / 2.0, remainder / 2.0]

        # G-statistic (power divergence with lambda=0)
        g_stat = 0.0
        for obs, exp in zip(sorted_counts, expected):
            if obs > 0 and exp > 0:
                g_stat += 2.0 * obs * math.log(obs / exp)

        # Conservative p-value: chi-squared with df=1
        return _chi2_sf(g_stat, df=1)

    @staticmethod
    def _classify_quartet(
        counts: list[int], alpha: float, beta: float
    ) -> dict[str, object]:
        """Classify a quartet following the NANUQ algorithm.

        1. Star test: if p_star > beta, classify as 'unresolved' (star tree).
        2. Tree test (T3 model): if p_T3 > alpha, classify as 'tree'.
        3. Otherwise: classify as 'hybrid' (4-cycle / reticulation).

        Parameters match MSCquartets defaults: alpha=0.05, beta=0.95.
        """
        total = sum(counts)
        if total == 0:
            return {
                "classification": "unresolved",
                "p_star": 1.0,
                "p_tree": 1.0,
            }

        p_star = QuartetNetwork._star_test(counts)
        p_tree = QuartetNetwork._tree_test(counts)

        if p_star > beta:
            classification = "unresolved"
        elif p_tree > alpha:
            classification = "tree"
        else:
            classification = "hybrid"

        return {
            "classification": classification,
            "p_star": p_star,
            "p_tree": p_tree,
        }

    # ------------------------------------------------------------------
    # NANUQ distance matrix
    # ------------------------------------------------------------------

    @staticmethod
    def _cross_pairs(topo_idx, ia, ib, ic, id_):
        """Return cross-split taxon index pairs for a given topology."""
        if topo_idx == 0:  # ab|cd
            return [(ia, ic), (ia, id_), (ib, ic), (ib, id_)]
        elif topo_idx == 1:  # ac|bd
            return [(ia, ib), (ia, id_), (ic, ib), (ic, id_)]
        else:  # ad|bc
            return [(ia, ib), (ia, ic), (id_, ib), (id_, ic)]

    @staticmethod
    def _dominant_topology_index(counts):
        """Return the first index of the largest of three topology counts."""
        c0, c1, c2 = counts
        if c1 > c0:
            if c2 > c1:
                return 2
            return 1
        if c2 > c0:
            return 2
        return 0

    @staticmethod
    def _top_two_topology_indices(counts):
        """Return the first two topology indices by descending count."""
        c0, c1, c2 = counts
        if c0 >= c1:
            if c1 >= c2:
                return 0, 1
            if c0 >= c2:
                return 0, 2
            return 2, 0
        if c0 >= c2:
            return 1, 0
        if c1 >= c2:
            return 1, 2
        return 2, 1

    @staticmethod
    def _compute_nanuq_distance(all_taxa, quartet_results):
        """Compute NANUQ distance matrix from quartet classifications.

        Following Allman, Banos & Rhodes (2019):
        - Star quartets: +1 to all 6 pairwise distances
        - Tree quartets: +1 to the 4 cross-split pairs
        - Hybrid quartets: +0.5 to cross-split pairs for the two elevated topologies
        """
        taxa_list = sorted(all_taxa)
        n = len(taxa_list)
        taxa_idx = {t: i for i, t in enumerate(taxa_list)}
        dist = [[0.0] * n for _ in range(n)]

        for quartet, result in quartet_results.items():
            a, b, c, d = quartet
            ia, ib, ic, id_ = taxa_idx[a], taxa_idx[b], taxa_idx[c], taxa_idx[d]
            classification = result["classification"]
            counts = result["counts"]

            if classification == "unresolved":
                for x, y in [(ia, ib), (ia, ic), (ia, id_), (ib, ic), (ib, id_), (ic, id_)]:
                    dist[x][y] += 1
                    dist[y][x] += 1
            elif classification == "tree":
                dominant_idx = QuartetNetwork._dominant_topology_index(counts)
                for x, y in QuartetNetwork._cross_pairs(dominant_idx, ia, ib, ic, id_):
                    dist[x][y] += 1
                    dist[y][x] += 1
            else:  # hybrid
                for topo_idx in QuartetNetwork._top_two_topology_indices(counts):
                    for x, y in QuartetNetwork._cross_pairs(topo_idx, ia, ib, ic, id_):
                        dist[x][y] += 0.5
                        dist[y][x] += 0.5

        return taxa_list, dist, taxa_idx

    # ------------------------------------------------------------------
    # Circular ordering via Neighbor-Joining
    # ------------------------------------------------------------------

    @staticmethod
    def _neighbor_joining_order(taxa_list, dist_matrix):
        """Compute a circular ordering of taxa using Neighbor-Joining.

        Standard NJ agglomeration; leaf order at each merged node is
        preserved to derive a circular ordering suitable for splits graphs.
        """
        n = len(taxa_list)
        if n <= 2:
            return taxa_list[:]

        node_taxa = {i: [taxa_list[i]] for i in range(n)}
        active = list(range(n))
        d = {}
        for i in range(n):
            for j in range(n):
                d[(i, j)] = dist_matrix[i][j]

        next_id = n
        while len(active) > 2:
            m = len(active)
            r = {}
            for i in active:
                r[i] = sum(d.get((i, j), 0) for j in active if j != i)

            best_q = float("inf")
            best_i, best_j = active[0], active[1]
            for ai in range(m):
                for bi in range(ai + 1, m):
                    i, j = active[ai], active[bi]
                    q = (m - 2) * d.get((i, j), 0) - r[i] - r[j]
                    if q < best_q:
                        best_q = q
                        best_i, best_j = i, j

            u = next_id
            next_id += 1
            node_taxa[u] = node_taxa[best_i] + node_taxa[best_j]

            for k in active:
                if k != best_i and k != best_j:
                    d_uk = (d.get((best_i, k), 0) + d.get((best_j, k), 0)
                            - d.get((best_i, best_j), 0)) / 2
                    d[(u, k)] = d_uk
                    d[(k, u)] = d_uk
            d[(u, u)] = 0

            active = [k for k in active if k != best_i and k != best_j]
            active.append(u)

        return node_taxa[active[0]] + node_taxa[active[1]]

    # ------------------------------------------------------------------
    # Circular split weights from distance matrix
    # ------------------------------------------------------------------

    @staticmethod
    def _compute_circular_split_weights(ordering, dist_matrix, taxa_idx):
        """Compute isolation-index weights for all non-trivial circular splits.

        For a circular split with boundary gaps between (x_{i-1}, x_i) and
        (x_j, x_{j+1}):
            w = [d(x_{i-1}, x_j) + d(x_i, x_{j+1})
                 - d(x_{i-1}, x_{j+1}) - d(x_i, x_j)] / 2
        """
        n = len(ordering)
        all_taxa_set = frozenset(ordering)
        seen = {}

        for arc_len in range(2, n - 1):
            for start in range(n):
                arc_taxa = frozenset(ordering[(start + k) % n] for k in range(arc_len))
                complement = all_taxa_set - arc_taxa

                if len(arc_taxa) < len(complement):
                    canonical = arc_taxa
                elif len(arc_taxa) > len(complement):
                    canonical = complement
                else:
                    canonical = arc_taxa if sorted(arc_taxa) < sorted(complement) else complement

                if canonical in seen:
                    continue

                x_before = ordering[(start - 1) % n]
                x_start = ordering[start % n]
                x_end = ordering[(start + arc_len - 1) % n]
                x_after = ordering[(start + arc_len) % n]

                ib = taxa_idx[x_before]
                is_ = taxa_idx[x_start]
                ie = taxa_idx[x_end]
                ia = taxa_idx[x_after]

                w = (dist_matrix[ib][ie] + dist_matrix[is_][ia]
                     - dist_matrix[ib][ia] - dist_matrix[is_][ie]) / 2.0
                seen[canonical] = max(w, 0.0)

        return [(s, w) for s, w in seen.items() if w > 1e-10]

    # ------------------------------------------------------------------
    # Splits graph construction and drawing
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
        return QuartetNetwork._circular_gap_positions(split, ordering) is not None

    @staticmethod
    def _compute_split_directions(ordering, circular_splits,
                                  gap_positions_by_split=None):
        """Compute 2D direction vectors for each circular split."""
        n = len(ordering)
        angles = {taxon: 2 * math.pi * i / n for i, taxon in enumerate(ordering)}
        directions = {}
        for split, count, freq in circular_splits:
            if gap_positions_by_split is None:
                gap_positions = QuartetNetwork._circular_gap_positions(
                    split,
                    ordering,
                )
            else:
                gap_positions = gap_positions_by_split.get(split)
            if gap_positions is None:
                continue
            g1 = math.pi * (2 * gap_positions[0] + 1) / n
            g2 = math.pi * (2 * gap_positions[1] + 1) / n
            cx = math.cos(g2) - math.cos(g1)
            cy = math.sin(g2) - math.sin(g1)
            dx = -cy
            dy = cx
            length = math.sqrt(dx * dx + dy * dy)
            if length > 1e-10:
                dx /= length
                dy /= length
            else:
                dx, dy = 1.0, 0.0
            cx_split = sum(math.cos(angles[t]) for t in split) / len(split)
            cy_split = sum(math.sin(angles[t]) for t in split) / len(split)
            if dx * cx_split + dy * cy_split < 0:
                dx = -dx
                dy = -dy
            directions[split] = (dx, dy)
        return directions

    @staticmethod
    def _build_splits_graph(circular_splits, all_taxa):
        """Build the Buneman splits graph: valid sign vectors and edges."""
        splits_list = [s[0] for s in circular_splits]
        weights = [s[2] for s in circular_splits]
        n_splits = len(splits_list)
        if n_splits == 0:
            return {}, set(), [], [], []
        taxon_signs = {}
        for taxon in all_taxa:
            signs = tuple(1 if taxon in sp else -1 for sp in splits_list)
            taxon_signs[taxon] = signs
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
        valid_nodes = set()
        constraints_by_split = [[] for _ in range(n_splits)]
        for (i, j), fb_pairs in forbidden.items():
            constraints_by_split[j].append((i, fb_pairs))

        signs = [-1] * n_splits

        def add_valid_signs(split_idx):
            if split_idx == n_splits:
                valid_nodes.add(tuple(signs))
                return
            for sign in (-1, 1):
                valid = True
                for prev_idx, fb_pairs in constraints_by_split[split_idx]:
                    if (signs[prev_idx], sign) in fb_pairs:
                        valid = False
                        break
                if valid:
                    signs[split_idx] = sign
                    add_valid_signs(split_idx + 1)

        add_valid_signs(0)
        for signs in taxon_signs.values():
            valid_nodes.add(signs)
        edges = []
        for node in valid_nodes:
            for split_idx in range(n_splits):
                flipped = (
                    node[:split_idx]
                    + (-node[split_idx],)
                    + node[split_idx + 1:]
                )
                if flipped in valid_nodes and node < flipped:
                    edges.append((node, flipped, split_idx))
        return taxon_signs, valid_nodes, edges, splits_list, weights

    def _draw_quartet_network(
        self,
        all_taxa: frozenset,
        quartet_results: dict[tuple, dict],
        output_path: str,
    ):
        """Draw a NANUQ-style NeighborNet splits graph."""
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.collections import LineCollection

        # Compute NANUQ distance matrix and circular ordering
        taxa_list, dist_matrix, taxa_idx = self._compute_nanuq_distance(
            all_taxa, quartet_results
        )
        ordering = self._neighbor_joining_order(taxa_list, dist_matrix)
        splits_with_weights = self._compute_circular_split_weights(
            ordering, dist_matrix, taxa_idx
        )

        # Normalize weights to [0, 1]
        if splits_with_weights:
            max_w = max(w for _, w in splits_with_weights)
            if max_w > 0:
                splits_with_weights = [(s, w / max_w) for s, w in splits_with_weights]

        # Limit to top 20 splits to keep computation tractable
        splits_with_weights.sort(key=lambda x: -x[1])
        splits_with_weights = splits_with_weights[:20]

        # Format as (split, count, freq) tuples and filter to circular
        circular_splits = []
        gap_positions_by_split = {}
        for split, weight in splits_with_weights:
            gap_positions = self._circular_gap_positions(split, ordering)
            if gap_positions is not None:
                circular_splits.append((split, 1, weight))
                gap_positions_by_split[split] = gap_positions

        n = len(ordering)
        angles = {taxon: 2 * math.pi * i / n for i, taxon in enumerate(ordering)}

        directions = self._compute_split_directions(
            ordering,
            circular_splits,
            gap_positions_by_split,
        )
        circular_splits = [
            (s, c, f) for s, c, f in circular_splits if s in directions
        ]

        taxon_signs, valid_nodes, edges, splits_list, weights = (
            self._build_splits_graph(circular_splits, all_taxa)
        )

        config = self.plot_config
        config.resolve(n_rows=len(all_taxa), n_cols=None)
        fig, ax = plt.subplots(1, 1, figsize=(config.fig_width, config.fig_height))
        ax.set_aspect("equal")

        if not splits_list:
            for taxon in ordering:
                angle = angles[taxon]
                x, y = math.cos(angle), math.sin(angle)
                deg = math.degrees(angle)
                ha = "left" if -90 < deg < 90 or deg > 270 else "right"
                ax.text(x, y, taxon, ha=ha, va="center", fontsize=10)
        else:
            node_positions = {}
            for node in valid_nodes:
                x, y = 0.0, 0.0
                for i, sp in enumerate(splits_list):
                    dx, dy = directions[sp]
                    x += node[i] * weights[i] * dx / 2
                    y += node[i] * weights[i] * dy / 2
                node_positions[node] = (x, y)

            extent = _position_extent(node_positions.values())
            pendant_len = max(0.15, extent * 0.3)

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

            pendant_segments = []
            for taxon in ordering:
                angle = angles[taxon]
                if taxon in taxon_signs and taxon_signs[taxon] in node_positions:
                    nx, ny = node_positions[taxon_signs[taxon]]
                else:
                    nx, ny = 0.0, 0.0
                tx = nx + pendant_len * math.cos(angle)
                ty = ny + pendant_len * math.sin(angle)
                pendant_segments.append(((nx, ny), (tx, ty)))
                lx = tx + 0.03 * math.cos(angle)
                ly = ty + 0.03 * math.sin(angle)
                deg = math.degrees(angle)
                ha = "left" if -90 < deg < 90 or deg > 270 else "right"
                ax.text(lx, ly, taxon, ha=ha, va="center", fontsize=10, zorder=4)
            if pendant_segments:
                ax.add_collection(
                    LineCollection(
                        pendant_segments,
                        colors="black",
                        linewidths=1.5,
                        zorder=2,
                    ),
                    autolim=True,
                )
            ax.autoscale_view()

        ax.axis("off")

        plt.tight_layout()
        if config.show_title and config.title:
            ax.set_title(config.title, fontsize=config.title_fontsize)
        plt.savefig(output_path, dpi=config.dpi, bbox_inches="tight")
        plt.close()

    # ------------------------------------------------------------------
    # Output helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _format_quartet(quartet: tuple[str, str, str, str], topo_idx: int) -> str:
        a, b, c, d = quartet
        pairs_by_topo = {
            0: (f"{{{a}, {b}}}", f"{{{c}, {d}}}"),
            1: (f"{{{a}, {c}}}", f"{{{b}, {d}}}"),
            2: (f"{{{a}, {d}}}", f"{{{b}, {c}}}"),
        }
        left, right = pairs_by_topo[topo_idx]
        return f"{left} | {right}"

    # ------------------------------------------------------------------
    # Main entry point
    # ------------------------------------------------------------------

    def run(self):
        trees = self._parse_trees_from_source(self.trees)
        if len(trees) == 0:
            raise PhykitUserError(["No trees were parsed from input."], code=2)

        trees, pruned, all_taxa_set = self._normalize_taxa(trees)
        all_taxa = frozenset(all_taxa_set)

        if len(all_taxa) < 4:
            raise PhykitUserError(
                [
                    "At least 4 taxa are required for quartet analysis.",
                    f"Found {len(all_taxa)} taxa.",
                ],
                code=2,
            )

        n_trees = len(trees)

        # Compute quartet concordance factors
        quartet_cfs = self._compute_quartet_cfs(trees, all_taxa)

        # Classify each quartet
        quartet_results = {}
        tree_count = 0
        hybrid_count = 0
        unresolved_count = 0

        for quartet, counts in quartet_cfs.items():
            result = self._classify_quartet(counts, self.alpha, self.beta)
            result["counts"] = counts
            total = sum(counts)
            result["cfs"] = [c / total if total > 0 else 0.0 for c in counts]
            quartet_results[quartet] = result

            if result["classification"] == "tree":
                tree_count += 1
            elif result["classification"] == "hybrid":
                hybrid_count += 1
            else:
                unresolved_count += 1

        total_quartets = len(quartet_results)

        if self.json_output:
            quartets_list = []
            for quartet, result in quartet_results.items():
                dominant_idx = self._dominant_topology_index(result["counts"])
                quartets_list.append({
                    "taxa": list(quartet),
                    "counts": result["counts"],
                    "cfs": [round(cf, 4) for cf in result["cfs"]],
                    "classification": result["classification"],
                    "p_star": round(result["p_star"], 6),
                    "p_tree": round(result["p_tree"], 6),
                    "dominant_topology": self._format_quartet(quartet, dominant_idx),
                })
            print_json(
                dict(
                    input_tree_count=n_trees,
                    taxa=sorted(all_taxa),
                    alpha=self.alpha,
                    beta=self.beta,
                    total_quartets=total_quartets,
                    tree_count=tree_count,
                    hybrid_count=hybrid_count,
                    unresolved_count=unresolved_count,
                    pruned_to_shared_taxa=pruned,
                    quartets=quartets_list,
                )
            )
        else:
            lines = [
                f"Number of input trees: {n_trees}",
                f"Number of taxa: {len(all_taxa)}",
                f"Significance level (alpha): {self.alpha}",
                f"Star tree threshold (beta): {self.beta}",
            ]
            if pruned:
                lines.append("Pruned to shared taxa: yes")
            tree_pct = 100 * tree_count / total_quartets if total_quartets > 0 else 0
            hybrid_pct = 100 * hybrid_count / total_quartets if total_quartets > 0 else 0
            unresolved_pct = 100 * unresolved_count / total_quartets if total_quartets > 0 else 0
            lines.extend(
                [
                    f"Total quartets: {total_quartets}",
                    f"Tree-like: {tree_count} ({tree_pct:.1f}%)",
                    f"Hybrid: {hybrid_count} ({hybrid_pct:.1f}%)",
                    f"Unresolved: {unresolved_count} ({unresolved_pct:.1f}%)",
                    "---",
                ]
            )
            for quartet, result in quartet_results.items():
                counts = result["counts"]
                cfs = result["cfs"]
                dominant_idx = self._dominant_topology_index(counts)
                topo_str = self._format_quartet(quartet, dominant_idx)
                cf_str = "\t".join(f"{cf:.4f}" for cf in cfs)
                classification = result["classification"]
                p_star = result["p_star"]
                p_tree = result["p_tree"]
                lines.append(
                    f"{topo_str}\t{cf_str}\t{classification}\t"
                    f"p_star={p_star:.4f}\tp_tree={p_tree:.4f}"
                )
            print("\n".join(lines))

        if self.plot_output:
            self._draw_quartet_network(all_taxa, quartet_results, self.plot_output)
