import itertools
import math
from collections import Counter
from io import StringIO
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

from Bio import Phylo

from .base import Tree
from ...errors import PhykitUserError
from ...helpers.json_output import print_json


def _chi2_cdf(x, df):
    """Compute chi-squared CDF using scipy (lazy import)."""
    from scipy.stats import chi2
    return chi2.cdf(x, df)


def chisquare(*args, **kwargs):
    from scipy.stats import chisquare as _chisquare
    return _chisquare(*args, **kwargs)


class QuartetNetwork(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(trees=parsed["trees"])
        self.alpha = parsed["alpha"]
        self.beta = parsed["beta"]
        self.missing_taxa = parsed["missing_taxa"]
        self.plot_output = parsed["plot_output"]
        self.json_output = parsed["json_output"]

    def process_args(self, args) -> Dict[str, str]:
        return dict(
            trees=args.trees,
            alpha=args.alpha,
            beta=getattr(args, "beta", 0.95),
            missing_taxa=args.missing_taxa,
            plot_output=getattr(args, "plot_output", None),
            json_output=getattr(args, "json", False),
        )

    # ------------------------------------------------------------------
    # Tree parsing (copied per codebase convention)
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
    def _extract_bipartitions(tree, all_taxa: frozenset) -> List[Tuple[frozenset, frozenset]]:
        bipartitions = []
        for clade in tree.get_nonterminals():
            tips = frozenset(tip.name for tip in clade.get_terminals())
            if len(tips) <= 1 or len(tips) >= len(all_taxa) - 1:
                continue
            if tips == all_taxa:
                continue
            complement = all_taxa - tips
            bipartitions.append((tips, complement))
        return bipartitions

    # ------------------------------------------------------------------
    # Quartet topology determination
    # ------------------------------------------------------------------

    @staticmethod
    def _determine_quartet_topology(
        quartet: Tuple[str, str, str, str],
        bipartitions: List[Tuple[frozenset, frozenset]],
    ) -> Optional[int]:
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

    # ------------------------------------------------------------------
    # Compute quartet concordance factors
    # ------------------------------------------------------------------

    @staticmethod
    def _compute_quartet_cfs(
        trees: List,
        all_taxa: frozenset,
    ) -> Dict[Tuple[str, str, str, str], List[int]]:
        taxa_sorted = sorted(all_taxa)
        all_bipartitions = []
        for tree in trees:
            bips = QuartetNetwork._extract_bipartitions(tree, all_taxa)
            all_bipartitions.append(bips)

        results = {}
        for quartet in itertools.combinations(taxa_sorted, 4):
            counts = [0, 0, 0]
            for bips in all_bipartitions:
                topo = QuartetNetwork._determine_quartet_topology(quartet, bips)
                if topo is not None:
                    counts[topo] += 1
            results[quartet] = counts
        return results

    # ------------------------------------------------------------------
    # Statistical tests matching MSCquartets/NANUQ
    # ------------------------------------------------------------------

    @staticmethod
    def _star_test(counts: List[int]) -> float:
        """Star tree test (Pearson chi-squared, H0: all three topologies
        equally likely). Matches MSCquartets quartetStarTest.

        Returns p-value.
        """
        total = sum(counts)
        if total == 0:
            return 1.0
        expected = [total / 3.0] * 3
        stat, p = chisquare(counts, f_exp=expected)
        return float(p)

    @staticmethod
    def _tree_test(counts: List[int]) -> float:
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
        p = 1.0 - _chi2_cdf(g_stat, df=1)
        return float(p)

    @staticmethod
    def _classify_quartet(
        counts: List[int], alpha: float, beta: float
    ) -> Dict[str, object]:
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
                dominant_idx = counts.index(max(counts))
                for x, y in QuartetNetwork._cross_pairs(dominant_idx, ia, ib, ic, id_):
                    dist[x][y] += 1
                    dist[y][x] += 1
            else:  # hybrid
                sorted_with_idx = sorted(enumerate(counts), key=lambda p: -p[1])
                for topo_idx, _ in sorted_with_idx[:2]:
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
        for combo in range(2 ** n_splits):
            signs = tuple(1 if (combo >> j) & 1 else -1 for j in range(n_splits))
            valid = True
            for (i, j), fb_pairs in forbidden.items():
                if (signs[i], signs[j]) in fb_pairs:
                    valid = False
                    break
            if valid:
                valid_nodes.add(signs)
        for signs in taxon_signs.values():
            valid_nodes.add(signs)
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

    def _draw_quartet_network(
        self,
        all_taxa: frozenset,
        quartet_results: Dict[Tuple, Dict],
        output_path: str,
    ):
        """Draw a NANUQ-style NeighborNet splits graph."""
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

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
        circular_splits = [
            (split, 1, weight)
            for split, weight in splits_with_weights
            if self._is_circular_split(split, ordering)
        ]

        n = len(ordering)
        angles = {taxon: 2 * math.pi * i / n for i, taxon in enumerate(ordering)}

        directions = self._compute_split_directions(ordering, circular_splits)
        circular_splits = [
            (s, c, f) for s, c, f in circular_splits if s in directions
        ]

        taxon_signs, valid_nodes, edges, splits_list, weights = (
            self._build_splits_graph(circular_splits, all_taxa)
        )

        fig, ax = plt.subplots(1, 1, figsize=(8, 8))
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

            all_x = [p[0] for p in node_positions.values()]
            all_y = [p[1] for p in node_positions.values()]
            extent = max(
                max(all_x) - min(all_x) if all_x else 0,
                max(all_y) - min(all_y) if all_y else 0,
            )
            pendant_len = max(0.15, extent * 0.3)

            for s1, s2, split_idx in edges:
                x1, y1 = node_positions[s1]
                x2, y2 = node_positions[s2]
                ax.plot(
                    [x1, x2], [y1, y2], "-",
                    color="black", linewidth=1.5, zorder=2,
                )

            for taxon in ordering:
                angle = angles[taxon]
                if taxon in taxon_signs and taxon_signs[taxon] in node_positions:
                    nx, ny = node_positions[taxon_signs[taxon]]
                else:
                    nx, ny = 0.0, 0.0
                tx = nx + pendant_len * math.cos(angle)
                ty = ny + pendant_len * math.sin(angle)
                ax.plot(
                    [nx, tx], [ny, ty], "-",
                    color="black", linewidth=1.5, zorder=2,
                )
                lx = tx + 0.03 * math.cos(angle)
                ly = ty + 0.03 * math.sin(angle)
                deg = math.degrees(angle)
                ha = "left" if -90 < deg < 90 or deg > 270 else "right"
                ax.text(lx, ly, taxon, ha=ha, va="center", fontsize=10, zorder=4)

        ax.axis("off")

        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close()

    # ------------------------------------------------------------------
    # Output helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _format_quartet(quartet: Tuple[str, str, str, str], topo_idx: int) -> str:
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
                dominant_idx = result["counts"].index(max(result["counts"]))
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
            print(f"Number of input trees: {n_trees}")
            print(f"Number of taxa: {len(all_taxa)}")
            print(f"Significance level (alpha): {self.alpha}")
            print(f"Star tree threshold (beta): {self.beta}")
            if pruned:
                print("Pruned to shared taxa: yes")
            print(f"Total quartets: {total_quartets}")
            tree_pct = 100 * tree_count / total_quartets if total_quartets > 0 else 0
            hybrid_pct = 100 * hybrid_count / total_quartets if total_quartets > 0 else 0
            unresolved_pct = 100 * unresolved_count / total_quartets if total_quartets > 0 else 0
            print(f"Tree-like: {tree_count} ({tree_pct:.1f}%)")
            print(f"Hybrid: {hybrid_count} ({hybrid_pct:.1f}%)")
            print(f"Unresolved: {unresolved_count} ({unresolved_pct:.1f}%)")
            print("---")
            for quartet, result in quartet_results.items():
                counts = result["counts"]
                cfs = result["cfs"]
                dominant_idx = counts.index(max(counts))
                topo_str = self._format_quartet(quartet, dominant_idx)
                cf_str = "\t".join(f"{cf:.4f}" for cf in cfs)
                classification = result["classification"]
                p_star = result["p_star"]
                p_tree = result["p_tree"]
                print(
                    f"{topo_str}\t{cf_str}\t{classification}\t"
                    f"p_star={p_star:.4f}\tp_tree={p_tree:.4f}"
                )

        if self.plot_output:
            self._draw_quartet_network(all_taxa, quartet_results, self.plot_output)
