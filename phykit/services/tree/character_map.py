"""
Character mapping on a phylogenetic tree.

Performs Fitch parsimony reconstruction on a discrete character matrix,
classifies changes as synapomorphies, convergences, or reversals, and
produces a phylogram/cladogram plot with annotated character changes.

Uses the generalized parsimony utilities in phykit.helpers.parsimony_utils.
"""
from __future__ import annotations

from collections import Counter, defaultdict

from .base import Tree
from ...errors import PhykitUserError


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


_PARSIMONY_UTILS = None


def _parsimony_utils():
    global _PARSIMONY_UTILS
    if _PARSIMONY_UTILS is not None:
        return _PARSIMONY_UTILS

    from ...helpers import parsimony_utils

    _PARSIMONY_UTILS = parsimony_utils
    return parsimony_utils


def build_parent_map(*args, **kwargs):
    return _parsimony_utils().build_parent_map(*args, **kwargs)


def resolve_polytomies(*args, **kwargs):
    return _parsimony_utils().resolve_polytomies(*args, **kwargs)


def fitch_downpass(*args, **kwargs):
    return _parsimony_utils().fitch_downpass(*args, **kwargs)


def fitch_uppass_acctran(*args, **kwargs):
    return _parsimony_utils().fitch_uppass_acctran(*args, **kwargs)


def fitch_uppass_deltran(*args, **kwargs):
    return _parsimony_utils().fitch_uppass_deltran(*args, **kwargs)


def detect_changes(*args, **kwargs):
    return _parsimony_utils().detect_changes(*args, **kwargs)


def classify_changes(*args, **kwargs):
    return _parsimony_utils().classify_changes(*args, **kwargs)


def consistency_index(*args, **kwargs):
    return _parsimony_utils().consistency_index(*args, **kwargs)


class CharacterMap(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"])
        self.data_path = parsed["data_path"]
        self.output_path = parsed["output_path"]
        self.optimization = parsed["optimization"]
        self.phylogram = parsed["phylogram"]
        self.characters_filter = parsed["characters_filter"]
        self.verbose = parsed["verbose"]
        self.json_output = parsed["json_output"]
        self.plot_config = parsed["plot_config"]

    def process_args(self, args) -> dict:
        from ...helpers.plot_config import PlotConfig

        # Parse characters filter: "0,1,3" -> [0, 1, 3]
        chars_str = getattr(args, "characters", None)
        characters_filter = None
        if chars_str is not None:
            characters_filter = [int(c.strip()) for c in chars_str.split(",")]

        return dict(
            tree_file_path=args.tree,
            data_path=args.data,
            output_path=args.output,
            optimization=getattr(args, "optimization", "acctran"),
            phylogram=getattr(args, "phylogram", False),
            characters_filter=characters_filter,
            verbose=getattr(args, "verbose", False),
            json_output=getattr(args, "json", False),
            plot_config=PlotConfig.from_args(args),
        )

    def run(self) -> None:
        tree = self.read_tree_file_unmodified()

        char_names, tip_states = self._parse_character_matrix(self.data_path)
        n_chars = len(char_names)

        # Compute shared taxa and prune
        tree_tip_names = self.get_tip_names_from_tree(tree)
        tree_tips = set(tree_tip_names)
        matrix_taxa = set(tip_states.keys())
        shared = tree_tips & matrix_taxa
        if len(shared) < 3:
            raise PhykitUserError(
                [
                    f"Only {len(shared)} shared taxa between tree and character matrix.",
                    "At least 3 shared taxa are required.",
                ],
                code=2,
            )

        preorder_clades = list(self._iter_preorder(tree.root))
        tips_to_prune = [tip for tip in tree_tip_names if tip not in shared]
        needs_polytomy_resolution = any(
            len(clade.clades) > 2 for clade in preorder_clades
        )
        needs_branch_length_fill = any(
            clade.branch_length is None for clade in preorder_clades
        )
        needs_working_copy = bool(tips_to_prune) or (
            needs_polytomy_resolution
            or needs_branch_length_fill
            or self.plot_config.ladderize
        )
        if needs_working_copy:
            tree = self._fast_copy(tree)
        if tips_to_prune:
            tree = self.prune_tree_using_taxa_list(tree, tips_to_prune)

        # Filter tip_states to shared taxa
        tip_states = {t: tip_states[t] for t in shared}

        # Resolve polytomies
        if needs_polytomy_resolution:
            resolve_polytomies(tree)

        # Ladderize if requested
        if self.plot_config.ladderize:
            tree.ladderize()

        # Ensure branch lengths
        if needs_branch_length_fill:
            for clade in self._iter_preorder(tree.root):
                if clade.branch_length is None:
                    clade.branch_length = 1e-8 if clade != tree.root else 0.0

        # Build parent map
        parent_map = build_parent_map(tree)

        # Fitch downpass
        node_state_sets, scores = fitch_downpass(tree, tip_states)

        # Uppass: ACCTRAN or DELTRAN
        if self.optimization == "deltran":
            node_states = fitch_uppass_deltran(tree, node_state_sets, parent_map)
        else:
            node_states = fitch_uppass_acctran(tree, node_state_sets, parent_map)

        # Detect and classify changes
        branch_changes = detect_changes(tree, node_states, parent_map)
        classified = classify_changes(tree, branch_changes, node_states, parent_map)

        # Compute CI and RI
        # n_states per character: count distinct non-wildcard states in tips
        n_states_per_char, tip_states_per_char, n_informative, counts_per_char = (
            self._summarize_character_states_with_counts(
                tip_states,
                include_states=False,
            )
        )

        # Observed changes per character from scores (downpass)
        # But we should count from actual branch changes for consistency
        observed_per_char = [0] * n_chars
        for cid, changes in classified.items():
            for char_idx, old, new, cls_type in changes:
                observed_per_char[char_idx] += 1

        ci_per_char, ci_overall = consistency_index(n_states_per_char, observed_per_char)
        ri_per_char, ri_overall = self._retention_index_from_counts(
            counts_per_char,
            observed_per_char,
        )

        # Tree length = sum of all parsimony steps
        tree_length = sum(observed_per_char)

        # Assign node labels for output
        node_labels = self._assign_node_labels(tree)

        # Plot
        self._plot_character_map(tree, classified, node_labels, parent_map)

        # Output
        if self.json_output:
            self._print_json(
                char_names, n_chars, n_informative, tree_length,
                ci_overall, ri_overall, ci_per_char, ri_per_char,
                observed_per_char, classified, node_labels,
            )
        elif self.verbose:
            self._print_verbose(
                char_names, n_chars, n_informative, tree_length,
                ci_overall, ri_overall, ci_per_char, ri_per_char,
                observed_per_char, classified, node_labels,
            )
        else:
            self._print_summary(
                n_chars, n_informative, tree_length,
                ci_overall, ri_overall,
            )

    @staticmethod
    def _summarize_character_states(
        tip_states: dict[str, list[str]],
    ) -> tuple[list[int], list[list[str]], int]:
        n_states_per_char, tip_states_per_char, n_informative, _ = (
            CharacterMap._summarize_character_states_with_counts(tip_states)
        )
        return n_states_per_char, tip_states_per_char, n_informative

    @staticmethod
    def _summarize_character_states_with_counts(
        tip_states: dict[str, list[str]],
        include_states: bool = True,
    ) -> tuple[list[int], list[list[str]], int, list[Counter]]:
        if not include_states:
            counts_only = CharacterMap._summarize_character_counts_only(tip_states)
            if counts_only is not None:
                return counts_only
        else:
            ascii_summary = CharacterMap._summarize_character_ascii(tip_states)
            if ascii_summary is not None:
                return ascii_summary

        n_states_per_char = []
        tip_states_per_char = []
        counts_per_char = []
        n_informative = 0

        for states_column in zip(*tip_states.values()):
            states_i = states_column
            if include_states:
                states_i = list(states_column)
                tip_states_per_char.append(states_i)
            counts = Counter(states_i)
            counts.pop("?", None)
            counts.pop("-", None)
            counts_per_char.append(counts)
            n_states_per_char.append(len(counts))
            repeated_states = 0
            for count in counts.values():
                if count >= 2:
                    repeated_states += 1
                    if repeated_states >= 2:
                        n_informative += 1
                        break

        return n_states_per_char, tip_states_per_char, n_informative, counts_per_char

    @staticmethod
    def _summarize_character_counts_only(
        tip_states: dict[str, list[str]],
    ) -> tuple[list[int], list[list[str]], int, list[Counter]] | None:
        ascii_summary = CharacterMap._summarize_character_ascii(
            tip_states,
            include_states=False,
        )
        if ascii_summary is not None:
            return ascii_summary

        return None

    @staticmethod
    def _summarize_character_ascii(
        tip_states: dict[str, list[str]],
        include_states: bool = True,
    ) -> tuple[list[int], list[list[str]], int, list[Counter]] | None:
        import numpy as np

        if not include_states:
            return CharacterMap._summarize_character_ascii_counts_only(tip_states)

        values = list(tip_states.values())
        if not values:
            return [], [], 0, []

        n_taxa = len(values)
        n_chars = len(values[0])
        if any(len(row) != n_chars for row in values):
            return None

        try:
            data = "".join("".join(row) for row in values).encode("ascii")
        except UnicodeEncodeError:
            return None
        if len(data) != n_taxa * n_chars:
            return None

        matrix = np.frombuffer(data, dtype=np.uint8).reshape(n_taxa, n_chars)
        states_by_char = (
            [list(bytes(row).decode("ascii")) for row in matrix.T]
            if include_states
            else []
        )
        return CharacterMap._summarize_ascii_matrix(matrix, states_by_char)

    @staticmethod
    def _summarize_character_ascii_counts_only(
        tip_states: dict[str, list[str]],
    ) -> tuple[list[int], list[list[str]], int, list[Counter]] | None:
        import numpy as np

        iterator = iter(tip_states.values())
        try:
            first = next(iterator)
        except StopIteration:
            return [], [], 0, []

        n_chars = len(first)
        chunks = []
        try:
            first_data = "".join(first).encode("ascii")
            if len(first_data) != n_chars:
                return None
            chunks.append(first_data)
            n_taxa = 1
            for row in iterator:
                if len(row) != n_chars:
                    return None
                row_data = "".join(row).encode("ascii")
                if len(row_data) != n_chars:
                    return None
                chunks.append(row_data)
                n_taxa += 1
        except UnicodeEncodeError:
            return None

        data = b"".join(chunks)
        matrix = np.frombuffer(data, dtype=np.uint8).reshape(n_taxa, n_chars)
        return CharacterMap._summarize_ascii_matrix(matrix, [])

    @staticmethod
    def _summarize_ascii_matrix(
        matrix,
        states_by_char: list[list[str]],
    ) -> tuple[list[int], list[list[str]], int, list[Counter]]:
        import numpy as np

        n_chars = matrix.shape[1]
        symbols = np.unique(matrix)
        symbols = symbols[(symbols != ord("?")) & (symbols != ord("-"))]
        if symbols.size == 0:
            return (
                [0] * n_chars,
                states_by_char,
                0,
                [Counter() for _ in range(n_chars)],
            )

        symbol_counts = CharacterMap._ascii_symbol_counts_by_char(
            matrix,
            symbols,
        )
        n_states_per_char = np.count_nonzero(symbol_counts, axis=0).tolist()
        n_informative = int(
            np.count_nonzero(np.count_nonzero(symbol_counts >= 2, axis=0) >= 2)
        )
        symbol_labels = [chr(int(symbol)) for symbol in symbols]
        counts_per_char = [
            Counter(
                {
                    label: int(count)
                    for label, count in zip(symbol_labels, symbol_counts[:, index])
                    if count
                }
            )
            for index in range(n_chars)
        ]
        return n_states_per_char, states_by_char, n_informative, counts_per_char

    @staticmethod
    def _ascii_symbol_counts_by_char(matrix, symbols):
        import numpy as np

        if symbols.size >= 16:
            n_chars = matrix.shape[1]
            max_code = int(matrix.max()) + 1
            encoded = matrix.astype(np.int64)
            encoded += np.arange(n_chars, dtype=np.int64) * max_code
            counts = np.bincount(
                encoded.ravel(),
                minlength=n_chars * max_code,
            ).reshape(n_chars, max_code)
            return counts[:, symbols].T

        return np.vstack(
            [np.count_nonzero(matrix == symbol, axis=0) for symbol in symbols]
        )

    @staticmethod
    def _retention_index_from_counts(
        counts_per_char: list[Counter],
        observed_per_char: list[int],
    ) -> tuple[list[float | None], float | None]:
        ri_per_char: list[float | None] = []
        sum_num = 0
        sum_den = 0

        for counts, observed in zip(counts_per_char, observed_per_char):
            n_taxa = sum(counts.values())
            if n_taxa == 0:
                ri_per_char.append(None)
                continue

            n_states = len(counts)
            f_max = max(counts.values())
            max_changes = n_taxa - f_max
            min_changes = n_states - 1

            if max_changes == min_changes:
                ri_per_char.append(None)
            else:
                ri = (max_changes - observed) / (max_changes - min_changes)
                ri_per_char.append(ri)
                sum_num += max_changes - observed
                sum_den += max_changes - min_changes

        ri_overall = sum_num / sum_den if sum_den > 0 else None
        return ri_per_char, ri_overall

    @staticmethod
    def _parse_character_matrix(path: str) -> tuple[list[str], dict[str, list[str]]]:
        """Parse a TSV character matrix.

        Expected format:
            taxon\\tchar0\\tchar1\\t...
            A\\t0\\t1\\t...
            B\\t1\\t0\\t...

        Returns:
            (char_names, tip_states) where char_names is a list of character
            names from the header, and tip_states maps taxon name to a list
            of state strings.
        """
        try:
            with open(path) as f:
                header = None
                for line in f:
                    if not line.strip():
                        continue
                    header = line.rstrip("\n\r").split("\t")
                    break

                if header is None:
                    raise PhykitUserError(
                        ["Character matrix must have a header row and at least one data row."],
                        code=2,
                    )

                char_names = header[1:]  # first column is taxon label
                n_expected = len(char_names)

                tip_states: dict[str, list[str]] = {}
                row_num = 2
                for line in f:
                    if not line.strip():
                        continue
                    fields = line.rstrip("\n\r").split("\t")
                    taxon = fields[0]
                    states = fields[1:]
                    if len(states) != n_expected:
                        raise PhykitUserError(
                            [
                                f"Row {row_num} (taxon '{taxon}') has {len(states)} values "
                                f"but header has {n_expected} characters."
                            ],
                            code=2,
                        )
                    tip_states[taxon] = states
                    row_num += 1
        except FileNotFoundError:
            raise PhykitUserError(
                [f"{path} corresponds to no such file or directory."],
                code=2,
            )

        if not tip_states:
            raise PhykitUserError(
                ["Character matrix must have a header row and at least one data row."],
                code=2,
            )

        return char_names, tip_states

    @staticmethod
    def _iter_preorder(root):
        stack = [root]
        pop = stack.pop
        append = stack.append
        while stack:
            clade = pop()
            yield clade
            children = clade.clades
            if children:
                append(children[-1])
                if len(children) == 2:
                    append(children[0])
                else:
                    for idx in range(len(children) - 2, -1, -1):
                        append(children[idx])

    @staticmethod
    def _assign_node_labels(tree) -> dict[int, str]:
        """Assign labels to all nodes: tip names for terminals, node_N for internals."""
        labels = {}
        internal_counter = 0
        for clade in CharacterMap._iter_preorder(tree.root):
            if clade.clades:
                labels[id(clade)] = f"node_{internal_counter}"
                internal_counter += 1
            else:
                labels[id(clade)] = clade.name
        return labels

    def _print_summary(
        self,
        n_chars: int,
        n_informative: int,
        tree_length: int,
        ci_overall: float | None,
        ri_overall: float | None,
    ) -> None:
        try:
            print(
                "\n".join(
                    [
                        f"Optimization: {self.optimization}",
                        f"Characters: {n_chars}",
                        f"Parsimony-informative: {n_informative}",
                        f"Tree length: {tree_length}",
                        (
                            f"CI: {ci_overall:.4f}"
                            if ci_overall is not None
                            else "CI: N/A"
                        ),
                        (
                            f"RI: {ri_overall:.4f}"
                            if ri_overall is not None
                            else "RI: N/A"
                        ),
                        f"Output: {self.output_path}",
                    ]
                )
            )
        except BrokenPipeError:
            pass

    def _print_verbose(
        self,
        char_names: list[str],
        n_chars: int,
        n_informative: int,
        tree_length: int,
        ci_overall: float | None,
        ri_overall: float | None,
        ci_per_char: list[float | None],
        ri_per_char: list[float | None],
        observed_per_char: list[int],
        classified: dict[int, list[tuple[int, str, str, str]]],
        node_labels: dict[int, str],
    ) -> None:
        try:
            # Summary header
            lines = [
                f"Optimization: {self.optimization}",
                f"Characters: {n_chars}",
                f"Parsimony-informative: {n_informative}",
                f"Tree length: {tree_length}",
                f"CI: {ci_overall:.4f}" if ci_overall is not None else "CI: N/A",
                f"RI: {ri_overall:.4f}" if ri_overall is not None else "RI: N/A",
                f"Output: {self.output_path}",
                "",
            ]

            changes_by_char = defaultdict(list)
            for cid, changes in classified.items():
                label = node_labels.get(cid, f"id_{cid}")
                for char_idx, old, new, cls_type in changes:
                    changes_by_char[char_idx].append((label, old, new, cls_type))

            # Per-character details
            for i, name in enumerate(char_names):
                ci_str = f"{ci_per_char[i]:.4f}" if ci_per_char[i] is not None else "N/A"
                ri_str = f"{ri_per_char[i]:.4f}" if ri_per_char[i] is not None else "N/A"
                lines.append(
                    f"Character {i} ({name}): steps={observed_per_char[i]}, "
                    f"CI={ci_str}, RI={ri_str}"
                )

                # List changes on branches for this character
                for label, old, new, cls_type in changes_by_char.get(i, []):
                    lines.append(f"  {label}: {old}->{new} ({cls_type})")
            lines.append("")
            print("\n".join(lines))
        except BrokenPipeError:
            pass

    def _print_json(
        self,
        char_names: list[str],
        n_chars: int,
        n_informative: int,
        tree_length: int,
        ci_overall: float | None,
        ri_overall: float | None,
        ci_per_char: list[float | None],
        ri_per_char: list[float | None],
        observed_per_char: list[int],
        classified: dict[int, list[tuple[int, str, str, str]]],
        node_labels: dict[int, str],
    ) -> None:
        characters = []
        for i, name in enumerate(char_names):
            changes = []
            for cid, branch_changes in classified.items():
                for char_idx, old, new, cls_type in branch_changes:
                    if char_idx == i:
                        label = node_labels.get(cid, f"id_{cid}")
                        changes.append({
                            "branch": label,
                            "from": old,
                            "to": new,
                            "type": cls_type,
                        })
            characters.append({
                "index": i,
                "name": name,
                "steps": observed_per_char[i],
                "ci": round(ci_per_char[i], 4) if ci_per_char[i] is not None else None,
                "ri": round(ri_per_char[i], 4) if ri_per_char[i] is not None else None,
                "changes": changes,
            })

        payload = {
            "n_characters": n_chars,
            "n_informative": n_informative,
            "tree_length": tree_length,
            "ci": round(ci_overall, 4) if ci_overall is not None else None,
            "ri": round(ri_overall, 4) if ri_overall is not None else None,
            "optimization": self.optimization,
            "output_file": self.output_path,
            "characters": characters,
        }
        print_json(payload)

    def _plot_character_map(
        self,
        tree,
        classified: dict[int, list[tuple[int, str, str, str]]],
        node_labels: dict[int, str],
        parent_map: dict[int, object],
    ) -> None:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            from matplotlib.collections import LineCollection
            from matplotlib.patches import Patch
        except ImportError:
            print("matplotlib is required for character_map. Install matplotlib and retry.")
            raise SystemExit(2)

        config = self.plot_config
        root = tree.root
        preorder_clades = list(self._iter_preorder(root))
        tips = [clade for clade in preorder_clades if not clade.clades]
        n_tips = len(tips)
        config.resolve(n_rows=n_tips, n_cols=None)

        # Default colors: synapomorphy, convergence, reversal
        default_colors = ["#2b8cbe", "#d62728", "#969696"]
        colors = config.merge_colors(default_colors)
        color_map = {
            "synapomorphy": colors[0],
            "convergence": colors[1],
            "reversal": colors[2],
        }

        from ...helpers.plot_config import compute_node_positions

        node_x, node_y = compute_node_positions(
            tree,
            parent_map,
            cladogram=(not self.phylogram or self.plot_config.cladogram),
            preorder_clades=preorder_clades,
        )

        fig, ax = plt.subplots(figsize=(config.fig_width, config.fig_height))

        if self.plot_config.circular:
            import math
            from ...helpers.circular_layout import (
                circular_branch_points,
                compute_circular_coords,
                draw_circular_branches,
                draw_circular_tip_labels,
                radial_offset,
            )

            # --- Circular mode ---
            coords = compute_circular_coords(
                tree,
                node_x,
                parent_map,
                preorder_clades=preorder_clades,
                terminal_clades=tips,
            )
            ax.set_aspect("equal")
            ax.axis("off")

            draw_circular_branches(ax, tree, coords, parent_map)
            label_fontsize = config.ylabel_fontsize if config.ylabel_fontsize and config.ylabel_fontsize > 0 else 9
            max_x = max(node_x.values()) if node_x else 1.0
            draw_circular_tip_labels(ax, tree, coords, fontsize=label_fontsize, offset=max_x * 0.03)

            # Apply color annotations
            if self.plot_config.color_file:
                from ...helpers.color_annotations import (
                    build_color_legend_handles,
                    apply_label_colors,
                    draw_range_wedge,
                    get_clade_branch_ids,
                    parse_color_file,
                    resolve_mrca,
                )

                color_data = parse_color_file(self.plot_config.color_file)
                for taxa_list, clr, lbl in color_data["ranges"]:
                    mrca = resolve_mrca(tree, taxa_list)
                    if mrca is not None:
                        draw_range_wedge(ax, tree, mrca, clr, coords)
                for taxa_list, clade_color, lbl in color_data["clades"]:
                    mrca = resolve_mrca(tree, taxa_list)
                    if mrca is not None:
                        clade_ids = get_clade_branch_ids(tree, mrca, parent_map)
                        radial_segments = []
                        for cl in preorder_clades:
                            if cl == tree.root:
                                continue
                            if id(cl) in clade_ids and id(cl) in parent_map:
                                parent_coords = coords[id(parent_map[id(cl)])]
                                child_coords = coords[id(cl)]
                                angle = child_coords["angle"]
                                cos_angle = math.cos(angle)
                                sin_angle = math.sin(angle)
                                r_parent = parent_coords["radius"]
                                r_child = child_coords["radius"]
                                radial_segments.append(
                                    [
                                        (r_parent * cos_angle, r_parent * sin_angle),
                                        (r_child * cos_angle, r_child * sin_angle),
                                    ]
                                )
                        if radial_segments:
                            ax.add_collection(
                                LineCollection(
                                    radial_segments,
                                    colors=clade_color,
                                    linewidths=1.5,
                                    capstyle="round",
                                    zorder=2,
                                )
                            )
                apply_label_colors(ax, color_data["labels"])

            # Character change circles on branches
            # Smaller markers in circular mode — radial branches are shorter
            marker_size = max(10, min(50, 300 / max(n_tips, 1)))
            change_fontsize = max(3.0, min(5.0, 6.0 - n_tips * 0.03))

            filter_set = set(self.characters_filter) if self.characters_filter else None
            change_x = []
            change_y = []
            change_colors = []

            for cid, changes in classified.items():
                if cid not in parent_map:
                    continue
                parent_obj = parent_map[cid]
                pid = id(parent_obj)
                if pid not in coords or cid not in coords:
                    continue

                # Filter changes to draw
                changes_to_draw = []
                for char_idx, old, new, cls_type in changes:
                    if filter_set is not None and char_idx not in filter_set:
                        continue
                    changes_to_draw.append((char_idx, old, new, cls_type))

                if not changes_to_draw:
                    continue

                n_changes = len(changes_to_draw)
                # Get evenly spaced points along the radial branch
                branch_pts = circular_branch_points(
                    coords[pid], coords[cid], n_changes + 2
                )
                # Use positions 1..n_changes (skip endpoints)
                for j, (char_idx, old, new, cls_type) in enumerate(changes_to_draw):
                    pt_x, pt_y, pt_angle = branch_pts[j + 1]

                    color = color_map.get(cls_type, "#999999")

                    change_x.append(pt_x)
                    change_y.append(pt_y)
                    change_colors.append(color)

                    # Character index: offset radially outward
                    dx_out, dy_out = radial_offset(pt_angle, 8)
                    ax.annotate(
                        str(char_idx), (pt_x, pt_y),
                        textcoords="offset points", xytext=(dx_out, dy_out),
                        ha="center", va="center",
                        fontsize=change_fontsize, zorder=6,
                    )
                    # State transition: offset radially inward
                    dx_in, dy_in = radial_offset(pt_angle, -8)
                    ax.annotate(
                        f"{old}\u2192{new}", (pt_x, pt_y),
                        textcoords="offset points", xytext=(dx_in, dy_in),
                        ha="center", va="center",
                        fontsize=change_fontsize, zorder=6,
                    )

            if change_x:
                ax.scatter(
                    change_x, change_y, s=marker_size, c=change_colors,
                    edgecolors="black", linewidths=0.5, zorder=5,
                )

            # Legend
            legend_handles = [
                Patch(facecolor=colors[0], edgecolor="black", linewidth=0.5,
                      label="Synapomorphy"),
                Patch(facecolor=colors[1], edgecolor="black", linewidth=0.5,
                      label="Convergence"),
                Patch(facecolor=colors[2], edgecolor="black", linewidth=0.5,
                      label="Reversal"),
            ]
            if self.plot_config.color_file:
                color_legend = build_color_legend_handles(color_data)
                legend_handles.extend(color_legend)
            legend_loc = config.legend_position or "upper right"
            if legend_loc != "none":
                ax.legend(handles=legend_handles, loc=legend_loc, fontsize=8, frameon=True)

            if config.show_title:
                ax.set_title(
                    config.title or "Character Map",
                    fontsize=config.title_fontsize,
                )

            fig.subplots_adjust(left=0.05, right=0.80, top=0.92, bottom=0.08)

            fig.savefig(self.output_path, dpi=config.dpi)
            plt.close(fig)

        else:
            # --- Rectangular mode ---
            from ...helpers.plot_config import draw_tree_branches

            draw_tree_branches(
                ax,
                tree,
                node_x,
                node_y,
                parent_map,
                lw=1.5,
                vertical_lw=1.5,
            )

            # Tip labels
            max_x = max(node_x.values()) if node_x else 1.0
            offset = max_x * 0.03
            label_fontsize = config.ylabel_fontsize if config.ylabel_fontsize and config.ylabel_fontsize > 0 else 9
            for tip in tips:
                ax.text(
                    node_x[id(tip)] + offset, node_y[id(tip)],
                    tip.name, va="center", fontsize=label_fontsize,
                )

            # Apply color annotations
            if self.plot_config.color_file:
                from ...helpers.color_annotations import (
                    build_color_legend_handles,
                    apply_label_colors,
                    draw_range_rect,
                    get_clade_branch_ids,
                    parse_color_file,
                    resolve_mrca,
                )

                color_data = parse_color_file(self.plot_config.color_file)
                for taxa_list, clr, lbl in color_data["ranges"]:
                    mrca = resolve_mrca(tree, taxa_list)
                    if mrca is not None:
                        draw_range_rect(ax, tree, mrca, clr, node_x, node_y)
                for taxa_list, clade_color, lbl in color_data["clades"]:
                    mrca = resolve_mrca(tree, taxa_list)
                    if mrca is not None:
                        clade_ids = get_clade_branch_ids(tree, mrca, parent_map)
                        horizontal_segments = []
                        vertical_segments = []
                        for cl in preorder_clades:
                            if cl == tree.root:
                                continue
                            if id(cl) in clade_ids and id(cl) in parent_map:
                                pid_val = id(parent_map[id(cl)])
                                cid_val = id(cl)
                                x0, x1 = node_x[pid_val], node_x[cid_val]
                                y0 = node_y.get(pid_val, 0)
                                y1 = node_y.get(cid_val, 0)
                                horizontal_segments.append([(x0, y1), (x1, y1)])
                                vertical_segments.append([(x0, y0), (x0, y1)])
                        if vertical_segments:
                            ax.add_collection(
                                LineCollection(
                                    vertical_segments,
                                    colors=clade_color,
                                    linewidths=1.5,
                                    zorder=2,
                                )
                            )
                        if horizontal_segments:
                            ax.add_collection(
                                LineCollection(
                                    horizontal_segments,
                                    colors=clade_color,
                                    linewidths=1.5,
                                    zorder=2,
                                )
                            )
                apply_label_colors(ax, color_data["labels"])

            # Character change circles on branches
            # Use scatter (marker size in points²) so circles stay round
            # regardless of axis aspect ratio.
            marker_size = max(15, min(80, 600 / max(n_tips, 1)))
            change_fontsize = max(3.0, min(6.0, 7.0 - n_tips * 0.03))

            # Filter characters if requested
            filter_set = set(self.characters_filter) if self.characters_filter else None
            change_x = []
            change_y = []
            change_colors = []

            for cid, changes in classified.items():
                if cid not in parent_map:
                    continue
                parent = parent_map[cid]
                pid = id(parent)
                if pid not in node_x or cid not in node_x:
                    continue

                # Filter changes to draw
                changes_to_draw = []
                for char_idx, old, new, cls_type in changes:
                    if filter_set is not None and char_idx not in filter_set:
                        continue
                    changes_to_draw.append((char_idx, old, new, cls_type))

                if not changes_to_draw:
                    continue

                x0 = node_x[pid]
                x1 = node_x[cid]
                y_branch = node_y.get(cid, 0)

                n_changes = len(changes_to_draw)
                for j, (char_idx, old, new, cls_type) in enumerate(changes_to_draw):
                    # Position along branch
                    frac = (j + 1) / (n_changes + 1)
                    cx = x0 + (x1 - x0) * frac
                    cy = y_branch

                    color = color_map.get(cls_type, "#999999")

                    change_x.append(cx)
                    change_y.append(cy)
                    change_colors.append(color)

                    # Character index above (offset in points via annotate)
                    ax.annotate(
                        str(char_idx), (cx, cy),
                        textcoords="offset points", xytext=(0, 6),
                        ha="center", va="bottom",
                        fontsize=change_fontsize, zorder=6,
                    )
                    # State transition below
                    ax.annotate(
                        f"{old}\u2192{new}", (cx, cy),
                        textcoords="offset points", xytext=(0, -6),
                        ha="center", va="top",
                        fontsize=change_fontsize, zorder=6,
                    )

            if change_x:
                ax.scatter(
                    change_x, change_y, s=marker_size, c=change_colors,
                    edgecolors="black", linewidths=0.5, zorder=5,
                )

            # Legend
            legend_handles = [
                Patch(facecolor=colors[0], edgecolor="black", linewidth=0.5,
                      label="Synapomorphy"),
                Patch(facecolor=colors[1], edgecolor="black", linewidth=0.5,
                      label="Convergence"),
                Patch(facecolor=colors[2], edgecolor="black", linewidth=0.5,
                      label="Reversal"),
            ]
            if self.plot_config.color_file:
                color_legend = build_color_legend_handles(color_data)
                legend_handles.extend(color_legend)
            legend_loc = config.legend_position or "upper right"
            if legend_loc != "none":
                ax.legend(handles=legend_handles, loc=legend_loc, fontsize=8, frameon=True)

            # Axes formatting
            ax.set_yticks([])
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.spines["left"].set_visible(False)

            if self.phylogram:
                ax.set_xlabel("Branch length (subs/site)")
            else:
                ax.set_xlabel("")
                ax.set_xticks([])
                ax.spines["bottom"].set_visible(False)

            if config.show_title:
                ax.set_title(
                    config.title or "Character Map",
                    fontsize=config.title_fontsize,
                )
            if config.axis_fontsize:
                ax.xaxis.label.set_fontsize(config.axis_fontsize)

            # Adjust margins to leave room for labels
            fig.subplots_adjust(left=0.05, right=0.80, top=0.92, bottom=0.08)

            fig.savefig(self.output_path, dpi=config.dpi)
            plt.close(fig)
