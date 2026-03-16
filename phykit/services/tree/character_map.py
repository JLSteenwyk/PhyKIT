"""
Character mapping on a phylogenetic tree.

Performs Fitch parsimony reconstruction on a discrete character matrix,
classifies changes as synapomorphies, convergences, or reversals, and
produces a phylogram/cladogram plot with annotated character changes.

Uses the generalized parsimony utilities in phykit.helpers.parsimony_utils.
"""
import copy
from collections import Counter
from typing import Dict, List, Optional, Set, Tuple

from .base import Tree
from ...helpers.json_output import print_json
from ...helpers.plot_config import PlotConfig
from ...helpers.parsimony_utils import (
    build_parent_map,
    resolve_polytomies,
    fitch_downpass,
    fitch_uppass_acctran,
    fitch_uppass_deltran,
    detect_changes,
    classify_changes,
    consistency_index,
    retention_index,
)
from ...errors import PhykitUserError


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

    def process_args(self, args) -> Dict:
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
        tree = self.read_tree_file()
        tree = copy.deepcopy(tree)

        char_names, tip_states = self._parse_character_matrix(self.data_path)
        n_chars = len(char_names)

        # Compute shared taxa and prune
        tree_tips = set(t.name for t in tree.get_terminals())
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

        tips_to_prune = list(tree_tips - shared)
        if tips_to_prune:
            tree = self.prune_tree_using_taxa_list(tree, tips_to_prune)

        # Filter tip_states to shared taxa
        tip_states = {t: tip_states[t] for t in shared}

        # Resolve polytomies
        resolve_polytomies(tree)

        # Ladderize if requested
        if self.plot_config.ladderize:
            tree.ladderize()

        # Ensure branch lengths
        for clade in tree.find_clades():
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
        wildcard = {"?", "-"}
        n_states_per_char = []
        tip_states_per_char = []
        for i in range(n_chars):
            states_i = [tip_states[t][i] for t in tip_states]
            tip_states_per_char.append(states_i)
            distinct = set(s for s in states_i if s not in wildcard)
            n_states_per_char.append(len(distinct))

        # Observed changes per character from scores (downpass)
        # But we should count from actual branch changes for consistency
        observed_per_char = [0] * n_chars
        for cid, changes in classified.items():
            for char_idx, old, new, cls_type in changes:
                observed_per_char[char_idx] += 1

        ci_per_char, ci_overall = consistency_index(n_states_per_char, observed_per_char)
        ri_per_char, ri_overall = retention_index(tip_states_per_char, observed_per_char)

        # Tree length = sum of all parsimony steps
        tree_length = sum(observed_per_char)

        # Count parsimony-informative characters
        # A character is parsimony-informative if at least 2 states each occur >= 2 times
        n_informative = 0
        for i in range(n_chars):
            states_i = [s for s in tip_states_per_char[i] if s not in wildcard]
            counts = Counter(states_i)
            states_gte2 = sum(1 for c in counts.values() if c >= 2)
            if states_gte2 >= 2:
                n_informative += 1

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
    def _parse_character_matrix(path: str) -> Tuple[List[str], Dict[str, List[str]]]:
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
                lines = [line.rstrip("\n\r") for line in f if line.strip()]
        except FileNotFoundError:
            raise PhykitUserError(
                [f"{path} corresponds to no such file or directory."],
                code=2,
            )

        if len(lines) < 2:
            raise PhykitUserError(
                ["Character matrix must have a header row and at least one data row."],
                code=2,
            )

        header = lines[0].split("\t")
        char_names = header[1:]  # first column is taxon label
        n_expected = len(char_names)

        tip_states: Dict[str, List[str]] = {}
        for line_num, line in enumerate(lines[1:], start=2):
            fields = line.split("\t")
            taxon = fields[0]
            states = fields[1:]
            if len(states) != n_expected:
                raise PhykitUserError(
                    [
                        f"Row {line_num} (taxon '{taxon}') has {len(states)} values "
                        f"but header has {n_expected} characters."
                    ],
                    code=2,
                )
            tip_states[taxon] = states

        return char_names, tip_states

    @staticmethod
    def _assign_node_labels(tree) -> Dict[int, str]:
        """Assign labels to all nodes: tip names for terminals, node_N for internals."""
        labels = {}
        internal_counter = 0
        for clade in tree.find_clades(order="preorder"):
            if clade.is_terminal():
                labels[id(clade)] = clade.name
            else:
                labels[id(clade)] = f"node_{internal_counter}"
                internal_counter += 1
        return labels

    def _print_summary(
        self,
        n_chars: int,
        n_informative: int,
        tree_length: int,
        ci_overall: Optional[float],
        ri_overall: Optional[float],
    ) -> None:
        try:
            print(f"Optimization: {self.optimization}")
            print(f"Characters: {n_chars}")
            print(f"Parsimony-informative: {n_informative}")
            print(f"Tree length: {tree_length}")
            print(f"CI: {ci_overall:.4f}" if ci_overall is not None else "CI: N/A")
            print(f"RI: {ri_overall:.4f}" if ri_overall is not None else "RI: N/A")
            print(f"Output: {self.output_path}")
        except BrokenPipeError:
            pass

    def _print_verbose(
        self,
        char_names: List[str],
        n_chars: int,
        n_informative: int,
        tree_length: int,
        ci_overall: Optional[float],
        ri_overall: Optional[float],
        ci_per_char: List[Optional[float]],
        ri_per_char: List[Optional[float]],
        observed_per_char: List[int],
        classified: Dict[int, List[Tuple[int, str, str, str]]],
        node_labels: Dict[int, str],
    ) -> None:
        try:
            # Summary header
            print(f"Optimization: {self.optimization}")
            print(f"Characters: {n_chars}")
            print(f"Parsimony-informative: {n_informative}")
            print(f"Tree length: {tree_length}")
            print(f"CI: {ci_overall:.4f}" if ci_overall is not None else "CI: N/A")
            print(f"RI: {ri_overall:.4f}" if ri_overall is not None else "RI: N/A")
            print(f"Output: {self.output_path}")
            print()

            # Per-character details
            for i, name in enumerate(char_names):
                ci_str = f"{ci_per_char[i]:.4f}" if ci_per_char[i] is not None else "N/A"
                ri_str = f"{ri_per_char[i]:.4f}" if ri_per_char[i] is not None else "N/A"
                print(f"Character {i} ({name}): steps={observed_per_char[i]}, CI={ci_str}, RI={ri_str}")

                # List changes on branches for this character
                for cid, changes in classified.items():
                    for char_idx, old, new, cls_type in changes:
                        if char_idx == i:
                            label = node_labels.get(cid, f"id_{cid}")
                            print(f"  {label}: {old}->{new} ({cls_type})")
            print()
        except BrokenPipeError:
            pass

    def _print_json(
        self,
        char_names: List[str],
        n_chars: int,
        n_informative: int,
        tree_length: int,
        ci_overall: Optional[float],
        ri_overall: Optional[float],
        ci_per_char: List[Optional[float]],
        ri_per_char: List[Optional[float]],
        observed_per_char: List[int],
        classified: Dict[int, List[Tuple[int, str, str, str]]],
        node_labels: Dict[int, str],
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
        classified: Dict[int, List[Tuple[int, str, str, str]]],
        node_labels: Dict[int, str],
        parent_map: Dict[int, object],
    ) -> None:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            from matplotlib.patches import Patch
        except ImportError:
            print("matplotlib is required for character_map. Install matplotlib and retry.")
            raise SystemExit(2)

        import numpy as np

        config = self.plot_config
        tips = list(tree.get_terminals())
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

        # Compute node positions
        node_x: Dict[int, float] = {}
        node_y: Dict[int, float] = {}

        # Assign y-positions: tips get integer indices
        for i, tip in enumerate(tips):
            node_y[id(tip)] = float(i)

        root = tree.root

        if self.phylogram:
            # Phylogram mode: node_x = parent_x + branch_length
            for clade in tree.find_clades(order="preorder"):
                if clade == root:
                    node_x[id(clade)] = 0.0
                elif id(clade) in parent_map:
                    parent = parent_map[id(clade)]
                    bl = clade.branch_length if clade.branch_length else 0.0
                    node_x[id(clade)] = node_x.get(id(parent), 0.0) + bl
        else:
            # Cladogram mode: all tips at max_x, internal at depth * step_size
            # Compute depth (number of edges from root) for each node
            node_depth: Dict[int, int] = {}
            for clade in tree.find_clades(order="preorder"):
                if clade == root:
                    node_depth[id(clade)] = 0
                elif id(clade) in parent_map:
                    parent = parent_map[id(clade)]
                    node_depth[id(clade)] = node_depth.get(id(parent), 0) + 1

            max_depth = max(node_depth.values()) if node_depth else 1
            step_size = 1.0 / max(max_depth, 1)

            for clade in tree.find_clades(order="preorder"):
                cid = id(clade)
                if clade.is_terminal():
                    node_x[cid] = float(max_depth) * step_size
                else:
                    node_x[cid] = float(node_depth.get(cid, 0)) * step_size

        # Assign y-positions for internal nodes (average of children)
        for clade in tree.find_clades(order="postorder"):
            if not clade.is_terminal() and id(clade) not in node_y:
                child_ys = [
                    node_y[id(c)] for c in clade.clades if id(c) in node_y
                ]
                if child_ys:
                    node_y[id(clade)] = np.mean(child_ys)
                else:
                    node_y[id(clade)] = 0.0

        fig, ax = plt.subplots(figsize=(config.fig_width, config.fig_height))

        # Draw branches (horizontal + vertical connectors)
        for clade in tree.find_clades(order="preorder"):
            if clade == root:
                continue
            cid = id(clade)
            if cid not in parent_map:
                continue
            parent = parent_map[cid]
            pid = id(parent)
            if pid not in node_x or cid not in node_x:
                continue

            x0 = node_x[pid]
            x1 = node_x[cid]
            y0 = node_y.get(pid, 0)
            y1 = node_y.get(cid, 0)

            # Horizontal branch
            ax.plot([x0, x1], [y1, y1], color="black", lw=1.5)
            # Vertical connector
            ax.plot([x0, x0], [y0, y1], color="black", lw=1.5)

        # Tip labels
        max_x = max(node_x.values()) if node_x else 1.0
        offset = max_x * 0.03
        label_fontsize = config.ylabel_fontsize if config.ylabel_fontsize and config.ylabel_fontsize > 0 else 9
        for tip in tips:
            ax.text(
                node_x[id(tip)] + offset, node_y[id(tip)],
                tip.name, va="center", fontsize=label_fontsize,
            )

        # Character change circles on branches
        # Use scatter (marker size in points²) so circles stay round
        # regardless of axis aspect ratio.
        marker_size = max(15, min(80, 600 / max(n_tips, 1)))
        change_fontsize = max(3.0, min(6.0, 7.0 - n_tips * 0.03))

        # Filter characters if requested
        filter_set = set(self.characters_filter) if self.characters_filter else None

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

                ax.scatter(
                    cx, cy, s=marker_size, c=color,
                    edgecolors="black", linewidths=0.5, zorder=5,
                )

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

        # Legend
        legend_handles = [
            Patch(facecolor=colors[0], edgecolor="black", linewidth=0.5,
                  label="Synapomorphy"),
            Patch(facecolor=colors[1], edgecolor="black", linewidth=0.5,
                  label="Convergence"),
            Patch(facecolor=colors[2], edgecolor="black", linewidth=0.5,
                  label="Reversal"),
        ]
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
