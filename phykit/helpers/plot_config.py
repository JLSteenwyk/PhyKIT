import sys


VALID_LEGEND_POSITIONS = frozenset([
    "best", "upper right", "upper left", "lower left", "lower right",
    "right", "center left", "center right", "lower center", "upper center",
    "center", "none",
])


class PlotConfig:
    def __init__(
        self,
        fig_width: float | None = None,
        fig_height: float | None = None,
        dpi: int = 300,
        show_title: bool = True,
        title: str | None = None,
        legend_position: str | None = None,
        ylabel_fontsize: float | None = None,
        xlabel_fontsize: float | None = None,
        title_fontsize: float | None = None,
        axis_fontsize: float | None = None,
        colors: list[str] | None = None,
        ladderize: bool = False,
        cladogram: bool = False,
        circular: bool = False,
        color_file: str | None = None,
    ):
        self.fig_width = fig_width
        self.fig_height = fig_height
        self.dpi = dpi
        self.show_title = show_title
        self.title = title
        self.legend_position = legend_position
        self.ylabel_fontsize = ylabel_fontsize
        self.xlabel_fontsize = xlabel_fontsize
        self.title_fontsize = title_fontsize
        self.axis_fontsize = axis_fontsize
        self.colors = colors
        self.ladderize = ladderize
        self.cladogram = cladogram
        self.circular = circular
        self.color_file = color_file

    def validate(self) -> None:
        if self.fig_width is not None and self.fig_width <= 0:
            print(f"Error: --fig-width must be positive, got {self.fig_width}")
            sys.exit(2)
        if self.fig_height is not None and self.fig_height <= 0:
            print(f"Error: --fig-height must be positive, got {self.fig_height}")
            sys.exit(2)
        if self.dpi <= 0:
            print(f"Error: --dpi must be positive, got {self.dpi}")
            sys.exit(2)
        for name, val in [
            ("--ylabel-fontsize", self.ylabel_fontsize),
            ("--xlabel-fontsize", self.xlabel_fontsize),
            ("--title-fontsize", self.title_fontsize),
            ("--axis-fontsize", self.axis_fontsize),
        ]:
            if val is not None and val < 0:
                print(f"Error: {name} must be non-negative, got {val}")
                sys.exit(2)
        if (
            self.legend_position is not None
            and self.legend_position not in VALID_LEGEND_POSITIONS
        ):
            print(
                f"Error: --legend-position '{self.legend_position}' is not valid. "
                f"Choose from: {', '.join(sorted(VALID_LEGEND_POSITIONS))}"
            )
            sys.exit(2)

    @classmethod
    def auto_scale(cls, n_rows=None, n_cols=None) -> "PlotConfig":
        # fig_height
        if n_rows is not None:
            raw_height = 3.0 + n_rows * 0.18
            fig_height = max(5.0, min(200.0, raw_height))
            if raw_height > 200.0:
                print(
                    f"Warning: computed figure height ({raw_height:.0f} in) "
                    f"exceeds 200 in; capping at 200. Use --fig-height to override.",
                    file=sys.stderr,
                )
        else:
            fig_height = 8.0

        # fig_width
        if n_cols is not None:
            fig_width = max(10.0, min(20.0, 8.0 + n_cols * 0.15))
        else:
            fig_width = 14.0

        # ylabel_fontsize
        if n_rows is not None:
            if n_rows > 800:
                ylabel_fontsize = 0.0
            else:
                ylabel_fontsize = max(3.0, min(7.0, 7.0 - (n_rows - 50) * 0.008))
        else:
            ylabel_fontsize = 7.0

        # xlabel_fontsize
        if n_cols is not None:
            if n_cols > 60:
                xlabel_fontsize = 0.0
            else:
                xlabel_fontsize = max(3.0, min(7.0, 7.0 - (n_cols - 20) * 0.1))
        else:
            xlabel_fontsize = 7.0

        return cls(
            fig_width=fig_width,
            fig_height=fig_height,
            ylabel_fontsize=ylabel_fontsize,
            xlabel_fontsize=xlabel_fontsize,
            title_fontsize=12.0,
            axis_fontsize=10.0,
        )

    def resolve(self, n_rows=None, n_cols=None) -> "PlotConfig":
        if self.circular and self.fig_width is None and self.fig_height is None:
            size = max(8.0, min(14.0, 3.0 + (n_rows or 0) * 0.08))
            self.fig_width = size
            self.fig_height = size
        defaults = PlotConfig.auto_scale(n_rows, n_cols)
        for fld in [
            "fig_width", "fig_height", "ylabel_fontsize", "xlabel_fontsize",
            "title_fontsize", "axis_fontsize",
        ]:
            if getattr(self, fld) is None:
                setattr(self, fld, getattr(defaults, fld))
        return self

    def apply_to_figure(self, fig, ax, default_title, default_colors):
        # Title
        if self.show_title:
            title_text = self.title if self.title is not None else default_title
            ax.set_title(title_text, fontsize=self.title_fontsize)
        else:
            ax.set_title("")

        # Legend
        if self.legend_position is not None:
            if self.legend_position == "none":
                legend = ax.get_legend()
                if legend is not None:
                    legend.set_visible(False)
            else:
                legend = ax.get_legend()
                if legend is not None:
                    handles = legend.legend_handles
                    labels = [t.get_text() for t in legend.get_texts()]
                    legend.remove()
                    ax.legend(handles=handles, labels=labels, loc=self.legend_position)

        # Axis label font sizes
        if self.axis_fontsize is not None:
            ax.xaxis.label.set_fontsize(self.axis_fontsize)
            ax.yaxis.label.set_fontsize(self.axis_fontsize)

        # Y-axis tick labels
        if self.ylabel_fontsize is not None:
            if self.ylabel_fontsize == 0.0:
                ax.set_yticklabels([])
            else:
                for label in ax.get_yticklabels():
                    label.set_fontsize(self.ylabel_fontsize)

        # X-axis tick labels
        if self.xlabel_fontsize is not None:
            if self.xlabel_fontsize == 0.0:
                ax.set_xticklabels([])
            else:
                for label in ax.get_xticklabels():
                    label.set_fontsize(self.xlabel_fontsize)

        # Colors
        return self.merge_colors(default_colors)

    def merge_colors(self, defaults: list[str]) -> list[str]:
        if self.colors is None:
            return list(defaults)
        result = list(defaults)
        for i, user_color in enumerate(self.colors[:len(defaults)]):
            if user_color:  # non-empty string overrides
                result[i] = user_color
        return result

    @classmethod
    def from_args(cls, args) -> "PlotConfig":
        colors_str = getattr(args, "colors", None)
        colors = None
        if colors_str is not None:
            colors = [c.strip() for c in colors_str.split(",")]

        no_title = getattr(args, "no_title", False)

        config = cls(
            fig_width=getattr(args, "fig_width", None),
            fig_height=getattr(args, "fig_height", None),
            dpi=getattr(args, "dpi", 300),
            show_title=not no_title,
            title=getattr(args, "title", None),
            legend_position=getattr(args, "legend_position", None),
            ylabel_fontsize=getattr(args, "ylabel_fontsize", None),
            xlabel_fontsize=getattr(args, "xlabel_fontsize", None),
            title_fontsize=getattr(args, "title_fontsize", None),
            axis_fontsize=getattr(args, "axis_fontsize", None),
            colors=colors,
            ladderize=getattr(args, "ladderize", False),
            cladogram=getattr(args, "cladogram", False),
            circular=getattr(args, "circular", False),
            color_file=getattr(args, "color_file", None),
        )
        config.validate()
        return config


def add_plot_arguments(parser) -> None:
    group = parser.add_argument_group("plot options")
    group.add_argument("--fig-width", type=float, default=None, help="Figure width in inches (auto-scaled if omitted)")
    group.add_argument("--fig-height", type=float, default=None, help="Figure height in inches (auto-scaled if omitted)")
    group.add_argument("--dpi", type=int, default=300, help="Resolution in DPI (default: 300)")
    group.add_argument("--no-title", action="store_true", default=False, help="Hide the plot title")
    group.add_argument("--title", type=str, default=None, help="Custom title text")
    group.add_argument("--legend-position", type=str, default=None, help="Legend location (e.g., 'upper right', 'none' to hide)")
    group.add_argument("--ylabel-fontsize", type=float, default=None, help="Font size for y-axis labels; 0 to hide")
    group.add_argument("--xlabel-fontsize", type=float, default=None, help="Font size for x-axis labels; 0 to hide")
    group.add_argument("--title-fontsize", type=float, default=None, help="Font size for the title")
    group.add_argument("--axis-fontsize", type=float, default=None, help="Font size for axis labels")
    group.add_argument("--colors", type=str, default=None, help="Comma-separated colors (hex or named)")
    group.add_argument("--ladderize", action="store_true", default=False, help="Ladderize (sort) the tree before plotting")
    group.add_argument("--cladogram", action="store_true", default=False, help="Draw cladogram (equal branch lengths, tips aligned) instead of phylogram")
    group.add_argument("--circular", action="store_true", default=False, help="Draw circular (radial/fan) phylogram instead of rectangular")
    group.add_argument("--color-file", type=str, default=None, help="Color annotation file for tip labels, clade ranges, and branch colors")


def compute_node_x_cladogram(tree, parent_map):
    """Compute cladogram x-coordinates: tips aligned at right, internal nodes at depth.

    Returns dict mapping node id -> x-coordinate (0.0 to 1.0 range).
    """
    root = tree.root
    node_depth = {}
    for clade in tree.find_clades(order="preorder"):
        if clade == root:
            node_depth[id(clade)] = 0
        elif id(clade) in parent_map:
            parent = parent_map[id(clade)]
            node_depth[id(clade)] = node_depth.get(id(parent), 0) + 1

    max_depth = max(node_depth.values()) if node_depth else 1
    step_size = 1.0 / max(max_depth, 1)

    node_x = {}
    for clade in tree.find_clades(order="preorder"):
        cid = id(clade)
        if clade.is_terminal():
            node_x[cid] = float(max_depth) * step_size
        else:
            node_x[cid] = float(node_depth.get(cid, 0)) * step_size
    return node_x


# ---- Shared rectangular tree plotting utilities ----


def build_parent_map(tree):
    """Build a dict mapping child node id -> parent node."""
    try:
        root = tree.root
        root.clades
    except AttributeError:
        return _build_parent_map_legacy(tree)

    parent_map = {}
    stack = [root]
    pop = stack.pop
    extend = stack.extend
    while stack:
        clade = pop()
        children = getattr(clade, "clades", None)
        if not isinstance(children, list):
            return _build_parent_map_legacy(tree)
        for child in children:
            parent_map[id(child)] = clade
        if children:
            extend(children)
    return parent_map


def _build_parent_map_legacy(tree):
    parent_map = {}
    for clade in tree.find_clades(order="preorder"):
        for child in clade.clades:
            parent_map[id(child)] = clade
    return parent_map


def _mean_child_y(children, node_y):
    total = 0.0
    count = 0
    for child in children:
        child_id = id(child)
        if child_id in node_y:
            total += node_y[child_id]
            count += 1
    return float(total / count) if count else 0.0


def _assign_internal_y_from_children(clade, node_y):
    children = clade.clades
    if len(children) == 2:
        first_y = node_y.get(id(children[0]))
        second_y = node_y.get(id(children[1]))
        if first_y is not None and second_y is not None:
            return (first_y + second_y) * 0.5
    return _mean_child_y(children, node_y)


def compute_node_positions(tree, parent_map, cladogram=False, preorder_clades=None):
    """Compute (node_x, node_y) for a rectangular tree layout.

    Parameters
    ----------
    tree : Bio.Phylo tree
    parent_map : dict from build_parent_map()
    cladogram : if True, use equal-depth x-positions (tips aligned)
    preorder_clades : optional precomputed preorder clade list

    Returns
    -------
    (node_x, node_y) : dicts mapping node id -> coordinate
    """
    if preorder_clades is not None:
        return _compute_node_positions_from_preorder(
            tree, parent_map, cladogram, preorder_clades
        )

    try:
        root = tree.root
        root.clades
    except AttributeError:
        return _compute_node_positions_legacy(tree, parent_map, cladogram=cladogram)

    node_x = {}
    node_y = {}
    node_depth = {} if cladogram else None
    preorder = []

    stack = [(root, 0.0, 0)]
    while stack:
        clade, x_value, depth = stack.pop()
        children = getattr(clade, "clades", None)
        if not isinstance(children, list):
            return _compute_node_positions_legacy(
                tree, parent_map, cladogram=cladogram
            )

        preorder.append(clade)
        cid = id(clade)
        if cladogram:
            node_depth[cid] = depth
        else:
            node_x[cid] = x_value

        for child in reversed(children):
            branch_length = child.branch_length if child.branch_length else 0.0
            stack.append((child, x_value + branch_length, depth + 1))

    terminal_count = 0
    for clade in preorder:
        if not clade.clades:
            node_y[id(clade)] = terminal_count
            terminal_count += 1

    if cladogram:
        max_depth = max(node_depth.values()) if node_depth else 1
        step_size = 1.0 / max(max_depth, 1)
        for clade in reversed(preorder):
            cid = id(clade)
            if clade.clades:
                node_x[cid] = float(node_depth.get(cid, 0)) * step_size
            else:
                node_x[cid] = float(max_depth) * step_size

    for clade in reversed(preorder):
        children = clade.clades
        cid = id(clade)
        if children and cid not in node_y:
            node_y[cid] = _assign_internal_y_from_children(clade, node_y)

    return node_x, node_y


def _compute_node_positions_from_preorder(
    tree, parent_map, cladogram, preorder_clades
):
    node_x = {}
    node_y = {}
    root = tree.root

    if cladogram:
        node_depth = {id(root): 0}
        max_depth = 0
        for clade in preorder_clades:
            cid = id(clade)
            if clade is not root:
                parent = parent_map.get(cid)
                if parent is None:
                    continue
                depth = node_depth.get(id(parent), 0) + 1
                node_depth[cid] = depth
                if depth > max_depth:
                    max_depth = depth
            elif cid not in node_depth:
                node_depth[cid] = 0

        step_size = 1.0 / max(max_depth, 1)
        for clade in preorder_clades:
            cid = id(clade)
            if clade.clades:
                node_x[cid] = float(node_depth.get(cid, 0)) * step_size
            else:
                node_x[cid] = float(max_depth) * step_size
    else:
        node_x[id(root)] = 0.0
        for clade in preorder_clades:
            if clade is root:
                continue
            cid = id(clade)
            parent = parent_map.get(cid)
            if parent is None:
                continue
            branch_length = clade.branch_length if clade.branch_length else 0.0
            node_x[cid] = node_x.get(id(parent), 0.0) + branch_length

    terminal_count = 0
    for clade in preorder_clades:
        if not clade.clades:
            node_y[id(clade)] = terminal_count
            terminal_count += 1

    for clade in reversed(preorder_clades):
        children = clade.clades
        cid = id(clade)
        if children and cid not in node_y:
            node_y[cid] = _assign_internal_y_from_children(clade, node_y)

    return node_x, node_y


def _compute_node_positions_legacy(tree, parent_map, cladogram=False):
    tips = list(tree.get_terminals())
    root = tree.root

    node_y = {}
    for i, tip in enumerate(tips):
        node_y[id(tip)] = i

    if cladogram:
        node_x = compute_node_x_cladogram(tree, parent_map)
    else:
        node_x = {}
        for clade in tree.find_clades(order="preorder"):
            if clade == root:
                node_x[id(clade)] = 0.0
            elif id(clade) in parent_map:
                parent = parent_map[id(clade)]
                t = clade.branch_length if clade.branch_length else 0.0
                node_x[id(clade)] = node_x.get(id(parent), 0.0) + t

    for clade in tree.find_clades(order="postorder"):
        if not clade.is_terminal() and id(clade) not in node_y:
            node_y[id(clade)] = _mean_child_y(clade.clades, node_y)

    return node_x, node_y


def _terminal_clades(tree):
    direct_terminals = _terminal_clades_direct(tree)
    if direct_terminals is not None:
        return direct_terminals
    return list(tree.get_terminals())


def _terminal_clades_direct(tree):
    try:
        root = tree.root
        root.clades
    except AttributeError:
        return None

    terminals = []
    stack = [root]
    pop = stack.pop
    append = stack.append
    append_terminal = terminals.append
    try:
        while stack:
            clade = pop()
            children = clade.clades
            if children:
                child_count = len(children)
                if child_count == 2:
                    append(children[1])
                    append(children[0])
                else:
                    for index in range(child_count - 1, -1, -1):
                        append(children[index])
            else:
                append_terminal(clade)
    except AttributeError:
        return None
    return terminals


def _preorder_clades(tree):
    direct_clades = _preorder_clades_direct(tree)
    if direct_clades is not None:
        return direct_clades
    return list(tree.find_clades(order="preorder"))


def _preorder_clades_direct(tree):
    try:
        root = tree.root
        root.clades
    except AttributeError:
        return None

    clades = []
    stack = [root]
    pop = stack.pop
    append = stack.append
    append_clade = clades.append
    try:
        while stack:
            clade = pop()
            append_clade(clade)
            children = clade.clades
            if children:
                child_count = len(children)
                if child_count == 2:
                    append(children[1])
                    append(children[0])
                else:
                    for index in range(child_count - 1, -1, -1):
                        append(children[index])
    except AttributeError:
        return None
    return clades


def draw_tree_branches(
    ax, tree, node_x, node_y, parent_map,
    color="black", lw=1.5, vertical_color="black", vertical_lw=0.8,
):
    """Draw rectangular tree branches (horizontal + vertical connectors).

    Override color per branch by passing a callable for `color`:
        color=lambda clade: "red" if ... else "black"
    """
    if hasattr(ax, "add_collection"):
        _draw_tree_branches_collections(
            ax,
            tree,
            node_x,
            node_y,
            parent_map,
            color=color,
            lw=lw,
            vertical_color=vertical_color,
            vertical_lw=vertical_lw,
        )
        return

    root = tree.root
    parent_get = parent_map.get
    node_x_get = node_x.get
    node_y_get = node_y.get
    plot = ax.plot
    color_is_callable = callable(color)
    missing = object()
    for clade in _preorder_clades(tree):
        if clade == root:
            continue
        cid = id(clade)
        parent = parent_get(cid, missing)
        if parent is missing:
            continue
        pid = id(parent)
        x0 = node_x_get(pid, missing)
        x1 = node_x_get(cid, missing)
        if x0 is missing or x1 is missing:
            continue

        y0 = node_y_get(pid, 0)
        y1 = node_y_get(cid, 0)

        branch_color = color(clade) if color_is_callable else color
        plot([x0, x1], [y1, y1], color=branch_color, lw=lw)
        plot([x0, x0], [y0, y1], color=vertical_color, lw=vertical_lw)


def _draw_tree_branches_collections(
    ax,
    tree,
    node_x,
    node_y,
    parent_map,
    color="black",
    lw=1.5,
    vertical_color="black",
    vertical_lw=0.8,
):
    from matplotlib.collections import LineCollection

    root = tree.root
    horizontal_segments = []
    color_is_callable = callable(color)
    horizontal_colors = [] if color_is_callable else color
    vertical_segments = []
    append_horizontal = horizontal_segments.append
    append_vertical = vertical_segments.append
    append_horizontal_color = (
        horizontal_colors.append if color_is_callable else None
    )
    parent_get = parent_map.get
    node_x_get = node_x.get
    node_y_get = node_y.get
    missing = object()

    for clade in _preorder_clades(tree):
        if clade == root:
            continue
        cid = id(clade)
        parent = parent_get(cid, missing)
        if parent is missing:
            continue
        pid = id(parent)
        x0 = node_x_get(pid, missing)
        x1 = node_x_get(cid, missing)
        if x0 is missing or x1 is missing:
            continue

        y0 = node_y_get(pid, 0)
        y1 = node_y_get(cid, 0)

        append_horizontal(((x0, y1), (x1, y1)))
        append_vertical(((x0, y0), (x0, y1)))
        if color_is_callable:
            append_horizontal_color(color(clade))

    if horizontal_segments:
        ax.add_collection(
            LineCollection(
                horizontal_segments,
                colors=horizontal_colors,
                linewidths=lw,
                capstyle="round",
            ),
            autolim=True,
        )
    if vertical_segments:
        ax.add_collection(
            LineCollection(
                vertical_segments,
                colors=vertical_color,
                linewidths=vertical_lw,
                capstyle="round",
            ),
            autolim=True,
        )
    ax.autoscale_view()


def draw_tip_labels(
    ax, tree, node_x, node_y, fontsize=9, offset_fraction=0.03,
):
    """Draw taxon name labels at tree tips."""
    max_x = max(node_x.values()) if node_x else 1.0
    offset = max_x * offset_fraction

    if fontsize <= 0:
        return

    for tip in _terminal_clades(tree):
        ax.text(
            node_x[id(tip)] + offset, node_y[id(tip)],
            tip.name, va="center", fontsize=fontsize,
        )


def cleanup_tree_axes(ax, show_xlabel=True):
    """Standard axis cleanup for rectangular tree plots."""
    ax.set_yticks([])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    if show_xlabel:
        ax.set_xlabel("Branch length")
