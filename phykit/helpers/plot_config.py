import sys
from dataclasses import dataclass, field
from typing import List, Optional


VALID_LEGEND_POSITIONS = frozenset([
    "best", "upper right", "upper left", "lower left", "lower right",
    "right", "center left", "center right", "lower center", "upper center",
    "center", "none",
])


@dataclass
class PlotConfig:
    fig_width: Optional[float] = None
    fig_height: Optional[float] = None
    dpi: int = 300
    show_title: bool = True
    title: Optional[str] = None
    legend_position: Optional[str] = None
    ylabel_fontsize: Optional[float] = None
    xlabel_fontsize: Optional[float] = None
    title_fontsize: Optional[float] = None
    axis_fontsize: Optional[float] = None
    colors: Optional[List[str]] = None
    ladderize: bool = False

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

    def merge_colors(self, defaults: List[str]) -> List[str]:
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
