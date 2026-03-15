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
