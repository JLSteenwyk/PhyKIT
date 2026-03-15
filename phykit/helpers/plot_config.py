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
