from collections import Counter
from typing import List, Dict
import numpy as np

from Bio.Align import MultipleSeqAlignment

from .base import Alignment
from ...helpers.json_output import print_json
from ...helpers.plot_config import PlotConfig


class EvolutionaryRatePerSite(Alignment):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(alignment_file_path=parsed["alignment_file_path"])
        self.json_output = parsed["json_output"]
        self.plot = parsed["plot"]
        self.plot_output = parsed["plot_output"]
        self.plot_config = parsed["plot_config"]

    def run(self):
        alignment, _, is_protein = self.get_alignment_and_format()
        pic_values = self.calculate_evolutionary_rate_per_site(alignment, is_protein)
        rows = [
            dict(site=idx + 1, evolutionary_rate=round(value, 4))
            for idx, value in enumerate(pic_values)
        ]

        if self.plot:
            self._plot_evolutionary_rate_per_site(rows)

        if self.json_output:
            payload = dict(rows=rows, sites=rows)
            if self.plot:
                payload["plot_output"] = self.plot_output
            print_json(payload)
            return

        for row in rows:
            print(f"{row['site']}\t{row['evolutionary_rate']}")

        if self.plot:
            print(f"Saved evolutionary-rate plot: {self.plot_output}")

    def process_args(self, args):
        return dict(
            alignment_file_path=args.alignment,
            json_output=getattr(args, "json", False),
            plot=getattr(args, "plot", False),
            plot_output=getattr(args, "plot_output", "evolutionary_rate_per_site_plot.png"),
            plot_config=PlotConfig.from_args(args),
        )

    def _plot_evolutionary_rate_per_site(self, rows):
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        except ImportError:
            print("matplotlib is required for --plot in evolutionary_rate_per_site. Install matplotlib and retry.")
            raise SystemExit(2)

        if not rows:
            return

        sites = np.array([row["site"] for row in rows], dtype=np.int32)
        rates = np.array([row["evolutionary_rate"] for row in rows], dtype=np.float64)

        config = self.plot_config
        config.resolve(n_rows=len(sites), n_cols=None)
        default_colors = ["#2b8cbe"]
        colors = config.merge_colors(default_colors)

        fig, ax = plt.subplots(figsize=(config.fig_width, config.fig_height))
        ax.plot(sites, rates, color=colors[0], linewidth=1.2, alpha=0.9)
        ax.scatter(sites, rates, s=8, color=colors[0], alpha=0.75, edgecolors="none")
        ax.set_xlabel("Alignment site")
        ax.set_ylabel("Evolutionary rate")
        ax.set_xlim(1, int(np.max(sites)))
        ax.set_ylim(0, max(1.0, float(np.max(rates) * 1.1)))

        if config.show_title:
            ax.set_title(config.title or "Evolutionary Rate Per Site", fontsize=config.title_fontsize)
        if config.axis_fontsize:
            ax.xaxis.label.set_fontsize(config.axis_fontsize)
            ax.yaxis.label.set_fontsize(config.axis_fontsize)

        fig.tight_layout()
        fig.savefig(self.plot_output, dpi=config.dpi, bbox_inches="tight")
        plt.close(fig)

    def remove_gap_characters(self, seq: str, gap_chars: List[str]) -> str:
        return ''.join([char for char in seq if char not in gap_chars]).upper()

    def get_number_of_occurrences_per_character(
        self,
        alignment: MultipleSeqAlignment,
        idx: int,
        gap_chars: List[str]
    ) -> Dict[str, int]:
        seq_at_position = alignment[:, idx]
        clean_seq = self.remove_gap_characters(seq_at_position, gap_chars)

        return Counter(clean_seq)

    def calculate_pic(
        self,
        num_occurrences: Dict[str, int],
    ) -> float:
        total_frequencies = sum(num_occurrences.values())
        sum_of_frequencies = sum(
            (frequency / total_frequencies) ** 2
            for frequency in num_occurrences.values()
        )
        return 1 - sum_of_frequencies

    def calculate_evolutionary_rate_per_site(
        self,
        alignment: MultipleSeqAlignment,
        is_protein: bool = False,
    ) -> List[float]:
        aln_len = alignment.get_alignment_length()
        gap_chars = set(self.get_gap_chars(is_protein))

        # Convert alignment to numpy array for vectorized operations
        alignment_array = np.array([
            [c.upper() for c in str(record.seq)]
            for record in alignment
        ], dtype='U1')

        pic_values = []

        # Process each column
        for col_idx in range(aln_len):
            column = alignment_array[:, col_idx]

            # Filter out gaps
            non_gap_mask = ~np.isin(column, list(gap_chars))
            filtered_column = column[non_gap_mask]

            if len(filtered_column) > 0:
                # Count occurrences using numpy
                unique_chars, counts = np.unique(filtered_column, return_counts=True)
                total_frequencies = len(filtered_column)

                # Calculate PIC (Probability of Identical Characters)
                sum_of_frequencies = np.sum((counts / total_frequencies) ** 2)
                pic = 1 - sum_of_frequencies
            else:
                pic = 0

            pic_values.append(pic)

        return pic_values
