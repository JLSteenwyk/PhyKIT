from typing import Dict, List

import numpy as np

from .base import Alignment
from ...helpers.json_output import print_json
from ...helpers.plot_config import PlotConfig


class AlignmentEntropy(Alignment):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(alignment_file_path=parsed["alignment_file_path"], verbose=parsed["verbose"])
        self.json_output = parsed["json_output"]
        self.plot = parsed["plot"]
        self.plot_output = parsed["plot_output"]
        self.plot_config = parsed["plot_config"]

    def run(self) -> None:
        alignment, _, is_protein = self.get_alignment_and_format()
        entropies = self.calculate_site_entropies(alignment, is_protein)
        rows = [
            dict(site=idx, entropy=round(entropy, 4))
            for idx, entropy in enumerate(entropies, start=1)
        ]

        if self.plot:
            self._plot_alignment_entropy(rows)

        if self.json_output:
            if self.verbose:
                payload = dict(
                    verbose=True,
                    rows=rows,
                    sites=rows,
                )
            else:
                payload = dict(
                    verbose=False,
                    mean_entropy=round(float(np.mean(entropies)), 4) if entropies else 0.0,
                )
            if self.plot:
                payload["plot_output"] = self.plot_output
            print_json(payload)
            return

        if self.verbose:
            for row in rows:
                print(f"{row['site']}\t{row['entropy']}")
        else:
            if entropies:
                print(round(float(np.mean(entropies)), 4))
            else:
                print(0.0)

        if self.plot:
            print(f"Saved alignment entropy plot: {self.plot_output}")

    def process_args(self, args) -> Dict[str, str]:
        return dict(
            alignment_file_path=args.alignment,
            verbose=args.verbose,
            json_output=getattr(args, "json", False),
            plot=getattr(args, "plot", False),
            plot_output=getattr(args, "plot_output", "alignment_entropy_plot.png"),
            plot_config=PlotConfig.from_args(args),
        )

    def _plot_alignment_entropy(self, rows: List[Dict[str, float]]) -> None:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        except ImportError:
            print("matplotlib is required for --plot in alignment_entropy. Install matplotlib and retry.")
            raise SystemExit(2)

        if not rows:
            return

        sites = np.array([int(row["site"]) for row in rows], dtype=np.int32)
        entropies = np.array([float(row["entropy"]) for row in rows], dtype=np.float64)

        config = self.plot_config
        config.resolve(n_rows=len(sites), n_cols=None)
        default_colors = ["#2b8cbe"]
        colors = config.merge_colors(default_colors)

        fig, ax = plt.subplots(figsize=(config.fig_width, config.fig_height))
        ax.plot(sites, entropies, color=colors[0], linewidth=1.2, alpha=0.9)
        ax.scatter(sites, entropies, s=8, color=colors[0], alpha=0.75, edgecolors="none")
        ax.set_xlabel("Alignment site")
        ax.set_ylabel("Shannon entropy")
        ax.set_xlim(1, int(np.max(sites)))
        max_entropy = float(np.max(entropies)) if entropies.size else 0.0
        ax.set_ylim(0, max(1.0, max_entropy * 1.1))

        if config.show_title:
            ax.set_title(config.title or "Alignment Entropy Per Site", fontsize=config.title_fontsize)
        if config.axis_fontsize:
            ax.xaxis.label.set_fontsize(config.axis_fontsize)
            ax.yaxis.label.set_fontsize(config.axis_fontsize)

        fig.tight_layout()
        fig.savefig(self.plot_output, dpi=config.dpi, bbox_inches="tight")
        plt.close(fig)

    def calculate_site_entropies(self, alignment, is_protein: bool) -> List[float]:
        if is_protein:
            invalid_chars = np.array(["-", "?", "*", "X"], dtype="U1")
        else:
            invalid_chars = np.array(["-", "?", "*", "X", "N"], dtype="U1")
        alignment_array = np.array(
            [[c.upper() for c in str(record.seq)] for record in alignment],
            dtype="U1",
        )

        site_entropies: List[float] = []
        for col_idx in range(alignment_array.shape[1]):
            column = alignment_array[:, col_idx]
            valid = column[~np.isin(column, invalid_chars)]
            if valid.size == 0:
                site_entropies.append(0.0)
                continue

            _, counts = np.unique(valid, return_counts=True)
            probs = counts / counts.sum()
            entropy = float(-np.sum(probs * np.log2(probs)))
            site_entropies.append(entropy)

        return site_entropies
