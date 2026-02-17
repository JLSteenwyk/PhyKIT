import numpy as np

from .base import Alignment
from ...helpers.json_output import print_json


class RelativeCompositionVariabilityTaxon(Alignment):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(alignment_file_path=parsed["alignment_file_path"])
        self.json_output = parsed["json_output"]
        self.plot = parsed["plot"]
        self.plot_output = parsed["plot_output"]

    def run(self):
        alignment, _, is_protein = self.get_alignment_and_format()
        rows = self.calculate_rows(alignment, is_protein)

        if self.json_output:
            payload = dict(rows=rows, taxa=rows)
            if self.plot:
                self._plot_rcvt(rows)
                payload["plot_output"] = self.plot_output
            print_json(payload)
            return

        for row in rows:
            print(f"{row['taxon']}\t{row['rcvt']}")

        if self.plot:
            self._plot_rcvt(rows)
            print(f"Saved RCVT plot: {self.plot_output}")

    def _plot_rcvt(self, rows):
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        except ImportError:
            print("matplotlib is required for --plot in rcvt. Install matplotlib and retry.")
            raise SystemExit(2)

        sorted_rows = sorted(rows, key=lambda row: row["rcvt"], reverse=True)
        taxa = [row["taxon"] for row in sorted_rows]
        rcvt_values = [row["rcvt"] for row in sorted_rows]

        fig, ax = plt.subplots(figsize=(10, 5))
        ax.bar(np.arange(len(taxa)), rcvt_values, color="#2b8cbe")
        ax.set_title("RCVT Per Taxon")
        ax.set_ylabel("RCVT")
        ax.set_xticks(np.arange(len(taxa)))
        ax.set_xticklabels(taxa, rotation=90, fontsize=8)
        fig.tight_layout()
        fig.savefig(self.plot_output, dpi=300, bbox_inches="tight")
        plt.close(fig)

    def calculate_rows(self, alignment, is_protein: bool):
        num_records = len(alignment)
        if num_records == 0:
            return []

        alignment_array = np.array([
            list(str(record.seq).upper()) for record in alignment
        ], dtype="U1")

        if is_protein:
            invalid_chars = np.array(["-", "?", "*", "X"], dtype="U1")
        else:
            invalid_chars = np.array(["-", "?", "*", "X", "N"], dtype="U1")
        valid_mask = ~np.isin(alignment_array, invalid_chars)
        valid_lengths = np.sum(valid_mask, axis=1).astype(np.float64)

        valid_chars = alignment_array[valid_mask]
        if valid_chars.size == 0:
            return [dict(taxon=record.id, rcvt=0.0) for record in alignment]

        unique_chars = np.unique(valid_chars)
        count_matrix = np.zeros((num_records, len(unique_chars)), dtype=np.float32)
        for seq_idx, seq in enumerate(alignment_array):
            seq_valid = valid_mask[seq_idx]
            for char_idx, char in enumerate(unique_chars):
                count_matrix[seq_idx, char_idx] = np.sum((seq == char) & seq_valid)

        average_counts = np.sum(count_matrix, axis=0) / num_records
        deviations = np.abs(count_matrix - average_counts)
        seq_sums = np.sum(deviations, axis=1)
        denom = num_records * valid_lengths
        rcv_values = np.divide(
            seq_sums,
            denom,
            out=np.zeros_like(seq_sums, dtype=np.float64),
            where=denom > 0,
        )

        return [
            dict(taxon=record.id, rcvt=round(float(rcv_values[i]), 4))
            for i, record in enumerate(alignment)
        ]

    def process_args(self, args):
        return dict(
            alignment_file_path=args.alignment,
            json_output=getattr(args, "json", False),
            plot=getattr(args, "plot", False),
            plot_output=getattr(args, "plot_output", "rcvt_plot.png"),
        )
