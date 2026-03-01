import sys
import os
from textwrap import dedent
from typing import Dict, List, Tuple
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import partial
import multiprocessing as mp
from collections import defaultdict
import numpy as np

from .base import Alignment
from ...helpers.files import read_single_column_file_to_list
from ...helpers.json_output import print_json


class CreateConcatenationMatrix(Alignment):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(
            alignment_list_path=parsed["alignment_list_path"],
            prefix=parsed["prefix"],
        )
        self.json_output = parsed["json_output"]
        self.plot_occupancy = parsed["plot_occupancy"]
        self.plot_output = parsed["plot_output"]
        self.threshold = parsed["threshold"]

    def run(self) -> None:
        self.create_concatenation_matrix(
            self.alignment_list_path,
            self.prefix
        )

    def process_args(self, args) -> Dict[str, str]:
        return dict(
            alignment_list_path=args.alignment_list,
            prefix=args.prefix,
            json_output=getattr(args, "json", False),
            plot_occupancy=getattr(args, "plot_occupancy", False),
            plot_output=getattr(args, "plot_output", None),
            threshold=getattr(args, "threshold", 0),
        )

    def _plot_concatenation_occupancy(
        self,
        taxa: List[str],
        alignment_paths: List[str],
        concatenated_seqs: Dict[str, List[str]],
        present_taxa_by_gene: List[set],
        gene_lengths: List[int],
        output_file: str,
    ) -> None:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            from matplotlib.colors import ListedColormap
            from matplotlib.patches import Patch
        except ImportError:
            print("matplotlib is required for --plot-occupancy. Install matplotlib and retry.")
            sys.exit(2)

        # 0: absent gene block, 1: present but gap/ambiguous character, 2: represented character
        total_len = int(sum(gene_lengths))
        state_matrix = np.zeros((len(taxa), total_len), dtype=np.uint8)
        invalid_chars = set(["-", "?", "*", "X", "x", "N", "n"])

        gene_boundaries = []
        cursor = 0
        for gene_idx, gene_len in enumerate(gene_lengths):
            start = cursor
            end = cursor + gene_len
            gene_boundaries.append(end)
            present_taxa = present_taxa_by_gene[gene_idx]

            for taxon_idx, taxon in enumerate(taxa):
                if taxon not in present_taxa:
                    state_matrix[taxon_idx, start:end] = 0
                    continue
                seq = concatenated_seqs[taxon][gene_idx]
                for pos_idx, char in enumerate(seq):
                    if char in invalid_chars:
                        state_matrix[taxon_idx, start + pos_idx] = 1
                    else:
                        state_matrix[taxon_idx, start + pos_idx] = 2
            cursor = end

        # Sort taxa by total represented occupancy (state == 2), descending
        represented_counts = np.sum(state_matrix == 2, axis=1)
        order = np.argsort(-represented_counts)
        state_matrix = state_matrix[order, :]
        taxa_sorted = [taxa[idx] for idx in order]

        fig_height = max(5.0, min(18.0, 3.0 + len(taxa_sorted) * 0.18))
        fig, ax = plt.subplots(figsize=(14, fig_height))
        cmap = ListedColormap(["#525252", "#d9d9d9", "#2b8cbe"])
        ax.imshow(state_matrix, aspect="auto", interpolation="nearest", cmap=cmap, vmin=0, vmax=2)

        for boundary in gene_boundaries[:-1]:
            ax.axvline(boundary - 0.5, color="black", linewidth=0.6, alpha=0.8)

        # Label genes at centers when feasible
        if len(alignment_paths) <= 40:
            starts = [0] + gene_boundaries[:-1]
            centers = [((s + e) / 2) - 0.5 for s, e in zip(starts, gene_boundaries)]
            labels = [os.path.basename(path) for path in alignment_paths]
            ax.set_xticks(centers)
            ax.set_xticklabels(labels, rotation=90, fontsize=7)
        else:
            ax.set_xticks([])
            ax.set_xlabel("Concatenated alignment positions (gene boundaries shown)")

        ax.set_yticks(np.arange(len(taxa_sorted)))
        ax.set_yticklabels(taxa_sorted, fontsize=7)
        ax.set_ylabel("Taxa (sorted by represented occupancy)")
        ax.set_title("Concatenation Occupancy Map")

        legend_handles = [
            Patch(facecolor="#2b8cbe", label="Represented character"),
            Patch(facecolor="#d9d9d9", label="Gap/Ambiguous in present gene"),
            Patch(facecolor="#525252", label="Gene absent (placeholder block)"),
        ]
        ax.legend(handles=legend_handles, loc="upper right", fontsize=8, frameon=True)

        fig.tight_layout()
        fig.savefig(output_file, dpi=300, bbox_inches="tight")
        plt.close(fig)

    def read_alignment_paths(self, alignment_list_path: str) -> List[str]:
        try:
            return read_single_column_file_to_list(alignment_list_path)
        except FileNotFoundError:
            print("Alignment list file (-a) is not found. Please check pathing.")
            sys.exit(2)

    @staticmethod
    def _get_taxa_from_alignment(alignment_path: str) -> set:
        """Extract taxa names from a single alignment file."""
        return {seq_record.id for seq_record in SeqIO.parse(alignment_path, "fasta")}

    def get_taxa_names(self, alignment_paths: List[str]) -> List[str]:
        """Get all unique taxa names from alignment files in parallel."""
        taxa = set()

        # Process files in parallel if there are many
        if len(alignment_paths) > 10:
            try:
                with ProcessPoolExecutor(max_workers=min(mp.cpu_count(), len(alignment_paths))) as executor:
                    futures = [executor.submit(self._get_taxa_from_alignment, path) for path in alignment_paths]
                    for future in as_completed(futures):
                        taxa.update(future.result())
            except (PermissionError, OSError, RuntimeError, NotImplementedError):
                for alignment_path in alignment_paths:
                    taxa.update(self._get_taxa_from_alignment(alignment_path))
        else:
            # Process sequentially for small datasets
            for alignment_path in alignment_paths:
                taxa.update(self._get_taxa_from_alignment(alignment_path))

        return sorted(taxa)

    def print_start_message(
        self,
        taxa: List[str],
        alignment_paths: List[str],
        file_partition: str,
        fasta_output: str,
        file_occupancy: str,
    ) -> None:
        start_message = dedent(f"""
            --------------------
            | General features |
            --------------------
            Total number of taxa: {len(taxa)}
            Total number of alignments: {len(alignment_paths)}


            ----------------
            | Output files |
            ----------------
            Partition file output: {file_partition}
            Concatenated fasta output: {fasta_output}
            Occupancy report: {file_occupancy}
        """)
        print(start_message)

    def get_list_of_taxa_and_records(
        self, alignment_path: str
    ) -> Tuple[set, List[SeqRecord]]:
        records = list(SeqIO.parse(alignment_path, "fasta"))
        og_taxa = {record.id for record in records}
        return og_taxa, records

    def create_missing_seq_str(self, records: List[SeqRecord]) -> Tuple[str, int]:
        """Create a placeholder string for sequences with missing taxa."""
        if not records:
            print("No sequence records found. Exiting...")
            sys.exit(2)

        og_len = len(records[0].seq)
        missing_seq = '?' * og_len
        return missing_seq, og_len

    def process_taxa_sequences(
        self,
        records: List[SeqRecord],
        taxa: List[str],
        concatenated_seqs: Dict[str, List[str]],
        missing_seq: str,
    ) -> None:
        present_taxa = {record.id for record in records}
        missing_taxa = set(taxa) - present_taxa

        # Add sequences for present taxa
        for record in records:
            concatenated_seqs[record.id].append(str(record.seq))

        # Add missing sequences for missing taxa
        for taxon in missing_taxa:
            concatenated_seqs[taxon].append(missing_seq)

    def add_to_partition_info(
        self,
        partition_info: List[str],
        og_len: int,
        field_one: str,
        fasta: str,
        first_len: int,
        second_len: int,
    ) -> Tuple[List[str], int, int]:
        second_len += og_len
        partition_info.append(f"{field_one}, {fasta}={first_len}-{second_len}\n")
        return partition_info, second_len + 1, second_len

    def add_to_occupancy_info(
        self,
        occupancy_info: List[str],
        present_taxa: set,
        taxa: List[str],
        fasta: str,
    ) -> List[str]:
        missing_taxa = sorted(set(taxa) - present_taxa)
        num_present = len(present_taxa)
        num_missing = len(missing_taxa)
        percent_occupancy = num_present / len(taxa)
        occupancy_info.append(f"{fasta}\t{num_present}\t{num_missing}\t{percent_occupancy:.4f}\t{';'.join(missing_taxa)}\n")
        return occupancy_info

    def fasta_file_write(self, fasta_output: str, concatenated_seqs: Dict[str, List[str]]) -> None:
        """Write concatenated sequences to FASTA file with buffered I/O."""
        # Use larger buffer for better I/O performance
        with open(fasta_output, "w", buffering=8192) as final_fasta_file:
            for taxon, sequences in concatenated_seqs.items():
                # Join sequences once instead of in the write statement
                concatenated = ''.join(sequences)
                final_fasta_file.write(f">{taxon}\n{concatenated}\n")

    def write_occupancy_or_partition_file(self, info: List[str], output_file_name: str) -> None:
        with open(output_file_name, "w") as f:
            f.writelines(info)

    @staticmethod
    def _process_alignment_file(alignment_path: str, taxa: List[str]) -> Tuple[str, Dict[str, str], set, int]:
        """Process a single alignment file and return its data."""
        records = list(SeqIO.parse(alignment_path, "fasta"))
        present_taxa = {record.id for record in records}

        if not records:
            return alignment_path, {}, present_taxa, 0

        og_len = len(records[0].seq)
        missing_seq = '?' * og_len

        # Create sequence dict for this alignment
        seq_dict = {}
        for taxon in taxa:
            if taxon in present_taxa:
                # Find the sequence for this taxon
                for record in records:
                    if record.id == taxon:
                        seq_dict[taxon] = str(record.seq)
                        break
            else:
                seq_dict[taxon] = missing_seq

        return alignment_path, seq_dict, present_taxa, og_len

    def _compute_effective_occupancy(
        self,
        concatenated_seqs: Dict[str, List[str]],
        threshold: float,
    ) -> Tuple[Dict[str, float], set]:
        """Compute per-taxon effective occupancy and return taxa to exclude.

        Effective occupancy is the fraction of informative (non-gap,
        non-ambiguous) characters across all concatenated positions.
        """
        invalid_chars = {"-", "?", "*", "X", "x", "N", "n"}
        occupancy = {}
        excluded = set()

        for taxon, seq_parts in concatenated_seqs.items():
            total_positions = 0
            informative = 0
            for part in seq_parts:
                total_positions += len(part)
                for ch in part:
                    if ch not in invalid_chars:
                        informative += 1
            eff = informative / total_positions if total_positions > 0 else 0.0
            occupancy[taxon] = eff
            if eff < threshold:
                excluded.add(taxon)

        return occupancy, excluded

    def create_concatenation_matrix(self, alignment_list_path: str, prefix: str) -> None:
        alignment_paths = self.read_alignment_paths(alignment_list_path)
        taxa = self.get_taxa_names(alignment_paths)

        # Create output directory if needed
        output_dir = os.path.dirname(prefix)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir, exist_ok=True)

        # Assign output file names
        file_partition = f"{prefix}.partition"
        fasta_output = f"{prefix}.fa"
        file_occupancy = f"{prefix}.occupancy"

        if not self.json_output:
            self.print_start_message(taxa, alignment_paths, file_partition, fasta_output, file_occupancy)

        # Initialize placeholders for partition info
        first_len, second_len = 1, 0
        partition_info, occupancy_info = [], []
        concatenated_seqs = defaultdict(list)
        present_taxa_by_gene = []
        gene_lengths = []

        # Process alignment files in parallel if there are many
        if len(alignment_paths) > 2:
            try:
                with ProcessPoolExecutor(max_workers=min(mp.cpu_count(), 8)) as executor:
                    process_func = partial(self._process_alignment_file, taxa=taxa)
                    # Keep results indexed by path to maintain order
                    futures = {executor.submit(process_func, path): path for path in alignment_paths}
                    results = {}

                    for future in as_completed(futures):
                        path = futures[future]
                        results[path] = future.result()

                    # Process results in original order
                    for alignment_path in alignment_paths:
                        _, seq_dict, present_taxa, og_len = results[alignment_path]
                        present_taxa_by_gene.append(present_taxa)
                        gene_lengths.append(og_len)

                        # Add sequences to concatenated dict
                        for taxon in taxa:
                            concatenated_seqs[taxon].append(seq_dict[taxon])

                        # Add to partition and occupancy info
                        partition_info, first_len, second_len = self.add_to_partition_info(
                            partition_info, og_len, "AUTO", alignment_path, first_len, second_len
                        )
                        occupancy_info = self.add_to_occupancy_info(occupancy_info, present_taxa, taxa, alignment_path)
            except (PermissionError, OSError, RuntimeError, NotImplementedError):
                for alignment_path in alignment_paths:
                    _, seq_dict, present_taxa, og_len = self._process_alignment_file(alignment_path, taxa)
                    present_taxa_by_gene.append(present_taxa)
                    gene_lengths.append(og_len)
                    for taxon in taxa:
                        concatenated_seqs[taxon].append(seq_dict[taxon])
                    partition_info, first_len, second_len = self.add_to_partition_info(
                        partition_info, og_len, "AUTO", alignment_path, first_len, second_len
                    )
                    occupancy_info = self.add_to_occupancy_info(occupancy_info, present_taxa, taxa, alignment_path)
        else:
            # Process sequentially for small datasets
            for alignment_path in alignment_paths:
                present_taxa, records = self.get_list_of_taxa_and_records(alignment_path)
                missing_seq, og_len = self.create_missing_seq_str(records)
                present_taxa_by_gene.append(present_taxa)
                gene_lengths.append(og_len)

                # Process taxa sequences and add to the concatenated sequences
                self.process_taxa_sequences(records, taxa, concatenated_seqs, missing_seq)

                # Add to partition and occupancy info
                partition_info, first_len, second_len = self.add_to_partition_info(
                    partition_info, og_len, "AUTO", alignment_path, first_len, second_len
                )
                occupancy_info = self.add_to_occupancy_info(occupancy_info, present_taxa, taxa, alignment_path)

        # Convert defaultdict to regular dict for writing
        if isinstance(concatenated_seqs, defaultdict):
            concatenated_seqs = dict(concatenated_seqs)

        # Apply threshold filtering
        excluded_taxa_info = []
        if self.threshold > 0:
            occupancy_scores, excluded = self._compute_effective_occupancy(
                concatenated_seqs, self.threshold
            )
            if excluded:
                excluded_taxa_info = sorted(
                    [
                        {"taxon": t, "effective_occupancy": round(occupancy_scores[t], 4)}
                        for t in excluded
                    ],
                    key=lambda x: x["taxon"],
                )
                for t in excluded:
                    del concatenated_seqs[t]
                taxa = [t for t in taxa if t not in excluded]
                # Rebuild occupancy info with filtered taxa
                occupancy_info = []
                for gene_idx, alignment_path in enumerate(alignment_paths):
                    present_taxa = present_taxa_by_gene[gene_idx]
                    # Only consider taxa that survived filtering
                    filtered_present = present_taxa - excluded
                    occupancy_info = self.add_to_occupancy_info(
                        occupancy_info, filtered_present, taxa, alignment_path
                    )
                if not self.json_output:
                    print(f"Threshold ({self.threshold}): excluded {len(excluded)} taxa")
                    for info in excluded_taxa_info:
                        print(f"  {info['taxon']}: effective occupancy = {info['effective_occupancy']}")

        # Write output files
        self.fasta_file_write(fasta_output, concatenated_seqs)
        self.write_occupancy_or_partition_file(occupancy_info, file_occupancy)
        self.write_occupancy_or_partition_file(partition_info, file_partition)

        plot_output = self.plot_output or f"{prefix}.occupancy.png"
        if self.plot_occupancy:
            self._plot_concatenation_occupancy(
                taxa=taxa,
                alignment_paths=alignment_paths,
                concatenated_seqs=concatenated_seqs,
                present_taxa_by_gene=present_taxa_by_gene,
                gene_lengths=gene_lengths,
                output_file=plot_output,
            )

        if self.json_output:
            payload = dict(
                input_alignment_list=alignment_list_path,
                total_taxa=len(taxa),
                total_alignments=len(alignment_paths),
                concatenated_length=second_len,
                threshold=self.threshold,
                excluded_taxa=excluded_taxa_info,
                output_files=dict(
                    fasta=fasta_output,
                    partition=file_partition,
                    occupancy=file_occupancy,
                ),
            )
            if self.plot_occupancy:
                payload["output_files"]["occupancy_plot"] = plot_output
            print_json(payload)
        else:
            if self.plot_occupancy:
                print(f"Occupancy plot output: {plot_output}")
            print("Complete!\n")
