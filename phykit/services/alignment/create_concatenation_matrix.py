from __future__ import annotations

import sys
import os
from collections import defaultdict

from ._fasta import (
    _clean_sequence,
    read_fasta_first_token_records,
    read_fasta_first_token_set,
)
from .base import Alignment


_OCCUPANCY_INVALID_CHARS = {"-", "?", "*", "X", "x", "N", "n"}
_OCCUPANCY_INVALID_BYTES = b"-?*XxNn"
_OCCUPANCY_STATE_LOOKUP = None
_FASTA_WRITE_CHUNK_ROWS = 4096
_OCCUPANCY_INVALID_SCAN_BYTES = 4096


class _LazyNumpy:
    _module = None

    def __getattr__(self, name):
        module = self._module
        if module is None:
            import numpy as _np

            module = _np
            self._module = module

        value = getattr(module, name)
        setattr(self, name, value)
        return value


np = _LazyNumpy()


def _dedent(*args, **kwargs):
    from textwrap import dedent

    return dedent(*args, **kwargs)


def _has_occupancy_invalid_bytes(sequence_bytes: bytes) -> bool:
    if len(sequence_bytes) <= _OCCUPANCY_INVALID_SCAN_BYTES:
        return any(code in sequence_bytes for code in _OCCUPANCY_INVALID_BYTES)

    if any(
        code in sequence_bytes[:_OCCUPANCY_INVALID_SCAN_BYTES]
        for code in _OCCUPANCY_INVALID_BYTES
    ):
        return True
    if any(
        code in sequence_bytes[-_OCCUPANCY_INVALID_SCAN_BYTES:]
        for code in _OCCUPANCY_INVALID_BYTES
    ):
        return True
    return any(code in sequence_bytes for code in _OCCUPANCY_INVALID_BYTES)


class _LazyMultiprocessing:
    def cpu_count(self):
        import multiprocessing as _mp

        return _mp.cpu_count()


mp = _LazyMultiprocessing()


def ProcessPoolExecutor(*args, **kwargs):
    from concurrent.futures import ProcessPoolExecutor as _ProcessPoolExecutor

    return _ProcessPoolExecutor(*args, **kwargs)


def as_completed(*args, **kwargs):
    from concurrent.futures import as_completed as _as_completed

    return _as_completed(*args, **kwargs)


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


def read_single_column_file_to_list(*args, **kwargs):
    from ...helpers.files import (
        read_single_column_file_to_list as _read_single_column_file_to_list,
    )

    return _read_single_column_file_to_list(*args, **kwargs)


def _occupancy_state_lookup():
    global _OCCUPANCY_STATE_LOOKUP
    if _OCCUPANCY_STATE_LOOKUP is None:
        lookup = np.full(256, 2, dtype=np.uint8)
        lookup[
            np.frombuffer(
                "".join(_OCCUPANCY_INVALID_CHARS).encode("ascii"),
                dtype=np.uint8,
            )
        ] = 1
        _OCCUPANCY_STATE_LOOKUP = lookup
    return _OCCUPANCY_STATE_LOOKUP


class _ParsedFastaRecord:
    __slots__ = ("id", "seq")

    def __init__(self, record_id: str, sequence: str) -> None:
        self.id = record_id
        self.seq = sequence


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
        self.plot_config = parsed["plot_config"]

    def run(self) -> None:
        self.create_concatenation_matrix(
            self.alignment_list_path,
            self.prefix
        )

    def process_args(self, args) -> dict[str, str]:
        from ...helpers.plot_config import PlotConfig

        return dict(
            alignment_list_path=args.alignment_list,
            prefix=args.prefix,
            json_output=getattr(args, "json", False),
            plot_occupancy=getattr(args, "plot_occupancy", False),
            plot_output=getattr(args, "plot_output", None),
            threshold=getattr(args, "threshold", 0),
            plot_config=PlotConfig.from_args(args),
        )

    @staticmethod
    def _sequence_occupancy_states(sequence: str) -> np.ndarray:
        try:
            sequence_array = np.frombuffer(sequence.encode("ascii"), dtype=np.uint8)
            return _occupancy_state_lookup()[sequence_array]
        except UnicodeEncodeError:
            return np.fromiter(
                (
                    1 if character in _OCCUPANCY_INVALID_CHARS else 2
                    for character in sequence
                ),
                dtype=np.uint8,
                count=len(sequence),
            )

    @classmethod
    def _build_occupancy_state_matrix(
        cls,
        taxa: list[str],
        concatenated_seqs: dict[str, list[str]],
        present_taxa_by_gene: list[set],
        gene_lengths: list[int],
    ) -> tuple[np.ndarray, list[int]]:
        taxa_count = len(taxa)
        if taxa_count and present_taxa_by_gene:
            taxa_set = set(taxa)
            if all(present_taxa == taxa_set for present_taxa in present_taxa_by_gene):
                return cls._build_complete_occupancy_state_matrix(
                    taxa,
                    concatenated_seqs,
                    gene_lengths,
                )

        # State values: 0 absent block, 1 gap/ambiguous, 2 represented.
        total_len = int(sum(gene_lengths))
        state_matrix = np.zeros((taxa_count, total_len), dtype=np.uint8)

        gene_boundaries = []
        taxon_to_idx = {taxon: idx for idx, taxon in enumerate(taxa)}
        cursor = 0
        for gene_idx, gene_len in enumerate(gene_lengths):
            start = cursor
            end = cursor + gene_len
            gene_boundaries.append(end)

            row_indices = []
            sequences = []
            for taxon in present_taxa_by_gene[gene_idx]:
                taxon_idx = taxon_to_idx.get(taxon)
                if taxon_idx is None:
                    continue
                row_indices.append(taxon_idx)
                sequences.append(concatenated_seqs[taxon][gene_idx])

            if sequences and gene_len and all(
                len(sequence) == gene_len for sequence in sequences
            ):
                try:
                    sequence_array = np.frombuffer(
                        "".join(sequences).encode("ascii"),
                        dtype=np.uint8,
                    )
                except UnicodeEncodeError:
                    sequence_array = None

                if sequence_array is not None:
                    state_matrix[
                        np.asarray(row_indices),
                        start:end,
                    ] = _occupancy_state_lookup()[sequence_array].reshape(
                        len(sequences),
                        gene_len,
                    )
                    cursor = end
                    continue

            for taxon_idx, sequence in zip(row_indices, sequences):
                seq_states = cls._sequence_occupancy_states(sequence)
                state_matrix[taxon_idx, start:start + len(seq_states)] = seq_states
            cursor = end

        return state_matrix, gene_boundaries

    @classmethod
    def _build_complete_occupancy_state_matrix(
        cls,
        taxa: list[str],
        concatenated_seqs: dict[str, list[str]],
        gene_lengths: list[int],
    ) -> tuple[np.ndarray, list[int]]:
        total_len = int(sum(gene_lengths))
        state_matrix = np.zeros((len(taxa), total_len), dtype=np.uint8)
        gene_boundaries = []
        cursor = 0
        lookup = _occupancy_state_lookup()
        for gene_idx, gene_len in enumerate(gene_lengths):
            start = cursor
            end = cursor + gene_len
            gene_boundaries.append(end)
            sequences = [concatenated_seqs[taxon][gene_idx] for taxon in taxa]
            if gene_len and all(len(sequence) == gene_len for sequence in sequences):
                try:
                    sequence_array = np.frombuffer(
                        "".join(sequences).encode("ascii"),
                        dtype=np.uint8,
                    )
                except UnicodeEncodeError:
                    sequence_array = None

                if sequence_array is not None:
                    state_matrix[:, start:end] = lookup[sequence_array].reshape(
                        len(sequences),
                        gene_len,
                    )
                    cursor = end
                    continue

            for taxon_idx, sequence in enumerate(sequences):
                seq_states = cls._sequence_occupancy_states(sequence)
                state_matrix[taxon_idx, start:start + len(seq_states)] = seq_states
            cursor = end

        return state_matrix, gene_boundaries

    def _plot_concatenation_occupancy(
        self,
        taxa: list[str],
        alignment_paths: list[str],
        concatenated_seqs: dict[str, list[str]],
        present_taxa_by_gene: list[set],
        gene_lengths: list[int],
        output_file: str,
    ) -> None:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            from matplotlib.collections import LineCollection
            from matplotlib.colors import ListedColormap
            from matplotlib.patches import Patch
        except ImportError:
            print("matplotlib is required for --plot-occupancy. Install matplotlib and retry.")
            sys.exit(2)

        config = self.plot_config
        config.resolve(n_rows=len(taxa), n_cols=len(alignment_paths))

        state_matrix, gene_boundaries = self._build_occupancy_state_matrix(
            taxa,
            concatenated_seqs,
            present_taxa_by_gene,
            gene_lengths,
        )

        # Sort taxa by total represented occupancy (state == 2), descending
        represented_counts = np.count_nonzero(state_matrix == 2, axis=1)
        order = np.argsort(-represented_counts)
        state_matrix = state_matrix[order, :]
        taxa_sorted = [taxa[idx] for idx in order]

        default_colors = ["#525252", "#d9d9d9", "#2b8cbe"]
        colors = config.merge_colors(default_colors)

        fig, ax = plt.subplots(figsize=(config.fig_width, config.fig_height))
        cmap = ListedColormap(colors)
        ax.imshow(state_matrix, aspect="auto", interpolation="nearest", cmap=cmap, vmin=0, vmax=2)

        boundary_segments = [
            ((boundary - 0.5, 0.0), (boundary - 0.5, 1.0))
            for boundary in gene_boundaries[:-1]
        ]
        if boundary_segments:
            ax.add_collection(
                LineCollection(
                    boundary_segments,
                    colors="black",
                    linewidths=0.6,
                    alpha=0.8,
                    transform=ax.get_xaxis_transform(),
                )
            )

        # Label genes at centers when feasible (controlled by xlabel_fontsize)
        if config.xlabel_fontsize and config.xlabel_fontsize > 0:
            starts = [0] + gene_boundaries[:-1]
            centers = [((s + e) / 2) - 0.5 for s, e in zip(starts, gene_boundaries)]
            labels = [os.path.basename(path) for path in alignment_paths]
            ax.set_xticks(centers)
            ax.set_xticklabels(labels, rotation=90, fontsize=config.xlabel_fontsize)
        else:
            ax.set_xticks([])
            ax.set_xlabel("Concatenated alignment positions (gene boundaries shown)")

        # Y-axis tick labels (controlled by ylabel_fontsize)
        if config.ylabel_fontsize and config.ylabel_fontsize > 0:
            ax.set_yticks(np.arange(len(taxa_sorted)))
            ax.set_yticklabels(taxa_sorted, fontsize=config.ylabel_fontsize)
        else:
            ax.set_yticks([])

        ax.set_ylabel("Taxa (sorted by represented occupancy)")

        legend_handles = [
            Patch(facecolor=colors[2], label="Represented character"),
            Patch(facecolor=colors[1], label="Gap/Ambiguous in present gene"),
            Patch(facecolor=colors[0], label="Gene absent (placeholder block)"),
        ]

        legend_loc = config.legend_position or "upper right"
        if legend_loc != "none":
            ax.legend(handles=legend_handles, loc=legend_loc, fontsize=8, frameon=True)

        # Apply title via config
        if config.show_title:
            title_text = config.title if config.title is not None else "Concatenation Occupancy Map"
            ax.set_title(title_text, fontsize=config.title_fontsize)

        if config.axis_fontsize:
            ax.xaxis.label.set_fontsize(config.axis_fontsize)
            ax.yaxis.label.set_fontsize(config.axis_fontsize)

        fig.tight_layout()
        fig.savefig(output_file, dpi=config.dpi, bbox_inches="tight")
        plt.close(fig)

    def read_alignment_paths(self, alignment_list_path: str) -> list[str]:
        try:
            return read_single_column_file_to_list(alignment_list_path)
        except FileNotFoundError:
            print("Alignment list file (-a) is not found. Please check pathing.")
            sys.exit(2)

    @staticmethod
    def _get_taxa_from_alignment(alignment_path: str) -> set:
        """Extract taxa names from a single alignment file."""
        return read_fasta_first_token_set(alignment_path)

    def get_taxa_names(self, alignment_paths: list[str]) -> list[str]:
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
        taxa: list[str],
        alignment_paths: list[str],
        file_partition: str,
        fasta_output: str,
        file_occupancy: str,
    ) -> None:
        start_message = _dedent(f"""
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
    ) -> tuple[set, list[object]]:
        return read_fasta_first_token_records(alignment_path, _ParsedFastaRecord)

    def create_missing_seq_str(self, records: list[object]) -> tuple[str, int]:
        """Create a placeholder string for sequences with missing taxa."""
        if not records:
            print("No sequence records found. Exiting...")
            sys.exit(2)

        og_len = len(records[0].seq)
        missing_seq = '?' * og_len
        return missing_seq, og_len

    def process_taxa_sequences(
        self,
        records: list[object],
        taxa: list[str],
        concatenated_seqs: dict[str, list[str]],
        missing_seq: str,
        taxa_set: set | None = None,
        present_taxa: set | None = None,
        missing_taxa: list[str] | None = None,
    ) -> None:
        if taxa_set is None:
            taxa_set = set(taxa)

        seqs_by_taxon = concatenated_seqs
        records_are_parsed_fasta = (
            bool(records) and type(records[0]) is _ParsedFastaRecord
        )
        if present_taxa is None:
            present_taxa = set()
            add_present_taxon = present_taxa.add
            if records_are_parsed_fasta:
                for record in records:
                    record_id = record.id
                    add_present_taxon(record_id)
                    seqs_by_taxon[record_id].append(record.seq)
            else:
                for record in records:
                    record_id = record.id
                    add_present_taxon(record_id)
                    seqs_by_taxon[record_id].append(str(record.seq))
        else:
            if records_are_parsed_fasta:
                for record in records:
                    seqs_by_taxon[record.id].append(record.seq)
            else:
                for record in records:
                    seqs_by_taxon[record.id].append(str(record.seq))

        if len(present_taxa) == len(taxa_set):
            return

        if missing_taxa is not None:
            for taxon in missing_taxa:
                seqs_by_taxon[taxon].append(missing_seq)
        elif len(taxa) == len(taxa_set):
            for taxon in taxa:
                if taxon not in present_taxa:
                    seqs_by_taxon[taxon].append(missing_seq)
        else:
            for taxon in taxa_set - present_taxa:
                seqs_by_taxon[taxon].append(missing_seq)

    def add_to_partition_info(
        self,
        partition_info: list[str],
        og_len: int,
        field_one: str,
        fasta: str,
        first_len: int,
        second_len: int,
    ) -> tuple[list[str], int, int]:
        second_len += og_len
        partition_info.append(f"{field_one}, {fasta}={first_len}-{second_len}\n")
        return partition_info, second_len + 1, second_len

    def add_to_occupancy_info(
        self,
        occupancy_info: list[str],
        present_taxa: set,
        taxa: list[str],
        fasta: str,
        sorted_taxa: list[str] | None = None,
        total_taxa_count: int | None = None,
        missing_taxa: list[str] | None = None,
    ) -> list[str]:
        if total_taxa_count is None:
            total_taxa_count = len(taxa)
        if missing_taxa is None:
            if sorted_taxa is None:
                sorted_taxa = sorted(set(taxa))
            missing_taxa = [taxon for taxon in sorted_taxa if taxon not in present_taxa]

        num_present = len(present_taxa)
        num_missing = len(missing_taxa)
        percent_occupancy = num_present / total_taxa_count
        occupancy_info.append(
            f"{fasta}\t{num_present}\t{num_missing}\t"
            f"{percent_occupancy:.4f}\t{';'.join(missing_taxa)}\n"
        )
        return occupancy_info

    def fasta_file_write(
        self,
        fasta_output: str,
        concatenated_seqs: dict[str, list[str]],
    ) -> None:
        """Write concatenated sequences to FASTA file with buffered I/O."""
        with open(fasta_output, "w", buffering=8192) as final_fasta_file:
            rows = []
            append_row = rows.append
            chunk_size = _FASTA_WRITE_CHUNK_ROWS
            for row_index, (taxon, sequences) in enumerate(
                concatenated_seqs.items(),
                start=1,
            ):
                append_row(f">{taxon}\n{''.join(sequences)}\n")
                if row_index % chunk_size == 0:
                    final_fasta_file.write("".join(rows))
                    rows.clear()
            if rows:
                final_fasta_file.write("".join(rows))

    def write_occupancy_or_partition_file(self, info: list[str], output_file_name: str) -> None:
        with open(output_file_name, "w") as f:
            f.writelines(info)

    @staticmethod
    def _append_ordered_sequences(
        taxa: list[str],
        seq_dict: dict[str, str],
        taxon_seq_lists: list[list[str]],
    ) -> None:
        for seqs, taxon in zip(taxon_seq_lists, taxa):
            seqs.append(seq_dict[taxon])

    @staticmethod
    def _process_alignment_file(alignment_path: str, taxa: list[str]) -> tuple[str, dict[str, str], set, int]:
        """Process a single alignment file and return its data."""
        record_sequences = {}
        present_taxa = set()
        og_len = 0
        with open(alignment_path) as handle:
            record_id = None
            sequence_parts = []
            for line in handle:
                if line[0] == ">":
                    if record_id is not None:
                        sequence = _clean_sequence(sequence_parts)
                        if not og_len:
                            og_len = len(sequence)
                        if record_id not in record_sequences:
                            record_sequences[record_id] = sequence
                    record_id = line[1:].split(None, 1)[0]
                    present_taxa.add(record_id)
                    sequence_parts = []
                elif record_id is not None:
                    sequence_parts.append(line.rstrip())

            if record_id is not None:
                sequence = _clean_sequence(sequence_parts)
                if not og_len:
                    og_len = len(sequence)
                if record_id not in record_sequences:
                    record_sequences[record_id] = sequence

        if not record_sequences:
            return alignment_path, {}, present_taxa, 0

        if len(present_taxa) == len(taxa):
            return alignment_path, record_sequences, present_taxa, og_len

        missing_seq = '?' * og_len
        seq_dict = {
            taxon: record_sequences.get(taxon, missing_seq)
            for taxon in taxa
        }

        return alignment_path, seq_dict, present_taxa, og_len

    def _compute_effective_occupancy(
        self,
        concatenated_seqs: dict[str, list[str]],
        threshold: float,
    ) -> tuple[dict[str, float], set]:
        """Compute per-taxon effective occupancy and return taxa to exclude.

        Effective occupancy is the fraction of informative (non-gap,
        non-ambiguous) characters across all concatenated positions.
        """
        occupancy = {}
        excluded = set()

        for taxon, seq_parts in concatenated_seqs.items():
            sequence = ''.join(seq_parts)
            total_positions = len(sequence)
            try:
                sequence_bytes = sequence.encode("ascii")
                if _has_occupancy_invalid_bytes(sequence_bytes):
                    informative = len(
                        sequence_bytes.translate(
                            None,
                            _OCCUPANCY_INVALID_BYTES,
                        )
                    )
                else:
                    informative = len(sequence_bytes)
            except UnicodeEncodeError:
                informative = sum(
                    character not in _OCCUPANCY_INVALID_CHARS
                    for character in sequence
                )
            eff = informative / total_positions if total_positions > 0 else 0.0
            occupancy[taxon] = eff
            if eff < threshold:
                excluded.add(taxon)

        return occupancy, excluded

    @staticmethod
    def _excluded_taxa_info(
        excluded: set,
        occupancy_scores: dict[str, float],
    ) -> list[dict[str, str | float]]:
        return [
            {"taxon": taxon, "effective_occupancy": round(occupancy_scores[taxon], 4)}
            for taxon in sorted(excluded)
        ]

    def create_concatenation_matrix(self, alignment_list_path: str, prefix: str) -> None:
        alignment_paths = self.read_alignment_paths(alignment_list_path)
        taxa = self.get_taxa_names(alignment_paths)
        sorted_taxa = taxa
        taxa_set = set(taxa)
        total_taxa_count = len(taxa)

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
        taxon_seq_lists = [concatenated_seqs[taxon] for taxon in taxa]
        present_taxa_by_gene = []
        gene_lengths = []

        # Process alignment files in parallel if there are many
        if len(alignment_paths) > 2:
            try:
                from functools import partial

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
                        missing_taxa = [
                            taxon for taxon in sorted_taxa if taxon not in present_taxa
                        ]
                        present_taxa_by_gene.append(present_taxa)
                        gene_lengths.append(og_len)

                        # Add sequences to concatenated dict
                        self._append_ordered_sequences(
                            taxa,
                            seq_dict,
                            taxon_seq_lists,
                        )

                        # Add to partition and occupancy info
                        partition_info, first_len, second_len = self.add_to_partition_info(
                            partition_info, og_len, "AUTO", alignment_path, first_len, second_len
                        )
                        occupancy_info = self.add_to_occupancy_info(
                            occupancy_info,
                            present_taxa,
                            taxa,
                            alignment_path,
                            sorted_taxa=sorted_taxa,
                            total_taxa_count=total_taxa_count,
                            missing_taxa=missing_taxa,
                        )
            except (PermissionError, OSError, RuntimeError, NotImplementedError):
                for alignment_path in alignment_paths:
                    _, seq_dict, present_taxa, og_len = self._process_alignment_file(
                        alignment_path,
                        taxa,
                    )
                    missing_taxa = [
                        taxon for taxon in sorted_taxa if taxon not in present_taxa
                    ]
                    present_taxa_by_gene.append(present_taxa)
                    gene_lengths.append(og_len)
                    self._append_ordered_sequences(
                        taxa,
                        seq_dict,
                        taxon_seq_lists,
                    )
                    partition_info, first_len, second_len = self.add_to_partition_info(
                        partition_info, og_len, "AUTO", alignment_path, first_len, second_len
                    )
                    occupancy_info = self.add_to_occupancy_info(
                        occupancy_info,
                        present_taxa,
                        taxa,
                        alignment_path,
                        sorted_taxa=sorted_taxa,
                        total_taxa_count=total_taxa_count,
                        missing_taxa=missing_taxa,
                    )
        else:
            # Process sequentially for small datasets
            for alignment_path in alignment_paths:
                present_taxa, records = self.get_list_of_taxa_and_records(alignment_path)
                missing_seq, og_len = self.create_missing_seq_str(records)
                missing_taxa = [
                    taxon for taxon in sorted_taxa if taxon not in present_taxa
                ]
                present_taxa_by_gene.append(present_taxa)
                gene_lengths.append(og_len)

                # Process taxa sequences and add to the concatenated sequences
                self.process_taxa_sequences(
                    records,
                    taxa,
                    concatenated_seqs,
                    missing_seq,
                    taxa_set=taxa_set,
                    present_taxa=present_taxa,
                    missing_taxa=missing_taxa,
                )

                # Add to partition and occupancy info
                partition_info, first_len, second_len = self.add_to_partition_info(
                    partition_info, og_len, "AUTO", alignment_path, first_len, second_len
                )
                occupancy_info = self.add_to_occupancy_info(
                    occupancy_info,
                    present_taxa,
                    taxa,
                    alignment_path,
                    sorted_taxa=sorted_taxa,
                    total_taxa_count=total_taxa_count,
                    missing_taxa=missing_taxa,
                )

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
                excluded_taxa_info = self._excluded_taxa_info(
                    excluded,
                    occupancy_scores,
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
                    missing_taxa = [
                        taxon for taxon in taxa if taxon not in filtered_present
                    ]
                    occupancy_info = self.add_to_occupancy_info(
                        occupancy_info,
                        filtered_present,
                        taxa,
                        alignment_path,
                        sorted_taxa=taxa,
                        total_taxa_count=len(taxa),
                        missing_taxa=missing_taxa,
                    )
                if not self.json_output:
                    threshold_lines = [
                        f"Threshold ({self.threshold}): excluded {len(excluded)} taxa"
                    ]
                    threshold_lines.extend(
                        "  "
                        f"{info['taxon']}: effective occupancy = "
                        f"{info['effective_occupancy']}"
                        for info in excluded_taxa_info
                    )
                    print("\n".join(threshold_lines))

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
