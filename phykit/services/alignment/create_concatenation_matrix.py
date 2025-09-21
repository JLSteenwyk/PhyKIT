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

from .base import Alignment
from ...helpers.files import read_single_column_file_to_list


class CreateConcatenationMatrix(Alignment):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self) -> None:
        self.create_concatenation_matrix(
            self.alignment_list_path,
            self.prefix
        )

    def process_args(self, args) -> Dict[str, str]:
        return dict(alignment_list_path=args.alignment_list, prefix=args.prefix)

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
            with ProcessPoolExecutor(max_workers=min(mp.cpu_count(), len(alignment_paths))) as executor:
                futures = [executor.submit(self._get_taxa_from_alignment, path) for path in alignment_paths]
                for future in as_completed(futures):
                    taxa.update(future.result())
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
            print(f"No sequence records found. Exiting...")
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

        self.print_start_message(taxa, alignment_paths, file_partition, fasta_output, file_occupancy)

        # Initialize placeholders for partition info
        first_len, second_len = 1, 0
        partition_info, occupancy_info = [], []
        concatenated_seqs = defaultdict(list)

        # Process alignment files in parallel if there are many
        if len(alignment_paths) > 2:
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

                    # Add sequences to concatenated dict
                    for taxon in taxa:
                        concatenated_seqs[taxon].append(seq_dict[taxon])

                    # Add to partition and occupancy info
                    partition_info, first_len, second_len = self.add_to_partition_info(
                        partition_info, og_len, "AUTO", alignment_path, first_len, second_len
                    )
                    occupancy_info = self.add_to_occupancy_info(occupancy_info, present_taxa, taxa, alignment_path)
        else:
            # Process sequentially for small datasets
            for alignment_path in alignment_paths:
                present_taxa, records = self.get_list_of_taxa_and_records(alignment_path)
                missing_seq, og_len = self.create_missing_seq_str(records)

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

        # Write output files
        self.fasta_file_write(fasta_output, concatenated_seqs)
        self.write_occupancy_or_partition_file(occupancy_info, file_occupancy)
        self.write_occupancy_or_partition_file(partition_info, file_partition)

        print("Complete!\n")
