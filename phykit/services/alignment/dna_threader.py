import sys
from typing import Dict, List
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from .base import Alignment


class DNAThreader(Alignment):
    """
    Threads DNA on top of protein alignment
    """

    def __init__(self, args) -> None:
        self.process_args(args)

    def process_args(self, args):
        self.remove_stop_codon = args.stop
        self.protein_file_path = args.protein
        self.nucleotide_file_path = args.nucleotide
        self.clipkit_log_file = args.clipkit_log_file

    @property
    def clipkit_log_data(self) -> List[List[str]]:
        if self.clipkit_log_file:
            with open(self.clipkit_log_file) as f:
                return [line.rstrip("\n").split(" ") for line in f.readlines()]
        return None

    def run(self) -> None:
        prot_records = SeqIO.parse(self.protein_file_path, "fasta")
        pal2nal = self.thread(prot_records)

        for gene_id, sequence in pal2nal.items():
            print(f">{gene_id}")
            print(f"{sequence}")

    def create_mask(self, length: int) -> List[bool]:
        if not self.clipkit_log_data:
            return [True] * length

        # create a mask that replicates the 'keep' and 'remove' status 3 times for each amino acid
        return [
            True if row[1] == "keep" else False
            for row in self.clipkit_log_data
            for _ in range(3)
        ]

    def normalize_p_seq(self, p_seq: Seq) -> str:
        # triplicate each amino acid
        return ''.join(c * 3 for c in p_seq)

    def normalize_n_seq(self, n_seq: Seq, p_seq: Seq) -> str:
        # Pre-split codons for faster access
        codons = [str(n_seq[i:i+3]) for i in range(0, len(n_seq), 3)]
        normalized_n_seq = []
        gap_chars = {'-', '?', '*', 'X', 'x'}

        codon_idx = 0
        for aa in p_seq:
            if aa in gap_chars:
                normalized_n_seq.append("---")
            else:
                if codon_idx < len(codons):
                    normalized_n_seq.append(codons[codon_idx])
                    codon_idx += 1
                else:
                    normalized_n_seq.append("---")  # fallback in case of misalignment

        return ''.join(normalized_n_seq)

    def thread(self, prot_records) -> Dict[str, str]:
        pal2nal = dict()
        prot_dict = SeqIO.to_dict(prot_records)

        if not prot_dict:
            print("Protein file is empty or incorrectly formatted.")
            sys.exit(2)

        # Pre-load nucleotide sequences only for proteins we have
        nucl_records = {}
        for record in SeqIO.parse(self.nucleotide_file_path, "fasta"):
            if record.id in prot_dict:
                nucl_records[record.id] = record

        length = len(next(iter(prot_dict.values())).seq)
        keep_mask = self.create_mask(length * 3)

        # Convert keep_mask to numpy array for faster operations
        keep_mask_arr = np.array(keep_mask)
        gap_chars = {'-', '?', '*', 'X', 'x'}

        for gene_id, protein_seq_record in prot_dict.items():
            try:
                if gene_id not in nucl_records:
                    print(f"Nucleotide sequence for {gene_id} not found.")
                    sys.exit(2)

                p_seq = protein_seq_record.seq
                n_seq = nucl_records[gene_id].seq

                # Get normalized sequences
                normalized_p_seq = self.normalize_p_seq(p_seq)
                normalized_n_seq = self.normalize_n_seq(n_seq, normalized_p_seq)

                # Convert to numpy arrays for faster operations
                p_arr = np.array(list(normalized_p_seq), dtype='U1')
                n_arr = np.array(list(normalized_n_seq), dtype='U1')

                # Create mask for non-gap positions in protein
                non_gap_mask_protein = ~np.isin(p_arr, list(gap_chars))

                # Expand protein mask to nucleotide positions (each AA = 3 nucleotides)
                non_gap_mask = np.repeat(non_gap_mask_protein, 3)

                # Ensure masks have same shape
                min_len = min(len(non_gap_mask), len(keep_mask_arr), len(n_arr))
                non_gap_mask = non_gap_mask[:min_len]
                keep_mask_arr_trimmed = keep_mask_arr[:min_len]
                n_arr = n_arr[:min_len]

                # Combine masks
                final_mask = keep_mask_arr_trimmed & non_gap_mask

                # Apply masks and build result
                result = np.where(final_mask, n_arr, '-')

                # Handle stop codon if needed
                if self.remove_stop_codon and p_seq[-1] == "*":
                    # Find the last 3 positions that were kept
                    kept_indices = np.where(keep_mask_arr_trimmed)[0]
                    if len(kept_indices) >= 3:
                        last_3_indices = kept_indices[-3:]
                        for idx in last_3_indices:
                            if idx < len(result) and idx < len(n_arr):
                                result[idx] = n_arr[idx]

                # Only keep positions marked in keep_mask
                pal2nal[gene_id] = ''.join(result[keep_mask_arr_trimmed[:len(result)]])

            except KeyError:
                print(f"Nucleotide sequence for {gene_id} not found.")
                sys.exit(2)

        return pal2nal
