import sys
from typing import Dict, List
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

    def normalize_n_seq(self, n_seq: Seq, normalized_p_seq: str) -> str:
        # insert gaps at corresponding positions in the nucleotide sequence based on the protein sequence
        normalized_n_seq = list(n_seq)
        for idx, c in enumerate(normalized_p_seq):
            if c in "-?*Xx":
                normalized_n_seq.insert(idx, "-")
        return ''.join(normalized_n_seq)

    def thread(self, prot_records) -> Dict[str, str]:
        pal2nal = dict()
        nucl_records = SeqIO.to_dict(SeqIO.parse(self.nucleotide_file_path, "fasta"))
        prot_dict = SeqIO.to_dict(prot_records)

        if not prot_dict:
            print("Protein file is empty or incorrectly formatted.")
            sys.exit()

        length = len(next(iter(prot_dict.values())).seq)
        keep_mask = self.create_mask(length * 3)

        for gene_id, protein_seq_record in prot_dict.items():
            try:
                p_seq = protein_seq_record.seq
                n_seq = nucl_records[gene_id].seq

                normalized_p_seq = self.normalize_p_seq(p_seq)
                normalized_n_seq = self.normalize_n_seq(n_seq, normalized_p_seq)

                sequence = [
                    normalized_n_seq[idx] if c not in "-?*Xx" and keep_mask[idx] else "-"
                    for idx, c in enumerate(normalized_p_seq)
                    if keep_mask[idx]
                ]

                if self.remove_stop_codon and normalized_p_seq[-1] == "*":
                    sequence[-3:] = normalized_n_seq[-3:]

                pal2nal[gene_id] = ''.join(sequence)

            except KeyError:
                print(f"Nucleotide sequence for {gene_id} not found.")
                sys.exit()

        return pal2nal
