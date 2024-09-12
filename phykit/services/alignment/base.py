from collections import Counter
import sys

from typing import List

from ..base import BaseService
from ...helpers.files import (
    get_alignment_and_format as get_alignment_and_format_helper
)


class Alignment(BaseService):
    def __init__(
        self,
        *args,
        alignment_file_path=None,
        code=None,
        fasta=None,
        output_file_path=None,
        protein_file_path=None,
        nucleotide_file_path=None,
        alignment_list_path=None,
        prefix=None,
        idmap=None,
        reference=None,
        verbose=None,
        entry=None,
        exclude_gaps=None,
    ):
        self.alignment_file_path = alignment_file_path
        self.code = code,
        self.output_file_path = output_file_path
        self.protein_file_path = (protein_file_path,)
        self.nucleotide_file_path = nucleotide_file_path
        self.alignment_list_path = alignment_list_path
        self.prefix = prefix
        self.fasta = fasta
        self.idmap = idmap
        self.reference = reference
        self.verbose = verbose
        self.entry = entry
        self.exclude_gaps = exclude_gaps

    def get_alignment_and_format(self):
        """
        automatic file type determination
        """
        try:
            return get_alignment_and_format_helper(self.alignment_file_path)
        except FileNotFoundError:
            print("Input corresponds to no such file or directory.")
            print("Please double check pathing and filenames")
            sys.exit()

    def calculate_rcv(self) -> float:
        alignment, _, _ = self.get_alignment_and_format()
        aln_len = alignment.get_alignment_length()

        concat_seq = []
        num_records = len(alignment)

        concat_seq = "".join(str(record.seq) for record in alignment)

        total_counts = Counter(concat_seq)

        average_d = {
            seq: total_counts[seq] / num_records for seq in total_counts
        }

        indiv_rcv_values = []

        for record in alignment:
            record_counts = Counter(record.seq)
            temp_rcv = sum(
                abs(
                    record_counts[seq_letter] - average_d[seq_letter]
                ) for seq_letter in total_counts
            )
            indiv_rcv_values.append(temp_rcv / (num_records * aln_len))

        relative_composition_variability = sum(indiv_rcv_values)

        return relative_composition_variability

    def get_gap_chars(is_protein: bool) -> List[str]:
        if is_protein:
            return ["-", "?", "*", "X", "x"]
        else:
            return ["-", "?", "*", "X", "x", "N", "n"]
