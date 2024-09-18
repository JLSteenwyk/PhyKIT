from collections import Counter
from multiprocessing import cpu_count
from multiprocessing import Pool
import sys
from typing import Dict, List

from Bio.Seq import Seq

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
        cpu=None,
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
        self.cpu = cpu,
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
        try:
            return get_alignment_and_format_helper(self.alignment_file_path)
        except FileNotFoundError:
            print("Input corresponds to no such file or directory.")
            print("Please double check pathing and filenames")
            sys.exit()

    def calculate_rcv(self) -> float:
        alignment, _, _ = self.get_alignment_and_format()
        aln_len = alignment.get_alignment_length()
        num_records = len(alignment)

        concat_seq = "".join(str(record.seq) for record in alignment)
        total_counts = Counter(concat_seq)

        average_d = {
            seq: total_counts[seq] / num_records for seq in total_counts
        }

        cpu = self.set_cpu()

        with Pool(cpu) as pool:
            indiv_rcv_values = pool.starmap(
                self.calculate_indiv_rcv,
                [
                    (record.seq, average_d, total_counts, aln_len, num_records)
                    for record in alignment
                ]
            )

        relative_composition_variability = sum(indiv_rcv_values)

        return relative_composition_variability

    def calculate_indiv_rcv(
        self,
        seq: Seq,
        average_d: Dict[str, float],
        total_counts: Counter,
        aln_len: int,
        num_records: int,
    ) -> float:
        record_counts = Counter(seq)
        temp_rcv = sum(
            abs(record_counts[seq_letter] - average_d[seq_letter])
            for seq_letter in total_counts
        )
        return temp_rcv / (num_records * aln_len)

    def get_gap_chars(is_protein: bool) -> List[str]:
        if is_protein:
            return ["-", "?", "*", "X", "x"]
        else:
            return ["-", "?", "*", "X", "x", "N", "n"]

    def set_cpu(self) -> int:
        cpu_available = cpu_count()

        cpu = int(self.cpu[0]) if self.cpu[0] else cpu_count()

        if cpu > cpu_available:
            return cpu_available
        else:
            return cpu
