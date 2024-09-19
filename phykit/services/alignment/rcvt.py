from collections import Counter
from multiprocessing import Pool
from typing import Dict, Tuple

from .base import Alignment


class RelativeCompositionVariabilityTaxon(Alignment):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        alignment, _, _ = self.get_alignment_and_format()
        aln_len = alignment.get_alignment_length()
        num_records = len(alignment)

        cpu = self.set_cpu()

        concat_seq = "".join(str(record.seq) for record in alignment)
        total_counts = Counter(concat_seq)

        average_d = {char: total_counts[char] / num_records for char in total_counts}

        with Pool(cpu) as pool:
            results = pool.starmap(
                self.calculate_indiv_rcv,
                [
                    (record.id, record.seq, average_d, total_counts, aln_len, num_records) 
                    for record in alignment
                ]
            )

        for record_id, rcv_value in results:
            print(f"{record_id}\t{round(rcv_value, 4)}")

    def calculate_indiv_rcv(
        self,
        record_id: str,
        sequence: str,
        average_d: Dict[str, float],
        total_counts: Counter,
        aln_len: int,
        num_records: int
    ) -> Tuple[str, float]:
        """
        Calculate the RCV for an individual sequence.
        """
        record_counts = Counter(sequence)
        temp_rcv = sum(
            abs(record_counts[seq_letter] - average_d[seq_letter])
            for seq_letter in total_counts
        )
        rcv_value = temp_rcv / (num_records * aln_len)
        return record_id, rcv_value

    def process_args(self, args) -> Dict[str, str]:
        return dict(alignment_file_path=args.alignment, cpu=args.cpu)
