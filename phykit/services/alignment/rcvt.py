from collections import Counter

from .base import Alignment


class RelativeCompositionVariabilityTaxon(Alignment):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        alignment, _, _ = self.get_alignment_and_format()
        aln_len = alignment.get_alignment_length()
        num_records = len(alignment)

        concat_seq = "".join(str(record.seq) for record in alignment)
        total_counts = Counter(concat_seq)

        average_d = {
            char: total_counts[char] / num_records for char in total_counts
        }

        for record in alignment:
            record_counts = Counter(record.seq)
            temp_rcv = \
                sum(
                    abs(
                        record_counts[seq_letter] - average_d[seq_letter]
                        ) for seq_letter in total_counts
                )
            rcv_value = temp_rcv / (num_records * aln_len)
            print(f"{record.id}\t{round(rcv_value, 4)}")

    def process_args(self, args):
        return dict(alignment_file_path=args.alignment)
