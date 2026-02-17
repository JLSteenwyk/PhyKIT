import sys
from typing import Dict
import numpy as np

from Bio.Align import MultipleSeqAlignment

from .base import Alignment
from ...helpers.files import get_alignment_and_format
from ...helpers.json_output import print_json


class GCContent(Alignment):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(fasta=parsed["fasta"], verbose=parsed["verbose"])
        self.json_output = parsed["json_output"]

    def run(self):
        records, _, is_protein = get_alignment_and_format(self.fasta)

        if is_protein:
            print("GC content can't be calculated for protein sequences")
            sys.exit(2)

        if self.json_output:
            if self.verbose:
                rows = self.calculate_gc_per_sequence_data(records, is_protein)
                row_payload = [
                    dict(taxon=taxon, gc_content=round(gc_content, 4))
                    for taxon, gc_content in rows
                ]
                print_json(
                    dict(
                        verbose=True,
                        rows=row_payload,
                        sequences=row_payload,
                    )
                )
            else:
                print_json(
                    dict(
                        verbose=False,
                        gc_content=self.calculate_gc_total_value(records, is_protein),
                    )
                )
            return

        if self.verbose:
            self.calculate_gc_per_sequence(records, is_protein)
        else:
            self.calculate_gc_total(records, is_protein)

    def process_args(self, args) -> Dict[str, str]:
        return dict(
            fasta=args.fasta,
            verbose=args.verbose,
            json_output=getattr(args, "json", False),
        )

    def calculate_gc_per_sequence_data(
        self, records: MultipleSeqAlignment, is_protein: bool = False
    ) -> list[tuple[str, float]]:
        gap_chars = set(self.get_gap_chars(is_protein))
        output = []
        for record in records:
            seq_arr = np.array(list(str(record.seq).upper()), dtype='U1')
            non_gap_mask = ~np.isin(seq_arr, list(gap_chars))
            cleaned_seq = seq_arr[non_gap_mask]

            if len(cleaned_seq) > 0:
                gc_count = np.sum((cleaned_seq == 'G') | (cleaned_seq == 'C'))
                gc_content = gc_count / len(cleaned_seq)
            else:
                gc_content = 0
            output.append((record.id, float(gc_content)))
        return output

    def calculate_gc_per_sequence(
        self, records: MultipleSeqAlignment, is_protein: bool = False
    ) -> None:
        for record_id, gc_content in self.calculate_gc_per_sequence_data(
            records, is_protein
        ):
            try:
                print(f"{record_id}\t{round(gc_content, 4)}")
            except BrokenPipeError:
                pass

    def calculate_gc_total_value(
        self, records: MultipleSeqAlignment, is_protein: bool = False
    ) -> float:
        gap_chars = set(self.get_gap_chars(is_protein))

        # Combine all sequences into one array
        all_seqs = [list(str(record.seq).upper()) for record in records]
        combined_arr = np.concatenate([np.array(seq, dtype='U1') for seq in all_seqs])

        # Filter out gaps
        non_gap_mask = ~np.isin(combined_arr, list(gap_chars))
        cleaned_seq = combined_arr[non_gap_mask]

        if len(cleaned_seq) > 0:
            gc_count = np.sum((cleaned_seq == 'G') | (cleaned_seq == 'C'))
            return round(gc_count / len(cleaned_seq), 4)

        print(
            "Input file has an unacceptable format. Please check input file argument."
        )
        sys.exit(2)

    def calculate_gc_total(
        self, records: MultipleSeqAlignment, is_protein: bool = False
    ) -> None:
        print(self.calculate_gc_total_value(records, is_protein))
