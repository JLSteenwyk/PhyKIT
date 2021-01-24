import itertools
import statistics as stat

from Bio.Align import MultipleSeqAlignment
import numpy as np

from .base import Alignment
from ...helpers.stats_summary import calculate_summary_statistics_from_dict, print_summary_statistics

class PairwiseIdentity(Alignment):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        # get aln
        alignment, alignment_format = self.get_alignment_and_format()
        
        # get entry indices
        entries = self.get_entry_indices(alignment)

        # determine pairwise combinations of entry indices
        combos = list(itertools.combinations(entries, 2))

        pairwise_identities, stats = self.calculate_pairwise_identities(alignment, combos)
        
        if self.verbose:
            try:
                for pair, identity in pairwise_identities.items():
                    print(f"{pair}\t{identity}")
            except BrokenPipeError:
                pass
        else:
            print_summary_statistics(stats)

    def process_args(self, args):
        return dict(alignment_file_path=args.alignment, verbose=args.verbose)

    def get_entry_indices(self, alignment) -> list:
        entries = []
        entries_count = 0
        for record in alignment:
            entries.append(entries_count)
            entries_count+=1

        return entries

    def calculate_pairwise_identities(self, alignment, combos):
        # get aln length
        aln_len = alignment.get_alignment_length()

        # determine pairwise identity of each combination
        pairwise_identities = {}
        for combo in combos:
            identities = 0
            seq_one = alignment[combo[0]].seq
            seq_two = alignment[combo[1]].seq
            for idx in range(0, aln_len):
                if seq_one[idx] == seq_two[idx]:
                    identities += 1
            ids = alignment[combo[0]].id + '-' + alignment[combo[1]].id
            pairwise_identities[ids] = identities / aln_len

        stats = calculate_summary_statistics_from_dict(pairwise_identities)

        return pairwise_identities, stats