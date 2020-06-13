import itertools
import statistics as stat

from Bio.Align import MultipleSeqAlignment
import numpy as np

from .base import Alignment

class PairwiseIdentity(Alignment):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        alignment, alignment_format = self.get_alignment_and_format()
        pairwise_identities, stats = self.calculate_pairwise_identities(alignment)
        if self.verbose:
            for pair, identity in pairwise_identities.items():
                print(f"{pair}\t{identity}")
        else:
            print(f"mean: {stats['mean']}")
            print(f"median: {stats['median']}")
            print(f"25th percentile: {stats['twenty_fifth']}")
            print(f"75th percentile: {stats['seventy_fifth']}")
            print(f"minimum: {stats['minimum']}")
            print(f"maximum: {stats['maximum']}")
            print(f"standard deviation: {stats['standard_deviation']}")
            print(f"variance: {stats['variance']}")

    def process_args(self, args):
        return dict(alignment_file_path=args.alignment, verbose=args.verbose)

    def calculate_pairwise_identities(self, alignment):
        # get aln length
        aln_len = alignment.get_alignment_length()
        
        # get entry indices
        entries = []
        entries_count = 0
        for record in alignment:
            entries.append(entries_count)
            entries_count += 1

        # determine pairwise combinations of entry indices
        combos = list(itertools.combinations(entries, 2))

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

        stats = dict(
            mean=stat.mean([*pairwise_identities.values()]),
            median=stat.median([*pairwise_identities.values()]),
            twenty_fifth=np.percentile([*pairwise_identities.values()], 25),
            seventy_fifth=np.percentile([*pairwise_identities.values()], 75),
            minimum=np.min([*pairwise_identities.values()]),
            maximum=np.max([*pairwise_identities.values()]),
            standard_deviation=stat.stdev([*pairwise_identities.values()]),
            variance=stat.variance([*pairwise_identities.values()])
        )

        return pairwise_identities, stats