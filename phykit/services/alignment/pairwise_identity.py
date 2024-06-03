import itertools

from .base import Alignment
from ...helpers.stats_summary import (
    calculate_summary_statistics_from_dict,
    print_summary_statistics,
)


class PairwiseIdentity(Alignment):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        # get aln
        alignment, _ = self.get_alignment_and_format()

        # get entry indices
        entries = self.get_entry_indices(alignment)

        # determine pairwise combinations of entry indices
        combos = list(itertools.combinations(entries, 2))

        pair_ids, pairwise_identities, stats = self.calculate_pairwise_identities(
            alignment, combos, self.exclude_gaps
        )

        if self.verbose:
            try:
                zipped = zip(pairwise_identities.items(), pair_ids)
                for (_, identity), pair in zipped:
                    print(f"{pair[0]}\t{pair[1]}\t{identity}")
            except BrokenPipeError:
                pass
        else:
            print_summary_statistics(stats)

    def process_args(self, args):
        return dict(alignment_file_path=args.alignment, verbose=args.verbose, exclude_gaps=args.exclude_gaps)

    def get_entry_indices(self, alignment) -> list:
        entries = []
        entries_count = 0
        for _ in alignment:
            entries.append(entries_count)
            entries_count += 1

        return entries

    def calculate_pairwise_identities(
        self,
        alignment,
        combos: list,
        exclude_gaps: bool
    ):
        # get aln length
        aln_len = alignment.get_alignment_length()

        # determine pairwise identity of each combination
        pairwise_identities = {}
        pair_ids = []
        for combo in combos:
            identities = 0
            seq_one = alignment[combo[0]].seq
            seq_two = alignment[combo[1]].seq
            for idx in range(0, aln_len):
                if seq_one[idx] == seq_two[idx]:
                    if exclude_gaps:
                        if seq_one[idx] != "-" and seq_two[idx] != "-":
                            identities += 1
                    else:
                        identities += 1
            ids = alignment[combo[0]].id + "-" + alignment[combo[1]].id
            pair_ids.append([alignment[combo[0]].id, alignment[combo[1]].id])
            pairwise_identities[ids] = identities / aln_len

        stats = calculate_summary_statistics_from_dict(pairwise_identities)

        return pair_ids, pairwise_identities, stats
