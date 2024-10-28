import itertools
from typing import Dict, List, Tuple

from Bio.Align import MultipleSeqAlignment

from .base import Alignment
from ...helpers.stats_summary import (
    calculate_summary_statistics_from_dict,
    print_summary_statistics,
)


class PairwiseIdentity(Alignment):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        alignment, _, is_protein = self.get_alignment_and_format()

        pair_ids, pairwise_identities, stats = \
            self.calculate_pairwise_identities(
                alignment, self.exclude_gaps
            )

        if self.verbose:
            try:
                for pair, identity in zip(
                    pair_ids, pairwise_identities.values()
                ):
                    print(f"{pair[0]}\t{pair[1]}\t{round(identity, 4)}")
            except BrokenPipeError:
                pass
        else:
            print_summary_statistics(stats)

    def process_args(self, args) -> Dict[str, str]:
        return dict(
            alignment_file_path=args.alignment,
            verbose=args.verbose,
            exclude_gaps=args.exclude_gaps,
        )

    def calculate_pairwise_identities(
        self,
        alignment: MultipleSeqAlignment,
        exclude_gaps: bool,
    ) -> Tuple[List[List[str]], Dict[str, float], Dict[str, float]]:
        gap_chars = self.get_gap_chars()

        pairwise_identities = {}
        pair_ids = []

        for idx1, idx2 in itertools.combinations(range(len(alignment)), 2):
            seq_one = alignment[idx1].seq
            seq_two = alignment[idx2].seq
            identities = 0
            total_compared = 0

            for res_one, res_two in zip(seq_one, seq_two):
                res_one = res_one.upper()
                res_two = res_two.upper()

                total_compared += 1

                if exclude_gaps:
                    if (res_one not in gap_chars or res_two not in gap_chars):
                        if res_one == res_two:
                            identities += 1
                else:
                    if res_one == res_two:
                        identities += 1

            if total_compared > 0:
                identity_score = identities / total_compared
            else:
                identity_score = 0

            pair_id = [alignment[idx1].id, alignment[idx2].id]
            pair_ids.append(pair_id)
            pairwise_identities["-".join(pair_id)] = identity_score

        stats = calculate_summary_statistics_from_dict(pairwise_identities)

        return pair_ids, pairwise_identities, stats
