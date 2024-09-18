import itertools
from multiprocessing import Pool
from typing import Dict, List, Tuple

from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord

from .base import Alignment
from ...helpers.stats_summary import (
    calculate_summary_statistics_from_dict,
    print_summary_statistics,
)


class PairwiseIdentity(Alignment):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        alignment, _, _ = self.get_alignment_and_format()

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
            cpu=args.cpu,
        )

    def calculate_pairwise_identities(
        self,
        alignment: MultipleSeqAlignment,
        exclude_gaps: bool,
    ) -> Tuple[
        List[List[str]],
        Dict[str, float],
        Dict[str, float],
    ]:

        cpu = self.set_cpu()

        gap_chars = self.get_gap_chars()
        pairs = list(itertools.combinations(range(len(alignment)), 2))

        with Pool(cpu) as pool:
            results = pool.starmap(
                self.calculate_identity_for_pair,
                [
                    (alignment[idx1], alignment[idx2], gap_chars, exclude_gaps)
                    for idx1, idx2 in pairs
                ]
            )

        pair_ids = [
            [alignment[idx1].id, alignment[idx2].id]
            for idx1, idx2 in pairs
        ]
        pairwise_identities = {
            f"{alignment[idx1].id}-{alignment[idx2].id}": score
            for ((idx1, idx2), score) in zip(pairs, results)
        }

        stats = calculate_summary_statistics_from_dict(pairwise_identities)

        return pair_ids, pairwise_identities, stats

    def calculate_identity_for_pair(
        self,
        seq_one_record: SeqRecord,
        seq_two_record: SeqRecord,
        gap_chars: List[str],
        exclude_gaps: bool,
    ) -> float:

        seq_one = seq_one_record.seq
        seq_two = seq_two_record.seq
        identities = 0
        total_compared = 0

        for res_one, res_two in zip(seq_one, seq_two):
            if exclude_gaps and (res_one in gap_chars or res_two in gap_chars):
                continue
            total_compared += 1
            if res_one == res_two:
                identities += 1

        return identities / total_compared if total_compared > 0 else 0
