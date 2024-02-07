from .base import Alignment


class RelativeCompositionVariabilityTaxon(Alignment):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        alignment, _ = self.get_alignment_and_format()
        aln_len = alignment.get_alignment_length()

        # string to hold all sequences
        concat_seq = ""
        # initialize a counter for the number of sequences in the input fasta file
        num_records = 0

        # for each record join concatSeq string and sequence as well as keeping track
        # of the number of records
        for record in alignment:
            concat_seq += record.seq
            num_records += 1

        # dictionary to hold the average occurence of each sequence letter
        average_d = {}
        # loop through the different sequences that appear in the fasta file
        # population dictionary with how many times that sequence appears
        for seq in set(concat_seq):
            average_d[seq] = concat_seq.count(seq) / num_records

        # intiailize list to hold the RCV values per ith taxa
        # that will later be summed
        indiv_rcv_values = []

        # loop through records again and calculate RCV for
        # each taxa and append to indivRCVvalues
        for record in alignment:
            # temp holds a temporary value of the numerator before appending
            # to numeratorRCVvalues and then is reassigned to 0 when it goes
            # through the loop again
            temp = 0
            # calculates the absolute value of the ith sequence letter minus the average
            for seq_letter in set(concat_seq):
                temp += abs(
                    record.seq.count(seq_letter) - average_d[seq_letter]
                )
            indiv_rcv_values.append(temp / (num_records * aln_len))
            print(f"{record.id}\t{round((temp / (num_records * aln_len)),4)}")

        # relative_composition_variability = sum(indiv_rcv_values)

        # print the sum of all RCV values
        # print(round(relative_composition_variability, 4))

    def process_args(self, args):
        return dict(alignment_file_path=args.alignment)
