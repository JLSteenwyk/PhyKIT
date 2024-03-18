from .base import Alignment


class EvolutionaryRatePerSite(Alignment):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        alignment, _ = self.get_alignment_and_format()
        pic_values = self.calculate_evolutionary_rate_per_site(
            alignment
        )
        for idx, value in enumerate(pic_values):
            print(f"{idx+1}\t{round(value,4)}")

    def process_args(self, args):
        return dict(alignment_file_path=args.alignment)

    def get_number_of_occurences_per_character(self, alignment, idx: int) -> dict:
        # obtain sequence at position, remove gaps, and make
        # all characters uppercase
        seq_at_position = ""
        seq_at_position += alignment[:, idx]
        seq_at_position = seq_at_position.upper().replace("-", "")
        num_occurences = {}
        for char in set(seq_at_position.replace("-", "")):
            num_occurences[char] = seq_at_position.count(char)

        return num_occurences

    def calculate_pic(self, num_occurences: dict):
        sum_of_frequencies = 0
        total_frequencies = sum(num_occurences.values())
        for _, frequency in num_occurences.items():
            sum_of_frequencies += ((frequency/total_frequencies) ** 2)
        return 1 - sum_of_frequencies

    def calculate_evolutionary_rate_per_site(self, alignment):
        # get aln length
        aln_len = alignment.get_alignment_length()

        pic_values = []

        # count number of parsimony informative sites
        for idx in range(0, aln_len, int(1)):
            # count occurneces of each character at site idx
            num_occurences = self.get_number_of_occurences_per_character(
                alignment, idx
            )

            pic_values.append(self.calculate_pic(num_occurences))

        return pic_values
