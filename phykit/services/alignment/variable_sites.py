from Bio.Align import MultipleSeqAlignment

from .base import Alignment

class VariableSites(Alignment):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        alignment, alignment_format = self.get_alignment_and_format()
        var_sites, aln_len, var_sites_per = self.calculate_variable_sites(alignment)
        if (var_sites, aln_len, var_sites_per):
            print(f"{var_sites}\t{aln_len}\t{round(var_sites_per, 4)}")


    def process_args(self, args):
        return dict(alignment_file_path=args.alignment)

    def calculate_variable_sites(self, alignment):
        aln_len = alignment.get_alignment_length()
        
        var_sites = 0
        # count number of variable sites
        for i in range(0, aln_len, int(1)):
            # obtain sequence at position, remove gaps, and make
            # all characters uppercase
            seq_at_position = ""
            seq_at_position += alignment[:, i]
            seq_at_position = seq_at_position.upper().replace("-", "")
            if len(set(seq_at_position)) > 1:
                var_sites += 1
        # calculate percent of variable sites
        var_sites_per = (var_sites / aln_len)*100
        
        return var_sites, aln_len, var_sites_per