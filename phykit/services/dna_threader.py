from Bio import SeqIO, SeqRecord

from phykit.services.base import BaseService


class DNAThreader(BaseService):
    """
    Threads DNA on top of protein alignment
    """

    def __init__(self, args) -> None:
        self.process_args(args)

    def process_args(self, args):
        self.include_stop_codon = args.stop
        self.protein_file_path = args.protein
        self.nucleotide_file_path = args.nucleotide

    def run(self):
        prot = self.read(self.protein_file_path)
        nucl = self.read(self.nucleotide_file_path)

        pal2nal = self.thread(prot, nucl)

        # print out threaded DNA alignment
        for entry in range(0, len(PROTa)):
            print(">{}\n{}".format(PROTa[entry].id, pal2nal[PROTa[entry].id]))

    def read_file(self, file_path: str, file_format: str = "fasta") -> SeqRecord:
        return SeqIO.parse(file_path, file_format)

    def print_threaded_alignment(self, pal2nal: dict, protein: SeqRecord) -> None:
        for protein_seq_record in protein:
            gene_id = protein_seq_record.id
            print(f">{gene_id}\n{pal2nal[gene_id]}")

    def thread(self, protein: SeqRecord, nucleotide: SeqRecord) -> dict:
        # protein alignment to nucleotide alignment
        pal2nal = {}

        for protein_seq_record, nucleotide_seq_record in zip(protein, nucleotide):
            gene_id = protein_seq_record.id

            # save protein sequence to p_seq
            p_seq = protein_seq_record.seq

            # save nucleotide sequence to n_seq
            n_seq = nucleotide_seq_record.seq

            pal2nal[gene_id] = ""
            gap_count = 0

            # loop through the sequence
            for AA in range(0, (int(len(p_seq)) + 1) - 1, 1):
                if self.include_stop_codon:
                    # if AA is a gap insert a codon of gaps
                    if p_seq[AA] == "-":
                        pal2nal[gene_id] += "---"
                        gap_count += 1
                    # if AA is not a gap, insert the corresponding codon
                    elif p_seq[AA] != "-":
                        NTwin = (AA - gap_count) * 3
                        pal2nal[gene_id] += n_seq[NTwin : NTwin + 3]
                else:
                    # if AA is a gap insert a codon of gaps
                    if p_seq[AA] == "-":
                        pal2nal[gene_id] += "---"
                        gap_count += 1
                    # if AA is not a gap, insert the corresponding codon
                    elif p_seq[AA] != "-":
                        # if AA is a stop or ambiguous insert a codon of gaps
                        if p_seq[AA] == "X" or p_seq[AA] == "*":
                            pal2nal[gene_id] += "---"
                        else:
                            NTwin = (AA - gap_count) * 3
                            pal2nal[gene_id] += n_seq[NTwin : NTwin + 3]

                # TODO: do we need to port this?
                ## this commented code will check if the nucleotide window
                ## translates to the corresponding codon
                # if Pseq[AA] != Nseq[NTwin:NTwin+3].translate():
                #     print("\nAmino acid position", AA, "(",Pseq[AA],")", "does not correspond to codon")
                #     print(Nseq[NTwin:NTwin+3], "in nucleotide window", NTwin,"-",NTwin+3)
                #     print(Nseq[NTwin:NTwin+3], "translates to", Nseq[NTwin:NTwin+3].translate())
                #     print("Nucleotides cannot be threaded ontop of the protein sequence.")
                #     print("Exiting now...\n")
                #     sys.exit()

        return pal2nal
