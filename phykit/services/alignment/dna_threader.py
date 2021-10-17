import sys

from Bio import SeqIO, SeqRecord

from .base import Alignment


class DNAThreader(Alignment):
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
        prot = self.read_file(self.protein_file_path)
        nucl = self.read_file(self.nucleotide_file_path)

        pal2nal = self.thread(prot, nucl)

        for record in pal2nal:
            sequence = ''.join(pal2nal[record])
            print(f">{record}")
            print(f"{sequence}")

    def read_file(
        self,
        file_path: str,
        file_format: str = "fasta"
    ) -> SeqRecord:
        return SeqIO.parse(file_path, file_format)

    def thread(
        self,
        protein: SeqRecord,
        nucleotide: SeqRecord
    ) -> dict:
        # protein alignment to nucleotide alignment
        pal2nal = {}

        try:
            for protein_seq_record, nucleotide_seq_record in zip(protein, nucleotide):
                # get gene id and sequences
                gene_id, p_seq, n_seq = self.get_id_and_seqs(protein_seq_record, nucleotide_seq_record)

                # initialize gap counter and gene in pal2nal dict
                gap_count = 0
                pal2nal[gene_id] = []

                # loop through the sequence
                for aa_idx in range(0, len(p_seq), 1):
                    # if AA is a gap insert a codon of gaps
                    if p_seq[aa_idx] == "-":
                        pal2nal = self.add_gap(pal2nal, gene_id)
                        gap_count += 1
                    else:
                        if self.include_stop_codon:
                            # if AA is not a gap, insert the corresponding codon
                            if p_seq[aa_idx] != "-":
                                pal2nal = self.add_codon(
                                    aa_idx,
                                    gap_count,
                                    n_seq,
                                    gene_id,
                                    pal2nal
                                )
                        else:
                            # if AA is a stop or ambiguous insert a codon of gaps
                            if p_seq[aa_idx] == "X" or p_seq[aa_idx] == "*":
                                pal2nal = self.add_gap(pal2nal, gene_id)
                            else:
                                pal2nal = self.add_codon(
                                    aa_idx,
                                    gap_count,
                                    n_seq,
                                    gene_id,
                                    pal2nal
                                )
            return pal2nal
        except FileNotFoundError:
            try:
                print("One input corresponds to no such file or directory.")
                print("Please double checking pathing and filenames")
                sys.exit()
            except BrokenPipeError:
                pass

    def get_id_and_seqs(
        self,
        protein_seq_record: SeqRecord,
        nucleotide_seq_record: SeqRecord
    ):
        """
        get gene id, protein sequence, and nucleotide sequence
        """

        # save gene id to gene_id
        gene_id = protein_seq_record.id

        # save protein sequence to p_seq
        p_seq = protein_seq_record.seq

        # save nucleotide sequence to n_seq
        n_seq = nucleotide_seq_record.seq

        return gene_id, p_seq, n_seq

    def add_gap(
        self,
        pal2nal: dict,
        gene_id: str
    ) -> dict:
        """
        add a gap to the growing sequence
        """

        pal2nal[gene_id].append("---")

        return pal2nal

    def add_codon(
        self,
        aa_idx: int,
        gap_count: int,
        n_seq: SeqRecord,
        gene_id: str,
        pal2nal: dict
    ) -> dict:
        """
        add a gap to the growing sequence
        """

        nt_window = (aa_idx - gap_count) * 3
        seq = n_seq[nt_window : nt_window + 3]._data
        seq = seq.decode("utf-8")
        pal2nal[gene_id].append(seq)

        return pal2nal





