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
        self.remove_stop_codon = args.stop
        self.protein_file_path = args.protein
        self.nucleotide_file_path = args.nucleotide
        self.clipkit_log_file = args.clipkit_log_file

    @property
    def clipkit_log_data(self):
        data = None
        if self.clipkit_log_file:
            with open(self.clipkit_log_file) as f:
                data = [
                    line.rstrip('\n').split(" ") for line in f.readlines()
                ]
        return data

    def run(self):
        prot = self.read_file(self.protein_file_path)

        pal2nal = self.thread(prot)

        for record in pal2nal:
            sequence = "".join(pal2nal[record])
            print(f">{record}")
            print(f"{sequence}")

    def read_file(self, file_path: str, file_format: str = "fasta") -> SeqRecord:
        return SeqIO.parse(file_path, file_format)

    def create_mask(self, length):
        keep_mask = [True] * length
        if self.clipkit_log_data:
            mask = []
            for row in self.clipkit_log_data:
                if row[1] == "keep":
                    mask.append(True)
                    mask.append(True)
                    mask.append(True)
                else:
                    mask.append(False)
                    mask.append(False)
                    mask.append(False)
            keep_mask = mask

        return keep_mask

    def thread(self, protein: SeqRecord) -> dict:
        # protein alignment to nucleotide alignment
        pal2nal = {}

        nucl_records = SeqIO.to_dict(SeqIO.parse(self.nucleotide_file_path, "fasta"))

        prot_dict = SeqIO.to_dict(SeqIO.parse(self.protein_file_path, "fasta"))
        length = len(prot_dict.get(next(iter(prot_dict))))
        keep_mask = self.create_mask(length * 3)

        try:
            for protein_seq_record in protein:
                gene_id = protein_seq_record.id
                p_seq = protein_seq_record.seq
                n_seq = nucl_records[gene_id].seq

                if self.clipkit_log_data:
                    untrimmed = []
                    offset = 0
                    for idx, value in enumerate(keep_mask[::3]):
                        if value is True:
                            untrimmed.append(p_seq[idx - offset])
                        else:
                            offset += 1
                            untrimmed.append("#")
                    p_seq = "".join(untrimmed)

                sequence = []
                transformed = "".join([c * 3 for c in p_seq])

                for idx, c in enumerate(transformed):
                    if not keep_mask[idx]:
                        continue
                    if c not in "-?*Xx":
                        try:
                            n_seq_c = n_seq[idx]
                            sequence.append(n_seq_c)
                        except IndexError:
                            break
                    else:
                        sequence.append("-")

                if self.remove_stop_codon is True:
                    if transformed[-1] == "*":
                        sequence[-1] = n_seq[-1]
                        sequence[-2] = n_seq[-2]
                        sequence[-3] = n_seq[-3]

                pal2nal[gene_id] = "".join(sequence)

            return pal2nal
        except FileNotFoundError:
            try:
                print("One input corresponds to no such file or directory.")
                print("Please double checking pathing and filenames")
                sys.exit()
            except BrokenPipeError:
                pass

    def get_id_and_seqs(
        self, protein_seq_record: SeqRecord, nucleotide_seq_record: SeqRecord
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

    def add_gap(self, pal2nal: dict, gene_id: str) -> dict:
        """
        add a gap to the growing sequence
        """

        pal2nal[gene_id].append("---")

        return pal2nal

    def add_codon(
        self, aa_idx: int, gap_count: int, n_seq: SeqRecord, gene_id: str, pal2nal: dict
    ) -> dict:
        """
        add a gap to the growing sequence
        """

        nt_window = (aa_idx - gap_count) * 3
        seq = n_seq[nt_window:nt_window + 3]._data
        seq = seq.decode("utf-8")
        if not len(seq):
            seq = "---"
        pal2nal[gene_id].append(seq)

        return pal2nal
