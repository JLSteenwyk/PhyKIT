from enum import Enum

from Bio import AlignIO
from ..base import BaseService

class FileFormat(Enum):
    fasta = "fasta"
    clustal = "clustal"
    maf = "maf"
    mauve = "mauve"
    phylip = "phylip"
    phylip_seq = "phylip-sequential"
    phylip_rel = "phylip-relaxed"
    stockholm = "stockholm"

class Alignment(BaseService):
    def __init__(
        self,
        *args,
        alignment_file_path=None,
        output_file_path=None,
        protein_file_path=None,
        nucleotide_file_path=None,
        alignment_list=None,
        prefix=None
    ):
        self.alignment_file_path = alignment_file_path
        self.output_file_path = output_file_path
        self.protein_file_path = protein_file_path,
        self.nucleotide_file_path = nucleotide_file_path 
        self.alignment_list = alignment_list
        self.prefix = prefix

    def get_alignment_and_format(self):
        """
        automatic file type determination
        """

        # if file format is provided, read the file according to the user's file format
        for fileFormat in FileFormat:
            try:
                alignment = AlignIO.read(open(self.alignment_file_path), fileFormat.value)
                return alignment, fileFormat.value
            # the following exceptions refer to skipping over errors
            # associated with reading the wrong input file
            except ValueError:
                continue
            except AssertionError:
                continue

        raise Exception("Input file could not be read")

    def calculate_alignment_length(self, alignment):
        """
        calculate alignment length
        """
        aln_len = alignment.get_alignment_length()
        return aln_len

