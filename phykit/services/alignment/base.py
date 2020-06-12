from ..base import BaseService
from ...helpers.files import get_alignment_and_format as get_alignment_and_format_helper

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
        return get_alignment_and_format_helper(self.alignment_file_path)

    def calculate_alignment_length(self, alignment):
        """
        calculate alignment length
        """
        aln_len = alignment.get_alignment_length()
        return aln_len

