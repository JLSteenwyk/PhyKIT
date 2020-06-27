from enum import Enum

from Bio import AlignIO


class FileFormat(Enum):
    fasta = "fasta"
    clustal = "clustal"
    maf = "maf"
    mauve = "mauve"
    phylip = "phylip"
    phylip_seq = "phylip-sequential"
    phylip_rel = "phylip-relaxed"
    stockholm = "stockholm"


def get_alignment_and_format(alignment_file_path: str):
    # if file format is provided, read the file according to the user's file format
    for fileFormat in FileFormat:
        try:
            alignment = AlignIO.read(open(alignment_file_path), fileFormat.value)
            return alignment, fileFormat.value
        # the following exceptions refer to skipping over errors
        # associated with reading the wrong input file
        except ValueError:
            continue
        except AssertionError:
            continue

    raise Exception("Input file could not be read")

def read_single_column_file_to_list(single_col_file_path: str) -> list:
    with open(single_col_file_path) as f:
        return [line.rstrip('\n') for line in f]
