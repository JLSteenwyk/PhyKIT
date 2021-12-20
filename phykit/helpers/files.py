from enum import Enum
import sys

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
        except FileNotFoundError:
            print(f"{alignment_file_path} corresponds to no such file.")
            print("Please check file name and pathing")
            sys.exit()

def read_single_column_file_to_list(single_col_file_path: str) -> list:
    try:
        with open(single_col_file_path) as f:
            return [line.rstrip('\n').strip() for line in f]
    except FileNotFoundError:
        print(f"{single_col_file_path} corresponds to no such file or directory.")
        print("Please check file name and pathing")
        sys.exit()

