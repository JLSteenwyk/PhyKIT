import statistics as stat

from Bio import SeqIO, SeqRecord
from Bio.Seq import Seq

from .base import Alignment


class CreateConcatenationMatrix(Alignment):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        alignment_list = self.alignment_list
        prefix = self.prefix
        
        self.create_concatenation_matrix(alignment_list, prefix)

    def process_args(self, args):
        return dict(
            alignment_list=args.alignment_list,
            prefix=args.prefix
        )

    def create_concatenation_matrix(self, alignment_list, prefix):
        # read lists into python lists
        alignments  = [line.rstrip('\n') for line in open(alignment_list)]

        # get taxa names
        taxa = []
        for alignment in alignments:
            alignment = SeqIO.parse(alignment, "fasta")
            for indiv in alignment:
                if indiv.name not in taxa:
                    taxa.append(indiv.name)

        # assignment the output file names
        file_partition = prefix + ".partition"
        fasta_output   = prefix + ".fa"
        file_occupancy = prefix + ".occupancy"

        # print log message
        print("\n"+"-"*(len("- General features -")))
        print("| General features |")
        print("-"*(len("- General features -")))
        print("Total number of taxa:", len(taxa))
        print("Total number of alignments:", len(alignments))
        print("\n\n"+"-"*(len("- Output files -")))
        print("| Output files |")
        print("-"*(len("- Output files -")))
        print("partition file output:", file_partition)
        print("concatenated fasta output:", fasta_output)
        print("occupancy report:", file_occupancy)

        # # assigning placeholders for lengths
        first_len    = 1
        second_len   = 0

        # create a dictionary with keys as taxa and empty list values
        concat = {key: [] for key in taxa}

        # string for first field of partition file
        field_one = 'AUTO'

        # create an empty file_partition
        open(file_partition, "w")
        # create an empty occupancy per OG file
        open(file_occupancy, "w")
        # loop through USCO trimal files
        for fasta in alignments:
            # list for keeping track of taxa per USCO
            og_taxa    = []
            # open alignment file
            records = list(SeqIO.parse(fasta, "fasta"))

            # create lists of USCOtaxa and missing_taxa
            for record in records:
                # append taxa per USCO to USCOtaxa
                og_taxa.append(record.id)
                
            # initialize list and determine missing taxa
            missing_taxa = []
            missing_taxa = list(set(taxa)-set(og_taxa))

            # create string for missing seq
            og_len     = ''
            missing_seq = ''
            og_len   = len(records[0].seq)
            missing_seq = og_len*'?'

            # create missing taxa sequence based on either it is protein or nucleotide
            for indiv in missing_taxa:
                indiv_record = SeqRecord(Seq(missing_seq), 
                    id = indiv, name = indiv, description = indiv)
                concat[indiv].append(indiv_record.seq)        

            # create record for present data
            for record in records:
                if record.id in og_taxa:
                    concat[record.id].append(record.seq)

            # create partition file
            with open(file_partition, "a") as f:
                # second value in partition file
                second_len += og_len
                entry = field_one+", "+str(fasta)+"="+str(first_len)+"-"+str(second_len)+"\n"
                f.write(str(entry))
                # add to first value for partition file
                first_len += og_len

            # create partition file
            with open(file_occupancy, "a") as f:
                # determine number of present taxa
                num_present  = len(og_taxa)
                # determine number of missing taxa
                num_missing  = len(missing_taxa)
                # determine percent occupied
                percent_occupancy = num_present/(num_present+num_missing)
                if num_missing == 0:
                    missing_taxa = ["None"]
                entry = str(fasta)+'\t'+str(num_present)+'\t'+str(num_missing)+'\t'+str(percent_occupancy)+'\t'+';'.join(missing_taxa)+'\t'+';'.join(og_taxa)+"\n"
                f.write(str(entry))

        # join seqs of genes in value (list) and write to concat_fa   
        with open(fasta_output, "w") as final_fasta_file: 
            for x in concat:
                concatenated = Seq("")
                for s in concat[x]:
                    concatenated += s
                concat[x] = concatenated
                entry = '>'+x+"\n"+concat[x]+"\n"
                final_fasta_file.write(str(entry))             

        # print log message
        print("Complete!\n")
            
            

