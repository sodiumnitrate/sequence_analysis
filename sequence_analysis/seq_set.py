"""
This file holds the sequence set (seq_set) class and attributed methods.
Includes I/O of .fasta files.
"""
from Bio import SeqIO
from sequence_analysis.sequence import sequence
import pdb

class seq_set:
    def __init__(self,list_of_sequences=[],file_name=None):
        self.records = list_of_sequences

        if file_name is not None:
            self.read_fasta(file_name)

    def write_fasta(self,file_name):
        # function to write all sequences within the list to a fasta file
        f = open(file_name,'w')
        for seq in self.records:
            s = seq.seq
            name = "> " + seq.name
            f.write(name)
            f.write('\n')
            for ct, char in enumerate(s):
                if ct % 79 == 0 and ct != 0:
                    f.write('\n')
                f.write(char)
            f.write('\n')
        f.close()


    def read_fasta(self,file_name):
        if len(self.records) > 0:
            print("Warning: overwriting existing data")
        for record in SeqIO.parse(file_name,"fasta"):
            seq = sequence(str(record.seq), record.name)
            self.records.append(seq)

    def add_sequence(self,sequence_obj):
        self.records.append(sequence_obj)

