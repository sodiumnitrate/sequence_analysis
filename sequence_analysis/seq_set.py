"""
This file holds the sequence set (seq_set) class and attributed methods.
Includes I/O of .fasta files.
"""
from Bio import SeqIO


class seq_set:
    def __init__(self,list_of_sequences=[]):
        self.records = list_of_sequences
        self.names = None

    def write_fasta(self,file_name):
        # function to write all sequences within the list to a fasta file
        f = open(file_name,'w')
        for i,seq in self.records:
            s = seq.seq
            if self.names is not None:
                record_name = self.names[i]
            else:
                record_name = f"record|{i}"
            name = "> " + record_name 
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

