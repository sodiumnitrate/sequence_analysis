"""
This file holds the sequence set (seq_set) class and attributed methods.
Includes I/O of .fasta files.
"""


class seq_set:
    def __init__(self,list_of_sequences=None):
        self.records = list_of_sequences

    def write_fasta(self,file_name):
        for seq in self.records:
            pass
