"""
Class for reading SAM files while filtering.
"""
import numpy as np
from sequence_analysis.seq_set import seq_set
import sam_reader_cpp

class SamReader:
    def __init__(self, file_name, start=0, end=-1, mapped_onto="", min_score=0):
        self.file_name = file_name
        self.start = start
        self.end = end
        self.mapped_onto = mapped_onto
        self.min_score = min_score

        self.seq_start = self.start
        self.seq_end = self.end

        # if you want to read full lengths of sequences from another file
        self.full_lengths = {}

        self.sam_string_list = []
        self.header = ""

    def get_header(self):
        """
        Read header based on reference name (self.mapped_onto)
        """
        headers = []
        with open(self.file_name, 'r') as f:
            for line in f:
                if line.startswith('@'):
                    headers.append(line)
                else:
                    break

        if self.mapped_onto == "":
            self.header = headers
        else:
            self.header = [a for a in headers if self.mapped_onto in a]

    def read_full_lengths(self, file_name):
        """
        Read from a .fasta file the full lengths of the sequences we will encounter in
        the .sam file.
        """

        sset = seq_set(file_name=file_name)
        for seq in sset:
            self.full_lengths[seq.name] = len(seq)


    def read(self):
        """
        Given all the options, filter and read into a list of strings.
        """

        self.get_header()

        self.sam_string_list = sam_reader_cpp.sam_reader(self.file_name,
                                                         self.start,
                                                         self.end,
                                                         self.mapped_onto,
                                                         self.min_score)

        r = self.sam_string_list.pop().split('-')
        self.seq_start = int(r[0])
        self.seq_end = int(r[1])