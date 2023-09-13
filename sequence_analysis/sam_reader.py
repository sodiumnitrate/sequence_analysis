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

    def get_unique_reference_names(self):
        """
        Give a list of all the references the queries were mapped to.
        Also return the minimum range that contains all mappings per reference.
        """
        if len(self.sam_string_list) == 0:
            print("ERROR: there are no sam entries.")
            raise ValueError

        ref_names = {}
        for read in self.sam_string_list:
            name = read.split()[2]
            pos = int(read.strip().split()[3])
            length = len(read.strip().split()[9])
            end = pos + length
            if name not in ref_names:
                ref_names[name] = (pos, end)
            else:
                ref_names[name] = (min(ref_names[name][0], pos), max(ref_names[name][1], end))

        return ref_names

    def normalize_AS_and_filter(self, min_score):
        """
        Sometimes, AS is given as a raw score. We need to convert these to percentages.
        """
        if len(self.full_lengths) == 0:
            print("ERROR: full sequence lengths are not set.")
            raise ValueError

        if len(self.sam_string_list) == 0:
            print("ERROR: there are no Sam entries read.")
            raise ValueError

        new_reads = []
        for read in self.sam_string_list:
            AS = float(read.strip().split("AS:i:")[1].split()[0])
            name = read.strip().split()[0]
            seqlen = self.full_lengths[name]
            norm_score = AS / (2 * seqlen) * 100 # percent
            if norm_score >= min_score:
                new_reads.append(read)

        self.sam_string_list = new_reads

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

        r = self.sam_string_list.pop()
        print(r)
        r = r.split('-')
        self.seq_start = int(r[0])
        self.seq_end = int(r[1])