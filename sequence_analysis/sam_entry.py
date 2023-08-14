"""
SamEntry class takes in a string from a .sam file, and determines various scores
and alignment properties.

I wrote this to deal with a standard STAR output. Not robust by any means.
"""
from sequence_analysis.sequence import sequence

class SamEntry:
    def __init__(self, bam_string):
        self.bam_string = bam_string

        self.q_name = None
        self.bitwise_flag = None
        self.r_name = None
        self.pos = None
        self.mapping_quality = None
        self.cigar_string = None
        self.r_next = None
        self.p_next = None
        self.t_len = None
        self.sequence_string = None
        self.quality_33 = None
        self.other_flags = None
        self.parse_data()

        self.sequence = None
        self.generate_sequence_object()

        self.probability_of_incorrect_read = None
        self.generate_probability_of_incorrect_read()

    def parse_data(self):
        """
        Given a BamEntry object with a bam_string, parse the data.
        """
        split = self.bam_string.split()
        self.q_name = split[0]
        self.bitwise_flag = int(split[1])
        self.r_name = split[2]
        # bring back to 0-index
        self.pos = int(split[3]) - 1
        self.mapping_quality = int(split[4])
        self.cigar_string = split[5]
        self.r_next = split[6]
        self.p_next = split[7]
        self.t_len = split[8]
        self.sequence_string = split[9]
        self.quality_33 = split[10]
        self.other_flags = [c for c in split[11:]]

    def generate_sequence_object(self):
        """
        Generate a sequence object from the sequence string.
        """
        self.sequence = sequence(self.sequence_string)
        # check if we need to reverse complement
        if len(bin(self.bitwise_flag)[2:]) > 5:
            if bin(self.bitwise_flag)[2:][::-1][4] == '1':
                self.sequence = self.sequence.reverse_complement()

        self.sequence.name = self.q_name

    def generate_probability_of_incorrect_read(self):
        """
        self.quality_33 gives the Illumina ASCII quality+33 values.

        The probability of an incorrect read is given by P = 10^(-Q/10),
        where Q is a quality value from 0 to 42(?).
        """
        Q_values = [ord(c)-33 for c in self.quality_33]
        self.probability_of_incorrect_read = [10**(-1*q/10) for q in Q_values]
