"""
This file holds the sequence class and related methods.
"""
import re
from Bio.Seq import Seq
import pdb
from sequence_analysis.utils import check_for_pattern

class sequence:
    # TODO: inherit from the Seq class?
    def __init__(self,seq,name=None,type=None):
        self.seq = seq.upper()
        self.name = name
        # TODO: assign type in seq_set creation
        self.type = type

    def frame_shift(self, frame=0):
        if self.type != 'rna' and self.type != 'dna':
            print(f"ERROR: type {self.type} is not meaningful for frame shifting. Please provide RNA or DNA sequences.")
            return None

        seq = self.seq

        # returns frame shifted sequence
        if frame == 0:
            pass

        if frame == 1:
            seq = seq[1:]

        if frame == 2:
            seq = seq[2:]

        N_to_be_added = 3 - len(seq) % 3
        seq = seq + "N" * N_to_be_added

        assert(len(seq) % 3 == 0)

        return seq

    def six_frame_check(self, regex):
        # does a 6-frame translation and checks for the existence of a given regex pattern
        selected = False
        true_sequence = ""
        for frame in range(3):
            for order in range(1):
                if order == 0:
                    s = Seq(self.frame_shift(frame=frame)).translate()
                    print(s)
                else:
                    s = Seq(self.frame_shift(frame=frame)).reverse_complement().translate()
                    print(s)

                # NOTE: the following implicitly assumes that only one reading frame is valid
                selected = check_for_pattern(s,regex)
                if selected:
                    true_sequence = str(s)

        return selected, true_sequence
